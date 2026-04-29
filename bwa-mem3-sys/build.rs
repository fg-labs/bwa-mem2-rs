use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

#[allow(clippy::too_many_lines)]
fn main() {
    // Match cc's deployment target to cargo's to avoid SIGBUS on macOS arm64.
    // Without this, cc builds with the host SDK's deployment version (e.g.,
    // macOS 26) while cargo links at a lower version (e.g., macOS 11),
    // causing the dynamic linker to fault on binary start.
    if env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos")
        && env::var("MACOSX_DEPLOYMENT_TARGET").is_err()
    {
        env::set_var("MACOSX_DEPLOYMENT_TARGET", "11.0");
    }

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=shim/bwa_shim.cpp");
    println!("cargo:rerun-if-changed=shim/bwa_shim_align.cpp");
    println!("cargo:rerun-if-changed=shim/bwa_shim.h");
    println!("cargo:rerun-if-changed=shim/bwa_shim_types.h");
    println!("cargo:rerun-if-changed=vendor/COMMIT");
    println!("cargo:rerun-if-changed=patches");

    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let build_dir = out.join("build");

    // 1. Verify MATE_SORT=0 default in vendored Makefile (shim semantics depend on it).
    let makefile_path = manifest.join("vendor/bwa-mem3/Makefile");
    let makefile = fs::read_to_string(&makefile_path)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", makefile_path.display(), e));
    assert!(
        makefile
            .lines()
            .any(|l| l.contains("CPPFLAGS") && l.contains("-DMATE_SORT=0")),
        "vendored bwa-mem3 Makefile must retain `-DMATE_SORT=0` default; shim pairing logic depends on it",
    );

    // 2. Copy vendor -> OUT_DIR/build (idempotent; vendor tree is never mutated in place).
    if build_dir.exists() {
        fs::remove_dir_all(&build_dir).unwrap();
    }
    fs::create_dir_all(&build_dir).unwrap();
    copy_dir(
        &manifest.join("vendor/bwa-mem3"),
        &build_dir.join("bwa-mem3"),
    );

    // 3. Apply any patches in patches/ lexicographic order. Expected empty in v1.
    let patches_dir = manifest.join("patches");
    if patches_dir.is_dir() {
        let mut patches: Vec<_> = fs::read_dir(&patches_dir)
            .unwrap()
            .filter_map(Result::ok)
            .filter(|e| {
                e.file_name()
                    .to_string_lossy()
                    .to_lowercase()
                    .ends_with(".patch")
            })
            .collect();
        patches.sort_by_key(std::fs::DirEntry::file_name);
        for p in patches {
            // Vendored source may have CRLF line endings; patch files we ship
            // are LF. `--ignore-whitespace` tolerates that mismatch.
            let status = Command::new("patch")
                .current_dir(build_dir.join("bwa-mem3"))
                .args(["-p1", "--ignore-whitespace", "-i"])
                .arg(p.path())
                .status()
                .expect("patch command not found");
            assert!(
                status.success(),
                "patch {} failed to apply",
                p.path().display()
            );
        }
    }

    let vendor_src = build_dir.join("bwa-mem3/src");
    let ssl = build_dir.join("bwa-mem3/ext/safestringlib");
    let s2n = build_dir.join("bwa-mem3/ext/sse2neon");

    // 4a. safestringlib is pure C — compile as C with its own cc::Build.
    if ssl.is_dir() {
        let safeclib = ssl.join("safeclib");
        if safeclib.is_dir() {
            let mut c_build = cc::Build::new();
            for entry in fs::read_dir(&safeclib).unwrap() {
                let e = entry.unwrap();
                if e.path().extension().is_some_and(|x| x == "c") {
                    c_build.file(e.path());
                }
            }
            c_build.include(ssl.join("include"));
            c_build.flag_if_supported("-Wno-unused-parameter");
            c_build.flag_if_supported("-Wno-unused-variable");
            c_build.flag_if_supported("-Wno-unused-function");
            c_build.flag_if_supported("-Wno-unused-but-set-variable");
            c_build.flag_if_supported("-Wno-format");
            // safestringlib relies on implicit declarations in places; newer
            // compilers error on these. Relax to a warning.
            c_build.flag_if_supported("-Wno-implicit-function-declaration");
            c_build.flag_if_supported("-Wno-incompatible-pointer-types");
            apply_simd_flags(&mut c_build);
            c_build.compile("safestring");
        }
    }

    // 4b. bwa-mem3 + shim: C++.
    let mut build = cc::Build::new();
    build.cpp(true);

    let skip: &[&str] = &[
        "main.cpp",     // CLI entry point
        "bwtindex.cpp", // index builder; out of scope (users run bwa-mem3 index)
        "runsimd.cpp",  // runtime SIMD-dispatch launcher; has unguarded main()
    ];
    // fastmap.cpp used to be excluded (CLI-side batch driver) but is now
    // built to expose worker_alloc/worker_free. Its entry point is
    // `main_mem`, not `main`, so no collision with the Rust test harness.
    for entry in fs::read_dir(&vendor_src).unwrap() {
        let e = entry.unwrap();
        let name = e.file_name().to_string_lossy().into_owned();
        if e.path().extension().is_some_and(|x| x == "cpp") && !skip.iter().any(|s| *s == name) {
            build.file(e.path());
        }
    }

    if ssl.is_dir() {
        build.include(ssl.join("include"));
    }
    if s2n.is_dir() {
        build.include(&s2n);
    }

    build.file(manifest.join("shim/bwa_shim.cpp"));
    build.file(manifest.join("shim/bwa_shim_align.cpp"));
    build.include(manifest.join("shim"));
    build.include(&vendor_src);

    apply_simd_flags(&mut build);

    build.std("c++17");
    build.define("ENABLE_PREFETCH", None);
    build.define("V17", Some("1"));
    build.define("MATE_SORT", Some("0"));
    build.flag_if_supported("-Wno-unused-parameter");
    build.flag_if_supported("-Wno-sign-compare");
    build.flag_if_supported("-Wno-unused-variable");
    build.flag_if_supported("-Wno-unused-function");
    build.flag_if_supported("-Wno-unused-but-set-variable");
    build.flag_if_supported("-Wno-deprecated-declarations");
    build.flag_if_supported("-Wno-format");
    build.flag_if_supported("-Wno-format-truncation");

    build.compile("bwa-mem3");

    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=pthread");
    let target = env::var("TARGET").unwrap();
    if target.contains("apple") {
        println!("cargo:rustc-link-lib=c++");
    } else {
        println!("cargo:rustc-link-lib=stdc++");
    }

    // 5. Generate Rust bindings for the shim header.
    generate_bindings(&manifest, &vendor_src, &out);
}

fn apply_simd_flags(build: &mut cc::Build) {
    let arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    match arch.as_str() {
        "x86_64" => {
            if cfg!(feature = "avx512") {
                build.flag_if_supported("-mavx512bw");
                build.flag_if_supported("-mavx512vl");
            } else if cfg!(feature = "sse42") {
                build.flag_if_supported("-msse4.2");
            } else if cfg!(feature = "native") {
                build.flag_if_supported("-march=native");
            } else {
                // default: AVX2
                build.flag_if_supported("-mavx2");
                build.flag_if_supported("-msse4.1");
            }
        }
        "aarch64" => {
            build.flag_if_supported("-march=armv8-a+simd");
            // Tell upstream's preprocessor to take the SSE2/SSE4.1 fallback path
            // (not the AVX-512 branch). sse2neon provides the intrinsics.
            build.define("__SSE2__", Some("1"));
            build.define("__SSE4_1__", Some("1"));
        }
        other => panic!("unsupported target arch: {other}"),
    }
}

fn copy_dir(src: &Path, dst: &Path) {
    fs::create_dir_all(dst).unwrap();
    for entry in fs::read_dir(src).unwrap() {
        let e = entry.unwrap();
        let ft = e.file_type().unwrap();
        let from = e.path();
        let to = dst.join(e.file_name());
        if ft.is_dir() {
            copy_dir(&from, &to);
        } else {
            fs::copy(&from, &to).unwrap();
        }
    }
}

fn generate_bindings(manifest: &Path, _vendor_src: &Path, out: &Path) {
    let shim_dir = manifest.join("shim");
    let bindings = bindgen::Builder::default()
        .header(shim_dir.join("bwa_shim.h").to_string_lossy())
        .clang_arg(format!("-I{}", shim_dir.display()))
        .allowlist_type("BwaReadPair")
        .allowlist_type("BwaIndex")
        .allowlist_type("BwaSeeds")
        .allowlist_type("BwaBatch")
        .allowlist_type("mem_opt_t")
        .allowlist_type("mem_pestat_t")
        .allowlist_function("bwa_shim_.*")
        .allowlist_var("MEM_F_.*")
        .derive_default(true)
        .generate()
        .expect("bindgen failed");
    bindings
        .write_to_file(out.join("bindings.rs"))
        .expect("failed to write bindings.rs");
}
