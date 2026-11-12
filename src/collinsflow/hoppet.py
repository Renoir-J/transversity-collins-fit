"""Acquire/build an unchanged HOPPET and the project adapter, outside numerics."""
from __future__ import annotations

import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

from .provenance import atomic_json, file_hash, read_json, utc_now


def prepare(root, destination, *, source=None, offline=False, cmake="cmake",
            fc="gfortran", fc_flags=(), jobs=2):
    root, destination = Path(root).resolve(), Path(destination).resolve()
    if os.name != "posix":
        raise ValueError("build HOPPET on Linux/WSL")
    if not destination.is_relative_to(root / ".runtime"):
        raise ValueError("HOPPET build must stay inside ignored .runtime")
    destination.mkdir(parents=True, exist_ok=True)
    contract = root / "config/hoppet-contract.json"
    upstream = read_json(contract)["upstream"]
    source = Path(source).resolve() if source else destination / "source"
    if not source.is_relative_to(root / ".runtime"):
        raise ValueError("use a complete pinned source under .runtime; external is unchanged")
    environment = os.environ.copy()
    environment.update(GIT_TERMINAL_PROMPT="0", OMP_NUM_THREADS="1")
    commands = []
    receipt_path = destination / "build.json"
    receipt_path.unlink(missing_ok=True)

    def run(command, cwd=destination):
        commands.append(list(map(str, command)))
        with (destination / "build.log").open("a") as log:
            subprocess.run(command, cwd=cwd, env=environment, check=True,
                           stdout=log, stderr=subprocess.STDOUT, timeout=1200)

    def capture(command):
        return subprocess.check_output(command, env=environment, text=True, timeout=30).strip()

    if not (source / ".git").is_dir():
        if offline:
            raise ValueError("complete pinned HOPPET source is unavailable offline")
        if source.exists() and any(source.iterdir()):
            raise ValueError("refuse to replace an existing source directory")
        source.mkdir(parents=True, exist_ok=True)
        run(["git", "init", str(source)])
        run(["git", "-C", str(source), "remote", "add", "origin", upstream["url"]])
        run(["git", "-C", str(source), "fetch", "--depth=1", "origin", upstream["commit"]])
        run(["git", "-C", str(source), "-c", "core.hooksPath=/dev/null",
             "-c", "core.autocrlf=false", "checkout", "--detach", upstream["commit"]])
    git = ["git", "-C", str(source)]
    if (capture([*git, "rev-parse", "HEAD"]) != upstream["commit"]
            or capture([*git, "rev-parse", "HEAD^{tree}"]) != upstream["tree"]
            or capture([*git, "status", "--porcelain"])):
        raise ValueError("HOPPET source does not match the clean locked commit/tree")
    if not (source / "cmake/Modules").is_dir():
        raise ValueError("a sparse reference snapshot cannot build HOPPET")
    project_sources = [contract, root / "config/evolution-contract.json",
        root / "config/matched-rows-contract.json",
        root / "fortran/hoppet_bridge.f90", root / "src/collinsflow/hoppet.py",
        root / "src/collinsflow/provenance.py", root / "scripts/prepare_hoppet.py"]
    upstream_paths = [source / p for p in capture([*git, "ls-files"]).splitlines()]
    source_links = {str(p): os.readlink(p) for p in upstream_paths if p.is_symlink()}
    if any(not Path(p).resolve().is_relative_to(source) for p in source_links):
        raise ValueError("upstream symlink escapes its source tree")
    source_hashes = {str(p): file_hash(p) for p in [*project_sources, *upstream_paths]
                     if not p.is_symlink()}
    binary = destination / "libcollins_hoppet.so"
    # Never attribute stale native objects to a newly checked source.
    with tempfile.TemporaryDirectory(prefix="compile-", dir=destination) as temporary:
        build = Path(temporary).resolve()
        if not build.is_relative_to(destination):
            raise ValueError("build escaped its runtime directory")
        run([cmake, "-S", str(source), "-B", str(build), "-G", "Unix Makefiles",
            "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
            "-DCMAKE_Fortran_COMPILER=" + str(Path(shutil.which(fc) or fc).resolve()),
            "-DCMAKE_Fortran_FLAGS=" + shlex.join(fc_flags),
            "-DCMAKE_BUILD_RPATH=" + ";".join(flag.removeprefix("-Wl,-rpath,")
                                            for flag in fc_flags if flag.startswith("-Wl,-rpath,")),
            "-DHOPPET_BUILD_EXAMPLES=ON", "-DHOPPET_ENABLE_TESTING=ON",
            "-DHOPPET_BUILD_BENCHMARK=OFF", "-DHOPPET_BUILD_PYINTERFACE=OFF",
            "-DHOPPET_GIT_WATCHER=OFF"])
        run([cmake, "--build", str(build), "--parallel", str(jobs)])
        ctest = str(Path(shutil.which(cmake) or cmake).resolve().with_name("ctest"))
        run([ctest, "--test-dir", str(build), "--output-on-failure", "--parallel", "1"])
        test_log = destination / "upstream-tests.log"
        shutil.copyfile(build / "Testing/Temporary/LastTest.log", test_log)
        run([fc, *fc_flags, "-O2", "-fPIC", "-shared", "-ffree-line-length-none",
            "-I" + str(build / "include/hoppet"), "-J" + str(build),
            str(root / "fortran/hoppet_bridge.f90"), str(build / "lib/libhoppet.a"),
            "-lstdc++", "-o", str(binary)])
    dependencies = {}
    for line in capture(["ldd", str(binary)]).splitlines():
        match = re.search(r"=>\s+(/\S+)", line)
        if match:
            dependency = Path(match[1]).resolve()
            dependencies[str(dependency)] = file_hash(dependency)
    if not any("libgfortran" in p for p in dependencies):
        raise ValueError("Fortran runtime dependency is unresolved")
    if capture([*git, "status", "--porcelain"]):
        raise ValueError("upstream source changed during its build")
    if any(file_hash(p) != digest for p, digest in source_hashes.items()):
        raise ValueError("source changed during HOPPET build")
    if any(os.readlink(p) != target for p, target in source_links.items()):
        raise ValueError("source link changed during HOPPET build")
    atomic_json(receipt_path, {"status": "BUILT_NOT_VALIDATED", "scope": "reconstructed_theory",
        "production_ready": False, "ABI_version": 1, "upstream": upstream,
        "source": str(source), "source_clean": True, "library": str(binary),
        "upstream_ctest_executed": True,
        "source_hashes": source_hashes, "source_links": source_links,
        "output_hashes": {str(binary): file_hash(binary), str(test_log): file_hash(test_log), **dependencies},
        "commands": commands, "versions": {"Fortran": capture([fc, "--version"]).splitlines()[0],
        "CMake": capture([cmake, "--version"]).splitlines()[0]}, "completed_at": utc_now()})
    return str(receipt_path)
