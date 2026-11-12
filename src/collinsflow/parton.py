"""Whole-job acquisition and compilation for reconstructed PDF/FF comparisons.

The coordinator does no PDF interpolation. Raw bytes, native libraries and legacy
reference executables live in an ignored runtime directory, never in external/.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

from .experimental import acquire
from .provenance import atomic_json, file_hash, read_json, utc_now


def extract_release(archive, destination):
    """Reject links, special files and paths outside a single release directory."""
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as stream:
        members = stream.getmembers()
        if len(members) > 5000 or sum(m.size for m in members) > 50_000_000:
            raise ValueError("LHAPDF archive exceeds extraction budget")
        for member in members:
            path = PurePosixPath(member.name)
            target = (destination / member.name).resolve()
            if (path.is_absolute() or ".." in path.parts or not path.parts
                    or path.parts[0] != "LHAPDF-6.5.6"
                    or not target.is_relative_to(destination)
                    or not (member.isfile() or member.isdir())):
                raise ValueError("unsafe LHAPDF archive member")
        stream.extractall(destination, members=members, filter="data")
    return destination / "LHAPDF-6.5.6"


def prepare(root, destination, *, offline=False, fc="gfortran", fc_flags=(),
            fortran_runtime=None, jobs=2):
    root, destination = Path(root).resolve(), Path(destination).resolve()
    # This is a local build, not an installer into external/ or the system.
    if not destination.is_relative_to(root / ".runtime"):
        raise ValueError("parton build must stay inside the ignored .runtime directory")
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "config/parton-sources.json"
    source_paths = {}
    for source in read_json(manifest_path)["sources"]:
        group = source["group"]
        if group not in {"CT10nlo", "DSSV08", "DSS14", "LHAPDF"}:
            raise ValueError("unexpected source group")
        source_paths[source["id"]] = acquire(source, destination / "inputs" / group,
            seed=root / source["seed_path"], offline=offline)
    receipt_path = destination / "build.json"
    # Remove admission evidence before a rebuild; an interrupted build has none.
    receipt_path.unlink(missing_ok=True)
    source_hashes = {str(p): file_hash(p) for p in source_paths.values()}
    project_sources = [manifest_path, root / "src/collinsflow/parton.py",
        root / "src/collinsflow/experimental.py", root / "src/collinsflow/provenance.py",
        root / "scripts/prepare_parton_inputs.py", root / "bridges/lhapdf.cpp",
        root / "bridges/check_lhapdf.cpp", root / "fortran/check_dssv08.f90",
        root / "fortran/check_dss14.f90"]
    source_hashes.update({str(p): file_hash(p) for p in project_sources})
    commands = []

    def run(command, cwd=destination):
        commands.append([str(x) for x in command])
        with (destination / "build.log").open("a") as log:
            subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=1200)

    prefix = destination / "install"
    # Always compile a freshly extracted tree, so stale native objects cannot be
    # attributed to a newly verified archive. Only the locked download is cached.
    with tempfile.TemporaryDirectory(prefix="source-", dir=destination) as source_directory:
        release = extract_release(source_paths["LHAPDF_source"], source_directory)
        run([str(release / "configure"), "--prefix=" + str(prefix),
             "--disable-python", "--disable-static"], cwd=release)
        run(["make", f"-j{int(jobs)}"], cwd=release)
        run(["make", "install"], cwd=release)
    output = destination / "bin"
    output.mkdir(exist_ok=True)
    library = output / "libcollins_lhapdf.so"
    reference = output / "check_lhapdf"
    cxx_flags = ["-std=c++11", "-O2", "-I" + str(prefix / "include")]
    link_flags = ["-L" + str(prefix / "lib"), "-Wl,-rpath," + str(prefix / "lib"), "-lLHAPDF"]
    run(["g++", *cxx_flags, "-shared", "-fPIC", str(root / "bridges/lhapdf.cpp"),
         *link_flags, "-o", str(library)])
    run(["g++", *cxx_flags, str(root / "bridges/check_lhapdf.cpp"),
         *link_flags, "-o", str(reference)])
    references = {}
    for name, source_id, driver, grid_id, grid_name in (
        ("DSSV08", "DSSV08_code", "check_dssv08", "DSSV08_grid",
         "analysis/qpdlib/sets/DSSV/Grids/DSSV-GRID.NLO"),
        ("DSS14", "DSS14_code", "check_dss14", "DSS14_grid", "PI_DSS14.GRID"),
    ):
        work = destination / "references" / name
        work.mkdir(parents=True, exist_ok=True)
        grid = work / grid_name
        grid.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_paths[grid_id], grid)
        objects = [work / "upstream.o", work / "driver.o"]
        for src, obj in zip((source_paths[source_id], root / "fortran" / (driver + ".f90")), objects):
            run([fc, *fc_flags, "-O2", "-std=legacy", "-ffixed-line-length-none",
                 "-c", str(src), "-o", str(obj)])
        executable = work / driver
        if fortran_runtime:
            # Isolated WSL compiler frontend can use the Julia-distributed runtime.
            runtime = Path(fortran_runtime).resolve()
            run(["gcc", *map(str, objects), str(runtime), "-lm",
                 "-Wl,-rpath," + str(runtime.parent), "-o", str(executable)])
        else:
            run([fc, *fc_flags, *map(str, objects), "-o", str(executable)])
        references[name] = {"executable": str(executable), "cwd": str(work)}
    outputs = [library, reference, prefix / "lib/libLHAPDF.so",
               prefix / "share/LHAPDF/lhapdf.conf"]
    outputs.extend(Path(ref["executable"]) for ref in references.values())
    outputs.extend((destination / "references/DSSV08/analysis/qpdlib/sets/DSSV/Grids/DSSV-GRID.NLO",
                    destination / "references/DSS14/PI_DSS14.GRID"))
    versions = {compiler: subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]
                for compiler in ("g++", fc)}
    if fortran_runtime:
        outputs.append(Path(fortran_runtime).resolve())
    for path, sha in source_hashes.items():
        if file_hash(path) != sha:
            raise ValueError("build inputs changed during compilation")
    receipt = {"scope": "reconstructed_theory", "production_ready": False,
        "status": "BUILT_NOT_VALIDATED", "created_at": utc_now(), "LHAPDF_version": "6.5.6",
        "input_root": str(destination / "inputs"), "ct10_library": str(library),
        "ct10_reference": str(reference), "references": references,
        "source_hashes": source_hashes, "output_hashes": {str(p): file_hash(p) for p in outputs},
        "compilers": versions, "commands": commands, "platform": os.uname().sysname}
    atomic_json(receipt_path, receipt)
    return receipt_path
