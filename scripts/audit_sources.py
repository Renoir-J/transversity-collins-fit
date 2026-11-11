"""Record source provenance and verify FeynFacet against canonical Git blobs.

The Git-blob check distinguishes upstream content from Windows checkout newline
conversion. Inventory is acquisition coverage, not a claim of semantic review.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPOS = ["FeynFacet", "jam3d", "jam3d_dev_lib", "jam3dlib", "tmd-parametrizations", "hoppet", "PhysMaster", "tParton"]


def git(folder, *args):
    return subprocess.check_output(["git", "-C", str(folder), *args])


def canonical_checksums():
    repo = ROOT / "external/FeynFacet"
    expected = dict(line.split("  ", 1)[::-1]
                    for line in git(repo, "show", "HEAD:SHA256SUMS").decode().splitlines())
    names = list(expected)
    requests = "".join(f"HEAD:{name}\n" for name in names).encode()
    output = subprocess.check_output(["git", "-C", str(repo), "cat-file", "--batch"], input=requests)
    position, mismatches, checkout_differences, newline_only = 0, [], [], []
    for name in names:
        stop = output.index(b"\n", position)
        header = output[position:stop].split()
        if len(header) != 3 or header[1] != b"blob":
            raise ValueError(f"Missing canonical blob: {name}")
        size = int(header[2])
        blob = output[stop + 1:stop + 1 + size]
        position = stop + 1 + size + 1
        if hashlib.sha256(blob).hexdigest() != expected[name]:
            mismatches.append(name)
        local = (repo / name).read_bytes()
        if local != blob:
            checkout_differences.append(name)
            if local.replace(b"\r\n", b"\n") == blob:
                newline_only.append(name)
    return {"canonical_blobs_checked": len(names), "canonical_mismatches": mismatches,
            "working_tree_differences": len(checkout_differences),
            "differences_explained_by_CRLF": len(newline_only),
            "unexplained_working_tree_differences": sorted(set(checkout_differences) - set(newline_only))}


def audit():
    repositories, files = [], []
    for name in REPOS:
        repo = ROOT / "external" / name
        head = git(repo, "rev-parse", "HEAD").decode().strip()
        url = git(repo, "remote", "get-url", "origin").decode().strip()
        tracked = git(repo, "ls-tree", "-r", "--name-only", "HEAD").decode().splitlines()
        present = [p for p in repo.rglob("*") if p.is_file() and ".git" not in p.parts
                   and "__pycache__" not in p.parts]
        for path in present:
            if path.suffix.lower() in {".py", ".f", ".f90", ".cc", ".cpp", ".h", ".wl", ".wls", ".md"}:
                payload = path.read_bytes()
                files.append({"repo": name, "path": path.relative_to(repo).as_posix(),
                              "sha256_worktree": hashlib.sha256(payload).hexdigest(),
                              "lines": len(payload.splitlines()), "bytes": len(payload)})
        repositories.append({"name": name, "origin": url, "commit": head,
                             "tracked_files": len(tracked), "materialized_files": len(present),
                             "snapshot": "full" if name in {"FeynFacet", "jam3dlib", "tmd-parametrizations", "tParton"}
                             else "sparse; source and selected data only",
                             "windows_reserved_source_omitted": [p for p in tracked if p.endswith("/aux.py")]})
    return {"audited_at": datetime.now(timezone.utc).isoformat(),
            "coverage_definition": "Inventory of acquired files; review scope is specified in code-review.json",
            "repositories": repositories, "source_files": files,
            "feynfacet_integrity": canonical_checksums()}


if __name__ == "__main__":
    result = audit()
    destination = ROOT / "references/source-inventory.json"
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"repositories": result["repositories"],
                      "source_files_inventoried": len(result["source_files"]),
                      "feynfacet_integrity": result["feynfacet_integrity"]}, indent=2))
    if result["feynfacet_integrity"]["canonical_mismatches"]:
        raise SystemExit(1)
