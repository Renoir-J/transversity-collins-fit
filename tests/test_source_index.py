from pathlib import Path

from collinsflow.provenance import read_json


def test_fortran_architecture_sources_exist_in_locked_index():
    root = Path(__file__).resolve().parents[1]
    inventory = read_json(root / "references/source-inventory.json")
    paths = {"external/" + entry["repo"] + "/" + entry["path"] for entry in inventory["source_files"]}
    numerics = read_json(root / "config/numerical-backends.json")
    for module in numerics["fortran_policy"]["module_actions"]:
        for source in module["sources"]:
            assert source in paths, source
