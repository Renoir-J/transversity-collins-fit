"""Build isolation and source/archive regressions, using synthetic files only."""
import io
import tarfile

import pytest

from collinsflow.parton import extract_release, prepare
from collinsflow.provenance import source_fingerprint


@pytest.mark.parametrize("name,kind", [
    ("../escape", "file"), ("/absolute", "file"),
    ("LHAPDF-6.5.6/link", "link"), ("wrong-version/configure", "file"),
])
def test_release_rejects_unsafe_members(tmp_path, name, kind):
    archive = tmp_path / "synthetic.tar.gz"
    with tarfile.open(archive, "w:gz") as output:
        entry = tarfile.TarInfo(name)
        if kind == "link":
            entry.type = tarfile.SYMTYPE
            entry.linkname = "../../outside"
            output.addfile(entry)
        else:
            entry.size = 1
            output.addfile(entry, io.BytesIO(b"x"))
    with pytest.raises(ValueError, match="unsafe"):
        extract_release(archive, tmp_path / "unpacked")
    assert not (tmp_path / "escape").exists()


def test_build_cannot_write_to_external_or_system(tmp_path):
    for destination in (tmp_path / "external", tmp_path / "install", tmp_path.parent / "other"):
        with pytest.raises(ValueError, match="runtime"):
            prepare(tmp_path, destination, offline=True)
        assert not destination.exists()


def test_native_bridge_changes_invalidate_receipts(tmp_path):
    directory = tmp_path / "bridges"
    directory.mkdir()
    path = directory / "synthetic.cpp"
    path.write_text("synthetic ABI version one")
    before = source_fingerprint(tmp_path)
    path.write_text("synthetic ABI version two")
    assert source_fingerprint(tmp_path) != before
