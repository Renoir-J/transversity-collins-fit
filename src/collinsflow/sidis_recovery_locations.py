"""Immutable location metadata scoped to one fully verified whole-job context."""
from collections import namedtuple
import copy
from pathlib import Path


class ScopedLocations(namedtuple("ScopedLocationsData", "root rows signature fields point_phases global_phases")):
    """No physical responses are cached; full filesystem checks stay at job boundaries."""

    __slots__ = ()

    def retained(self, phase, index):
        if phase not in self.point_phases or type(index) is not int or not 1 <= index <= 132:
            raise ValueError("declared phase and exact index1..132 required")
        for row in self.rows:
            if row[0] == phase and row[1] == index:
                return dict(zip(self.fields, row, strict=True))
        return None

    def current_path(self, *parts):
        path = self.root.joinpath(*parts)
        if not path.is_relative_to(self.root) or any(part in ("..", ".") for part in path.parts):
            raise ValueError("current-run paths cannot escape their immutable root")
        return path

    def folder(self, phase, index=0):
        if phase in self.point_phases:
            row = self.retained(phase, index)
            return Path(row["folder"]) if row else self.current_path(phase, f"p{index:04d}")
        if phase not in self.global_phases or type(index) is not int or index != 0:
            raise ValueError("declared fresh global phase/index0 required")
        return self.current_path(phase)

    def fresh_folder(self, phase, index=0):
        path = self.folder(phase, index)
        # Check the actual target again immediately before a new write.
        if not path.is_relative_to(self.root) or path.exists() or path.is_symlink() or path.resolve() != path:
            raise ValueError("never dispatch through retained, aliased or occupied folders")
        return path

    def verify_boundary(self, paths):
        restored = paths.RecoveryRunPaths(self.root, [dict(zip(self.fields, row, strict=True)) for row in self.rows])
        if restored.signature != self.signature:
            raise ValueError("changed immutable location snapshot")
        return restored


def snapshot(validated_view, paths):
    validated_view.valid()
    fields = paths.FIELDS
    rows = tuple(tuple(copy.deepcopy(row[key]) for key in fields) for row in validated_view.entries)
    return ScopedLocations(validated_view.root, rows, validated_view.signature, fields, paths.POINT_PHASES, paths.GLOBAL_PHASES)
