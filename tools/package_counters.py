#!/usr/bin/env python3
"""Reproducibly package the checked-in standalone counters, without regeneration.

Unlike package_engine.py, this command never imports or refreshes source from
the parent research repository. Run it after editing either downloadable package.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import zipfile


DOWNLOADS = Path(__file__).resolve().parent.parent / "downloads"
PACKAGES = ("direct-ruling-counter", "tropical-tree-counter")


def package_counter(name: str) -> tuple[int, str]:
    source = DOWNLOADS / name
    if not (source / "README.md").is_file():
        raise FileNotFoundError(f"missing package README: {source}")
    members = sorted(
        (path for path in source.rglob("*")
         if path.is_file()
         and (path.suffix in {".py", ".js", ".md"}
              or path.name in {"example_polygon.json", "examples.json"})
         and not any(part.startswith(".") or part == "__pycache__"
                     for part in path.relative_to(source).parts)),
        key=lambda path: path.relative_to(source).as_posix(),
    )
    archive_path = DOWNLOADS / f"{name}.zip"
    temporary = archive_path.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in members:
            info = zipfile.ZipInfo(
                f"{name}/{path.relative_to(source).as_posix()}",
                date_time=(1980, 1, 1, 0, 0, 0),
            )
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, path.read_bytes(),
                             compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    with zipfile.ZipFile(temporary) as archive:
        if archive.testzip() is not None:
            raise RuntimeError(f"invalid ZIP: {temporary}")
        for path in members:
            name_in_zip = f"{name}/{path.relative_to(source).as_posix()}"
            if archive.read(name_in_zip) != path.read_bytes():
                raise RuntimeError(f"ZIP content mismatch: {name_in_zip}")
    temporary.replace(archive_path)
    return len(members), hashlib.sha256(archive_path.read_bytes()).hexdigest()


def main() -> None:
    for name in PACKAGES:
        count, digest = package_counter(name)
        archive = DOWNLOADS / f"{name}.zip"
        print(f"{archive.name}: {count} files, {archive.stat().st_size} bytes")
        print(f"sha256={digest}")


if __name__ == "__main__":
    main()
