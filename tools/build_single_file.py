#!/usr/bin/env python3
"""Build a self-contained private prototype for browser-editor upload."""

from __future__ import annotations

import base64
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DIST = ROOT / "dist"


def svg_data_url(relative_path: str) -> str:
    payload = (ROOT / relative_path).read_bytes()
    encoded = base64.b64encode(payload).decode("ascii")
    return f"data:image/svg+xml;base64,{encoded}"


def main() -> int:
    html = (ROOT / "index.html").read_text(encoding="utf-8")
    css = (ROOT / "styles.css").read_text(encoding="utf-8")
    core = (ROOT / "core.js").read_text(encoding="utf-8")
    cases = (ROOT / "data" / "cases.js").read_text(encoding="utf-8")
    app = (ROOT / "app.js").read_text(encoding="utf-8")

    asset_paths = sorted(
        set(re.findall(r"assets/rulings/[A-Za-z0-9._-]+\.svg", cases))
    )
    for asset_path in asset_paths:
        cases = cases.replace(asset_path, svg_data_url(asset_path))

    html = html.replace(
        '    <link rel="stylesheet" href="styles.css">\n',
        f"    <style>\n{css}\n    </style>\n",
    )
    html = html.replace(
        '    <script src="core.js"></script>\n'
        '    <script src="data/cases.js"></script>\n'
        '    <script src="app.js"></script>\n',
        f"    <script>\n{core}\n    </script>\n"
        f"    <script>\n{cases}\n    </script>\n"
        f"    <script>\n{app}\n    </script>\n",
    )
    if "assets/rulings/" in html:
        raise AssertionError("the single-file build still references SVG assets")
    if 'src="core.js"' in html or 'href="styles.css"' in html:
        raise AssertionError("the single-file build still references source files")

    DIST.mkdir(parents=True, exist_ok=True)
    destination = DIST / "index.html"
    destination.write_text(html, encoding="utf-8")
    print(f"{destination} ({destination.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
