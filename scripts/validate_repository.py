#!/usr/bin/env python3
"""Rebuild the public displays and verify manuscript-facing artifacts."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]

EXPECTED = {
    "rebuild/source_data/figure1_source_data.csv": "ce799bc0351439d300647aa135fd662a412baa2b839f2c617144afa9fc5ca64e",
    "rebuild/source_data/figure2_source_data.csv": "aa9f9f4c30b40f62be682f21caf0098b87b0d788f6b0a2af5bf6b48f1a9a2671",
    "rebuild/source_data/figure3_source_data.csv": "370bc7275f640fbd3fc434a5248b9334d4970ae409b7ae3af6f89d79406d7b64",
    "rebuild/source_data/table1_source_data.csv": "2172d7cdaf5dfc57191a89a250b4892959ec81f9d81af860640210c157b3f4c4",
    "rebuild/source_data/table2_source_data.csv": "710ca9ad7c986b6f1fded4dcc34d5e042a301d5b9b7b8cc068ea5d25f7dfd0b4",
    "rebuild/source_data/table3_source_data.csv": "6d4b94b67b3b0202b1dcbd57f21cee6c8a2dc721ffe1b22bbeb4af98905f9d92",
    "rebuild/tables/table1_baseline.tex": "dee6968c47224f6ff07b4f7c600959cadafa70932771a6ecdcc0c7680eea72a5",
    "rebuild/tables/table2_primary_associations.tex": "f62925ad434ccaff1e8406b0f09c7bf1dc69f9776ed1e897998cf523ccdbbf8a",
    "rebuild/tables/table3_sensitivity.tex": "c3fcf286c59be719a2a993c940935cacbb9799d187bd1c86bd08c54be188b82f",
}

IMAGE_EXPECTED = {
    "rebuild/figures/figure1_cohort_flow.tiff": (4014, 3810),
    "rebuild/figures/figure2_adjustment_ladder.tiff": (4014, 2670),
    "rebuild/figures/figure3_measurement_robustness.tiff": (4014, 3720),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "build_displays.py")],
        cwd=ROOT,
        check=True,
    )
    failures: list[str] = []
    for relative, expected in EXPECTED.items():
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"missing: {relative}")
            continue
        observed = sha256(path)
        if observed != expected:
            failures.append(
                f"hash mismatch: {relative}; expected {expected}; observed {observed}"
            )
    for relative, expected_size in IMAGE_EXPECTED.items():
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"missing: {relative}")
            continue
        with Image.open(path) as image:
            dpi = tuple(round(float(value)) for value in image.info.get("dpi", (0, 0)))
            if image.size != expected_size:
                failures.append(
                    f"size mismatch: {relative}; expected {expected_size}; "
                    f"observed {image.size}"
                )
            if image.mode != "RGB":
                failures.append(
                    f"color-mode mismatch: {relative}; expected RGB; "
                    f"observed {image.mode}"
                )
            if dpi != (600, 600):
                failures.append(
                    f"resolution mismatch: {relative}; expected 600 dpi; "
                    f"observed {dpi}"
                )
            if image.info.get("compression") != "tiff_lzw":
                failures.append(
                    f"compression mismatch: {relative}; expected tiff_lzw; "
                    f"observed {image.info.get('compression')}"
                )
    if failures:
        raise RuntimeError("\n".join(failures))
    print(
        "REPOSITORY VALIDATION PASS: 9 source/table artifacts match the "
        "submission hashes; 3 TIFF figures pass dimension, RGB, LZW, and "
        "600-dpi checks."
    )


if __name__ == "__main__":
    main()
