#!/usr/bin/env python3
"""Build an aggregate MIMIC-IV download completion monitor.

This script reports only file-level integrity counts and relative filenames.
It does not read or export row-level clinical records.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import subprocess
from datetime import date
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MIMIC_ROOT = Path.home() / "secure_data" / "mimiciv" / "3.1"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"


def read_sha_manifest() -> list[tuple[str, str]]:
    sha_path = MIMIC_ROOT / "SHA256SUMS.txt"
    rows: list[tuple[str, str]] = []
    for line in sha_path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            rows.append((parts[0], parts[1]))
    return rows


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gzip_status(path: Path) -> str:
    try:
        with gzip.open(path, "rb") as handle:
            for _ in iter(lambda: handle.read(1024 * 1024), b""):
                pass
        return "pass"
    except Exception:
        return "fail"


def active_transfer_count() -> int:
    cmd = "ps aux | rg -i '(wget|curl).*physionet\\.org/files/mimiciv/3\\.1' | rg -v 'rg -i' || true"
    out = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True, check=False).stdout
    return len([line for line in out.splitlines() if line.strip()])


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    manifest = read_sha_manifest()
    rows: list[dict[str, object]] = []
    for expected_sha, rel_file in manifest:
        path = MIMIC_ROOT / rel_file
        exists = path.exists()
        gz_status = "not_applicable"
        sha_status = "missing"
        size_bytes = ""
        if exists:
            size_bytes = path.stat().st_size
            if rel_file.endswith(".gz"):
                gz_status = gzip_status(path)
            observed_sha = sha256_file(path)
            sha_status = "pass" if observed_sha == expected_sha else "fail"
        rows.append(
            {
                "file": rel_file,
                "exists": "yes" if exists else "no",
                "size_bytes": size_bytes,
                "gzip_status": gz_status,
                "sha256_status": sha_status,
            }
        )

    ok_files = [
        row
        for row in rows
        if row["exists"] == "yes"
        and row["sha256_status"] == "pass"
        and row["gzip_status"] in {"pass", "not_applicable"}
    ]
    missing = [row["file"] for row in rows if row["exists"] != "yes"]
    bad_gzip = [row["file"] for row in rows if row["gzip_status"] == "fail"]
    bad_sha = [row["file"] for row in rows if row["sha256_status"] == "fail"]

    table_path = TABLE_DIR / "mimiciv_download_completion_monitor_v19_0.csv"
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "exists", "size_bytes", "gzip_status", "sha256_status"])
        writer.writeheader()
        writer.writerows(rows)

    status_path = TABLE_DIR / "mimiciv_download_completion_status_v19_0.csv"
    total_size_gb = sum(int(row["size_bytes"] or 0) for row in rows) / 1024**3
    status_rows = [
        ("expected_files", len(rows)),
        ("ok_files", len(ok_files)),
        ("missing_files", len(missing)),
        ("bad_gzip_files", len(bad_gzip)),
        ("sha_fail_files", len(bad_sha)),
        ("total_present_size_gb", f"{total_size_gb:.3f}"),
        ("active_physionet_transfer_count", active_transfer_count()),
        ("download_gate", "pass" if len(ok_files) == len(rows) else "fail"),
    ]
    with status_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(status_rows)

    log_path = LOG_DIR / "mimiciv_download_completion_monitor_v19_0.md"
    lines = [
        "# MIMIC-IV Download Completion Monitor V19.0",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- Scope: aggregate file-level integrity monitor only.",
        "- Data governance: no row-level clinical records were exported; only relative filenames, sizes, gzip status and SHA256 status were written.",
        f"- Expected files: {len(rows)}.",
        f"- Integrity-passing files: {len(ok_files)}.",
        f"- Missing files: {len(missing)}.",
        f"- Gzip-failed files: {len(bad_gzip)}.",
        f"- SHA256-failed files: {len(bad_sha)}.",
        f"- Total present size: {total_size_gb:.3f} GB.",
        f"- Active PhysioNet transfer count: {active_transfer_count()}.",
        f"- Download gate: {'PASS' if len(ok_files) == len(rows) else 'FAIL'}.",
        "",
        "## Boundary",
        "",
        "This monitor updates only the download-integrity gate. It does not establish a MIMIC cohort, endpoint definition, model result, or manuscript evidence claim.",
        "",
    ]
    if missing or bad_gzip or bad_sha:
        lines.extend(
            [
                "## Non-pass Files",
                "",
                f"- Missing: {', '.join(missing) if missing else 'none'}.",
                f"- Gzip failed: {', '.join(bad_gzip) if bad_gzip else 'none'}.",
                f"- SHA256 failed: {', '.join(bad_sha) if bad_sha else 'none'}.",
                "",
            ]
        )
    lines.append("PASS: MIMIC-IV download completion monitor V19.0 built.")
    log_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {log_path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
