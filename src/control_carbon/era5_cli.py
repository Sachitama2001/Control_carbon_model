"""Command-line ERA5 point acquisition for VISIT site climate files."""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path
from typing import Sequence

from .era5 import (
    ERA5PointRequest,
    aggregate_era5_hourly_for_visit,
    fetch_cds_era5_point,
    write_visit_ncep_site_climate,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Retrieve one ERA5 point through the official Copernicus CDS and "
            "write the VISIT "
            "WMODE=1 default NCEP/NCAR climate format."
        )
    )
    parser.add_argument("--latitude", required=True, type=float)
    parser.add_argument("--longitude", required=True, type=float)
    parser.add_argument("--start-date", required=True, type=date.fromisoformat)
    parser.add_argument("--end-date", required=True, type=date.fromisoformat)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--archive-directory",
        type=Path,
        help="CDS ZIP directory (default: OUTPUT.parent/era5-cds)",
    )
    parser.add_argument(
        "--refresh-cds",
        action="store_true",
        help="replace valid cached CDS archives instead of reusing them",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="metadata path (default: OUTPUT.manifest.json)",
    )
    parser.add_argument(
        "--allow-partial-years",
        action="store_true",
        help="allow a diagnostic file that VISIT cannot use as a complete run input",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    request = ERA5PointRequest(
        latitude=args.latitude,
        longitude=args.longitude,
        start_date=args.start_date,
        end_date=args.end_date,
    )
    archive_directory = args.archive_directory or args.output.parent / "era5-cds"
    hourly = fetch_cds_era5_point(
        request, archive_directory, use_cache=not args.refresh_cds
    )
    daily = aggregate_era5_hourly_for_visit(hourly)
    manifest_path = args.manifest or Path(f"{args.output}.manifest.json")
    write_visit_ncep_site_climate(
        daily,
        args.output,
        manifest_path=manifest_path,
        require_complete_years=not args.allow_partial_years,
    )
    print(
        f"Retrieved {len(hourly.times)} ERA5 hours from {hourly.provider} at grid "
        f"({hourly.resolved_latitude}, {hourly.resolved_longitude})"
    )
    print(f"Wrote {len(daily.dates)} VISIT daily rows to {args.output}")
    print(f"Wrote provenance manifest to {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())