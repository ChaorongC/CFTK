#!/usr/bin/env python3
"""Backward-compatible source-checkout wrapper for CFTK evidence reporting."""

from pathlib import Path
import sys


SOURCE_ROOT = Path(__file__).resolve().parents[2] / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from validation_reports import EVIDENCE_FILENAMES, main, summarize


__all__ = ["EVIDENCE_FILENAMES", "main", "summarize"]


if __name__ == "__main__":
    raise SystemExit(main())
