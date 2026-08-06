"""Ensure every user-guide page exposes visible, concrete output guidance."""

from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[1]
GUIDE_DIR = REPO_ROOT / "docs" / "user_guide"
FIGURE_RE = re.compile(r"^\s*\.\.\s+(?:figure|image)::\s+(\S+)", re.MULTILINE)
OUTPUT_HEADING_RE = re.compile(
    r"^Expected Outputs(?: By Step| Location)?\s*\n[-=]{3,}\s*$",
    re.MULTILINE,
)


def test_every_user_guide_page_has_visible_output_and_contract():
    pages = sorted(GUIDE_DIR.glob("*.rst"))
    assert pages, "user-guide inventory is empty"

    failures = []
    for page in pages:
        text = page.read_text(encoding="utf-8")
        refs = FIGURE_RE.findall(text)
        local_refs = [ref for ref in refs if not re.match(r"^[a-z]+://", ref)]
        missing = [
            str((page.parent / ref).resolve())
            for ref in local_refs
            if not (page.parent / ref).is_file()
        ]
        if not local_refs:
            failures.append(f"{page.name}: no local figure/image directive")
        if missing:
            failures.append(f"{page.name}: missing assets: {', '.join(missing)}")
        if not OUTPUT_HEADING_RE.search(text):
            failures.append(f"{page.name}: no Expected Outputs section")

    assert not failures, "\n".join(failures)
