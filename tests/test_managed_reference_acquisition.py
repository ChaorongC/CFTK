import gzip
import hashlib
import io
import json
import threading
import time
from pathlib import Path
from types import SimpleNamespace
from urllib.error import URLError

import pytest

from test_init_modes import make_sample_sheet


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def profiles_module(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import reference_profiles

    return reference_profiles


@pytest.fixture
def init_module(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import init

    return init


class Response(io.BytesIO):
    def __init__(self, payload, url):
        super().__init__(payload)
        self._url = url

    def geturl(self):
        return self._url


class FixtureOpener:
    def __init__(self, payloads, *, delay=0):
        self.payloads = payloads
        self.delay = delay
        self.calls = []
        self._lock = threading.Lock()

    def __call__(self, request, timeout=None):
        url = request.full_url
        with self._lock:
            self.calls.append(url)
        if self.delay:
            time.sleep(self.delay)
        return Response(self.payloads[url], url)


def _digest(payload):
    return hashlib.sha256(payload).hexdigest()


def make_registry(
    *, corrupt_download_hash=False, escaped_path=False, derived_chrom_sizes=False
):
    installed = {
        "genome_fa": b">chr1\nACGT\n",
        "genome_2bit": b"2bit",
        "chrom_sizes": b"chr1\t4\n",
        "target_bed": b"chr1\t0\t4\n",
    }
    paths = {
        "genome_fa": "genome/hg38.fa",
        "genome_2bit": "genome/hg38.2bit",
        "chrom_sizes": "genome/hg38.chrom.sizes",
        "target_bed": "assay/covered_targets.bed",
    }
    downloads = dict(installed)
    downloads["genome_fa"] = gzip.compress(installed["genome_fa"], mtime=0)
    if derived_chrom_sizes:
        downloads["chrom_sizes"] = b"chr1\t4\t6\t4\t5\n"
    payloads = {}
    components = {}
    for name, payload in downloads.items():
        url = f"https://downloads.example.test/releases/1.0.0/{name}"
        payloads[url] = payload
        artifact_hash = _digest(payload)
        if corrupt_download_hash and name == "genome_fa":
            artifact_hash = "0" * 64
        components[name] = {
            "path": "../../outside.fa" if escaped_path and name == "genome_fa" else paths[name],
            "size": len(installed[name]),
            "sha256": _digest(installed[name]),
            "artifact": {
                "urls": [url],
                "size": len(payload),
                "sha256": artifact_hash,
                "compression": "gzip" if name == "genome_fa" else "none",
                "immutable": True,
            },
            "license": {
                "name": "CC0-1.0",
                "url": "https://creativecommons.org/publicdomain/zero/1.0/",
            },
            "source": {
                "name": "CFTK test fixture",
                "url": "https://example.test/reference-provenance",
            },
        }
        if derived_chrom_sizes and name == "chrom_sizes":
            components[name]["artifact"]["transform"] = "fai_to_chrom_sizes"
    registry = {
        "registry_version": 1,
        "profiles": {
            "twist_human_methylome_hg38": {
                "1.0.0": {
                    "profile_id": "twist_human_methylome_hg38",
                    "version": "1.0.0",
                    "assay": "twist_human_methylome",
                    "genome": "hg38",
                    "components": components,
                }
            }
        },
    }
    return registry, payloads


def acquire(profiles_module, root, registry, opener):
    return profiles_module.acquire_reference_profile(
        mode="managed",
        reference_root=root,
        profile_id="twist_human_methylome_hg38",
        version="1.0.0",
        registry=registry,
        opener=opener,
        verify_checksums=True,
        validate_compatibility=True,
    )


def test_registry_requires_license_source_size_hash_and_immutable_https(
    profiles_module,
):
    registry, _ = make_registry()
    component = registry["profiles"]["twist_human_methylome_hg38"]["1.0.0"][
        "components"
    ]["target_bed"]
    del component["license"]

    with pytest.raises(SystemExit, match="target_bed.*license"):
        profiles_module.validate_reference_registry(registry)

    registry, _ = make_registry()
    artifact = registry["profiles"]["twist_human_methylome_hg38"]["1.0.0"][
        "components"
    ]["target_bed"]["artifact"]
    artifact["urls"] = ["http://downloads.example.test/latest/targets.bed"]
    artifact["immutable"] = False

    with pytest.raises(SystemExit, match="target_bed.*immutable HTTPS"):
        profiles_module.validate_reference_registry(registry)


def test_managed_acquisition_verifies_decompresses_and_publishes_atomically(
    profiles_module, tmp_path
):
    registry, payloads = make_registry()
    opener = FixtureOpener(payloads)

    profile = acquire(profiles_module, tmp_path, registry, opener)

    profile_dir = tmp_path / "twist_human_methylome_hg38" / "1.0.0"
    assert Path(profile["profile_dir"]) == profile_dir
    assert (profile_dir / "genome/hg38.fa").read_bytes() == b">chr1\nACGT\n"
    assert len(opener.calls) == 4
    manifest = json.loads((profile_dir / "manifest.json").read_text())
    assert manifest["acquisition"]["registry_version"] == 1
    assert manifest["components"]["genome_fa"]["sha256"] == _digest(
        b">chr1\nACGT\n"
    )
    assert not list(tmp_path.glob(".cftk-install-*"))


def test_managed_acquisition_derives_two_column_chrom_sizes_from_fai(
    profiles_module, tmp_path
):
    registry, payloads = make_registry(derived_chrom_sizes=True)

    profile = acquire(
        profiles_module, tmp_path, registry, FixtureOpener(payloads)
    )

    assert Path(profile["components"]["chrom_sizes"]).read_bytes() == b"chr1\t4\n"


def test_registry_rejects_unknown_artifact_transform(profiles_module):
    registry, _ = make_registry(derived_chrom_sizes=True)
    registry["profiles"]["twist_human_methylome_hg38"]["1.0.0"][
        "components"
    ]["chrom_sizes"]["artifact"]["transform"] = "run_arbitrary_command"

    with pytest.raises(SystemExit, match="chrom_sizes.*artifact.transform"):
        profiles_module.validate_reference_registry(registry)


def test_fai_projection_rejects_malformed_input_and_cleans_stage(
    profiles_module, tmp_path
):
    registry, payloads = make_registry(derived_chrom_sizes=True)
    component = registry["profiles"]["twist_human_methylome_hg38"]["1.0.0"][
        "components"
    ]["chrom_sizes"]
    url = component["artifact"]["urls"][0]
    payloads[url] = b"chr1\tnot-an-integer\t6\t4\t5\n"
    component["artifact"]["size"] = len(payloads[url])
    component["artifact"]["sha256"] = _digest(payloads[url])

    with pytest.raises(SystemExit, match="invalid FAI sequence length.*chrom_sizes"):
        acquire(profiles_module, tmp_path, registry, FixtureOpener(payloads))
    assert not list(tmp_path.glob(".cftk-install-*"))
    assert not (tmp_path / "twist_human_methylome_hg38/1.0.0").exists()


def test_managed_acquisition_cleans_up_after_corruption(profiles_module, tmp_path):
    registry, payloads = make_registry(corrupt_download_hash=True)

    with pytest.raises(SystemExit, match="download checksum mismatch.*genome_fa"):
        acquire(profiles_module, tmp_path, registry, FixtureOpener(payloads))

    assert not (tmp_path / "twist_human_methylome_hg38" / "1.0.0").exists()
    assert not list(tmp_path.glob(".cftk-install-*"))


def test_managed_acquisition_is_idempotent_and_refuses_corrupt_installed_profile(
    profiles_module, tmp_path
):
    registry, payloads = make_registry()
    opener = FixtureOpener(payloads)
    acquire(profiles_module, tmp_path, registry, opener)
    calls_after_install = len(opener.calls)

    acquire(profiles_module, tmp_path, registry, opener)
    assert len(opener.calls) == calls_after_install

    target = tmp_path / "twist_human_methylome_hg38/1.0.0/genome/hg38.fa"
    target.write_text(">chr1\nAAAA\n")
    with pytest.raises(SystemExit, match="installed immutable profile.*invalid"):
        acquire(profiles_module, tmp_path, registry, opener)
    assert target.read_text() == ">chr1\nAAAA\n"
    assert len(opener.calls) == calls_after_install


def test_managed_acquisition_rejects_path_escape_before_network(
    profiles_module, tmp_path
):
    registry, payloads = make_registry(escaped_path=True)
    opener = FixtureOpener(payloads)

    with pytest.raises(SystemExit, match="genome_fa.*path.*relative"):
        acquire(profiles_module, tmp_path, registry, opener)

    assert opener.calls == []


def test_managed_acquisition_reports_offline_failure_and_cleans_stage(
    profiles_module, tmp_path
):
    registry, _ = make_registry()

    def offline(request, timeout=None):
        raise URLError("network unavailable")

    with pytest.raises(SystemExit, match="could not download.*genome_fa.*network unavailable"):
        acquire(profiles_module, tmp_path, registry, offline)
    assert not list(tmp_path.glob(".cftk-install-*"))


def test_managed_acquisition_cleans_stage_when_interrupted(profiles_module, tmp_path):
    registry, _ = make_registry()

    def interrupted(request, timeout=None):
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        acquire(profiles_module, tmp_path, registry, interrupted)
    assert not list(tmp_path.glob(".cftk-install-*"))
    assert not (tmp_path / "twist_human_methylome_hg38/1.0.0").exists()


def test_managed_acquisition_refuses_registry_redefinition(
    profiles_module, tmp_path
):
    registry, payloads = make_registry()
    opener = FixtureOpener(payloads)
    acquire(profiles_module, tmp_path, registry, opener)
    registry["profiles"]["twist_human_methylome_hg38"]["1.0.0"][
        "components"
    ]["target_bed"]["source"]["name"] = "Redefined source"

    with pytest.raises(SystemExit, match="differs from the registry.*refusing"):
        acquire(profiles_module, tmp_path, registry, opener)
    assert len(opener.calls) == 4


def test_concurrent_managed_acquisition_downloads_once(profiles_module, tmp_path):
    registry, payloads = make_registry()
    opener = FixtureOpener(payloads, delay=0.02)
    results = []
    errors = []

    def install():
        try:
            results.append(acquire(profiles_module, tmp_path, registry, opener))
        except BaseException as exc:
            errors.append(exc)

    threads = [threading.Thread(target=install) for _ in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=5)

    assert errors == []
    assert len(results) == 2
    assert len(opener.calls) == 4
    assert results[0]["component_hashes"] == results[1]["component_hashes"]


def test_noninteractive_managed_init_uses_default_root_and_records_provenance(
    init_module, tmp_path, monkeypatch
):
    registry, payloads = make_registry()
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(registry))
    reference_root = tmp_path / "managed-references"
    monkeypatch.setenv("CFTK_REFERENCE_REGISTRY", str(registry_path))
    monkeypatch.setenv("CFTK_REFERENCE_ROOT", str(reference_root))
    monkeypatch.setattr("reference_profiles.urlopen", FixtureOpener(payloads))
    sheet = make_sample_sheet(tmp_path)
    config_path = tmp_path / "cftk_init.json"
    args = SimpleNamespace(
        config=str(config_path),
        non_interactive=True,
        sample_sheet=str(sheet),
        reference_root=None,
        reference_mode=None,
        profile="twist_human_methylome_hg38",
        profile_version="1.0.0",
        project_name=None,
        output_dir=None,
        assay="twist_human_methylome",
        genome="hg38",
        skip_reference_prep=True,
        ref_index=False,
        ref_dict=False,
    )

    init_module.init(args)

    compact = json.loads(config_path.read_text())
    lock = json.loads((tmp_path / "cftk.lock.json").read_text())
    assert compact["reference_mode"] == "managed"
    assert compact["reference_root"] == str(reference_root)
    assert lock["reference_profile"]["acquisition"]["mode"] == "managed"
    assert len(lock["reference_profile"]["acquisition"]["registry_entry_sha256"]) == 64
