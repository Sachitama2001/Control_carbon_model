import json

import pytest

from control_carbon.provenance import SourceRef, provenance_manifest
from control_carbon.visit_source_map import (
    VISIT_LAND_COVER_NAMES,
    VISIT_PARAMETER_WORKBOOK_SHA256,
    VISIT_PARAMETER_WORKBOOK_SOURCE,
    VISIT_SOURCE_COMMIT,
    VISIT_SOURCE_REPOSITORY,
    VISIT_SOURCES,
    visit_provenance_manifest,
)


def test_source_ref_and_manifest_are_json_compatible():
    source = SourceRef(
        path="model.c",
        symbol="step",
        role="state update",
        repository="example/model",
        commit="abc123",
        function="step",
        assumptions=("fixed regime",),
        approximation_level="native-source",
    )

    manifest = provenance_manifest(
        model="example", approximation_level="reduced", sources=(source,)
    )
    assert json.loads(json.dumps(manifest)) == manifest
    assert manifest["sources"][0]["assumptions"] == ["fixed regime"]


def test_source_ref_rejects_missing_required_metadata():
    with pytest.raises(ValueError, match="path"):
        SourceRef(path="", symbol="step", role="state update")


def test_visit_manifest_pins_authoritative_revision():
    manifest = visit_provenance_manifest()

    assert manifest["model"] == "VISIT"
    assert len(manifest["sources"]) == len(VISIT_SOURCES)
    assert all(
        source["repository"] == VISIT_SOURCE_REPOSITORY
        and source["commit"] == VISIT_SOURCE_COMMIT
        for source in manifest["sources"]
    )


def test_supplementary_parameter_workbook_is_hash_pinned():
    assert len(VISIT_LAND_COVER_NAMES) == 16
    assert VISIT_LAND_COVER_NAMES[4] == "Deciduous Broadleaf Forest"
    assert len(VISIT_PARAMETER_WORKBOOK_SHA256) == 64
    assert VISIT_PARAMETER_WORKBOOK_SHA256 in VISIT_PARAMETER_WORKBOOK_SOURCE.assumptions[0]
    assert VISIT_PARAMETER_WORKBOOK_SOURCE.commit is None