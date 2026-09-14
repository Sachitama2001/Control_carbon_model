import json

from control_carbon.visitc_source_map import (
    VISITC_HYDROLOGY_SOURCE,
    VISITC_SOURCE_COMMIT,
    VISITC_SOURCE_REPOSITORY,
    VISITC_SOURCES,
    visitc_provenance_manifest,
)


def test_visitc_manifest_is_separate_and_commit_pinned():
    manifest = visitc_provenance_manifest()

    assert manifest["model"] == "VISITc"
    assert len(VISITC_SOURCE_COMMIT) == 40
    assert all(
        source.repository == VISITC_SOURCE_REPOSITORY
        and source.commit == VISITC_SOURCE_COMMIT
        for source in VISITC_SOURCES
    )
    assert VISITC_HYDROLOGY_SOURCE.path == "point/hydro_balance.c"
    assert json.loads(json.dumps(manifest)) == manifest

