"""Content-contract checks for the read-only explorer (no new model equations)."""
from pathlib import Path
import runpy

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def explorer():
    source = ROOT.parent / "VISIT-matrix/visit_local"
    if not source.is_dir():
        pytest.skip("pinned source checkout required to build explanatory inventory")
    build = runpy.run_path(str(ROOT / "examples/build_model_explorer.py"))
    return build["build_data"](source)


def test_inventory_keeps_physical_groups_and_source_identity(explorer):
    states = [n["state"] for n in explorer["nodes"].values() if "state" in n]
    assert len(states) == 37
    assert {g: sum(s["group"] == g for s in states) for g in "CNWD"} == dict(C=21, N=11, W=3, D=2)
    assert sum(s["name"].endswith("_nsc") for s in states) == 3
    assert "x_n_mic" in {s["name"] for s in states}
    assert all(n["state"]["source"] for n in explorer["nodes"].values() if "state" in n)


def test_every_state_is_reachable_from_x_and_breadcrumbs_terminate(explorer):
    nodes = explorer["nodes"]
    def visit(key, seen):
        if key in seen:
            return
        seen.add(key)
        for child in nodes[key]["children"]:
            visit(child, seen)
    seen = set()
    visit("x", seen)
    assert {k for k, n in nodes.items() if "state" in n} <= seen
    for key in nodes:
        path = set()
        while key is not None:
            assert key not in path, "cyclic breadcrumb"
            path.add(key)
            key = nodes[key]["parent"]
        assert "home" in path
    assert nodes["fixed"]["badge"] == "説明用の仮定"
    assert nodes["compare"]["badge"] == "仮定を比べる"


def test_offline_provenance_and_no_simulated_values(explorer):
    for source in explorer["sources"].values():
        assert source["excerpt"]
        assert len(source["sha256"]) == 64
        assert source["order"] and source["level"]
    for node in explorer["nodes"].values():
        if "state" in node:
            assert "value" not in node  # Seed values are not current forest values.
    assert "nsch_storage" in explorer["sources"]["flush"]["excerpt"]
    assert "opt_lai" in explorer["sources"]["target"]["excerpt"]


def test_nested_html_contains_adopted_equations_without_script(explorer):
    from html import unescape
    from html.parser import HTMLParser
    render = runpy.run_path(str(ROOT / 'examples/render_nested_explorer.py'))['render_nested']
    html = render(explorer)
    class Inspector(HTMLParser):
        def __init__(self):
            super().__init__()
            self.ids = []
            self.tags = []
        def handle_starttag(self, tag, attrs):
            self.tags.append(tag)
            if 'id' in dict(attrs):
                self.ids.append(dict(attrs)['id'])
    parser = Inspector()
    parser.feed(html)
    assert len(parser.ids) == len(set(parser.ids))
    assert 'script' not in parser.tags
    assert 'input' not in parser.tags
    assert parser.tags.count('details') > 37
    plain = unescape(html)
    for key, node in explorer['nodes'].items():
        if 'state' in node:
            assert explorer['equations']['m_rhs_'+key[2:]]['expression'] in plain
    assert 'node-B-allocation-dynamic-target-target_formula' in parser.ids
    assert 'eq-p_t_alloc_ass' in parser.ids
