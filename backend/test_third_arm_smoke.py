"""
Local structural smoke test for the F.A.D.E third arm (Boltz-2 affinity screening).

This test does NOT touch the HPC cluster or run Boltz-2. It:
  1. Stubs out heavy / network dependencies (paramiko, dotenv, etc.),
  2. Imports and compiles the LangGraph workflow,
  3. Verifies the new screening/analysis nodes and the rewired edges are present,
  4. Unit-tests the pure-Python affinity helpers (SMILES atom count, JSON parse,
     IC50->kcal conversion) and the screener/analyzer nodes with a MOCKED
     Boltz-2 client (no SLURM, no network).

Run: python test_third_arm_smoke.py
Exit code 0 = all checks passed.
"""

import sys
import os
import types
from pathlib import Path

# Ensure the backend package is importable.
BACKEND = Path(__file__).resolve().parent
sys.path.insert(0, str(BACKEND))

FAILURES = []


def check(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


# --------------------------------------------------------------------------- #
# 1. Stub heavy / network dependencies so imports succeed offline.
# --------------------------------------------------------------------------- #
def _install_stubs():
    # paramiko (SSH) - only needs to be importable.
    if "paramiko" not in sys.modules:
        paramiko = types.ModuleType("paramiko")
        paramiko.SSHClient = object
        paramiko.AutoAddPolicy = object
        paramiko.RSAKey = object
        sys.modules["paramiko"] = paramiko

    # python-dotenv is installed, but guard anyway.
    if "dotenv" not in sys.modules:
        try:
            import dotenv  # noqa
        except Exception:
            dotenv = types.ModuleType("dotenv")
            dotenv.load_dotenv = lambda *a, **k: None
            sys.modules["dotenv"] = dotenv

    # LLM provider packages are not needed for a structural test — stub them with
    # a dummy chat class so `from langchain_openai import ChatOpenAI` etc. succeed.
    class _DummyChat:
        def __init__(self, *a, **k):
            pass

        def invoke(self, *a, **k):
            raise RuntimeError("LLM disabled in smoke test")

    _llm_stub_specs = {
        "langchain_openai": ["ChatOpenAI", "OpenAIEmbeddings"],
        "langchain_anthropic": ["ChatAnthropic"],
        "langchain_google_genai": ["ChatGoogleGenerativeAI", "GoogleGenerativeAIEmbeddings"],
        "langchain_google_vertexai": ["ChatVertexAI"],
        "langchain_ollama": ["ChatOllama"],
    }
    for mod_name, attrs in _llm_stub_specs.items():
        if mod_name not in sys.modules:
            try:
                __import__(mod_name)
            except Exception:
                mod = types.ModuleType(mod_name)
                for attr in attrs:
                    setattr(mod, attr, _DummyChat)
                sys.modules[mod_name] = mod


# --------------------------------------------------------------------------- #
# Test 1: Pure-Python affinity helpers (no graph import needed).
# --------------------------------------------------------------------------- #
def test_affinity_helpers():
    print("\n[Test 1] Boltz-2 affinity helper functions")
    _install_stubs()
    os.environ.setdefault("AI_MODEL", "gpt-4o")  # config needs a model name
    from fade.tools.hpc_boltz2 import HPCBoltz2Client

    client = HPCBoltz2Client.__new__(HPCBoltz2Client)  # skip __init__ (no SSH creds)

    # IC50 -> kcal/mol: 1 uM should give a negative dG around -8.2 kcal/mol.
    kcal = HPCBoltz2Client._ic50_to_kcal(1.0)
    check("ic50_to_kcal(1uM) is negative & plausible", kcal is not None and -12 < kcal < -6,
          f"got {kcal}")

    # Parse a representative Boltz-2 affinity JSON payload.
    raw = '{"affinity_pred_value": -1.5, "affinity_probability_binary": 0.87}'
    parsed = client._parse_affinity_json(raw)
    check("parse affinity_pred_value", parsed["affinity_pred_value"] == -1.5)
    check("parse affinity_probability_binary", parsed["affinity_probability_binary"] == 0.87)
    check("derived predicted_ic50_uM present", parsed["predicted_ic50_uM"] is not None)
    check("derived binding_affinity_kcal_mol present", parsed["binding_affinity_kcal_mol"] is not None)

    # Malformed JSON should degrade gracefully.
    bad = client._parse_affinity_json("not json{")
    check("malformed JSON handled", bad["affinity_pred_value"] is None and "parse_error" in bad)


# --------------------------------------------------------------------------- #
# Test 2: Compile the graph and inspect structure.
# --------------------------------------------------------------------------- #
def test_graph_structure():
    print("\n[Test 2] Compile LangGraph workflow & inspect third-arm wiring")
    _install_stubs()
    os.environ.setdefault("AI_MODEL", "gpt-4o")

    import asyncio
    from fade.workflows.drug_discovery import create_drug_discovery_graph

    app = asyncio.get_event_loop().run_until_complete(create_drug_discovery_graph())
    graph = app.get_graph()
    nodes = set(graph.nodes.keys())
    edges = [(e.source, e.target) for e in graph.edges]

    # New nodes registered?
    check("node 'boltz2_screener' present", "boltz2_screener" in nodes)
    check("node 'results_analyzer' present", "results_analyzer" in nodes)
    check("node 'molecule_generator' present", "molecule_generator" in nodes)
    check("node 'medicinal_filter' present", "medicinal_filter" in nodes)

    # Rewired edges present? (conditional edges appear as source->target pairs)
    def has_edge(src, dst):
        return any(s == src and t == dst for s, t in edges)

    check("edge pocket_ranker -> molecule_generator", has_edge("pocket_ranker", "molecule_generator"))
    check("edge pocket_from_ligand -> molecule_generator", has_edge("pocket_from_ligand", "molecule_generator"))
    check("edge molecule_generator -> medicinal_filter", has_edge("molecule_generator", "medicinal_filter"))
    check("edge medicinal_filter -> boltz2_screener", has_edge("medicinal_filter", "boltz2_screener"))
    check("edge boltz2_screener -> results_analyzer", has_edge("boltz2_screener", "results_analyzer"))
    check("edge results_analyzer -> END", has_edge("results_analyzer", "__end__"))

    # Third-arm reachability: structure_wait (Boltz-2 predicted) must connect through
    # to the shared screening tail.
    check("edge structure_wait -> pocket_mapper (de-novo arm)", has_edge("structure_wait", "pocket_mapper"))

    # Full third-arm path check via simple reachability from structure_wait.
    adj = {}
    for s, t in edges:
        adj.setdefault(s, set()).add(t)
    seen, stack = set(), ["structure_wait"]
    while stack:
        n = stack.pop()
        if n in seen:
            continue
        seen.add(n)
        stack.extend(adj.get(n, ()))
    for required in ["pocket_mapper", "pocket_ranker", "molecule_generator",
                     "medicinal_filter", "boltz2_screener", "results_analyzer"]:
        check(f"third arm reaches '{required}'", required in seen)


# --------------------------------------------------------------------------- #
# Test 3: Run screener + analyzer nodes with a MOCKED Boltz-2 client.
# --------------------------------------------------------------------------- #
def test_nodes_mocked():
    print("\n[Test 3] Screener + analyzer nodes with mocked Boltz-2 (no HPC)")
    _install_stubs()
    os.environ.setdefault("AI_MODEL", "gpt-4o")

    import fade.tools as tools_pkg
    from fade.agents.screening import langgraph_nodes as screening_mod
    from fade.agents.analysis import langgraph_nodes as analysis_mod

    # Build a fake Boltz-2 client: deterministic affinity per ligand.
    class FakeBoltz2:
        def __init__(self):
            self._n = 0

        def predict_affinity(self, sequence, ligand_smiles, job_name=None, **kw):
            self._n += 1
            return {"job_id": f"fake_{self._n}", "status": "submitted"}

        def wait_for_affinity(self, job_id, max_wait=None, poll_interval=None):
            # Stronger (more negative) for the first ligand, weaker for later ones.
            idx = int(job_id.split("_")[1])
            pred_value = -2.0 + idx * 0.5
            prob = max(0.1, 0.95 - idx * 0.2)
            return {
                "affinity_pred_value": pred_value,
                "affinity_probability_binary": prob,
                "binding_affinity_kcal_mol": -9.0 + idx,
                "predicted_ic50_uM": 10 ** pred_value,
                "pose_path": None,
                "raw": {},
            }

    fake = FakeBoltz2()
    # Patch the getter used inside the screening node.
    screening_mod.get_boltz2_client = lambda: fake

    state = {
        "run_id": "smoke",
        "query": "Find inhibitors for TESTKINASE",
        "target_info": {"protein_name": "TESTKINASE", "sequence": "MKVLAAGGSTKPLRDEF" * 3},
        "structure": {"source": "Boltz2", "structure_path": "/tmp/fake.pdb"},
        "selected_pocket": {"pocket_id": "pocket_1", "druggability_score": 0.8},
        "filtered_molecules": [
            {"molecule_id": "mol_A", "smiles": "CCO", "qed_score": 0.7},
            {"molecule_id": "mol_B", "smiles": "c1ccccc1O", "qed_score": 0.6},
            {"molecule_id": "mol_C", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "qed_score": 0.55},
        ],
        "should_continue": True,
    }

    screen_out = screening_mod.boltz2_screener_node(state)
    results = screen_out.get("screening_results", [])
    check("screener produced results", len(results) == 3, f"got {len(results)}")
    check("results ranked best-first (rank 1 lowest pred_value)",
          results and results[0]["rank"] == 1 and
          results[0]["affinity_pred_value"] <= results[-1]["affinity_pred_value"])
    check("ScreeningResult has required schema fields",
          all(k in results[0] for k in ("molecule_id", "smiles", "binding_affinity",
                                         "affinity_confidence", "pose_path")))
    check("likely_binder flag computed", "likely_binder" in results[0])
    check("routing helper returns continue",
          screening_mod.should_continue_screening({**state, **screen_out}) == "continue")

    # Feed into analyzer.
    merged = {**state, **screen_out}
    analysis_out = analysis_mod.results_analyzer_node(merged)
    analysis = analysis_out.get("analysis", {})
    check("analyzer set top_molecule", analysis.get("top_molecule") is not None)
    check("analyzer top_molecule is rank 1",
          analysis.get("top_molecule", {}).get("molecule_id") == results[0]["molecule_id"])
    check("analyzer produced final_report markdown",
          isinstance(analysis_out.get("final_report"), str) and
          "F.A.D.E Drug Discovery Report" in analysis_out["final_report"])
    check("analyzer top_10 non-empty", len(analysis.get("top_10_molecules", [])) == 3)


def main():
    print("=" * 70)
    print("F.A.D.E Third Arm — Local Structural Smoke Test (no HPC, mocked)")
    print("=" * 70)
    test_affinity_helpers()
    test_graph_structure()
    test_nodes_mocked()

    print("\n" + "=" * 70)
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} check(s) FAILED: {FAILURES}")
        sys.exit(1)
    print("RESULT: ALL CHECKS PASSED")
    print("=" * 70)
    sys.exit(0)


if __name__ == "__main__":
    main()
