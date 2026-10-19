"""Regression tests for sequence-driven production workflow handoffs (no HPC required)."""
import asyncio
from datetime import datetime
from pathlib import Path

import pytest

from fade.tools.ssh_client import HPCSSHClient
from fade.tools.hpc_diffsbdd import HPCDiffSBDDClient
from fade.agents.invention.langgraph_nodes import _format_pocket_residues, _apply_simple_filters
from fade.agents.targeting.langgraph_nodes import determine_pocket_detection_strategy
from fade.workflows.drug_discovery import entry_router_node


def test_local_hpc_transport(monkeypatch, tmp_path):
    monkeypatch.setenv("FADE_LOCAL_HPC", "1")
    client = HPCSSHClient()
    assert client.connect()
    assert client.upload_content("hello", str(tmp_path / "source"))
    assert client.upload_file(str(tmp_path / "source"), str(tmp_path / "dest"))
    assert client.download_file(str(tmp_path / "dest"), str(tmp_path / "download"))
    assert (tmp_path / "download").read_text() == "hello"
    assert client.execute_command("printf worked; exit 7") == ("worked", "", 7)


def test_explicit_sequence_routes_without_llm():
    state = {"target_info": {"sequence": "ACDE"}, "query": "CRBP1", "should_continue": True}
    assert entry_router_node(state)["routing_decision"] == "structure_resolver"


def test_fpocket_override_and_source_case():
    state = {"structure": {"source": "Boltz2"}, "pocket_detection_mode": "fpocket"}
    assert determine_pocket_detection_strategy(state)["tools"] == ["fpocket"]
    assert not determine_pocket_detection_strategy(state)["use_consensus"]
    state = {"structure": {"source": "PDB", "has_ligand": True}}
    assert not determine_pocket_detection_strategy(state)["use_consensus"]


def test_residues_preserved_and_no_kras_fallback():
    pocket = {"residues": ["B:17", "B:52", "B:17"]}
    assert _format_pocket_residues(pocket) == ["B:17", "B:52"]
    assert HPCDiffSBDDClient()._format_residues_for_diffsbdd(pocket) == ["B:17", "B:52"]
    for invalid in ({}, {"residues": ["GLY12"]}):
        with pytest.raises(ValueError):
            _format_pocket_residues(invalid)
        with pytest.raises(ValueError):
            HPCDiffSBDDClient()._format_residues_for_diffsbdd(invalid)


def test_rdkit_filter():
    assert not _apply_simple_filters({"smiles": "invalid"})["passed"]
    assert not _apply_simple_filters({"smiles": "CCO"})["passed"]
    result = _apply_simple_filters({"smiles": "CNc1nc2c(OC)cccc2s1"})
    assert result["passed"]
    assert result["properties"]["mol_weight"] == pytest.approx(194.26, abs=0.02)
    assert result["properties"]["sa_score"] is None


def test_boltz2_wait_race_and_failure(monkeypatch, tmp_path):
    """The Boltz2 SLURM script prints its exit-code marker BEFORE copying the CIF
    (observed gap ~22 s). wait_for_structure must outlast that grace window and
    recover the CIF, yet still return None promptly on a genuine failure.
    Pure in-process: fake clock, patched client I/O — no SLURM, no polling."""
    import time as wall
    from fade.tools import hpc_boltz2 as hb

    client = hb.get_hpc_boltz2_client()
    monkeypatch.setattr(client, "_connect", lambda: None)
    monkeypatch.setattr(client, "_disconnect", lambda: None)
    monkeypatch.setattr(client, "remote_base", str(tmp_path))

    remote_dir = f"{tmp_path}/race"
    out_file = f"{remote_dir}/boltz2.out"
    cif = f"{remote_dir}/output/predicted_structure.cif"
    marker = "Boltz2 exit code: 0\n"
    cif_content = "data_test\n_atom_site.Cartn_x\n"

    clock = {"t": 0.0}
    monkeypatch.setattr(hb.time, "sleep", lambda s: clock.__setitem__("t", clock["t"] + s))
    monkeypatch.setattr(hb.time, "time", lambda: clock["t"])

    def run(appear_at):
        clock["t"] = 0.0

        def exists(path):
            if path == out_file:
                return True
            if path == cif:
                return clock["t"] >= appear_at
            return False

        monkeypatch.setattr(client, "_file_exists", exists)
        monkeypatch.setattr(client, "_download_file",
                            lambda p: marker if p == out_file else (cif_content if p == cif else ""))
        return client.wait_for_structure(job_id="race", max_wait=600, poll_interval=5)

    # Race regression: CIF appears 25 s after the marker -> must be recovered.
    assert run(25) == cif_content
    # Genuine failure: CIF never appears -> None after the ~60 s grace window,
    # measured instantly on the real wall clock.
    start = wall.monotonic()
    assert run(10**9) is None
    assert wall.monotonic() - start < 30, "genuine failure must not wait out the full timeout"


def test_cif_conversion_in_structure_wait(monkeypatch, tmp_path):
    from Bio.PDB import PDBParser, MMCIFIO
    from io import StringIO
    from fade.agents.structure import langgraph_nodes as nodes
    pdb = "ATOM      1  CA  ALA A   1      10.000  11.000  12.000  1.00 90.00           C  \nEND\n"
    writer = MMCIFIO()
    writer.set_structure(PDBParser(QUIET=True).get_structure("x", StringIO(pdb)))
    output = StringIO()
    writer.save(output)
    class Client:
        def wait_for_structure(self, **kwargs):
            return output.getvalue()
    monkeypatch.setattr(nodes, "get_boltz2_client", Client)
    monkeypatch.setattr(nodes, "_save_structure", lambda content, name: (tmp_path / name).write_text(content) and str(tmp_path / name))
    result = nodes.structure_wait_node({"structure": {"source": "Boltz2", "status": "pending", "boltz2_job_id": "test"}})
    assert Path(result["structure"]["structure_path"]).read_text().startswith("ATOM")
    assert result["structure"]["confidence_score"] is None


def test_graph_retains_stage_metadata(monkeypatch, tmp_path):
    from fade.workflows import drug_discovery as workflow
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(workflow, "structure_resolver_node", lambda s: {"structure": {"source": "Boltz2", "status": "pending"}})
    monkeypatch.setattr(workflow, "structure_wait_node", lambda s: {"structure": {"source": "Boltz2", "status": "completed"}})
    monkeypatch.setattr(workflow, "pocket_mapper_node", lambda s: {"pockets": [{"pocket_id": "p1"}], "pocket_detection_metadata": {"tools": ["fpocket"]}})
    monkeypatch.setattr(workflow, "pocket_ranker_node", lambda s: {"selected_pocket": {"pocket_id": "p1"}})
    monkeypatch.setattr(workflow, "molecule_generator_node", lambda s: {"generated_molecules": [{"smiles": "CCO"}], "generation_parameters": {"n_requested": 5}})
    monkeypatch.setattr(workflow, "medicinal_filter_node", lambda s: {"filtered_molecules": [{"smiles": "CCO"}], "filter_statistics": {"passed": 1}})
    monkeypatch.setattr(workflow, "boltz2_screener_node", lambda s: {"screening_results": [{"smiles": "CCO"}], "screening_metadata": {"n_screened": 1}})
    monkeypatch.setattr(workflow, "results_analyzer_node", lambda s: {"final_report": "test report"})
    async def run():
        app = await workflow.create_drug_discovery_graph()
        try:
            return await workflow.run_drug_discovery_pipeline(app, "test", target_info={"sequence": "ACDE"}, force_prediction=True, pocket_detection_mode="fpocket")
        finally:
            await app.checkpointer.conn.close()
    result = asyncio.run(run())
    assert not result.get("error")
    assert result["generation_parameters"]["n_requested"] == 5
    assert result["screening_metadata"]["n_screened"] == 1
    assert result["force_prediction"]
