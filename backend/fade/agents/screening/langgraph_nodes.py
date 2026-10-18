"""
Screening module LangGraph nodes.

This module contains the SHARED binding-affinity screening engine used by all three
branches of the F.A.D.E pipeline (holo / apo / de-novo). Regardless of how a branch
obtained its protein structure and selected pocket, by the time the pipeline reaches
screening the relevant state is identical:

    - target_info["sequence"]      : the protein sequence
    - structure                     : the resolved protein structure
    - selected_pocket               : the chosen binding pocket
    - filtered_molecules            : drug-like candidate ligands (with SMILES)

So a single ``boltz2_screener_node`` serves every branch; there is no per-arm
screening code. It uses Boltz-2's protein-ligand affinity mode to rank the
candidates by predicted binding.

Boltz-2 scores exactly one ligand per job, so screening a batch means submitting
one HPC job per molecule.
"""

from typing import Dict, Any, List, Optional

from langchain_core.messages import HumanMessage

from fade.state.langgraph_state import DrugDiscoveryState
from fade.tools import get_boltz2_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("nodes.screening")

__all__ = ["boltz2_screener_node", "should_continue_screening"]


def _resolve_protein_sequence(state: DrugDiscoveryState) -> Optional[str]:
    """
    Get the protein sequence for affinity screening.

    Prefers the sequence already stored on target_info; falls back to fetching it
    from UniProt if only an accession is available (same resolution the structure
    resolver uses).
    """
    target_info = state.get("target_info") or {}
    sequence = target_info.get("sequence")
    if sequence:
        return sequence

    uniprot_id = target_info.get("uniprot_id")
    if uniprot_id:
        try:
            from fade.tools.uniprot_api import get_uniprot_client
            sequence = get_uniprot_client().get_protein_sequence(uniprot_id)
            if sequence:
                logger.info(f"[Screener] Fetched sequence for {uniprot_id} ({len(sequence)} aa)")
                return sequence
        except Exception as e:
            logger.warning(f"[Screener] Could not fetch sequence for {uniprot_id}: {e}")

    return None


def _estimate_heavy_atoms(smiles: str) -> int:
    """
    Rough heavy-atom count from a SMILES string, used only to skip pathologically
    large ligands before spending an HPC job on them.

    Prefers RDKit when available; otherwise counts non-hydrogen atomic symbols
    approximately. This is a guard, not a precise cheminformatics measurement.
    """
    if not smiles:
        return 0
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return mol.GetNumHeavyAtoms()
    except Exception:
        pass

    # Approximate: count upper-case organic subset letters + bracketed atoms,
    # excluding explicit hydrogens.
    import re
    tokens = re.findall(r'Br|Cl|[BCNOPSFIbcnops]|\[[^\]]+\]', smiles)
    heavy = 0
    for t in tokens:
        if t.startswith('['):
            # bracket atom, e.g. [nH], [O-]; skip pure hydrogen
            inner = t.strip('[]')
            if inner.upper().startswith('H') and inner not in ('Hf', 'Hg', 'He'):
                continue
            heavy += 1
        else:
            heavy += 1
    return heavy


def boltz2_screener_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    Shared binding-affinity screening node (all branches).

    Runs Boltz-2 protein-ligand affinity prediction for each filtered candidate
    molecule and ranks them by predicted binding. One Boltz-2 job is submitted per
    ligand (Boltz-2 scores a single ligand per job).

    Args:
        state: Current pipeline state (converged from any of the three branches).

    Returns:
        State updates with ``screening_results`` (list of ScreeningResult dicts,
        ranked best-first) and ``screening_metadata``.
    """
    if hasattr(config, "BOLTZ2_ENABLE_AFFINITY") and not config.BOLTZ2_ENABLE_AFFINITY:
        logger.info("[Screener] Boltz-2 affinity screening disabled via config; skipping.")
        return {
            "screening_results": [],
            "screening_metadata": {"skipped": True, "reason": "BOLTZ2_ENABLE_AFFINITY is false"},
            "messages": [HumanMessage(content="Binding-affinity screening skipped (disabled in config).")],
            "current_step": "screening_skipped",
            "should_continue": True,
        }

    molecules = state.get("filtered_molecules") or []
    if not molecules:
        return {
            "error": "No filtered molecules available for screening",
            "should_continue": False,
            "current_step": "screening_failed",
        }

    sequence = _resolve_protein_sequence(state)
    if not sequence:
        return {
            "error": "No protein sequence available for binding-affinity screening",
            "should_continue": False,
            "current_step": "screening_failed",
        }

    # Cap how many ligands we screen (one HPC job each).
    max_ligands = getattr(config, "BOLTZ2_AFFINITY_MAX_LIGANDS", 25)
    max_atoms = getattr(config, "BOLTZ2_AFFINITY_MAX_LIGAND_ATOMS", 150)
    binder_threshold = getattr(config, "BOLTZ2_AFFINITY_BINDER_PROBABILITY_THRESHOLD", 0.5)

    # Select candidates with usable SMILES, respecting the size guard and the cap.
    candidates: List[Dict[str, Any]] = []
    skipped_no_smiles = 0
    skipped_too_large = 0
    for mol in molecules:
        smiles = mol.get("smiles")
        if not smiles:
            skipped_no_smiles += 1
            continue
        if _estimate_heavy_atoms(smiles) > max_atoms:
            skipped_too_large += 1
            continue
        candidates.append(mol)
        if len(candidates) >= max_ligands:
            break

    if not candidates:
        return {
            "error": "No screenable ligands (all lacked SMILES or exceeded the size guard)",
            "should_continue": False,
            "current_step": "screening_failed",
        }

    logger.info(
        f"[Screener] Screening {len(candidates)} ligand(s) against target "
        f"({len(sequence)} aa) with Boltz-2 affinity"
    )

    boltz2_client = get_boltz2_client()
    selected_pocket = state.get("selected_pocket") or {}
    pocket_id = selected_pocket.get("pocket_id")

    screening_results: List[Dict[str, Any]] = []
    n_failed = 0

    for idx, mol in enumerate(candidates):
        smiles = mol["smiles"]
        molecule_id = mol.get("molecule_id") or mol.get("id") or f"mol_{idx + 1}"

        try:
            submit = boltz2_client.predict_affinity(
                sequence=sequence,
                ligand_smiles=smiles,
                job_name=f"fade_screen_{molecule_id}",
            )
            job_id = submit.get("job_id")
            if not job_id:
                logger.warning(f"[Screener] No job_id returned for {molecule_id}; skipping.")
                n_failed += 1
                continue

            prediction = boltz2_client.wait_for_affinity(job_id=job_id)
            if not prediction:
                logger.warning(f"[Screener] Affinity prediction failed/timed out for {molecule_id}.")
                n_failed += 1
                continue

            prob = prediction.get("affinity_probability_binary")
            result = {
                "molecule_id": molecule_id,
                "smiles": smiles,
                # ScreeningResult.binding_affinity is the headline number the rest of
                # the codebase reads; we store the derived kcal/mol convenience value.
                "binding_affinity": prediction.get("binding_affinity_kcal_mol"),
                "affinity_confidence": prob,
                "pose_path": prediction.get("pose_path"),
                "interactions": None,
                "rmsd_from_pocket": None,
                "strain_energy": None,
                # Boltz-2 native outputs (authoritative for ranking):
                "affinity_pred_value": prediction.get("affinity_pred_value"),
                "affinity_probability_binary": prob,
                "predicted_ic50_uM": prediction.get("predicted_ic50_uM"),
                "likely_binder": (prob is not None and prob >= binder_threshold),
                "pocket_id": pocket_id,
                "boltz2_job_id": job_id,
                # carry forward useful medicinal-chem context if present
                "qed_score": mol.get("qed_score"),
                "mol_weight": mol.get("mol_weight"),
            }
            screening_results.append(result)
            logger.info(
                f"[Screener] {molecule_id}: pred_value={result['affinity_pred_value']}, "
                f"P(binder)={prob}"
            )

        except Exception as e:
            logger.error(f"[Screener] Error screening {molecule_id}: {e}")
            n_failed += 1
            continue

    if not screening_results:
        return {
            "error": "Binding-affinity screening produced no results",
            "should_continue": False,
            "current_step": "screening_failed",
            "screening_metadata": {
                "n_candidates": len(candidates),
                "n_failed": n_failed,
                "skipped_no_smiles": skipped_no_smiles,
                "skipped_too_large": skipped_too_large,
            },
        }

    # Rank best-first. Lower affinity_pred_value = stronger binding; break ties (and
    # missing values) by higher binder probability.
    def _rank_key(r: Dict[str, Any]):
        pred = r.get("affinity_pred_value")
        prob = r.get("affinity_probability_binary")
        pred_sort = pred if pred is not None else float("inf")
        prob_sort = -(prob if prob is not None else -1.0)
        return (pred_sort, prob_sort)

    screening_results.sort(key=_rank_key)
    for rank, r in enumerate(screening_results, 1):
        r["rank"] = rank

    n_binders = sum(1 for r in screening_results if r.get("likely_binder"))
    best = screening_results[0]

    message = (
        f"Screened {len(screening_results)} molecule(s) with Boltz-2 affinity "
        f"({n_binders} flagged as likely binders).\n"
        f"Top candidate: {best['molecule_id']} "
        f"(pred_value={best.get('affinity_pred_value')}, "
        f"P(binder)={best.get('affinity_probability_binary')})"
    )
    logger.info(message)

    screening_metadata = {
        "method": "Boltz-2 affinity",
        "n_candidates": len(candidates),
        "n_screened": len(screening_results),
        "n_failed": n_failed,
        "n_likely_binders": n_binders,
        "skipped_no_smiles": skipped_no_smiles,
        "skipped_too_large": skipped_too_large,
        "binder_probability_threshold": binder_threshold,
        "pocket_id": pocket_id,
    }

    return {
        "screening_results": screening_results,
        "screening_metadata": screening_metadata,
        "messages": [HumanMessage(content=message)],
        "current_step": "screening_completed",
        "should_continue": True,
    }


def should_continue_screening(state: DrugDiscoveryState) -> str:
    """Determine if the pipeline should continue after screening."""
    if state.get("error") or not state.get("should_continue", True):
        if state.get("error"):
            logger.warning(f"Screening failed: {state.get('error')}")
        return "end"

    if not state.get("screening_results"):
        logger.warning("No screening results produced")
        return "end"

    return "continue"
