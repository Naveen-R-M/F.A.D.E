"""
Analysis module LangGraph nodes.

Compiles the ranked Boltz-2 screening results into a final drug-discovery report.
This is the terminal stage shared by all three F.A.D.E branches.
"""

import json
from typing import Dict, Any, List, Optional
from datetime import datetime

from langchain_core.messages import HumanMessage

from fade.state.langgraph_state import DrugDiscoveryState
from fade.config import config
from fade.utils import get_logger

logger = get_logger("nodes.analysis")

__all__ = ["results_analyzer_node"]


def _fmt(value: Any, spec: str = "") -> str:
    """Format a possibly-None numeric value for the report."""
    if value is None:
        return "N/A"
    try:
        return format(value, spec) if spec else str(value)
    except Exception:
        return str(value)


def results_analyzer_node(state: DrugDiscoveryState) -> Dict[str, Any]:
    """
    Analyze Boltz-2 screening results and produce the final report.

    Builds an ``analysis`` dict (top candidate, top-10, key insights) plus a
    human-readable ``final_report`` string, and writes a JSON + markdown report to
    the run's output directory.

    Args:
        state: Current pipeline state with ``screening_results``.

    Returns:
        State updates with ``analysis`` and ``final_report``.
    """
    screening_results = state.get("screening_results") or []

    if not screening_results:
        return {
            "error": "No screening results available for analysis",
            "should_continue": False,
            "current_step": "analysis_failed",
        }

    logger.info(f"[Analyzer] Analyzing {len(screening_results)} screening result(s)")

    # Results arrive ranked from the screener; guard by re-sorting defensively.
    def _rank_key(r: Dict[str, Any]):
        pred = r.get("affinity_pred_value")
        prob = r.get("affinity_probability_binary")
        return (pred if pred is not None else float("inf"),
                -(prob if prob is not None else -1.0))

    ranked = sorted(screening_results, key=_rank_key)

    top_molecule = ranked[0]
    top_10 = ranked[:10]
    likely_binders = [r for r in ranked if r.get("likely_binder")]

    target_info = state.get("target_info") or {}
    target_name = (
        target_info.get("protein_name")
        or target_info.get("gene_name")
        or target_info.get("uniprot_id")
        or "the target"
    )
    structure = state.get("structure") or {}
    structure_source = structure.get("source", "unknown")
    selected_pocket = state.get("selected_pocket") or {}

    # Novel scaffolds = the SMILES of the top candidates (de-novo generated).
    novel_scaffolds = [r.get("smiles") for r in top_10 if r.get("smiles")]

    key_insights = [
        f"Screened {len(screening_results)} candidate ligand(s) against {target_name} "
        f"using Boltz-2 binding-affinity prediction.",
        f"{len(likely_binders)} candidate(s) exceeded the binder-probability threshold "
        f"({getattr(config, 'BOLTZ2_AFFINITY_BINDER_PROBABILITY_THRESHOLD', 0.5)}).",
        f"Structure source for this run: {structure_source} "
        f"(branch determines how the pocket was found, screening is identical).",
    ]
    if top_molecule.get("affinity_pred_value") is not None:
        key_insights.append(
            f"Top candidate {top_molecule.get('molecule_id')} has predicted "
            f"affinity value {_fmt(top_molecule.get('affinity_pred_value'), '.3f')} "
            f"(log10 IC50, lower is stronger) and P(binder)="
            f"{_fmt(top_molecule.get('affinity_probability_binary'), '.3f')}."
        )

    analysis: Dict[str, Any] = {
        "top_molecule": top_molecule,
        "top_10_molecules": top_10,
        "novel_scaffolds": novel_scaffolds,
        "key_insights": key_insights,
        "comparison_with_known": {},
        "synthesis_feasibility": {},
        "n_screened": len(screening_results),
        "n_likely_binders": len(likely_binders),
    }

    final_report = _build_markdown_report(
        state=state,
        target_name=target_name,
        structure_source=structure_source,
        selected_pocket=selected_pocket,
        ranked=ranked,
        likely_binders=likely_binders,
        key_insights=key_insights,
    )

    # Persist report artifacts to the run's output directory.
    report_path = None
    try:
        run_id = state.get("run_id", "run")
        reports_dir = config.OUTPUT_DIR / "reports" / str(run_id)
        reports_dir.mkdir(parents=True, exist_ok=True)

        md_path = reports_dir / "final_report.md"
        with open(md_path, "w") as f:
            f.write(final_report)

        json_path = reports_dir / "screening_results.json"
        with open(json_path, "w") as f:
            json.dump(
                {
                    "target": target_name,
                    "structure_source": structure_source,
                    "pocket_id": selected_pocket.get("pocket_id"),
                    "generated_at": datetime.now().isoformat(),
                    "screening_metadata": state.get("screening_metadata"),
                    "results": ranked,
                },
                f,
                indent=2,
                default=str,
            )

        report_path = str(md_path)
        analysis["report_path"] = report_path
        logger.info(f"[Analyzer] Wrote report to {md_path}")
    except Exception as e:
        logger.warning(f"[Analyzer] Could not write report files: {e}")

    message = (
        f"Analysis complete. Top candidate: {top_molecule.get('molecule_id')} "
        f"(P(binder)={_fmt(top_molecule.get('affinity_probability_binary'), '.3f')}). "
        f"{len(likely_binders)} likely binder(s) of {len(screening_results)} screened."
    )
    if report_path:
        message += f"\nReport: {report_path}"

    return {
        "analysis": analysis,
        "final_report": final_report,
        "messages": [HumanMessage(content=message)],
        "current_step": "analysis_completed",
        "should_continue": True,
    }


def _build_markdown_report(state: DrugDiscoveryState,
                           target_name: str,
                           structure_source: str,
                           selected_pocket: Dict[str, Any],
                           ranked: List[Dict[str, Any]],
                           likely_binders: List[Dict[str, Any]],
                           key_insights: List[str]) -> str:
    """Assemble the human-readable markdown report."""
    lines: List[str] = []
    lines.append(f"# F.A.D.E Drug Discovery Report: {target_name}")
    lines.append("")
    lines.append(f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*")
    lines.append("")
    lines.append(f"**Query:** {state.get('query', 'N/A')}")
    lines.append("")

    lines.append("## Pipeline Summary")
    lines.append("")
    lines.append(f"- **Target:** {target_name}")
    lines.append(f"- **Structure source:** {structure_source}")
    if selected_pocket:
        lines.append(
            f"- **Selected pocket:** {selected_pocket.get('pocket_id', 'N/A')} "
            f"(druggability {_fmt(selected_pocket.get('druggability_score'), '.2f')})"
        )
    meta = state.get("screening_metadata") or {}
    lines.append(f"- **Molecules screened:** {meta.get('n_screened', len(ranked))}")
    lines.append(f"- **Likely binders:** {len(likely_binders)}")
    lines.append(f"- **Screening method:** Boltz-2 binding-affinity prediction")
    lines.append("")

    lines.append("## Key Insights")
    lines.append("")
    for insight in key_insights:
        lines.append(f"- {insight}")
    lines.append("")

    lines.append("## Top Candidates")
    lines.append("")
    lines.append("Ranked by Boltz-2 predicted affinity value (lower log10 IC50 = stronger); "
                 "binder probability breaks ties.")
    lines.append("")
    lines.append("| Rank | Molecule | P(binder) | Affinity value (log10 IC50) | Pred. IC50 (µM) | QED |")
    lines.append("|------|----------|-----------|------------------------------|-----------------|-----|")
    for r in ranked[:10]:
        lines.append(
            f"| {r.get('rank', '?')} "
            f"| {r.get('molecule_id', 'N/A')} "
            f"| {_fmt(r.get('affinity_probability_binary'), '.3f')} "
            f"| {_fmt(r.get('affinity_pred_value'), '.3f')} "
            f"| {_fmt(r.get('predicted_ic50_uM'), '.3g')} "
            f"| {_fmt(r.get('qed_score'), '.2f')} |"
        )
    lines.append("")

    # Show the SMILES of the very top candidates for convenience.
    lines.append("### Top candidate structures (SMILES)")
    lines.append("")
    for r in ranked[:5]:
        lines.append(f"- **{r.get('molecule_id', 'N/A')}**: `{r.get('smiles', 'N/A')}`")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("*Affinity values are Boltz-2 predictions (log10 IC50); the derived kcal/mol "
                 "and IC50 columns are informational conversions and are assay-dependent. "
                 "Use predicted affinity and binder probability for ranking, and validate "
                 "top candidates experimentally.*")

    return "\n".join(lines)
