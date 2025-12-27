from datetime import datetime
from typing import List, Any

# Imports from other modules
from app.agents.planner import PlannerAgent
from app.agents.generator import GeneratorAgent
from app.agents.ranker import RankerAgent
from app.chemistry.rdkit_tools import MoleculeTools
from app.screening.filters import ScreeningFilter

# --- In-Memory State Store (Simulating DB) ---
# Structure: { run_id: { "status": str, "results": list, "trace": list, "config": dict } }
RUN_STORE = {}


def update_trace(run_id: str, agent: str, action: str, inputs: Any, outputs: Any, round_num: int = 0):
    entry = {
        "timestamp": datetime.now().isoformat(),
        "round": round_num,
        "agent": agent,
        "action": action,
        "inputs_summary": str(inputs)[:200],  # Truncate for readability
        "outputs_summary": str(outputs)[:200]
    }
    if run_id in RUN_STORE:
        RUN_STORE[run_id]["trace"].append(entry)


def run_pipeline(run_id: str, initial_smiles: List[str], num_rounds: int, generations_per_seed: int, top_k: int,
                 max_violations: int):
    # Initialize State
    RUN_STORE[run_id] = {
        "status": "RUNNING",
        "results": [],
        "trace": [],
        "config": {"max_violations": max_violations}
    }

    try:
        # 1. PLANNER: Init
        current_seeds = PlannerAgent.initialize(initial_smiles)
        update_trace(run_id, "Planner", "Initialize Seeds", initial_smiles, current_seeds, 0)

        for r in range(1, num_rounds + 1):
            round_candidates = []

            # 2. GENERATOR: Mutate seeds
            raw_generated = GeneratorAgent.mutate_batch(current_seeds, n_variants=generations_per_seed)
            update_trace(run_id, "Generator", "Mutate Batch", f"Seeds: {len(current_seeds)}",
                         f"Generated: {len(raw_generated)}", r)

            # 3. CHEMISTRY & SCREENING
            processed_mols = []
            for smi in raw_generated:
                # Validation & Descriptors
                props = MoleculeTools.compute_descriptors(smi)
                if not props["valid"]:
                    continue

                # Scoring & Filtering
                score_data = ScreeningFilter.evaluate(props, max_violations)

                # Merge data
                mol_data = {**props, **score_data}
                processed_mols.append(mol_data)

            update_trace(run_id, "RDKit+Screening", "Compute & Filter", f"Raw: {len(raw_generated)}",
                         f"Valid: {len(processed_mols)}", r)

            # 4. RANKER: Select Top K
            top_mols = RankerAgent.rank_and_select(processed_mols, top_k)
            update_trace(run_id, "Ranker", "Select Top K", f"Candidates: {len(processed_mols)}",
                         f"Selected: {len(top_mols)}", r)

            # Update seeds for next round
            if top_mols:
                current_seeds = [m["smiles"] for m in top_mols]
                # Update global results with the best seen so far
                RUN_STORE[run_id]["results"] = top_mols
            else:
                update_trace(run_id, "System", "Early Stop", "No valid molecules found", "Stop", r)
                break

        RUN_STORE[run_id]["status"] = "COMPLETED"

    except Exception as e:
        RUN_STORE[run_id]["status"] = "FAILED"
        update_trace(run_id, "System", "Error", str(e), "Pipeline Crashed", -1)
        print(f"Pipeline Error: {e}")


# --- Data Accessors ---
def get_run_status(run_id: str):
    data = RUN_STORE.get(run_id)
    if data:
        return {"run_id": run_id, "status": data["status"], "timestamp": datetime.now().isoformat()}
    return None


def get_run_results(run_id: str):
    data = RUN_STORE.get(run_id)
    return data["results"] if data else None


def get_run_trace(run_id: str):
    data = RUN_STORE.get(run_id)
    return data["trace"] if data else None
