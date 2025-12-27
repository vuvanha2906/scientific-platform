import uuid
from typing import List
from fastapi import FastAPI, BackgroundTasks, HTTPException
from pydantic import BaseModel

# Import runner
from app.workers.runner import run_pipeline, get_run_status, get_run_results, get_run_trace

app = FastAPI(title="AI Scientific Workflow Platform")


# --- Pydantic Models ---
class RunRequest(BaseModel):
    initial_smiles: List[str]
    num_rounds: int = 3
    num_generations_per_seed: int = 5
    top_k: int = 3
    max_violations: int = 1


class RunResponse(BaseModel):
    run_id: str
    status: str
    message: str


# --- Endpoints ---

@app.post("/runs", response_model=RunResponse)
async def create_run(request: RunRequest, background_tasks: BackgroundTasks):
    """Initiate a new drug discovery process."""
    run_id = str(uuid.uuid4())

    # Trigger background worker
    background_tasks.add_task(
        run_pipeline,
        run_id=run_id,
        initial_smiles=request.initial_smiles,
        num_rounds=request.num_rounds,
        generations_per_seed=request.num_generations_per_seed,
        top_k=request.top_k,
        max_violations=request.max_violations
    )

    return {
        "run_id": run_id,
        "status": "QUEUED",
        "message": "Workflow started in background."
    }


@app.get("/runs/{run_id}/status")
async def get_status(run_id: str):
    """Check the current state of Run."""
    status = get_run_status(run_id)
    if not status:
        raise HTTPException(status_code=404, detail="Run ID not found")
    return status


@app.get("/runs/{run_id}/results")
async def get_results(run_id: str):
    """Get the best molecular results."""
    results = get_run_results(run_id)
    if results is None:
        raise HTTPException(status_code=404, detail="Run ID not found")
    return {"run_id": run_id, "top_molecules": results}


@app.get("/runs/{run_id}/trace")
async def get_trace(run_id: str):
    """Retrieve the entire activity history (Trace logs)."""
    trace = get_run_trace(run_id)
    if trace is None:
        raise HTTPException(status_code=404, detail="Run ID not found")
    return {"run_id": run_id, "trace_log": trace}


@app.get("/")
async def root():
    return {"message": "AI Scientific Platform is running!", "docs_url": "/docs"}
