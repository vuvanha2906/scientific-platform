# AI Scientific Workflow Platform

An agentic, asynchronous backend service designed to simulate an AI-driven drug discovery pipeline. This platform orchestrates a multi-step workflow involving molecule generation, property screening, and evolutionary ranking using **FastAPI** and **RDKit**.

## 🚀 Key Features

* **Agentic Architecture:** Decoupled logic for Planning, Generation, and Ranking agents.
* **Asynchronous Execution:** Non-blocking API using FastAPI BackgroundTasks for long-running scientific computations.
* **Chemical Intelligence:** Integrated **RDKit** for creating valid SMILES, calculating descriptors (MW, LogP, TPSA), and evaluating drug-likeness (QED).
* **Traceability:** Full structured logging of every agent action, input, and output for scientific reproducibility.
* **Deterministic Logic:** Rule-based mutation and scoring for predictable, testable results (no black-box ML models).

## 📂 Project Structure

```text
scientific-platform/
├── app/
│   ├── main.py                 # API Entry point & Routes
│   ├── agents/                 # Intelligent Agents
│   │   ├── planner.py          # Setup & Config
│   │   ├── generator.py        # Molecule Mutation Logic
│   │   └── ranker.py           # Selection Logic
│   ├── chemistry/              # Domain Layer
│   │   └── rdkit_tools.py      # RDKit wrappers
│   ├── screening/              # Evaluation Layer
│   │   └── filters.py          # Scoring & Lipinski Rules
│   └── workers/                # Orchestration
│       └── runner.py           # Async Loop Manager
└── README.md
```

## 🛠️ Installation & Setup
**Prerequisites:**
* Python 3.10+ or higher


**Steps:**
1. Clone the repository:
```bash
git clone [https://github.com/vuvanha2906/scientific-platform.git](https://github.com/vuvanha2906/scientific-platform.git)
cd scientific-platform
```
2. Create a Virtual Environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate 
```
3. Install Dependencies:
```bash
pip install -r requirements.txt
```
4. Run the sever:
```bash
uvicorn app.main:app --reload
```
The API will be available at http://127.0.0.1:8000/docs.

## 🧪 Scientific Logic

**1. The Agentic Loop** 

The system runs an evolutionary loop for a specified number of rounds.

    a. Input: Seed molecules (SMILES).
    
    b. Generator: Creates N variants per seed using deterministic atom addition (e.g., Methylation).
    
    c. Screener: Calculates physicochemical properties and filters invalid molecules.
    
    d. Ranker: Selects the Top-K molecules based on the Scoring Function to serve as seeds for the next round.

**2. Scoring Function**

Molecules are evaluated based on QED (Quantitative Estimate of Drug-likeness) penalized by Lipinski Rule Violations.

$$ Score = QED - (0.1 \times Violations) $$

**Filters:**

* Molecular Weight (MW) ≤ 500

* LogP ≤ 5

* H-Bond Donors (HBD) ≤ 5

* H-Bond Acceptors (HBA) ≤ 10

* TPSA ≤ 140

* max_violations = 1

## 📡 API Usage (Example)

**1. Start a New Run (POST)**

Example (Using Aspirin, Ethoxybenzene,Acetanilide, and Ethylamine as a seed. You can try more molecules as SMILES from one or more seed SMILES) 
```
curl -X POST "[http://127.0.0.1:8000/runs](http://127.0.0.1:8000/runs)" \
     -H "Content-Type: application/json" \
     -d '{
           "initial_smiles": ["CCOc1ccc(CCNC(=O)C)cc1",
                             "c1ccccc1C(=O)NC",                
                             "CCN(CC)CCOC1=CC=CC=C1",
                             "CC(=O)OC1=CC=CC=C1C(=O)O"],     
           "num_rounds": 3,
           "num_generations_per_seed": 5,
           "top_k": 2,
           "max_violations": 1
         }'
```

Response: 

```JSON
{
  "run_id": "{run_id}",
  "status": "QUEUED",
  "message": "Workflow started in background."
}
```

**2. Check Status (GET)**

Check if the background worker is still processing or finished.
```
curl "[http://127.0.0.1:8000/runs/](http://127.0.0.1:8000/runs/){run_id}/status"
```

**3. Get Results (GET)**

Retrieve the top molecules from the final generation.
```
curl "[http://127.0.0.1:8000/runs/](http://127.0.0.1:8000/runs/){run_id}/results"
```
**4. Get Execution Trace (GET)**

Retrieve the full audit log of agent actions for debugging or analysis.

```
curl "[http://127.0.0.1:8000/runs/](http://127.0.0.1:8000/runs/){run_id}/trace"
```

## 📄 License
**MIT**




