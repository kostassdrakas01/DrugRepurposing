# KEGG-ID: User Guide & Examples

This guide provides detailed workflows for using the KEGG-ID Systemic Discovery Engine.

## 🏁 Quick Start

Always ensure your virtual environment is activated and dependencies are installed:
```bash
source venv/bin/activate
```

---

## 🔍 1. Traditional Discovery (Drug-to-Pathway)
Used when you have a specific molecule and want to predict its systemic effects across all biological pathways.

**Command:**
```bash
python main.py "Aspirin"
```

**With Tissue Context (Recommended):**
```bash
python main.py "Sirolimus" --tissue "Brain_Cortex"
```
*Tip: Using tissue context (from GTEx) filters out pathways not active in that specific organ, significantly increasing the signal-to-noise ratio.*

---

## 🔄 2. Reverse Discovery (Target-to-Drug)
Used when you have a specific protein target and want to analyze all approved drugs that interact with it.

**Command:**
```bash
python main.py --target "EGFR" --tissue "Lung"
```
*This will find drugs like Osimertinib and Gefitinib, then run the full V2 engine for each.*

---

## 🏥 3. Reverse Discovery (Disease-to-Drug)
Used when you have a clinical indication and want to find and analyze existing therapies.

**Command:**
```bash
python main.py --disease "Lung Cancer" --tissue "Lung"
```
*This resolves the disease name via OpenTargets, identifies associated drugs, and performs systemic discovery.*

---

## 🧬 4. Repurposing Discovery Mode (Network-Based)
The most advanced mode. It scans the entire molecular neighborhood of a disease to find non-obvious drug candidates (drugs not currently approved for the disease but targeting linked proteins).

**Command:**
```bash
python main.py --disease "Alzheimer" --repurpose --tissue "Brain_Cortex"
```

**What it does:**
1. **Driver Expansion**: Identifies the top 20 genetic drivers of the disease.
2. **Network Expansion**: Uses String-DB to find high-confidence interactors (neighbors) of these drivers.
3. **Pathway Filtering**: Cross-references neighbors with KEGG pathways to ensure functional relevance.
4. **Novelty Filtering**: Automatically excludes standard-of-care drugs to find true repurposing leads.
5. **Systemic Analysis**: Runs the full predictive engine for the top 10 discovered candidates.

---

## 📊 5. Post-Analysis: Analytical Visualization (R)
After the Python engine generates the `.csv` result file, use the R script to create production-grade analytical plots.

**Command:**
```bash
Rscript plot_results.R results/Osimertinib.csv
```

**Generated Visuals:**
1. **Target-Pathway Heatmap**: Shows which targets are driving which pathway hits.
2. **Discovery Quadrant Analysis**: Plots **Surprise Score** vs. **Z-Score** to identify high-confidence, non-obvious discoveries.

---

## ⌨️ Full Command-Line Reference

The `main.py` script supports several arguments to customize the discovery process.

| Flag | Parameter | Purpose |
| :--- | :--- | :--- |
| `[drug]` | `str` | **(Positional)** The primary molecule name. If provided, the engine runs a Forward Discovery (Drug -> Pathway). |
| `--disease` | `str` | Triggers **Reverse Discovery**. Resolves the disease to an EFO ID and finds associated drugs. |
| `--target` | `str` | Triggers **Reverse Discovery**. Finds all drugs targeting this gene symbol (e.g., `MTOR`, `EGFR`). |
| `--tissue` | `list` | Filters pathways by GTEx tissue expression. You can provide multiple tissues: `--tissue Lung Liver`. |
| `--repurpose` | *None* | **Boolean Flag**. Must be used with `--disease`. Switches to the network-expansion mode (String-DB + Drivers). |
| `--export` | `str` | Overrides the default report filename. Default: `results/[DrugName].md`. |

### Advanced Combinations
- **Strict Tissue Filtering**: `python main.py "Aspirin" --tissue "Blood" "Liver"` (Ensures targets are active in both/either tissue).
- **Deep Repurposing**: `python main.py --disease "Type 2 Diabetes" --repurpose --tissue "Pancreas"` (Searches the diabetes interactome for novel metabolic modulators).

---

## 📁 Understanding the Results
All results are stored in the `results/` directory:

- **`[DrugName].md`**: The primary report. Read this to understand the "Mechanism Narrative."
- **`[DrugName]_sankey.jpg`**: Visualizes the flow from drug to systemic outcome.
- **`[DrugName]_network.jpg`**: Shows the protein-protein interaction network.
- **`[DrugName].csv`**: Raw data used for R-plotting and further research.

---

## ⚠️ Common Issues
- **ModuleNotFoundError**: Ensure you are using the virtual environment (`source venv/bin/activate`).
- **Read Timeout**: Large disease categories (like "Cancer") can time out. Try more specific terms like "Lung Adenocarcinoma."
- **No Candidates Found**: Check that the Target Symbol (e.g., "EGFR") or Disease Name is correctly spelled.
- **Venv issues**: If `source venv/bin/activate` fails, ensure you created it with `python -m venv venv`.

