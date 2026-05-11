# DrugRepurposing: Bio-Network Discovery Engine

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

DrugRepurposing is a high-fidelity predictive engine designed to uncover novel therapeutic indications for existing drugs. By mapping drug-target interactions across complex biological networks, it identifies hidden crosstalk between pathways, predicts functional polarity, and ranks discoveries using a statistically robust multi-omic scoring system.

---

## Key Features

KEGG-ID is a high-fidelity drug discovery and repurposing engine that maps molecular entities to systemic biological pathways. It integrates **OpenTargets**, **String-DB**, **KEGG**, **MyGene**, **PubChem**, and **GTEx** to provide a multi-dimensional view of drug-target-pathway interactions.

## 🚀 Discovery Engine V2 Features

- **Multi-Directional Discovery**: Support for traditional Drug-to-Pathway analysis and new **Reverse Discovery** (Target-to-Drug and Disease-to-Drug).
- **Statistical Rigor**: Implemented Monte Carlo permutation testing (1,000 iterations) to calculate **Z-Scores**, neutralizing "hub bias" in pathway hits.
- **Interactome Expansion**: Automatic retrieval of high-affinity interactors via String-DB to deepen the biological context of single-target drugs.
- **Tissue Contextualization**: GTEx integration allows filtering pathways by tissue-specific expression (TPM > 10.0).
- **Mechanism Synthesis**: Automated generation of "Mechanism Narratives" explaining systemic flow from targets to pathways.
- **Analytical Plotting**: Integrated R-based visualization for discovery quadrant analysis (Surprise Score vs. Z-Score).

---

## 🛠️ Installation

### Prerequisites
- Python 3.9+ (Fully compatible with Python 3.13)
- R (for analytical plotting)
- Virtual Environment (recommended)

```bash
git clone https://github.com/kostassdrakas01/DrugRepurposing.git
cd DrugRepurposing
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

---

## 📖 Usage

### 1. Traditional Discovery (Drug → Pathways)
Analyze a specific drug and its systemic impact.
```bash
python main.py "Osimertinib" --tissue "Lung"
```

### 2. Reverse Discovery (Target → Drugs)
Find all approved drugs for a target and analyze their networks.
```bash
python main.py --target "EGFR" --tissue "Lung"
```

### 3. Reverse Discovery (Disease → Drugs)
Batch analyze drugs associated with a specific clinical indication.
```bash
python main.py --disease "Alzheimer" --tissue "Brain_Cortex"
```

---

## ⌨️ Command-Line Interface (CLI) Reference

| Argument | Description | Example |
| :--- | :--- | :--- |
| `drug` | **(Positional)** Name of the drug to analyze. | `python main.py "Aspirin"` |
| `--disease` | Name of a disease to find drugs for (Reverse Discovery). | `--disease "Alzheimer"` |
| `--target` | Name of a protein target to find drugs for (Gene Symbol). | `--target "EGFR"` |
| `--tissue` | GTEx Tissue IDs to filter for expressed genes (TPM > 10). | `--tissue "Lung" "Brain_Cortex"` |
| `--repurpose` | **(Flag)** Enables Network-Based Discovery for non-obvious leads. | `--repurpose` |
| `--export` | Custom filename for the generated Markdown report. | `--export "my_report.md"` |

---

### 4. Analytical Visualization
After running the Python discovery, generate the analytical plots using R:
```bash
Rscript plot_results.R results/Osimertinib.csv
```

### Supported Tissue Contexts
For a full list of supported GTEx tissue IDs, see [GTEX_TISSUE_LIST.md](GTEX_TISSUE_LIST.md).

---

## 📂 Output

All discoveries are saved in the results/ directory:

1.  **Discovery Reports (.md)**: Comprehensive summaries including Discovery Scores (Surprise + Z-Score), Polarity, and Mechanism Narratives.
2.  **Raw Discovery Data (.csv)**: Structured data for statistical analysis.
3.  **Systemic Visuals (.jpg)**: Sankey diagrams, Hub analysis, and Heatmaps.
4.  **Analytical Plots (.png)**: Quadrant analysis and statistical validation.

---

## Examples

This repository includes pre-generated analysis results demonstrating the engine's capabilities:
*   **Sirolimus**: Focusing on its role in Neurodegeneration (Brain Cortex).
*   **Semagacestat**: Exploring its failure in Alzheimer's vs. potential repurposing.
*   **Thalidomide**: Mapping its multi-target immunomodulatory effects.

---

## Security and Privacy
This tool uses public REST APIs (KEGG, EBI ChEMBL, OpenTargets). No private keys are required for basic usage. Ensure you comply with the terms of service for the respective data providers.

---

## License
This project is licensed under the MIT License - see the LICENSE file for details.

---
*Developed for advanced pharmacological research and bio-network exploration.*
