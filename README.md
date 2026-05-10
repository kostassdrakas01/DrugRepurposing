# DrugRepurposing: Bio-Network Discovery Engine

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

DrugRepurposing is a high-fidelity predictive engine designed to uncover novel therapeutic indications for existing drugs. By mapping drug-target interactions across complex biological networks, it identifies hidden crosstalk between pathways, predicts functional polarity, and ranks discoveries by biological "novelty."

---

## Key Features

*   **Multi-Omic Mapping**: Integrates data from KEGG, ChEMBL, OpenTargets, PubChem, and MyGene.
*   **Tissue-Specific Context (GTEx)**: Apply GTEx expression masks to filter results based on real-world tissue biology (e.g., analyzing drug effects specifically in the Brain Cortex or Liver).
*   **Functional Polarity Prediction**: Predicts whether a drug likely activates (via disinhibition) or inhibits a specific downstream pathway.
*   **Surprise Scoring**: A heuristic ranking system that prioritizes "cross-category" biological links (e.g., a Cardiovascular drug affecting Neurodegeneration).
*   **Interactive Visualizations**: Generates systemic ripple flows (Sankey diagrams) and hub-and-spoke network graphs automatically.

---

## Installation

1.  **Clone the Repository**:
    ```bash
    git clone https://github.com/kostassdrakas01/DrugRepurposing.git
    cd DrugRepurposing
    ```

2.  **Setup Environment**:
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate
    pip install -r requirements.txt
    ```

---

## Usage

Analyze any drug to discover its extended bio-network:

```bash
# Basic analysis
python main.py "Aspirin"

# Tissue-specific discovery (e.g., Sirolimus in the Brain)
python main.py "Sirolimus" --tissue "Brain_Cortex"

# Export a custom report
python main.py "Metformin" --export metformin_analysis.md
```

### Supported Tissue Contexts
For a full list of supported GTEx tissue IDs (e.g., Adipose_Subcutaneous, Heart_Left_Ventricle), see [GTEX_TISSUE_LIST.md](GTEX_TISSUE_LIST.md).

---

## Outputs and Insights

All discoveries are saved in the results/ directory:

1.  **Discovery Reports (.md)**: Comprehensive summaries including novelty scores and functional logic.
2.  **Systemic Ripple Flows (.html)**: Interactive Sankey diagrams showing the flow from Drug -> Target -> Functional Bridge -> Outcome.
3.  **Visual Diagnostics (.jpg)**: Heatmaps and Network graphs illustrating protein hubs and pathway convergence.

---

## Examples

This repository includes pre-generated analysis results for:
*   **Sirolimus**: Focusing on its role in Neurodegeneration (Brain Cortex).
*   **Semagacestat**: Exploring its failure in Alzheimer's vs. potential repurposing.
*   **Thalidomide**: Mapping its multi-target immunomodulatory effects.
*   **And more...** (Check the results/ folder)

---

## Security and Privacy
This tool uses public REST APIs (KEGG, EBI ChEMBL, OpenTargets). No private keys are required for basic usage. Ensure you comply with the terms of service for the respective data providers.

---

## License
This project is licensed under the MIT License - see the LICENSE file for details.

---
*Developed for advanced pharmacological research and bio-network exploration.*
