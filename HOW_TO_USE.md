# Bio-Network Discovery Tool: Instructions

This professional tool maps drug names to their molecular targets and biological pathways to discover hidden crosstalk, bottlenecks, and clinical insights (such as viral-cancer intersections).

## Installation

1. **Environment**: Ensure you have Python 3.9+ installed.
2. **Setup**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r requirements.txt
   ```

## Usage

Run a predictive discovery analysis for any drug:

```bash
python main.py "Osimertinib" --export report.pdf --tissue "Lung"
```
python main.py "Ibuprofen" --export ibuprofen_report.pdf

### Advanced Options
- **GTEx Tissue Masking**: Filter out pathways where targets are not expressed in your target tissue(s).
  ```bash
  # Single tissue
  python main.py "Amitriptyline" --tissue "Brain_Cortex"
  
  # Multiple tissues (Gene must be expressed in AT LEAST ONE to be kept)
  python main.py "Ibuprofen" --tissue "Liver" "Muscle_Skeletal"
  ```
- **Hide/Show Details**: Use `--detailed` to include full pathway deep dives in the PDF.

### Valid Tissue Types
For a comprehensive list of all 54+ valid tissue IDs, please refer to the [GTEX_TISSUE_LIST.md](GTEX_TISSUE_LIST.md) file.

## Predictive Outputs

All results are automatically organized in the `/results` folder:

1. **`<drug_name>.md`**: A professional report focusing on:
   - **Functional Polarity**: Predicts if the drug **activates** (disinhibits) or **inhibits** the pathway.
   - **Surprise Score**: Ranks connections by biological "novelty"—high scores for cross-category links (e.g., Metabolism → Neurodegeneration).
   - **Appendix Logic**: Common "In-Category" pathways are moved to the appendix to reduce noise in the results.
2. **`Systemic Ripple Flow (.html)`**: An **interactive Sankey diagram** mapping the flow from Drug → Target Protein → Functional Bridge → Outcome.
3. **`Functional Bridges`**: Grouped clusters of proteins that act as a single signaling unit across multiple pathways.

## Terminology Guide
- **Novelty Score**: 0.9+ suggests a rare cross-category link that warrants research.
- **Polarity**: "Activated (Disinhibition)" means the drug inhibits an inhibitor, leading to net activation of the downstream pathway.
- **Human-Readable IDs**: All database jargon (e.g., `hsa:1956`) is automatically translated to names (e.g., **EGFR**).
