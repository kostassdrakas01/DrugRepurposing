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
python main.py "Osimertinib" --tissue "Lung"
```

### Advanced Options
- **GTEx Tissue Masking**: Filter out pathways where targets are not expressed in your target tissue(s). Genes are considered expressed if TPM > 10.0 in the selected context.
  ```bash
  # Single tissue
  python main.py "Amitriptyline" --tissue "Brain_Cortex"
  
  # Multiple tissues (Gene must be expressed in AT LEAST ONE to be kept)
  python main.py "Ibuprofen" --tissue "Liver" "Muscle_Skeletal"
  ```
- **High-Affinity Expansion**: For drugs with a single primary target, the engine automatically integrates top 5 high-affinity interactors (via String-DB) to ensure discovery depth.

## Predictive Outputs

All results are automatically organized in the `/results` folder:

1. **`<drug_name>.md`**: A professional report focusing on:
   - **Functional Polarity**: Predicts if the drug **activates** (via disinhibition) or **inhibits** the pathway.
   - **Discovery Score**: A composite metric combining **Surprise Score** (biological novelty) and **Z-Score** (statistical significance).
   - **Mechanism Narratives**: Synthesized biological summaries explaining exactly how the drug influences each discovery pathway through specific "bridge" proteins.
   - **Appendix Logic**: Pathways with low statistical significance are moved to the appendix to reduce noise.
2. **`Systemic Ripple Flow (.html)`**: An **interactive Sankey diagram** mapping the flow from Drug -> Target Protein -> Functional Bridge -> Outcome.
3. **`Functional Bridges`**: Grouped clusters of proteins that act as a single signaling unit across multiple pathways.

## Terminology Guide

- **Discovery Score**: Normalized 0.0 to 1.0. Scores >= 0.80 are considered high-priority discoveries.
- **Statistical Z-Score**: Computed via 1,000 Monte Carlo permutations to neutralize "hub bias" (proteins that appear in many pathways by random chance).
- **Polarity**: "Activated (Disinhibition)" means the drug inhibits an inhibitor, leading to net activation of the downstream pathway.
- **Confidence Badges**: 
    - **High Confidence**: Supported by both statistical significance and tissue expression data.
    - **Tissue-Supported**: Primary bridge proteins are verified as expressed in the selected GTEx context.
