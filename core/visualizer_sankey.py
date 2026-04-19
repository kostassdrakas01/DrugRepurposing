import plotly.graph_objects as go
import os
from typing import Dict

class SankeyVisualizer:
    def __init__(self, output_dir: str):
        self.output_dir = output_dir

    def generate_sankey(self, analysis):
        """Generates an interactive Sankey diagram of the discovery flow."""
        data = analysis.sankey_data
        if not data or not data["links"]:
            return

        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=[node["name"] for node in data["nodes"]],
                color="royalblue"
            ),
            link=dict(
                source=[link["source"] for link in data["links"]],
                target=[link["target"] for link in data["links"]],
                value=[link["value"] for link in data["links"]],
                color="rgba(65, 105, 225, 0.4)" # Semi-transparent blue
            )
        )])

        fig.update_layout(title_text=f"Systemic Discovery Flow: {analysis.drug_name}", font_size=12)
        
        output_path = os.path.join(self.output_dir, f"{analysis.drug_name.lower()}_sankey_flow.html")
        fig.write_html(output_path)
        return output_path
