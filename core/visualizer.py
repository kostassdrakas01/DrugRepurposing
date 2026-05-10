from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from core.models import DrugAnalysis
import datetime
import os

class Visualizer:
    def __init__(self):
        self.console = Console()

    def display_analysis(self, analysis: DrugAnalysis):
        self.console.print(Panel.fit(
            f"[bold cyan]SYSTEMIC PREDICTIVE ANALYSIS:[/bold cyan] [bold white]{analysis.drug_name}[/bold white]\n"
            f"[bold cyan]Molecular Profile:[/bold cyan] {len(analysis.targets)} Targets | {len(analysis.pathways)} Novel Hits",
            subtitle=f"Generated {datetime.date.today()}",
            border_style="cyan"
        ))

        # 1. Discovery Table
        h_table = Table(title="Bio-Network Discovery Highlights", show_header=True, header_style="bold magenta")
        h_table.add_column("Pathway Discovery", style="white")
        h_table.add_column("System", style="dim")
        h_table.add_column("Predicted Polarity", justify="center")
        h_table.add_column("Discovery Score", justify="center")
        h_table.add_column("Z-Score", justify="center")

        discoveries = sorted(analysis.pathways.values(), key=lambda x: x.discovery_score, reverse=True)
        for p in discoveries[:10]:
            pol_color = "green" if "Activated" in p.polarity else "red" if "Inhibited" in p.polarity else "white"
            h_table.add_row(p.name, p.category, f"[{pol_color}]{p.polarity}[/{pol_color}]", f"{p.discovery_score:.2f}", f"{p.z_score:.2f}")
        
        if discoveries: self.console.print(h_table)

    def export_report(self, analysis: DrugAnalysis, filename: str, include_details: bool = False):
        """Saves enhanced Markdown discovery report."""
        # Ensure the filename ends with .md
        if not filename.lower().endswith(".md"):
            filename += ".md"
        
        self._export_markdown(analysis, filename)
        self.console.print(f"[bold green]Systemic Discovery Report exported to {filename}[/bold green]")

    def _export_markdown(self, analysis: DrugAnalysis, filename: str):
        """Generates a high-fidelity Markdown report with simplified human-readable terminology."""
        hub_image = f"{analysis.drug_name.lower()}_hub_and_spoke.jpg"
        sunburst_image = f"{analysis.drug_name.lower()}_discovery_sunburst.jpg"
        heatmap_image = f"{analysis.drug_name.lower()}_discovery_heatmap.jpg"
        
        with open(filename, "w") as f:
            f.write(f"# Systemic Discovery & Predictive Report: {analysis.drug_name}\n\n")
            
            f.write("## EXECUTIVE SUMMARY\n")
            f.write(f"**Target Analyzed:** {analysis.drug_name} (CID: {analysis.cid})\n")
            f.write(f"**Discovery Scope:** Identified {len(analysis.pathways)} novel disease links.\n")
            if analysis.target_tissue:
                f.write(f"**Tissue Context:** {analysis.target_tissue} (GTEx Verified)\n")
            f.write("\n")
            
            f.write("### 🗝️ Hub-and-Spoke Quick-Reference Map\n")
            f.write("The following table maps the numeric identifiers (1-15) displayed on the blue outcome nodes in the Meta-Network visual below to their assigned biological pathways.\n\n")
            f.write("| Node # | Pathway Discovery | Discovery Score |\n")
            f.write("|---|---|---|\n")
            if hasattr(analysis, 'network_legend'):
                for pid, (pname, score) in analysis.network_legend.items():
                    clean_name = pname.split(" - ")[0]
                    f.write(f"| {pid} | {clean_name} | {score:.2f} |\n")
            f.write("\n")

            f.write("### Visual Discovery Portfolio\n")
            f.write(f"![Discovery Meta-Network]({hub_image})\n\n")
            f.write(f"![Hierarchial Discovery Sunburst]({sunburst_image})\n\n")
            f.write(f"![Genetic Crosstalk Heatmap]({heatmap_image})\n\n")
            
            f.write("## I. NEW POTENTIAL DISEASE TARGETS\n")
            f.write("| Discovery Pathway | System Category | Predicted Effect | Discovery Score | Z-Score (Specificity) | Biological Mechanism Narrative |\n")
            f.write("|---|---|---|---|---|---|\n")
            
            primary_discoveries = sorted(analysis.pathways.values(), key=lambda x: x.discovery_score, reverse=True)
            for p in primary_discoveries:
                clean_name = p.name.split(" - ")[0]
                insight = p.discovery_insight or "Downstream systemic modulation detected."
                f.write(f"| [{clean_name}]({p.url}) | {p.category} | **{p.polarity}** | {p.discovery_score:.2f} | {p.z_score:.2f} | {insight} |\n")
            f.write("\n")

            f.write("## II. THE MOLECULAR CONNECTORS\n")
            f.write("| Connector Protein (Bridge) | Pathway Count | Discovery Context |\n")
            f.write("|---|---|---|\n")
            # Cluster summaries
            for group in sorted(analysis.convergence_groups, key=lambda x: x.pathway_count, reverse=True)[:10]:
                symbols = [n for n in group.nodes if n and n[0].isupper()]
                if not symbols:
                    if not group.nodes: continue
                    symbol = group.nodes[0]
                else:
                    symbol = symbols[0]
                
                cluster_name = f"**{symbol}**"
                if len(group.nodes) > 1: cluster_name += f" (+ {len(group.nodes)-1} others)"
                
                context = ", ".join([p.split(" - ")[0] for p in group.related_pathways[:3]])
                f.write(f"| {cluster_name} | {group.pathway_count} | {context}... |\n")
            f.write("\n")

            f.write("## III. DOWNSTREAM IMPACT ON CELLS\n")
            f.write("| Distal Pathway | System Branch | Discovery Score |\n")
            f.write("|---|---|---|\n")
            # Focus on higher connectivity distal nodes
            for pid, p in list(analysis.appendix_pathways.items())[:15]:
                clean_name = p.name.split(" - ")[0]
                f.write(f"| {clean_name} | {p.category} | {p.discovery_score:.2f} |\n")
            f.write("\n")

            f.write("--- \n")
            f.write("## IV. KNOWN & EXPECTED EFFECTS (APPENDIX)\n")
            f.write("| Known Mechanism | Logic | Evidence |\n")
            f.write("|---|---|---|\n")
            # Logic: Pathways with low discovery score or exact matches
            # Also include pathways that were explicitly marked as indications in surprise score
            for pid, p in {**analysis.pathways, **analysis.appendix_pathways}.items():
                if p.surprise_score <= 0.15:
                    clean_name = p.name.split(" - ")[0]
                    f.write(f"| {clean_name} | Primary Indication | High Confidence |\n")
            f.write("\n")

