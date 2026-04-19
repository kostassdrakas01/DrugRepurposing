import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import pandas as pd
import os
import plotly.graph_objects as go
from core.models import DrugAnalysis

class NetworkVisualizer:
    def __init__(self, output_dir="results"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def generate_all_visuals(self, analysis: DrugAnalysis):
        """Generates Hub-and-Spoke, Sunburst, and Heatmap."""
        self.plot_hub_and_spoke(analysis)
        self.plot_discovery_sunburst(analysis)
        self.plot_discovery_heatmap(analysis)

    def plot_hub_and_spoke(self, analysis: DrugAnalysis):
        """
        Refined Hub-and-Spoke (Directive: Visual Hierarchy Fix).
        Center: Drug | Inner: Targets (Red) | Middle: Bridges (Green) | Outer: Pathways (Blue).
        """
        G = nx.DiGraph()
        drug_node = analysis.drug_name.upper()
        G.add_node(drug_node, type='drug')
        
        # Targets
        primary_targets = [t.name for t in analysis.targets if t.is_primary]
        for t in primary_targets:
            G.add_node(t, type='target')
            G.add_edge(drug_node, t)
            
        # Bridges
        convergence_nodes = []
        bridge_labels = {}
        # Get top 5 bridges for clarity
        active_groups = sorted(analysis.convergence_groups, key=lambda x: x.pathway_count, reverse=True)[:5]
        for group in active_groups:
            if not group.nodes: continue
            # Find the actual gene symbol (robust check)
            symbols = [n for n in group.nodes if n and n[0].isupper()]
            if not symbols:
                if not group.nodes: continue
                symbol = group.nodes[0]
            else:
                symbol = symbols[0]
            
            label = symbol
            if len(group.nodes) > 1:
                label = f"{symbol} + {len(group.nodes)-1} others"
            
            node_id = f"Bridge: {symbol}"
            G.add_node(node_id, type='bridge')
            convergence_nodes.append(node_id)
            bridge_labels[node_id] = label
            
            # Connect targets to this bridge
            for pt in primary_targets:
                if pt in group.nodes or any(pt == s for s in group.nodes):
                    G.add_edge(pt, node_id)
                else:
                    # Weak fallback: if target shares any pathway with the bridge group
                    t_obj = [to for to in analysis.targets if to.name == pt][0]
                    if any(p in t_obj.pathways for p in group.related_pathways):
                        G.add_edge(pt, node_id)
        
        # Pathways
        pathway_ids = []
        pathway_legend_map = {}
        # Top 15 pathways to avoid hairball
        pathway_list = sorted(analysis.pathways.values(), key=lambda x: x.discovery_score, reverse=True)[:15]
        
        for i, p in enumerate(pathway_list, 1):
            pid = str(i)
            G.add_node(pid, type='pathway')
            pathway_ids.append(pid)
            pathway_legend_map[pid] = (p.name, p.discovery_score)
            
            # Connect bridges to pathways
            connected = False
            for bid, group in zip(convergence_nodes, active_groups):
                if p.name in group.related_pathways:
                    G.add_edge(bid, pid)
                    connected = True
            
            if not connected:
                # Direct target link if no bridge
                for pt in primary_targets:
                    t_obj = [to for to in analysis.targets if to.name == pt][0]
                    if p.id in t_obj.pathways:
                        G.add_edge(pt, pid)

        if not G.edges(): return

        plt.figure(figsize=(22, 22))
        shells = [[drug_node], primary_targets, convergence_nodes, pathway_ids]
        shells = [s for s in shells if s]
        pos = nx.shell_layout(G, shells)
        
        # Draw 1: Center Drug
        nx.draw_networkx_nodes(G, pos, nodelist=[drug_node], node_color='#FFD700', node_size=8000, node_shape='h')
        # Draw 2: Targets (Red)
        nx.draw_networkx_nodes(G, pos, nodelist=primary_targets, node_color='#FF0000', node_size=3500)
        # Draw 3: Bridges (Green Diamonds, 50% larger than all other nodes per mandate) -> 1.5 * 8000 = 12000
        if convergence_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=convergence_nodes, node_color='#00FF00', 
                                   node_size=12000, node_shape='d', edgecolors='#1B4332', linewidths=2)
        # Draw 4: Pathways (Blue Circles)
        nx.draw_networkx_nodes(G, pos, nodelist=pathway_ids, node_color='#2196F3', node_size=1600)
        
        # Labels
        labels = {n: n for n in G.nodes}
        for n in convergence_nodes: labels[n] = bridge_labels[n]
        
        # Draw Labels
        nx.draw_networkx_labels(G, pos, labels={drug_node: labels[drug_node]}, font_size=12, font_weight='bold')
        # Directive: Smaller names for target nodes (Red circles)
        nx.draw_networkx_labels(G, pos, labels={n: labels[n] for n in primary_targets}, font_size=7, font_weight='bold')
        nx.draw_networkx_labels(G, pos, labels={n: labels[n] for n in convergence_nodes}, font_size=10, font_weight='bold')
        nx.draw_networkx_labels(G, pos, labels={n: n for n in pathway_ids}, font_size=10, font_weight='bold', font_color='white')
        
        nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='#455A64', width=1.5, arrowsize=25, connectionstyle="arc3,rad=0.1")
        
        plt.title(f"Discovery Meta-Network: {drug_node}\n(Radial Discovery Logic | Functional Bridges Highlighted)", fontsize=26, fontweight='bold', pad=60)
        plt.axis('off')
        
        p_name = f"{analysis.drug_name.lower()}_hub_and_spoke.jpg"
        plt.savefig(os.path.join(self.output_dir, p_name), dpi=300, bbox_inches='tight')
        plt.close()

        # Save legend metadata for the MD report (Restructured: specifically for the 15-node plot)
        analysis.network_legend = pathway_legend_map

    def plot_discovery_sunburst(self, analysis: DrugAnalysis):
        """
        Creates a Hierarchical Sunburst Chart.
        Center: Drug | Inner: System Category | Outer: Pathways.
        Slices are weighted by Discovery Score.
        """
        import plotly.express as px
        
        # Prepare Data
        data = []
        # Center
        drug_name = analysis.drug_name.upper()
        
        # Get path entries: [System, Pathway, Score]
        # Include all discoveries > 0.81
        all_pathways = sorted(analysis.pathways.values(), key=lambda x: x.discovery_score, reverse=True)
        
        for p in all_pathways:
            data.append({
                "Drug": drug_name,
                "System": p.category,
                "Pathway": p.name.split(" - ")[0], # Human readable
                "Discovery Score": p.discovery_score
            })
        
        if not data: return
        
        df = pd.DataFrame(data)
        
        fig = px.sunburst(
            df,
            path=['Drug', 'System', 'Pathway'],
            values='Discovery Score',
            color='Discovery Score',
            color_continuous_scale='RdBu_r',
            title=f"Systemic Discovery Hierarchy: {drug_name}"
        )
        
        fig.update_layout(
            margin=dict(t=50, l=10, r=10, b=10),
            font=dict(size=14, family="Arial Black")
        )
        
        # Export as static image
        s_name = f"{analysis.drug_name.lower()}_discovery_sunburst.jpg"
        fig.write_image(os.path.join(self.output_dir, s_name), scale=2)
        
        # (Removed legacy legend overwrite from sunburst to preserve network map accuracy)

    def plot_discovery_heatmap(self, analysis: DrugAnalysis):
        """Heatmap limited to Top 25 pathways for maximum readability."""
        all_p = list(analysis.pathways.values()) + list(analysis.appendix_pathways.values())
        if len(all_p) < 2: return
        
        # Top 25 strict limit
        pathways = sorted(all_p, key=lambda x: x.discovery_score, reverse=True)[:25]
        
        names = [p.name.split(" - ")[0][:35] for p in pathways]
        data = [[len(set(p1.genes).intersection(set(p2.genes))) for p2 in pathways] for p1 in pathways]
        df = pd.DataFrame(data, index=names, columns=names)
        
        plt.figure(figsize=(18, 14))
        sns.heatmap(df, annot=True, fmt="d", cmap="YlGnBu", cbar=True, square=True, annot_kws={"size": 10})
        plt.title(f"Discovery Genetic Crosstalk: {analysis.drug_name}\n(Crosstalk Mapping Map)", fontsize=20, fontweight='bold', pad=25)
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)
        plt.tight_layout()
        
        h_name = f"{analysis.drug_name.lower()}_discovery_heatmap.jpg"
        plt.savefig(os.path.join(self.output_dir, h_name), dpi=300)
        plt.close()
