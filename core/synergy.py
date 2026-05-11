import networkx as nx
from typing import List, Dict
from core.models import DrugAnalysis

class SynergyCalculator:
    @staticmethod
    def calculate_synergy(analysis_a: DrugAnalysis, analysis_b: DrugAnalysis, global_graph: nx.DiGraph) -> Dict:
        """
        Calculates topology-based synergy scores between two drugs.
        Uses a simplified model:
        1. Pathway Overlap: Are they targeting the same pathway?
        2. Parallel Pathways: Are they targeting different pathways that converge?
        3. Network Proximity: How close are their targets in the global interactome?
        """
        pathways_a = set(analysis_a.pathways.keys())
        pathways_b = set(analysis_b.pathways.keys())

        shared_pathways = pathways_a.intersection(pathways_b)

        # Bliss-like score for shared pathways
        pathway_synergy = []
        for pid in shared_pathways:
            p_a = analysis_a.pathways[pid]
            p_b = analysis_b.pathways[pid]
            # If both have the same polarity, might be additive/synergistic
            # If opposite, might be antagonistic
            score = (p_a.discovery_score + p_b.discovery_score) / 2.0
            if p_a.polarity == p_b.polarity and p_a.polarity != "Neutral":
                score *= 1.2 # Synergy bonus
            pathway_synergy.append({"id": pid, "name": p_a.name, "score": min(0.99, score)})

        # Network Proximity
        targets_a = [t.kegg_id for t in analysis_a.targets if t.kegg_id]
        targets_b = [t.kegg_id for t in analysis_b.targets if t.kegg_id]

        avg_dist = 0
        pairs = 0
        for ta in targets_a:
            for tb in targets_b:
                try:
                    if ta in global_graph and tb in global_graph:
                        dist = nx.shortest_path_length(global_graph.to_undirected(), source=ta, target=tb)
                        avg_dist += dist
                        pairs += 1
                except nx.NetworkXNoPath:
                    pass

        proximity_score = 0
        if pairs > 0:
            avg_dist /= pairs
            # Closer targets (lower avg_dist) -> Higher synergy potential
            proximity_score = max(0, 1.0 - (avg_dist / 5.0))

        total_score = (proximity_score * 0.4) + (len(shared_pathways) / 10.0 * 0.6)

        return {
            "combined_score": min(0.95, total_score),
            "shared_pathways": pathway_synergy,
            "proximity_score": proximity_score,
            "description": f"Predicted synergy based on {len(shared_pathways)} shared pathways and network proximity."
        }
