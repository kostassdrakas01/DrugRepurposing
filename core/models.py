from dataclasses import dataclass, field
from typing import List, Dict, Optional

@dataclass
class Target:
    name: str
    kegg_id: Optional[str] = None
    pathways: List[str] = field(default_factory=list)
    is_primary: bool = True # True for Primary, False for Predicted Sub-Target
    interaction_label: str = "" # e.g., '+p', '-i'

@dataclass
class Pathway:
    id: str
    name: str
    description: str = ""
    summary: str = ""
    url: str = ""
    category: str = "Unknown"
    related_molecules: List[str] = field(default_factory=list)
    compounds: List[str] = field(default_factory=list)
    genes: List[str] = field(default_factory=list)
    downstream_effects: List[str] = field(default_factory=list)
    polarity: str = "Neutral" # 'Activated', 'Inhibited', 'Neutral'
    surprise_score: float = 0.0
    z_score: float = 0.0
    discovery_score: float = 0.0
    discovery_insight: str = ""

@dataclass
class Bottleneck:
    node_name: str
    pathway_count: int
    pathways: List[str]

@dataclass
class DiscoveryInsight:
    type: str  # 'Convergence', 'Oncomodulatory', 'Bottleneck', 'Cross-Category'
    description: str
    related_nodes: List[str]
    evidence: str
    novelty_score: float = 0.0

@dataclass
class CentralityNode:
    name: str
    score: float
    role: str  # 'Bridge', 'Traffic Controller'
    connections: List[str]

@dataclass
class PerturbationResult:
    target_name: str
    impacted_node: str
    change_direction: str  # 'Upstack Trigger', 'Downstream Suppression'
    estimated_impact: float  # Percentage
    evidence: str

@dataclass
class ConvergenceGroup:
    nodes: List[str]
    pathway_count: int
    related_pathways: List[str]
    description: str

@dataclass
class DrugAnalysis:
    drug_name: str
    cid: int
    target_tissue: Optional[str] = None # GTEx tissue ID
    targets: List[Target] = field(default_factory=list)
    pathways: Dict[str, Pathway] = field(default_factory=dict)
    appendix_pathways: Dict[str, Pathway] = field(default_factory=dict) # For "In-Category" noise
    connections: List[str] = field(default_factory=list)
    insights: List[DiscoveryInsight] = field(default_factory=list)
    bottlenecks: List[Bottleneck] = field(default_factory=list)
    centrality: List[CentralityNode] = field(default_factory=list)
    perturbations: List[PerturbationResult] = field(default_factory=list)
    convergence_groups: List[ConvergenceGroup] = field(default_factory=list)
    sankey_data: Dict = field(default_factory=dict)
    confidence_badge: str = "" # e.g. 'High Confidence'
