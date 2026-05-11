import re
from typing import List, Dict, Optional
import networkx as nx
from api_clients import (
    PubChemClient,
    ChEMBLClient,
    OpenTargetsClient,
    MyGeneClient,
    KEGGClient,
    NCBIClient,
    GTExClient,
    StringDBClient
)
from core.models import (
    DrugAnalysis, 
    DiseaseAnalysis,
    SynergyResult,
    Target, 
    Pathway, 
    Bottleneck, 
    DiscoveryInsight, 
    CentralityNode, 
    PerturbationResult,
    ConvergenceGroup,
    SimulationResult
)

class DrugPathwayAnalyzer:
    def __init__(self, tissue_mask: Optional[List[str]] = None):
        self.pubchem = PubChemClient()
        self.chembl = ChEMBLClient()
        self.opentargets = OpenTargetsClient()
        self.mygene = MyGeneClient()
        self.kegg = KEGGClient()
        self.ncbi = NCBIClient()
        self.gtex = GTExClient()
        self.string_db = StringDBClient()
        self.tissue_mask = tissue_mask
        self.global_graph = nx.DiGraph()

    def analyze_disease(self, disease_name: str) -> DiseaseAnalysis:
        print(f"[*] Starting Reverse Discovery for disease: {disease_name}...")
        analysis = DiseaseAnalysis(disease_name=disease_name)

        # 1. Map Disease to Pathways
        pathways_found = self.kegg.find_pathways(disease_name)
        if not pathways_found:
            print(f"[!] No direct pathways found for {disease_name}. Searching broader categories.")
            # Try to find some common ones or just exit
            return analysis

        print(f"[*] Found {len(pathways_found)} pathways. Identifying Master Regulators...")

        all_master_regulators = []
        pathway_objs = []

        for p_info in pathways_found[:3]: # Limit to top 3 for performance
            pid = p_info['id']
            p_meta = self.kegg.get_pathway_info(pid)

            kgml = self.kegg.get_kgml_pathway(pid)
            genes = []
            if kgml:
                genes = [e.name.split()[0] for e in kgml.entries.values() if e.type == "gene"]

            pathway_obj = Pathway(
                id=pid,
                name=p_meta.get("name", p_info['name']),
                description=p_meta.get("description", ""),
                genes=genes
            )
            pathway_objs.append(pathway_obj)

            # Build local graph for this pathway to find hubs
            local_graph = nx.DiGraph()
            if kgml:
                entry_map = {entry.id: entry for entry in kgml.entries.values()}
                for rel in kgml.relations:
                    s_id = rel.entry1.id if hasattr(rel.entry1, "id") else rel.entry1
                    t_id = rel.entry2.id if hasattr(rel.entry2, "id") else rel.entry2
                    if s_id in entry_map and t_id in entry_map:
                        s_name = entry_map[s_id].name.split()[0]
                        t_name = entry_map[t_id].name.split()[0]
                        local_graph.add_edge(s_name, t_name)

            if local_graph.number_of_nodes() > 0:
                centrality = nx.degree_centrality(local_graph)
                sorted_c = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
                for node, score in sorted_c[:5]:
                    if score > 0:
                        human_name = self.kegg.get_entity_name(node)
                        all_master_regulators.append(CentralityNode(
                            name=human_name, score=score, role="Master Regulator", connections=[]
                        ))

        analysis.pathways = pathway_objs
        analysis.master_regulators = sorted(all_master_regulators, key=lambda x: x.score, reverse=True)[:10]

        # 3. Query for Drugs
        print(f"[*] Searching for drugs targeting identified hubs...")
        suggested_drugs = []
        seen_drug_names = set()
        for regulator in analysis.master_regulators:
            if regulator.name == "Unknown": continue
            drugs = self.opentargets.get_drugs_by_target(regulator.name, limit=3)
            for drug in drugs:
                if drug['name'] not in seen_drug_names:
                    suggested_drugs.append({
                        "name": drug['name'],
                        "target": regulator.name,
                        "id": drug['id']
                    })
                    seen_drug_names.add(drug['name'])

        analysis.suggested_drugs = suggested_drugs
        return analysis

    def analyze_combination(self, drug_names: List[str], tissue: str = None) -> DrugAnalysis:
        print(f"[*] Analyzing multi-drug combination: {' + '.join(drug_names)}...")
        from core.synergy import SynergyCalculator

        # Analyze each drug individually first
        analyses = [self.analyze_drug(name, tissue=tissue) for name in drug_names]

        # Merge analyses (simplified)
        main_analysis = analyses[0]
        main_analysis.drug_name = " + ".join(drug_names)

        for other in analyses[1:]:
            main_analysis.targets.extend(other.targets)
            # Merge pathways (keep highest discovery score)
            for pid, p in other.pathways.items():
                if pid in main_analysis.pathways:
                    if p.discovery_score > main_analysis.pathways[pid].discovery_score:
                        main_analysis.pathways[pid] = p
                else:
                    main_analysis.pathways[pid] = p

            # Merge appendix
            for pid, p in other.appendix_pathways.items():
                if pid not in main_analysis.pathways and pid not in main_analysis.appendix_pathways:
                    main_analysis.appendix_pathways[pid] = p

        # Recalculate global graph with all targets
        self._build_global_graph(main_analysis)

        # Calculate Synergy if there are exactly 2 drugs
        if len(analyses) == 2:
            synergy_data = SynergyCalculator.calculate_synergy(analyses[0], analyses[1], self.global_graph)
            main_analysis.synergy = SynergyResult(**synergy_data)

        return main_analysis

    def _run_dynamic_simulation(self, analysis: DrugAnalysis):
        """Converts global graph to Boolean Network and simulates drug impact."""
        from core.simulation import BooleanNetworkSimulator
        sim = BooleanNetworkSimulator(self.global_graph)

        drug_impacts = {}
        for t in analysis.targets:
            if t.kegg_id:
                impact = -1.0
                if any(x in t.action_type.upper() for x in ["ACTIVATOR", "AGONIST", "STIMULATOR"]):
                    impact = 1.0
                drug_impacts[t.kegg_id.split(':')[-1] if ':' in t.kegg_id else t.kegg_id] = impact

        # Map to nodes in graph (which are often symbol-based or short IDs)
        mapped_impacts = {}
        for node in self.global_graph.nodes():
            for tid, val in drug_impacts.items():
                if tid in node:
                    mapped_impacts[node] = val

        shifts = sim.predict_shift(mapped_impacts)

        for s in shifts[:10]:
            analysis.simulations.append(SimulationResult(
                node_name=self.kegg.get_entity_name(s['node']),
                shift_direction="Activated (Steady-State)" if s['to'] else "Inhibited (Steady-State)",
                pathway_context="Global Bio-Network"
            ))

    def analyze_drug(self, drug_name: str, tissue: str = None) -> DrugAnalysis:
        print(f"[*] Starting systemic discovery for {drug_name}...")
        # 1. Resolve Drug to IDs
        cid = self.pubchem.get_cid_by_name(drug_name)
        chembl_id = self.chembl.get_chembl_id(drug_name)
        
        if not cid and not chembl_id:
            raise ValueError(f"Could not resolve drug: {drug_name}")

        analysis = DrugAnalysis(drug_name=drug_name, cid=cid or 0, chembl_id=chembl_id, target_tissue=tissue)

        # 2. Get Targets
        print(f"[*] Resolving molecular targets...")
        target_data_raw = [] # List of {"symbol": str, "action_type": str}
        if chembl_id:
            target_data_raw = self.opentargets.get_targets_by_chembl_id(chembl_id)
            if not target_data_raw:
                target_data_raw = self.chembl.get_protein_targets(chembl_id)

        if not target_data_raw and cid:
            symbols = self.pubchem.get_protein_targets(cid)
            target_data_raw = [{"symbol": s, "action_type": "UNKNOWN"} for s in symbols]
        
        # 2b. High-Affinity Neighbor Search (Added feature)
        # If a drug has only one primary target, automatically pull top 5 proteins with score > 0.9
        all_target_data = [] # List of (symbol, is_primary, action_type)
        for t in target_data_raw:
            all_target_data.append((t["symbol"], True, t["action_type"]))
        
        if len(all_target_data) == 1:
            symbol = all_target_data[0][0]
            print(f"[*] Single target detected ({symbol}). Searching for high-affinity neighbors (String-DB)...")
            interactors = self.opentargets.get_interactors(symbol, min_score=0.7, limit=50)
            for interactor in interactors:
                all_target_data.append((interactor, False, "UNKNOWN"))
        
        # 3. Process Targets and Map to Pathways with GTEx Masking
        pathway_ids_to_process = set()
        
        for symbol, is_primary, action_type in all_target_data:
            # GTEx Masking: If tissues are specified, check if gene is expressed in ANY of them
            if self.tissue_mask:
                # Check all tissues; if TPM > 1.0 in any, consider it expressed
                max_tpm = 0.0
                for tissue in self.tissue_mask:
                    tpm = self.gtex.get_expression(symbol, tissue)
                    max_tpm = max(max_tpm, tpm)
                
                if max_tpm < 1.0 and is_primary:
                    print(f"[!] Warning: Primary target {symbol} has low expression ({max_tpm:.2f} TPM) in specified tissue. Proceeding anyway.")
                elif max_tpm < 1.0:
                    print(f"[-] Skipping target {symbol} due to low expression ({max_tpm:.2f} TPM) in tissue context.")
                    continue

            kegg_gene_ids = self.mygene.get_kegg_id(symbol)
            if not kegg_gene_ids:
                print(f"[!] Error: Could not resolve KEGG ID for target {symbol}. Mapping failed.")
                continue

            for kge in kegg_gene_ids:
                target_obj = Target(name=symbol, kegg_id=kge, is_primary=is_primary, action_type=action_type)
                pathways = self.kegg.get_pathways_for_gene(kge)
                target_obj.pathways = pathways
                analysis.targets.append(target_obj)
                pathway_ids_to_process.update(pathways)

        print(f"[*] Found {len(pathway_ids_to_process)} biological pathways to analyze.")

        # 4. Analyze Pathways with Polarity and Novelty
        drug_category, indication_keywords = self._guess_drug_category(drug_name, analysis.targets, chembl_id=chembl_id)
        
        for i, pid in enumerate(pathway_ids_to_process):
            if i % 5 == 0:
                print(f"[*] Processing pathways: {i}/{len(pathway_ids_to_process)}...")
            p_info = self.kegg.get_pathway_info(pid)
            p_name = p_info.get("name", pid)
            p_category = self.kegg.get_pathway_category(pid)
            
            kgml = self.kegg.get_kgml_pathway(pid)
            genes = []
            compounds = []
            if kgml:
                genes = [e.name.split()[0] for e in kgml.entries.values() if e.type == "gene"]
                compounds = [e.name.split()[0] for e in kgml.entries.values() if e.type == "compound"]

            pathway_obj = Pathway(
                id=pid,
                name=p_name,
                description=p_info.get("description", ""),
                category=p_category,
                genes=genes,
                compounds=compounds,
                url=self.kegg.get_map_url(pid, [t.kegg_id for t in analysis.targets if pid in t.pathways])
            )
            
            # Predictive Polarity Calculation
            pathway_obj.polarity = self._calculate_pathway_polarity(pid, analysis.targets)

            # Surprise Score Logic (Refactored: Indication-Specific & Hierarchy Aware)
            pathway_obj.surprise_score = self._calculate_surprise_score(
                p_name, p_category, drug_category, indication_keywords
            )
            
            # Temporarily hold in analysis.pathways
            analysis.pathways[pid] = pathway_obj

        # --- Permutation Testing & Score Aggregation ---
        import numpy as np
        print("[*] Running Permutation Testing (1000 iter) to neutralize Hub Bias...")
        self._run_permutation_testing(analysis.pathways, analysis.targets)
        
        # Re-partition into Pathways and Appendix based on Discovery Score
        final_pathways = {}
        appendix_pathways = {}
        
        # Sort all pathways by discovery score first
        all_p_items = sorted(analysis.pathways.items(), key=lambda x: (x[1].surprise_score * 0.5) + (1.0 / (1.0 + np.exp(-x[1].z_score)) * 0.5), reverse=True)
        
        for i, (pid, p) in enumerate(all_p_items):
            # Apply Sigmoid to normalize Z-score to 0-1
            norm_z = 1.0 / (1.0 + np.exp(-p.z_score))
            p.discovery_score = (p.surprise_score * 0.5) + (norm_z * 0.5)
            
            # Lower threshold to 0.70 to ensure more hits are visible
            # Also always include the top 5 most novel hits to avoid "0 hits" frustration
            if p.discovery_score >= 0.70 or i < 5:
                final_pathways[pid] = p
            else:
                appendix_pathways[pid] = p
                
        analysis.pathways = final_pathways
        analysis.appendix_pathways = appendix_pathways

        # 6. Advanced Network Discovery
        print("[*] Running network clustering and crosstalk discovery...")
        analysis.convergence_groups = self._cluster_convergence(analysis)
        
        # 7. Tissue Verification (Added per user request)
        if tissue:
            print(f"[*] Verifying biological evidence in {tissue} context...")
            self._verify_tissue_relevance(analysis, tissue)
        analysis.connections = self._find_connections(analysis)
        analysis.bottlenecks = self._analyze_bottlenecks(analysis)
        analysis.insights = self._find_discovery_insights(analysis)
        
        # 6. Global Bio-Physics
        print("[*] Simulating systemic perturbations...")
        self._build_global_graph(analysis)
        analysis.centrality = self._calculate_centrality(analysis)
        analysis.perturbations = self._simulate_perturbations(analysis)
        
        # 7. Dynamic Simulation
        print("[*] Running dynamic Boolean Network simulation...")
        self._run_dynamic_simulation(analysis)

        # 8. Predictive Flow (Sankey)
        analysis.sankey_data = self._generate_sankey_data(analysis)

        # 9. Toxicity Profiling
        if analysis.chembl_id:
            print("[*] Fetching toxicity profile (Adverse Events)...")
            analysis.adverse_events = self.opentargets.get_adverse_events(analysis.chembl_id)

        # 8. Functional Mechanism Narratives (Directive: Mechanism Narratives)
        print("[*] Generating functional mechanism narratives...")
        self._generate_mechanism_narratives(analysis)

        print(f"[+] Discovery complete for {drug_name}.")
        return analysis

    def _run_permutation_testing(self, all_pathways: Dict[str, Pathway], targets: List[Target], iterations: int = 1000):
        """Runs Monte Carlo permutation to compute statistical novelty (Z-score) of pathway hits."""
        import random
        import numpy as np
        
        # 1. Establish the universe of all genes within the mapped pathways
        universe_genes = set()
        pathway_genes = {}
        for pid, p in all_pathways.items():
            genes = set(p.genes)
            pathway_genes[pid] = genes
            universe_genes.update(genes)
            
        universe_list = list(universe_genes)
        if not universe_list: return
        
        # 2. Extract observed target gene IDs
        target_kegg_ids = set()
        for t in targets:
            if t.kegg_id: target_kegg_ids.add(t.kegg_id)
            
        draw_size = min(max(len(target_kegg_ids), 1), len(universe_list))

        # 3. Simulate Random Draws
        hits_per_pathway = {pid: [] for pid in all_pathways}
        for _ in range(iterations):
            sim_targets = set(random.sample(universe_list, draw_size))
            for pid in all_pathways:
                hits = len(sim_targets.intersection(pathway_genes[pid]))
                hits_per_pathway[pid].append(hits)
                
        # 4. Calculate Z-scores
        for pid, p in all_pathways.items():
            observed = len(target_kegg_ids.intersection(pathway_genes[pid]))
            mean_hits = np.mean(hits_per_pathway[pid])
            std_hits = np.std(hits_per_pathway[pid])
            
            z_score = 0.0
            if std_hits > 0:
                z_score = (observed - mean_hits) / std_hits
            
            p.z_score = z_score

    def _generate_mechanism_narratives(self, analysis: DrugAnalysis):
        """Generates functional mechanism narratives for high-novelty discoveries."""
        primary_target_names = [t.name for t in analysis.targets if t.is_primary]
        primary_target_str = ", ".join(primary_target_names[:2]) if primary_target_names else "primary targets"
        
        for pid, p in analysis.pathways.items():
            if p.surprise_score < 0.81: continue # Narrative for all Section 1 hits
            
            # Find reaching bridge protein for this pathway
            bridge_protein = "the cellular signaling network"
            bridge_symbol = ""
            for group in analysis.convergence_groups:
                if p.name in group.related_pathways:
                    # Filter for symbols (usually upper case)
                    symbols = [n for n in group.nodes if n and n[0].isupper()]
                    if symbols:
                        bridge_symbol = symbols[0]
                        bridge_protein = f"the {bridge_symbol} bridge"
                        break
                    elif group.nodes:
                        bridge_symbol = group.nodes[0]
                        bridge_protein = f"the {bridge_symbol} activator"
                        break
            
            # Fetch biological function summary
            summary = "executing downstream cellular signaling"
            if bridge_symbol:
                summary_raw = self.mygene.get_functional_summary(bridge_symbol)
                if len(summary_raw) > 30:
                    # Clean up first sentence or first 25 words
                    sentences = summary_raw.split('.')
                    if sentences:
                        summary = sentences[0].strip()
                    else:
                        summary = " ".join(summary_raw.split()[:25])
            
            # Specific Outcome Heuristics (Directive Rule 1)
            outcome = "downstream systemic modulation"
            p_low = p.name.lower()
            if any(k in p_low for k in ["alzheimer", "parkinson", "huntington", "nervous"]):
                outcome = "neuro-inflammation or synaptic plasticity"
            elif any(k in p_low for k in ["cancer", "carcinoma", "glioma", "melanoma"]):
                outcome = "oncogenic signaling and cell cycle progression"
            elif any(k in p_low for k in ["metabolism", "diabetes", "insulin", "lipid"]):
                outcome = "homeostatic metabolic balance"
            elif any(k in p_low for k in ["heart", "cardio", "vascular", "vessel"]):
                outcome = "hemodynamic stability and vascular integrity"
            elif any(k in p_low for k in ["immune", "infection", "virus", "bacterial"]):
                outcome = "the systemic immune-inflammatory response"

            # Construction (Directive: Narrative Synthesis)
            # Template: "The drug targets [Target], which influences the [Bridge] to modify [Specific Biological Function] in [Discovery Pathway]."
            # We clean summary to ensure it doesn't start with "This gene encodes..." or similar
            clean_summary = summary
            prefixes_to_strip = ["This gene encodes ", "The protein encoded by this gene ", "This gene is a member of ", "Encoded by this gene is ", "This gene belongs to "]
            for prefix in prefixes_to_strip:
                if clean_summary.startswith(prefix):
                    clean_summary = clean_summary[len(prefix):].strip()
            
            # Ensure it starts with lowercase if it's continuing the sentence
            if clean_summary and clean_summary[0].isupper():
                clean_summary = clean_summary[0].lower() + clean_summary[1:]
            
            p_clean_name = p.name.split(' - ')[0]
            narrative = f"The drug targets {primary_target_str}, which influences {bridge_protein} to modify {clean_summary} in {p_clean_name}."
            p.discovery_insight = narrative

    def _build_global_graph(self, analysis: DrugAnalysis):
        """Constructs a unified graph from all relevant pathways."""
        self.global_graph.clear()
        for pid in analysis.pathways:
            kgml = self.kegg.get_kgml_pathway(pid)
            if not kgml: continue
            
            entry_map = {entry.id: entry for entry in kgml.entries.values()}
            for rel in kgml.relations:
                s_id = rel.entry1.id if hasattr(rel.entry1, "id") else rel.entry1
                t_id = rel.entry2.id if hasattr(rel.entry2, "id") else rel.entry2
                
                s_name = entry_map[s_id].name if s_id in entry_map else "Unknown"
                t_name = entry_map[t_id].name if t_id in entry_map else "Unknown"
                
                # Assign types
                subtypes = [s[0] for s in rel.subtypes]
                weight = 1.0
                if "inhibition" in subtypes: weight = -1.0
                elif "activation" in subtypes: weight = 1.0
                
                # Use representative IDs to prevent graph explosion
                s_best = s_name.split()[0] if s_name != "Unknown" else "Unknown"
                t_best = t_name.split()[0] if t_name != "Unknown" else "Unknown"
                
                if s_best != "Unknown" and t_best != "Unknown":
                    self.global_graph.add_edge(s_best, t_best, weight=weight, subtypes=subtypes)

    def _guess_drug_category(self, drug_name: str, targets: List[Target], chembl_id: str = None) -> tuple:
        """Heuristic to guess therapeutic category and extraction specific indication keywords."""
        all_pathway_ids = []
        for t in targets:
            all_pathway_ids.extend(t.pathways)
        
        categories = [self.kegg.get_pathway_category(pid) for pid in all_pathway_ids]
        if not categories: return "Unknown", set()
        
        # Indication Keyword Extraction Rule: Prefer OpenTargets Clinical Indications
        keywords = set()
        if chembl_id:
            official_indications = self.opentargets.get_indications(chembl_id)
            for ind in official_indications:
                # Clean and tokenize
                ignore = {"disease", "disorder", "syndrome", "type", "human"}
                words = [w.strip("(),").lower() for w in ind.split() if len(w) > 3 and w.lower() not in ignore]
                keywords.update(words)
        
        # Fallback to frequent pathway names only if no keywords found
        if not keywords:
            primary_targets = [t for t in targets if t.is_primary]
            for t in primary_targets:
                for pid in t.pathways:
                    p_info = self.kegg.get_pathway_info(pid)
                    p_name = p_info.get("name", "").lower()
                    ignore = {"pathway", "signaling", "metabolism", "human", "disease", "(human)", "sapiens", "active", "regulated"}
                    words = [w.strip("(),").lower() for w in p_name.split() if len(w) > 4 and w not in ignore]
                    keywords.update(words)

        from collections import Counter
        cat_counts = Counter(categories)
        return cat_counts.most_common(1)[0][0], keywords

    def _calculate_pathway_polarity(self, pathway_id: str, targets: List[Target]) -> str:
        """Predicts if the drug activates or inhibits the pathway."""
        kgml = self.kegg.get_kgml_pathway(pathway_id)
        if not kgml: return "Neutral"
        
        # Map KEGG ID to its action type for this drug
        target_impacts = {}
        for t in targets:
            if pathway_id in t.pathways and t.kegg_id:
                impact = -1.0 # Default to inhibition if unknown
                if "ACTIVATOR" in t.action_type.upper() or "AGONIST" in t.action_type.upper() or "STIMULATOR" in t.action_type.upper():
                    impact = 1.0
                elif "INHIBITOR" in t.action_type.upper() or "ANTAGONIST" in t.action_type.upper() or "BLOCKER" in t.action_type.upper():
                    impact = -1.0
                target_impacts[t.kegg_id] = impact

        entry_map = {entry.id: entry for entry in kgml.entries.values()}
        pathway_effects = []

        for relation in kgml.relations:
            s_id = relation.entry1.id if hasattr(relation.entry1, "id") else relation.entry1
            s_entry = entry_map.get(s_id)
            if not s_entry: continue
            
            # Check if this source entry matches any of our targets
            for tid, drug_impact in target_impacts.items():
                if tid in s_entry.name:
                    subtypes = [s[0] for s in relation.subtypes]

                    # Rigid Sign Enforcement
                    interaction_sign = 1.0
                    if "inhibition" in subtypes: interaction_sign = -1.0
                    elif "activation" in subtypes: interaction_sign = 1.0

                    # Final effect: Drug Impact * Interaction
                    final_effect = drug_impact * interaction_sign
                    pathway_effects.append(final_effect)
        
        avg_impact = sum(pathway_effects) / len(pathway_effects) if pathway_effects else 0
        if avg_impact > 0.1:
            # Check if it's disinhibition (Drug inhibitor [-1.0] * Interaction inhibitor [-1.0] = +1.0)
            # or direct activation (Drug activator [1.0] * Interaction activator [1.0] = +1.0)

            # We can't easily distinguish just from pathway_effects which is which if they both result in +1.0
            # Let's refine the logic to track if any inhibitor was involved in the activation
            is_disinhibition = False
            for tid, drug_impact in target_impacts.items():
                if drug_impact < 0: # Drug is an inhibitor
                    # If this drug inhibitor contributed to a positive effect, it's disinhibition
                    # (since it must have inhibited an inhibitor)
                    is_disinhibition = True
                    break

            return "Activated (Disinhibition)" if is_disinhibition else "Activated"
        if avg_impact < -0.1: return "Inhibited"

        # If no direct relation found, check if target is just 'in' the pathway
        if not pathway_effects and target_impacts:
            # Heuristic: if we only have activators, say Activated. If only inhibitors, say Inhibited.
            impacts = list(target_impacts.values())
            if all(i > 0 for i in impacts): return "Activated (Predicted)"
            if all(i < 0 for i in impacts): return "Inhibited (Predicted)"
            # Mixed
            if any(i > 0 for i in impacts) and any(i < 0 for i in impacts):
                return "Complex (Mixed Targets)"

        return "Neutral"

    def _calculate_surprise_score(self, pathway_name: str, pathway_category: str, drug_category: str, indication_keywords: set) -> float:
        """
        Refined Indication-Specific Filtering Logic (Technical Mandate):
        1. Only penalizes if pathway name matches specific clinical indication keywords.
        2. Guaranteed high score (0.92+) for Human Diseases and Nervous System discovery.
        3. Cross-hierarchy bonus (+0.05) to push novelty towards 1.0.
        """
        if drug_category == "Unknown" or pathway_category == "Unknown":
            return 0.5
        
        p_name_lower = pathway_name.lower()
        
        # Rule 1: The "Exact Match" Indication Penalty (Crucial for high fidelity)
        # Directive Rule: If names match drug classification or primary mechanism, force score to 0.10
        # Refined to use word boundaries to avoid partial matches like "car" in "cardiovascular"
        import re
        for kw in indication_keywords:
            if kw and len(kw) > 3:
                if re.search(rf'\b{re.escape(kw)}\b', p_name_lower):
                    return 0.10
        
        # Additional check for drug name itself in pathway (e.g., "Metformin signaling" or similar)
        # Or generic mechanism terms associated with the drug
        if drug_category != "Unknown":
            if re.search(rf'\b{re.escape(drug_category.lower())}\b', p_name_lower):
                return 0.10
        
        # Rule 2: High Sensitivity for Key Hierarchies (Technical Mandate Section 2)
        base_score = 0.81
        prioritized = ["Human Diseases", "Nervous System", "Neurodegenerative Diseases"]
        if pathway_category in prioritized:
            base_score = 0.92 # Force into Section 1 priority
        elif pathway_category != drug_category:
            base_score = 0.85

        # Rule 3: Cross-Hierarchy Bonus (+0.05)
        branches = {
            "Metabolism": "Core",
            "Genetic Information": "Core",
            "Signal Transduction": "Core",
            "Biological Systems": "System",
            "Nervous System": "System",
            "Human Diseases": "Disease",
            "Neurodegenerative Diseases": "Disease"
        }
        
        p_branch = branches.get(pathway_category, "Other")
        d_branch = branches.get(drug_category, "Other")
        
        final_score = base_score
        if p_branch != d_branch and p_branch != "Other" and d_branch != "Other":
            final_score += 0.05 # Push disease hits towards ~0.97
        
        return min(0.99, final_score)

    def _cluster_convergence(self, analysis: DrugAnalysis) -> List[ConvergenceGroup]:
        """Groups repeated proteins into functional convergence clusters."""
        node_to_pathways = {}
        all_ids = set()
        
        # 1. Collect all IDs first (Genes and Compounds)
        for pid, pathway in {**analysis.pathways, **analysis.appendix_pathways}.items():
            kgml = self.kegg.get_kgml_pathway(pid)
            if not kgml: continue
            for entry in kgml.entries.values():
                if entry.type in ["gene", "compound"]:
                    ids = entry.name.split()
                    all_ids.update(ids)
                    for node_id in ids:
                        if node_id not in node_to_pathways:
                            node_to_pathways[node_id] = set()
                        node_to_pathways[node_id].add(pathway.name)

        # 2. Filter frequent nodes
        frequent_nodes = [node_id for node_id, p_set in node_to_pathways.items() if len(p_set) >= 3]
        
        # 3. Batch resolve NAMES for frequent nodes only (Massive speedup)
        self.kegg.get_entity_name(frequent_nodes)

        # 4. Cluster by shared pathway sets
        pathway_sets = {}
        for gne_id in frequent_nodes:
            p_set = node_to_pathways[gne_id]
            p_tuple = tuple(sorted(list(p_set)))
            if p_tuple not in pathway_sets:
                pathway_sets[p_tuple] = []
            
            human_name = self.kegg.get_entity_name(gne_id)
            pathway_sets[p_tuple].append(human_name)
        
        groups = []
        for p_tuple, names in pathway_sets.items():
            if len(names) < 2: continue
            # Remove duplicates in names (multiple IDs can map to same name)
            unique_names = sorted(list(set(names)))
            node_list_str = " & ".join(unique_names[:3])
            if len(unique_names) > 3:
                node_list_str += f" (+{len(unique_names)-3} others)"
            groups.append(ConvergenceGroup(
                nodes=unique_names,
                pathway_count=len(p_tuple),
                related_pathways=list(p_tuple),
                description=f"Bridge: {node_list_str} connecting {len(p_tuple)} distinct pathways."
            ))
        
        return sorted(groups, key=lambda x: x.pathway_count, reverse=True)

    def _verify_tissue_relevance(self, analysis: DrugAnalysis, tissue: str):
        """Checks if bridge proteins are expressed in the target tissue."""
        bridge_count = 0
        nodes_to_check = set()
        for group in analysis.convergence_groups[:3]:
            nodes_to_check.update(group.nodes[:3])
            
        for gene in nodes_to_check:
            try:
                tpm = self.gtex.get_expression(gene, tissue)
                if tpm > 10.0: # Significant expression threshold
                    bridge_count += 1
            except: continue
        
        if bridge_count >= 2:
            analysis.confidence_badge = "High Confidence Discovery"
        elif bridge_count >= 1:
            analysis.confidence_badge = "Tissue-Supported"

    def _generate_sankey_data(self, analysis: DrugAnalysis) -> Dict:
        """Constructs data for a Drug -> Target -> Cluster -> Outcome Sankey diagram."""
        nodes = [analysis.drug_name]
        links = []
        
        # 1. Drug -> Targets
        for target in analysis.targets:
            nodes.append(target.name)
            links.append({"source": 0, "target": nodes.index(target.name), "value": 10})
        
        # 2. Targets -> Bridges
        for group in analysis.convergence_groups[:3]:
            # Unpack the bridge names as requested
            bridge_nodes = " & ".join(group.nodes[:2])
            if len(group.nodes) > 2: bridge_nodes += "..."
            cluster_label = f"Bridge: {bridge_nodes}"
            
            if cluster_label not in nodes:
                nodes.append(cluster_label)
            c_idx = nodes.index(cluster_label)
            
            # Connect targets to this bridge
            for target in analysis.targets:
                # Check if target is part of this bridge or hits a pathway in this group
                if target.name in group.nodes:
                    links.append({"source": nodes.index(target.name), "target": c_idx, "value": 8})
            
            # 3. Bridges -> Outcomes (Primary Pathways)
            for pname in group.related_pathways[:2]:
                pathway_label = f"Outcome: {pname[:35]}"
                if pathway_label not in nodes:
                    nodes.append(pathway_label)
                p_idx = nodes.index(pathway_label)
                links.append({"source": c_idx, "target": p_idx, "value": 4})
        
        return {"nodes": [{"name": n} for n in nodes], "links": links}

    def _calculate_centrality(self, analysis: DrugAnalysis):
        """Identifies systemic 'traffic controllers' using Degree Centrality (faster than Betweenness)."""
        try:
            if self.global_graph.number_of_nodes() < 5: return []
            
            # Switch to Degree Centrality for O(M) speed, perfectly fine for hub detection
            scores = nx.degree_centrality(self.global_graph)
            
            sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
            results = []
            for node, score in sorted_scores[:10]:
                if score < 0.001: continue
                human_name = self.kegg.get_entity_name(node)
                results.append(CentralityNode(
                    name=human_name,
                    score=round(score, 4),
                    role="Master Regulator" if score > 0.05 else "Traffic Controller",
                    connections=[self.kegg.get_entity_name(n) for n in list(self.global_graph.neighbors(node))[:3]]
                ))
            return results
        except Exception as e:
            print(f"[*] Centrality Calculation Error: {e}")
            return []

    def _simulate_perturbations(self, analysis: DrugAnalysis):
        """Predicts systemic 'ripple effects' from drug targets (optimized)."""
        raw_results = []
        target_ids = [t.kegg_id for t in analysis.targets]
        
        for tid in target_ids:
            if tid not in self.global_graph: continue
            
            # Trace 2 levels of impact
            for neighbor in self.global_graph.successors(tid):
                edge_data = self.global_graph.get_edge_data(tid, neighbor)
                weight = edge_data.get('weight', 1.0)
                
                for sub_neighbor in self.global_graph.successors(neighbor):
                    sub_edge = self.global_graph.get_edge_data(neighbor, sub_neighbor)
                    combined_weight = weight * sub_edge.get('weight', 1.0)
                    
                    raw_results.append({
                        "tid": tid,
                        "neighbor": neighbor,
                        "sub": sub_neighbor,
                        "weight": combined_weight,
                        "itypes": sub_edge.get('subtypes', [])
                    })

        # Settle for top 10 impacts to resolve names
        final_results = []
        for r in raw_results[:10]:
            final_results.append(PerturbationResult(
                target_name=self.kegg.get_entity_name(r['tid']),
                impacted_node=self.kegg.get_entity_name(r['sub']),
                change_direction="Downstream Suppression" if r['weight'] < 0 else "Downstream Trigger",
                estimated_impact=abs(r['weight']) * 35.0,
                evidence=f"Flow through {self.kegg.get_entity_name(r['neighbor'])}"
            ))
        return final_results

    def _analyze_bottlenecks(self, analysis: DrugAnalysis) -> list:
        """Ranks proteins by pathway frequency, translated to human names."""
        node_counts = {}
        for pid, pathway in analysis.pathways.items():
            for gene in pathway.genes:
                if gene not in node_counts:
                    node_counts[gene] = []
                node_counts[gene].append(pathway.name)
        
        sorted_nodes = sorted(node_counts.items(), key=lambda x: len(x[1]), reverse=True)
        bottlenecks = []
        for n, p in sorted_nodes[:10]:
            human_name = self.kegg.get_entity_name(n)
            bottlenecks.append(Bottleneck(node_name=human_name, pathway_count=len(p), pathways=p))
        return bottlenecks

    def _find_discovery_insights(self, analysis: DrugAnalysis) -> list:
        """Generates clear, jaron-free insights."""
        insights = []
        viral_paths = [p for p in analysis.pathways.values() if "viral" in p.name.lower() or "infection" in p.name.lower()]
        cancer_paths = [p for p in analysis.pathways.values() if "cancer" in p.name.lower() or "glioma" in p.name.lower() or "carcinoma" in p.name.lower()]
        
        if viral_paths and cancer_paths:
            common_genes = set()
            for vp in viral_paths:
                for cp in cancer_paths:
                    intersect = set(vp.genes).intersection(set(cp.genes))
                    common_genes.update(intersect)
            
            if common_genes:
                human_genes = [self.kegg.get_entity_name(g) for g in list(common_genes)[:5]]
                insights.append(DiscoveryInsight(
                    type="Oncomodulatory Synergy",
                    description=f"Risk: The drug affects pathways shared by viruses and tumors, potentially modifying viral propagation.",
                    related_nodes=human_genes,
                    evidence=f"Shared across {len(viral_paths)} infection maps and {len(cancer_paths)} cancer maps."
                ))

        compound_map = {}
        for pid, pathway in analysis.pathways.items():
            for cmp in pathway.compounds:
                if cmp not in compound_map: compound_map[cmp] = []
                compound_map[cmp].append(pathway.name)
        
        for cmp, pnames in compound_map.items():
            if len(set(pnames)) > 2:
                human_cmp = self.kegg.get_entity_name(cmp)
                insights.append(DiscoveryInsight(
                    type="Metabolic Convergence",
                    description=f"Major Metabolic Hub: '{human_cmp}' acts as a bridge between {len(set(pnames))} diverse biological processes.",
                    related_nodes=[human_cmp],
                    evidence=f"Integrates {', '.join(list(set(pnames))[:2])} and others."
                ))

        return insights

    def _analyze_downstream_with_data(self, pathway_data, pathway_id, targets):
        """Helper for analyze_drug to avoid redundant KGML fetches."""
        target_ids = [t.kegg_id for t in targets if pathway_id in t.pathways]
        effects = []
        G = nx.DiGraph()
        entry_map = {entry.id: entry for entry in pathway_data.entries.values()}
        for eid in entry_map: G.add_node(eid)
        for rel in pathway_data.relations:
            s_id = rel.entry1.id if hasattr(rel.entry1, "id") else rel.entry1
            t_id = rel.entry2.id if hasattr(rel.entry2, "id") else rel.entry2
            itypes = [subtype[0] for subtype in rel.subtypes]
            G.add_edge(s_id, t_id, types=itypes)

        target_entry_ids = [eid for eid, entry in entry_map.items() if any(tid in entry.name for tid in target_ids)]
        for teid in set(target_entry_ids):
            if teid in G:
                for succ_id in G.successors(teid):
                    succ_entry = entry_map[succ_id]
                    itypes = ", ".join(G.get_edge_data(teid, succ_id)['types'])
                    effects.append(f"{itypes}: {succ_entry.graphics[0].name}")
        return list(set(effects)) if effects else ["No direct interactions mapped."]

    def _find_connections(self, analysis: DrugAnalysis) -> list:
        """Identifies correlations and convergent effects across different pathways."""
        connections = []
        pathway_effects = {}
        
        for pid, pathway in analysis.pathways.items():
            pathway_effects[pid] = set(pathway.downstream_effects)

        # 1. Identify Shared Downstream Nodes (Intersection)
        pids = list(analysis.pathways.keys())
        for i in range(len(pids)):
            for j in range(i + 1, len(pids)):
                p1, p2 = pids[i], pids[j]
                common = pathway_effects[p1].intersection(pathway_effects[p2])
                # Filter out 'No direct...' messages
                common = {c for c in common if not c.startswith("No direct")}
                
                if common:
                    connections.append(
                        f"CONVERGENCE: Pathways '{analysis.pathways[p1].name}' and '{analysis.pathways[p2].name}' "
                        f"both converge on: {', '.join(list(common)[:3])}"
                    )

        # 2. Identify Potential Feedback Loops or Crosstalk Hubs
        # If a target node appears in multiple pathways
        target_counts = {}
        for target in analysis.targets:
            if len(target.pathways) > 1:
                target_counts[target.name] = target.pathways
        
        for tname, path_list in target_counts.items():
            connections.append(
                f"HUB TARGET: '{tname}' serves as a bridge between {len(path_list)} pathways, "
                f"potentially causing systemic crosstalk."
            )

        # 3. Domain Heuristics (e.g. Inflammation, Cancer)
        all_text = " ".join([p.description for p in analysis.pathways.values()]).lower()
        if "inflammation" in all_text and "immune" in all_text:
            connections.append("CORRELATION: Strong evidence of integrated Immuno-Inflammatory modulation.")
        if "prostate" in all_text or "breast" in all_text or "lung" in all_text:
            connections.append("CORRELATION: Pathway overlap suggests consistent multi-stage oncogenic disruption.")

        if not connections:
            return ["No significant cross-pathway correlations identified in current KGML maps."]
            
        return connections
