import argparse
import sys
import os
from core import DrugPathwayAnalyzer, Visualizer, NetworkVisualizer, SankeyVisualizer
from rich.console import Console
from rich.panel import Panel

def main():
    parser = argparse.ArgumentParser(description="Drug-to-Pathway Predictive Discovery Engine")
    parser.add_argument("drug", nargs="?", help="Name of the drug to analyze (e.g., Aspirin, Thalidomide)")
    parser.add_argument("--export", help="Filename for the Discovery Report (default: <drug_name>.md)")
    parser.add_argument("--tissue", nargs="+", help="GTEx Tissue mask(s) (e.g., Lung Brain_Cortex) to filter irrelevant pathways")
    parser.add_argument("--disease", help="Name of the disease to find drugs for (Reverse Discovery)")
    parser.add_argument("--target", help="Name of the protein target to find drugs for (e.g., EGFR, MTOR)")
    parser.add_argument("--repurpose", action="store_true", help="Discovery Mode: Find novel drugs by excluding existing therapies")
    
    args = parser.parse_args()
    if not args.drug and not args.disease and not args.target:
        parser.error("At least one of [drug], [--disease], or [--target] must be provided.")
    console = Console()
    
    # Ensure results directory exists
    os.makedirs("results", exist_ok=True)
    
    analyzer = DrugPathwayAnalyzer(tissue_mask=args.tissue)
    visualizer = Visualizer()
    net_visuals = NetworkVisualizer("results")
    sankey_visuals = SankeyVisualizer("results")
    
    try:
        if args.disease:
            console.print(f"[bold yellow]Reverse Discovery Mode:[/bold yellow] Indication: [bold white]{args.disease}[/bold white]")
            efo_id = analyzer.opentargets.get_efo_id(args.disease)
            if not efo_id:
                console.print(f"[bold red]Error:[/bold red] Could not resolve disease: {args.disease}")
                sys.exit(1)
            
            if args.repurpose:
                console.print(f"[bold cyan]Network-Based Discovery Mode:[/bold cyan] Scanning the {args.disease} molecular interactome...")
                
                # 1. Expand Driver Set (Top 20)
                console.print(" - Fetching top 20 primary disease drivers (Open Targets)...")
                drivers = analyzer.opentargets.get_targets_by_disease(efo_id, limit=20)
                if not drivers:
                    console.print("[bold red]Error:[/bold red] No drivers found for this disease.")
                    sys.exit(1)
                
                driver_symbols = {d['symbol'] for d in drivers}
                
                # 2. Interactome Expansion (String-DB neighbors for drivers)
                console.print(" - Expanding network via Protein-Protein Interactions (String-DB)...")
                expanded_nodes = [] # List of {"symbol": str, "source_driver": str, "score": float}
                seen_neighbors = set(driver_symbols)
                
                for d in drivers[:15]: # Expanded from 10 to 15
                    neighbors = analyzer.string_db.get_interactors(d['symbol'], limit=15, required_score=700)
                    for n in neighbors:
                        if n['symbol'] not in seen_neighbors:
                            expanded_nodes.append({
                                "symbol": n['symbol'],
                                "source_driver": d['symbol'],
                                "score": n['score']
                            })
                            seen_neighbors.add(n['symbol'])
                
                # 3. Functional Pathway Context & Semantic Prioritization (Graphwise AI inspired)
                console.print(" - Stratifying top 30 targets using Semantic AI metrics...")
                discovery_targets = []
                
                # Sort neighbors by PPI score first to pick the most likely functional partners
                sorted_expanded = sorted(expanded_nodes, key=lambda x: x['score'], reverse=True)
                
                for node in sorted_expanded[:30]: # Limit to top 30 for performance
                    k_ids = analyzer.mygene.get_kegg_id(node['symbol'])
                    if k_ids:
                        pathway_ids = analyzer.kegg.get_pathways_for_gene(k_ids[0])
                        is_relevant = False
                        for pid in pathway_ids:
                            cat = analyzer.kegg.get_pathway_category(pid)
                            if any(kw in cat.lower() for kw in ["signal", "disease", "biological", "nervous", "metabolism"]):
                                is_relevant = True; break
                        
                        if is_relevant:
                            # GRAPHWISE IDEA: Semantic Scoring (PPI + Literature + Tractability)
                            assoc = analyzer.opentargets.get_target_association_details(node['symbol'], efo_id)
                            tract = analyzer.opentargets.get_target_tractability(node['symbol'])
                            
                            # Semantic Priority Score (SPS)
                            semantic_score = (node['score'] / 1000.0 * 0.4) + \
                                             (assoc.get('literature', 0) * 0.4) + \
                                             ((1.0 if tract.get('small_molecule') else 0.5) * 0.2)
                            
                            node['semantic_score'] = semantic_score
                            discovery_targets.append(node)

                # 4. Tiered Drug Harvesting
                console.print(" - Harvesting drug candidates with Semantic Ranking...")
                # FUTURE IDEA: Side Effect Similarity Discovery (Sildenafil/Minoxidil pattern)
                # To be implemented: Fetch adverse events for standard drugs and find matches in candidates.
                
                exclude_keywords = set(args.disease.lower().split())
                candidate_pool = []
                seen_drug_ids = set()
                
                # Priority 1: High Semantic Score neighbors
                sorted_discovery = sorted(discovery_targets, key=lambda x: x.get('semantic_score', 0), reverse=True)
                # Priority 2: Hidden drivers (Targets 4-20)
                hidden_drivers = [{"symbol": d['symbol'], "score": d['score'], "is_driver": True} for d in drivers[3:]]
                
                all_search_nodes = sorted_discovery + hidden_drivers

                def harvest(nodes, check_novelty=True):
                    found = []
                    for node in nodes:
                        t_drugs = analyzer.opentargets.get_drugs_by_target(node['symbol'], limit=5)
                        for d in t_drugs:
                            if d['id'] in seen_drug_ids: continue
                            
                            is_standard = False
                            if check_novelty:
                                for ind in d.get("indication_names", []):
                                    if all(kw in ind.lower() for kw in exclude_keywords):
                                        is_standard = True; break
                                if efo_id in d.get("indications", []):
                                    is_standard = True
                            
                            if not is_standard or not check_novelty:
                                d['target_symbol'] = node['symbol']
                                d['is_neighbor'] = not node.get('is_driver', False)
                                d['semantic_rank'] = node.get('semantic_score', 0)
                                found.append(d)
                                seen_drug_ids.add(d['id'])
                                if len(found) >= 15: return found
                    return found

                candidate_pool = harvest(all_search_nodes, check_novelty=True)

                # 5. Fallbacks if empty
                if not candidate_pool:
                    console.print("[bold yellow]No novel neighbors found. Broadening search to Pathway Hubs...[/bold yellow]")
                    # Search for any drug targeting proteins in the top 3 pathways
                    pathway_ids = set()
                    for d in drivers[:5]:
                        k_ids = analyzer.mygene.get_kegg_id(d['symbol'])
                        if k_ids: pathway_ids.update(analyzer.kegg.get_pathways_for_gene(k_ids[0]))
                    
                    hub_genes = []
                    for pid in list(pathway_ids)[:3]:
                        hub_genes.extend([{"symbol": g} for g in analyzer.kegg.get_genes_by_pathway(pid)[:20]])
                    candidate_pool = harvest(hub_genes, check_novelty=True)

                if not candidate_pool:
                    console.print("[bold red]No novel candidates found. Reverting to Standard-of-Care analysis...[/bold red]")
                    # Final fallback: just get any drugs for the disease/drivers without novelty filter
                    candidate_pool = harvest([{"symbol": d['symbol']} for d in drivers[:10]], check_novelty=False)
                
                # Rank: Neighbors > Semantic Score > Phase
                candidate_pool.sort(key=lambda x: (x.get('is_neighbor', False), x.get('semantic_rank', 0)), reverse=True)
                process_drugs(analyzer, visualizer, net_visuals, sankey_visuals, candidate_pool[:10], args.tissue, console, context_name=args.disease)
            else:
                drugs = analyzer.opentargets.get_drugs_by_disease(efo_id, limit=10)
                process_drugs(analyzer, visualizer, net_visuals, sankey_visuals, drugs, args.tissue, console, context_name=args.disease)
            
        elif args.target:
            console.print(f"[bold yellow]Reverse Discovery Mode:[/bold yellow] Target: [bold white]{args.target}[/bold white]")
            drugs = analyzer.opentargets.get_drugs_by_target(args.target, limit=10)
            process_drugs(analyzer, visualizer, net_visuals, sankey_visuals, drugs, args.tissue, console, context_name=args.target)
            
        else:
            run_analysis(analyzer, visualizer, net_visuals, sankey_visuals, args.drug, args.tissue, console)
            
    except Exception as e:
        import traceback
        console.print(f"[bold red]Error during systemic discovery:[/bold red] {e}")
        traceback.print_exc()
        sys.exit(1)

def process_drugs(analyzer, visualizer, net_visuals, sankey_visuals, drugs, tissue, console, context_name=None):
    if not drugs:
        console.print("[bold red]No candidate drugs found for the given criteria.[/bold red]")
        sys.exit(1)
    
    # Create subdirectory for the context (Disease/Target)
    base_dir = "results"
    if context_name:
        base_dir = os.path.join("results", context_name.replace(" ", "_"))
        os.makedirs(base_dir, exist_ok=True)
        # Update visualizer export paths for the visuals
        net_visuals.output_dir = base_dir
        sankey_visuals.output_dir = base_dir

    console.print(f"[bold green]Found {len(drugs)} candidates. Running systemic discovery...[/bold green]")
    all_analyses = []
    
    # Process all drugs (Lightweight for batch)
    for drug_info in drugs:
        d_name = drug_info["name"]
        try:
            analysis = analyzer.analyze_drug(d_name, tissue=tissue[0] if tissue else None)
            all_analyses.append(analysis)
            
            # Show summary in terminal
            visualizer.display_analysis(analysis)
            
            # Export individual report to the subdirectory
            export_path = os.path.join(base_dir, f"{d_name}.md")
            visualizer.export_report(analysis, export_path)
            
        except Exception as e:
            console.print(f"[bold red]Error analyzing {d_name}:[/bold red] {e}")

    # Generate consolidated report in the base_dir
    if context_name:
        summary_path = os.path.join(base_dir, "Discovery_Summary.md")
        visualizer.export_consolidated_report(all_analyses, context_name, summary_path)

    # Generate visuals only for candidates with "Novel Hits"
    console.print("[bold cyan]Generating discovery visuals for high-probability leads...[/bold cyan]")
    for analysis in all_analyses:
        if len(analysis.pathways) > 0:
            with console.status(f"[bold blue]Rendering Discovery Portfolio for {analysis.drug_name}..."):
                net_visuals.generate_all_visuals(analysis)
                sankey_visuals.generate_sankey(analysis)
        else:
            console.print(f"[dim]Skipping visuals for {analysis.drug_name} (No novel hits detected).[/dim]")

def run_analysis(analyzer, visualizer, net_visuals, sankey_visuals, drug_name, tissue, console):
    with console.status(f"[bold green]Predicting Bio-Network for {drug_name}...") as status:
        if tissue:
            status.update(f"[bold green]Applying GTEx Context: {tissue} (TPM > 10.0)...")
        analysis = analyzer.analyze_drug(drug_name, tissue=tissue[0] if tissue else None)
    
    visualizer.display_analysis(analysis)
    
    # Only generate visuals if hits exist
    if len(analysis.pathways) > 0:
        with console.status(f"[bold blue]Generating Discovery visuals for {drug_name} (JPG)..."):
            net_visuals.generate_all_visuals(analysis)
            sankey_visuals.generate_sankey(analysis)
    
    export_path = os.path.join("results", f"{drug_name}.md")
    visualizer.export_report(analysis, export_path)
    csv_path = export_path.replace(".md", ".csv")
    visualizer.export_csv(analysis, csv_path)

if __name__ == "__main__":
    main()
