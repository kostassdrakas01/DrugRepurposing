import argparse
import sys
import os
from core import DrugPathwayAnalyzer, Visualizer, NetworkVisualizer, SankeyVisualizer
from rich.console import Console
from rich.panel import Panel

def main():
    parser = argparse.ArgumentParser(description="Drug-to-Pathway Predictive Discovery Engine")
    parser.add_argument("drug", nargs="*", help="Name(s) of the drug(s) to analyze (e.g., Aspirin, Thalidomide)")
    parser.add_argument("--disease", help="Name of a disease for Reverse Discovery (e.g., 'Alzheimer')")
    parser.add_argument("--export", help="Filename for the Discovery Report (default: <drug_name>.md)")
    parser.add_argument("--tissue", nargs="+", help="GTEx Tissue mask(s) (e.g., Lung Brain_Cortex) to filter irrelevant pathways")
    
    args = parser.parse_args()
    console = Console()
    
    # Ensure results directory exists
    os.makedirs("results", exist_ok=True)
    
    analyzer = DrugPathwayAnalyzer(tissue_mask=args.tissue)
    visualizer = Visualizer()
    net_visuals = NetworkVisualizer("results")
    sankey_visuals = SankeyVisualizer("results")
    
    try:
        if args.disease:
            with console.status(f"[bold magenta]Running Reverse Discovery for {args.disease}...") as status:
                analysis = analyzer.analyze_disease(args.disease)
            visualizer.display_disease_analysis(analysis)
            return

        if not args.drug:
            parser.print_help()
            return

        # Handle multiple drugs
        drug_name = " + ".join(args.drug)
        with console.status(f"[bold green]Predicting Bio-Network for {drug_name}...") as status:
            if args.tissue:
                status.update(f"[bold green]Applying GTEx Context: {args.tissue} (TPM > 10.0)...")

            if len(args.drug) == 1:
                analysis = analyzer.analyze_drug(args.drug[0], tissue=args.tissue[0] if args.tissue else None)
            else:
                analysis = analyzer.analyze_combination(args.drug, tissue=args.tissue[0] if args.tissue else None)

        # 1. Console Display
        visualizer.display_analysis(analysis)
        
        # 2. Advanced Visualizations (Heatmaps, Hubs, and Sankey)
        with console.status("[bold blue]Generating Discovery visuals (JPG)..."):
            net_visuals.generate_all_visuals(analysis)
            sankey_visuals.generate_sankey(analysis)
        
        # 3. Export Unified Discovery Report (Markdown + Embedded Visuals)
        export_name = args.export or ("_".join(args.drug) + ".md")
        export_path = export_name
        if not os.path.isabs(export_path) and not export_path.startswith("results/"):
            export_path = os.path.join("results", export_path)
        
        visualizer.export_report(analysis, export_path)
        
        # 4. Export Raw Data for R
        csv_path = export_path.replace(".md", ".csv")
        visualizer.export_csv(analysis, csv_path)
            
    except Exception as e:
        import traceback
        console.print(f"[bold red]Error during systemic discovery:[/bold red] {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
