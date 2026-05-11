from core import DrugPathwayAnalyzer
from rich.console import Console

def test_reverse_discovery():
    console = Console()
    analyzer = DrugPathwayAnalyzer()
    
    disease = "Alzheimer"
    console.print(f"[bold yellow]Testing Reverse Discovery for: {disease}[/bold yellow]")
    
    efo_id = analyzer.opentargets.get_efo_id(disease)
    if not efo_id:
        console.print("Could not resolve disease")
        return
    
    console.print(f"EFO ID: {efo_id}")
    
    # Simulate the logic in main.py
    drivers = analyzer.opentargets.get_targets_by_disease(efo_id, limit=5)
    console.print(f"Top 5 Drivers: {[d['symbol'] for d in drivers]}")
    
    # Expand from first driver
    if drivers:
        symbol = drivers[0]['symbol']
        console.print(f"Expanding from {symbol}...")
        neighbors = analyzer.string_db.get_interactors(symbol, limit=5, required_score=700)
        console.print(f"Neighbors: {[n['symbol'] for n in neighbors]}")
        
        if neighbors:
            n_symbol = neighbors[0]['symbol']
            console.print(f"Fetching drugs for neighbor {n_symbol}...")
            drugs = analyzer.opentargets.get_drugs_by_target(n_symbol, limit=3)
            for d in drugs:
                console.print(f" - Drug: {d['name']} (Indications: {d.get('indication_names', [])[:2]})")

if __name__ == "__main__":
    test_reverse_discovery()
