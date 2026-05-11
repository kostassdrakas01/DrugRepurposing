from api_clients import OpenTargetsClient
from rich.console import Console

def debug_alzheimer_drivers():
    console = Console()
    client = OpenTargetsClient()
    
    disease = "Alzheimer"
    efo_id = client.get_efo_id(disease)
    drivers = client.get_targets_by_disease(efo_id, limit=10)
    
    for d in drivers:
        symbol = d['symbol']
        drugs = client.get_drugs_by_target(symbol, limit=5)
        console.print(f"Driver: {symbol}, Drugs Found: {len(drugs)}")
        for dr in drugs:
            console.print(f" - {dr['name']} (Indications: {dr.get('indication_names', [])[:2]})")

if __name__ == "__main__":
    debug_alzheimer_drivers()
