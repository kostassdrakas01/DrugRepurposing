from api_clients import OpenTargetsClient
from rich.console import Console

def debug_ids():
    console = Console()
    client = OpenTargetsClient()
    
    symbols = ["APP", "PSEN1", "ACE", "EGFR"]
    for sym in symbols:
        # Simulate internal search
        q = 'query { search(queryString: "' + sym + '", entityNames: ["target"]) { hits { id } } }'
        r = client.session.post(client.URL, json={'query': q})
        hits = r.json().get("data", {}).get("search", {}).get("hits", [])
        eid = hits[0]["id"] if hits else "NOT_FOUND"
        
        # Now get drugs
        drugs = client.get_drugs_by_target(sym)
        console.print(f"Symbol: {sym}, Ensembl: {eid}, Drugs: {len(drugs)}")

if __name__ == "__main__":
    debug_ids()
