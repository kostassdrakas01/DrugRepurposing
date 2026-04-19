import requests

class MyGeneClient:
    BASE_URL = "https://mygene.info/v3/query"

    def __init__(self):
        self.session = requests.Session()

    def get_kegg_id(self, query):
        """
        Translates a gene symbol, UniProt ID, or name into a KEGG Gene ID.
        Returns a list of KEGG IDs (e.g., ['hsa:5742']).
        """
        params = {
            "q": query,
            "fields": "kegg,symbol",
            "species": "human", # Focusing on human for now as requested by 'hsa' examples
            "limit": 1
        }
        response = self.session.get(self.BASE_URL, params=params, timeout=15)
        kegg_ids = []
        if response.status_code == 200:
            data = response.json()
            if data.get("hits"):
                hit = data["hits"][0]
                kegg_data = hit.get("kegg")
                
                # Option 1: Explicit KEGG data
                if kegg_data:
                    if isinstance(kegg_data, dict):
                        if "id" in kegg_data:
                            kegg_ids.append(kegg_data["id"])
                        else:
                            for val in kegg_data.values():
                                if isinstance(val, str) and ":" in val:
                                    kegg_ids.append(val)
                    elif isinstance(kegg_data, str):
                        kegg_ids.append(kegg_data)
                    elif isinstance(kegg_data, list):
                        kegg_ids.extend(kegg_data)
                
                # Option 2: Fallback to Entrez Gene ID (which matches KEGG ID number)
                if not kegg_ids and hit.get("_id"):
                    # Map to 'hsa' for human as it's the most common request
                    kegg_ids.append(f"hsa:{hit['_id']}")
        
        # Format check: ensure it starts with a prefix
        return [kid for kid in kegg_ids if isinstance(kid, str) and ":" in kid]

    def get_functional_summary(self, symbol: str) -> str:
        """Fetches the RefSeq functional summary for a gene symbol."""
        params = {
            "q": symbol,
            "fields": "summary",
            "species": "human",
            "limit": 1
        }
        try:
            response = self.session.get(self.BASE_URL, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if data.get("hits"):
                    return data["hits"][0].get("summary", "Biological function details unavailable.")
        except Exception:
            pass
        return "Functional summary not found."

if __name__ == "__main__":
    client = MyGeneClient()
    print(client.get_kegg_id("PTGS1"))
