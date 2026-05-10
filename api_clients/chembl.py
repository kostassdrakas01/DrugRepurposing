import requests

class ChEMBLClient:
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})

    def get_chembl_id(self, name):
        """Resolves a drug name to a ChEMBL ID."""
        url = f"{self.BASE_URL}/molecule?pref_name__iexact={name}"
        response = self.session.get(url, timeout=15)
        if response.status_code == 200:
            data = response.json()
            if data.get("molecules"):
                return data["molecules"][0]["molecule_chembl_id"]
        return None

    def get_protein_targets(self, chembl_id):
        """Fetches biological targets and their action types for a ChEMBL ID."""
        url = f"{self.BASE_URL}/mechanism?molecule_chembl_id={chembl_id}&format=json"
        response = self.session.get(url, timeout=15)
        targets = []
        if response.status_code == 200:
            data = response.json()
            for mech in data.get("mechanisms", []):
                target_name = mech.get("target_name")
                action_type = mech.get("action_type", "UNKNOWN")
                target_chembl_id = mech.get("target_chembl_id")

                # Try to get gene symbol from target details if possible, otherwise use target_name
                symbol = target_name
                if target_chembl_id:
                    t_url = f"{self.BASE_URL}/target/{target_chembl_id}.json"
                    t_resp = self.session.get(t_url, timeout=10)
                    if t_resp.status_code == 200:
                        t_data = t_resp.json()
                        components = t_data.get("target_components", [])
                        for comp in components:
                            synonyms = comp.get("target_component_synonyms", [])
                            for syn in synonyms:
                                if syn.get("syn_type") == "GENE_SYMBOL":
                                    symbol = syn.get("component_synonym")
                                    break
                            if symbol != target_name: break

                if symbol:
                    targets.append({"symbol": symbol, "action_type": action_type})

        # Deduplicate while preserving action types
        unique_targets = {}
        for t in targets:
            sym = t["symbol"]
            if sym not in unique_targets:
                unique_targets[sym] = t["action_type"]

        return [{"symbol": sym, "action_type": act} for sym, act in unique_targets.items()]

if __name__ == "__main__":
    client = ChEMBLClient()
    chid = client.get_chembl_id("aspirin")
    print(f"ChEMBL ID: {chid}")
    if chid:
        targets = client.get_protein_targets(chid)
        print(f"Targets: {targets}")
