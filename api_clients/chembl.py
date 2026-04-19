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
        """Fetches biological targets for a ChEMBL ID."""
        url = f"{self.BASE_URL}/mechanism?molecule_chembl_id={chembl_id}"
        response = self.session.get(url, timeout=15)
        targets = []
        if response.status_code == 200:
            data = response.json()
            for tech in data.get("mechanisms", []):
                target_name = tech.get("target_name")
                if target_name:
                    targets.append(target_name)
        return list(set(targets))

if __name__ == "__main__":
    client = ChEMBLClient()
    chid = client.get_chembl_id("aspirin")
    print(f"ChEMBL ID: {chid}")
    if chid:
        targets = client.get_protein_targets(chid)
        print(f"Targets: {targets}")
