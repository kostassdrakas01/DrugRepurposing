import requests
import time

class PubChemClient:
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    VIEW_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data"

    def __init__(self):
        self.session = requests.Session()

    def get_cid_by_name(self, name):
        """Resolves a drug name to a CID."""
        url = f"{self.BASE_URL}/compound/name/{name}/cids/json"
        response = self.session.get(url)
        if response.status_code == 200:
            data = response.json()
            return data["IdentifierList"]["CID"][0]
        return None

    def get_protein_targets(self, cid):
        """
        Fetches protein targets for a given CID using the Annotation service (PUG View).
        Looks for 'Drug Targets' or 'Protein Targets'.
        """
        url = f"{self.VIEW_URL}/compound/{cid}/JSON?heading=Drug+Targets"
        response = self.session.get(url)
        targets = []
        
        if response.status_code == 200:
            data = response.json()
            targets.extend(self._parse_annotations(data))
        
        # Try another heading if empty
        if not targets:
            url = f"{self.VIEW_URL}/compound/{cid}/JSON?heading=Targets"
            response = self.session.get(url)
            if response.status_code == 200:
                data = response.json()
                targets.extend(self._parse_annotations(data))
                
        return list(set(targets))

    def _parse_annotations(self, data):
        """Helper to parse PUG View JSON for gene symbols or protein names."""
        symbols = []
        try:
            # Navigate the nested JSON structure of PUG View
            # This is complex because PUG View is deeply nested
            section = data.get("Record", {}).get("Section", [])
            for s in section:
                for info in s.get("Information", []):
                    value = info.get("Value", {})
                    # Look for StringWithMarkup
                    for swm in value.get("StringWithMarkup", []):
                        string = swm.get("String")
                        if string:
                            # Heuristic: PubChem often lists gene symbols in parentheses or as standalone strings
                            # We'll collect strings and let the analyst filter them later
                            symbols.append(string)
        except Exception:
            pass
        return symbols

if __name__ == "__main__":
    client = PubChemClient()
    cid = client.get_cid_by_name("aspirin")
    print(f"CID: {cid}")
    if cid:
        targets = client.get_protein_targets(cid)
        print(f"Targets: {targets}")
