import requests

class OpenTargetsClient:
    URL = "https://api.platform.opentargets.org/api/v4/graphql"

    def __init__(self):
        self.session = requests.Session()

    def get_targets_by_chembl_id(self, chembl_id):
        """Fetches targets for a ChEMBL ID using Open Targets GraphQL API."""
        query = """
        query drug($chemblId: String!) {
          drug(chemblId: $chemblId) {
            id
            name
            mechanismsOfAction {
              rows {
                targetName
                targets {
                  id
                  approvedSymbol
                }
              }
            }
          }
        }
        """
        variables = {"chemblId": chembl_id}
        response = self.session.post(self.URL, json={"query": query, "variables": variables}, timeout=15)
        
        targets = []
        if response.status_code == 200:
            data = response.json()
            drug_data = data.get("data", {}).get("drug")
            if drug_data and drug_data.get("mechanismsOfAction"):
                for row in drug_data["mechanismsOfAction"].get("rows", []):
                    for target in row.get("targets", []):
                        symbol = target.get("approvedSymbol")
                        if symbol:
                            targets.append(symbol)
        return list(set(targets))

    def get_indications(self, chembl_id):
        """Fetches therapeutic indications for a ChEMBL ID."""
        query = """
        query drug($chemblId: String!) {
          drug(chemblId: $chemblId) {
            indications {
              rows {
                disease {
                  name
                }
              }
            }
          }
        }
        """
        variables = {"chemblId": chembl_id}
        try:
            response = self.session.post(self.URL, json={"query": query, "variables": variables}, timeout=15)
            if response.status_code == 200:
                data = response.json()
                rows = data.get("data", {}).get("drug", {}).get("indications", {}).get("rows", [])
                return [row["disease"]["name"] for row in rows if "disease" in row and "name" in row["disease"]]
        except Exception:
            pass
        return []

    def get_interactors(self, symbol: str, min_score: float = 0.9, limit: int = 5):
        """Fetches high-affinity interactors for a gene symbol."""
        # 1. Resolve Symbol to Ensembl ID
        search_query = """
        query search($symbol: String!) {
          search(queryString: $symbol, entityNames: ["target"]) {
            hits {
              id
            }
          }
        }
        """
        try:
            r = self.session.post(self.URL, json={"query": search_query, "variables": {"symbol": symbol}}, timeout=10)
            hits = r.json().get("data", {}).get("search", {}).get("hits", [])
            if not hits: return []
            ensembl_id = hits[0]["id"]
            
            # 2. Fetch Interactions
            query = """
            query targetInteractions($ensemblId: String!) {
              target(ensemblId: $ensemblId) {
                interactions {
                  rows {
                    targetB {
                      approvedSymbol
                    }
                    score
                  }
                }
              }
            }
            """
            variables = {"ensemblId": ensembl_id}
            response = self.session.post(self.URL, json={"query": query, "variables": variables}, timeout=15)
            if response.status_code == 200:
                data = response.json()
                rows = data.get("data", {}).get("target", {}).get("interactions", {}).get("rows", [])
                
                # Filter by score and unique symbols
                interactors = []
                seen = {symbol} # Don't include self
                for row in rows:
                    target_b = row.get("targetB")
                    score = row.get("score", 0)
                    if target_b and score >= min_score:
                        b_symbol = target_b.get("approvedSymbol")
                        if b_symbol and b_symbol not in seen:
                            interactors.append(b_symbol)
                            seen.add(b_symbol)
                            if len(interactors) >= limit:
                                break
                return interactors
        except Exception as e:
            print(f"Error fetching interactors for {symbol}: {e}")
        return []

if __name__ == "__main__":
    client = OpenTargetsClient()
    # Aspirin: CHEMBL25
    print(f"Targets for Aspirin (CHEMBL25): {client.get_targets_by_chembl_id('CHEMBL25')}")
    # Osimertinib: CHEMBL3353410
    print(f"Targets for Osimertinib (CHEMBL3353410): {client.get_targets_by_chembl_id('CHEMBL3353410')}")
