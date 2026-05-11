import requests
from typing import List, Dict, Optional

class OpenTargetsClient:
    URL = "https://api.platform.opentargets.org/api/v4/graphql"

    def __init__(self):
        self.session = requests.Session()

    def get_targets_by_chembl_id(self, chembl_id):
        """Fetches targets and action types for a ChEMBL ID using Open Targets GraphQL API."""
        query = """
        query drug($chemblId: String!) {
          drug(chemblId: $chemblId) {
            id
            name
            mechanismsOfAction {
              rows {
                actionType
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
                    action_type = row.get("actionType", "UNKNOWN")
                    for target in row.get("targets", []):
                        symbol = target.get("approvedSymbol")
                        if symbol:
                            targets.append({"symbol": symbol, "action_type": action_type})

        # Deduplicate
        unique_targets = {}
        for t in targets:
            sym = t["symbol"]
            if sym not in unique_targets:
                unique_targets[sym] = t["action_type"]

        return [{"symbol": sym, "action_type": act} for sym, act in unique_targets.items()]

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

    def get_drugs_by_target(self, symbol: str, limit: int = 5):
        """Fetches drugs that target a specific gene symbol."""
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

            query = """
            query target($ensemblId: String!) {
              target(ensemblId: $ensemblId) {
                drugAndClinicalCandidates {
                    rows {
                      drug {
                        id
                        name
                      }
                    }
                }
              }
            }
            """
            variables = {"ensemblId": ensembl_id}
            response = self.session.post(self.URL, json={"query": query, "variables": variables}, timeout=15)
            if response.status_code == 200:
                data = response.json()
                rows = data.get("data", {}).get("target", {}).get("drugAndClinicalCandidates", {}).get("rows", [])
                drugs = []
                seen = set()
                for row in rows:
                    drug_info = row.get("drug")
                    if drug_info:
                        name = drug_info.get("name")
                        if name and name not in seen:
                            drugs.append({"name": name, "id": drug_info.get("id")})
                            seen.add(name)
                            if len(drugs) >= limit:
                                break
                return drugs
        except Exception as e:
            print(f"Error fetching drugs for {symbol}: {e}")
        return []

    def get_adverse_events(self, chembl_id: str, limit: int = 10):
        """Fetches adverse events for a drug."""
        query = """
        query drug($chemblId: String!) {
          drug(chemblId: $chemblId) {
            adverseEvents {
              count
              rows {
                name
                count
                logLR
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
                rows = data.get("data", {}).get("drug", {}).get("adverseEvents", {}).get("rows", [])
                return rows[:limit]
        except Exception as e:
            print(f"Error fetching adverse events for {chembl_id}: {e}")
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
            r = self.session.post(self.URL, json={"query": search_query, "variables": {"symbol": symbol}}, timeout=30)
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

    def get_efo_id(self, disease_name: str) -> Optional[str]:
        """Resolves a disease name to an EFO ID."""
        query = """
        query search($queryString: String!) {
          search(queryString: $queryString, entityNames: ["disease"]) {
            hits {
              id
              name
            }
          }
        }
        """
        try:
            r = self.session.post(self.URL, json={"query": query, "variables": {"queryString": disease_name}}, timeout=30)
            hits = r.json().get("data", {}).get("search", {}).get("hits", [])
            if hits:
                return hits[0]["id"]
        except Exception as e:
            print(f"Error resolving disease {disease_name}: {e}")
        return None

    def get_targets_by_disease(self, efo_id: str, limit: int = 10) -> List[Dict]:
        """Fetches the top protein targets associated with a disease."""
        query = """
        query associatedTargets($efoId: String!, $size: Int!) {
          disease(efoId: $efoId) {
            associatedTargets(page: {index: 0, size: $size}) {
              rows {
                target {
                  id
                  approvedSymbol
                }
                score
              }
            }
          }
        }
        """
        try:
            # Increased timeout to 30s for large target lists
            r = self.session.post(self.URL, json={"query": query, "variables": {"efoId": efo_id, "size": limit}}, timeout=30)
            if r.status_code == 200:
                rows = r.json().get("data", {}).get("disease", {}).get("associatedTargets", {}).get("rows", [])
                targets = []
                for row in rows:
                    t = row.get("target")
                    if t:
                        targets.append({
                            "symbol": t.get("approvedSymbol"),
                            "id": t.get("id"),
                            "score": row.get("score")
                        })
                return targets
        except Exception as e:
            print(f"Error fetching targets for disease {efo_id}: {e}")
        return []

    def get_drugs_by_disease(self, efo_id: str, limit: int = 10) -> List[Dict]:
        """Fetches approved drugs associated with a disease EFO ID."""
        query = """
        query drugsForDisease($efoId: String!) {
          disease(efoId: $efoId) {
            id
            name
            drugAndClinicalCandidates {
              count
              rows {
                drug {
                  id
                  name
                }
              }
            }
          }
        }
        """
        try:
            # Increased timeout to 30s for large disease results
            r = self.session.post(self.URL, json={"query": query, "variables": {"efoId": efo_id}}, timeout=30)
            if r.status_code == 200:
                data = r.json()
                disease_data = data.get("data", {}).get("disease")
                if not disease_data: return []
                rows = disease_data.get("drugAndClinicalCandidates", {}).get("rows", [])
                
                drugs = []
                seen_drugs = set()
                for row in rows:
                    drug_info = row.get("drug")
                    if drug_info:
                        d_name = drug_info.get("name")
                        d_id = drug_info.get("id")
                        if d_name and d_id not in seen_drugs:
                            drugs.append({
                                "name": d_name,
                                "id": d_id
                            })
                            seen_drugs.add(d_id)
                            if len(drugs) >= limit:
                                break
                return drugs
        except Exception as e:
            print(f"Error fetching drugs for disease {efo_id}: {e}")
        return []

    def get_target_tractability(self, symbol: str) -> Dict:
        """Fetches druggability (tractability) data for a gene symbol."""
        query = """
        query target($symbol: String!) {
          search(queryString: $symbol, entityNames: ["target"]) {
            hits {
              id
              target {
                tractability {
                  id
                  modality
                  value
                }
              }
            }
          }
        }
        """
        try:
            r = self.session.post(self.URL, json={"query": query, "variables": {"symbol": symbol}}, timeout=30)
            hits = r.json().get("data", {}).get("search", {}).get("hits", [])
            if not hits: return {}
            
            tract = hits[0].get("target", {}).get("tractability", [])
            # Simplify to high-level modalities
            results = {"small_molecule": False, "antibody": False, "other": False}
            for t in tract:
                if "SM" in t['id'] or "Small Molecule" in t['modality']:
                    if t['value']: results["small_molecule"] = True
                elif "AB" in t['id'] or "Antibody" in t['modality']:
                    if t['value']: results["antibody"] = True
            return results
        except Exception as e:
            print(f"Error fetching tractability for {symbol}: {e}")
        return {}

    def get_target_association_details(self, symbol: str, efo_id: str) -> Dict:
        """Fetches evidence-specific association scores between a target and disease."""
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
            r = self.session.post(self.URL, json={"query": search_query, "variables": {"symbol": symbol}}, timeout=30)
            hits = r.json().get("data", {}).get("search", {}).get("hits", [])
            if not hits: return {"overall": 0.0, "literature": 0.0, "genetics": 0.0}
            ensembl_id = hits[0]["id"]

            query = """
            query association($ensemblId: String!, $efoId: String!) {
              target(ensemblId: $ensemblId) {
                associatedDiseases(page: {index: 0, size: 100}) {
                  rows {
                    disease {
                      id
                    }
                    score
                    datasourceScores {
                      id
                      score
                    }
                  }
                }
              }
            }
            """
            r = self.session.post(self.URL, json={"query": query, "variables": {"ensemblId": ensembl_id, "efoId": efo_id}}, timeout=30)
            rows = r.json().get("data", {}).get("target", {}).get("associatedDiseases", {}).get("rows", [])
            
            for row in rows:
                if row["disease"]["id"] == efo_id:
                    results = {"overall": row["score"], "literature": 0.0, "genetics": 0.0}
                    for ds in row.get("datasourceScores", []):
                        if ds["id"] == "europepmc": results["literature"] = ds["score"]
                        elif ds["id"] == "ot_genetics_portal": results["genetics"] = ds["score"]
                    return results
        except Exception as e:
            print(f"Error fetching association details for {symbol}-{efo_id}: {e}")
        return {"overall": 0.0, "literature": 0.0, "genetics": 0.0}

    def get_drugs_by_target(self, symbol: str, limit: int = 10) -> List[Dict]:
        """Fetches approved drugs that target a specific protein symbol."""
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
            r = self.session.post(self.URL, json={"query": search_query, "variables": {"symbol": symbol}}, timeout=30)
            hits = r.json().get("data", {}).get("search", {}).get("hits", [])
            if not hits: return []
            ensembl_id = hits[0]["id"]

            # 2. Fetch Drugs
            query = """
            query drugsByTarget($ensemblId: String!) {
              target(ensemblId: $ensemblId) {
                id
                approvedSymbol
                drugAndClinicalCandidates {
                  rows {
                    drug {
                      id
                      name
                      indications {
                        rows {
                          disease {
                            id
                            name
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            """
            r = self.session.post(self.URL, json={"query": query, "variables": {"ensemblId": ensembl_id}}, timeout=30)
            if r.status_code == 200:
                data = r.json()
                target_data = data.get("data", {}).get("target")
                if not target_data: return []
                rows = target_data.get("drugAndClinicalCandidates", {}).get("rows", [])
                
                drugs = []
                seen_drugs = set()
                for row in rows:
                    drug_info = row.get("drug")
                    if drug_info:
                        d_name = drug_info.get("name")
                        d_id = drug_info.get("id")
                        if d_name and d_id not in seen_drugs:
                            # Extract the ID strings and Name strings from indication objects
                            raw_inds = drug_info.get("indications", {}).get("rows", [])
                            indication_ids = [ind.get("disease", {}).get("id") for ind in raw_inds if ind.get("disease")]
                            indication_names = [ind.get("disease", {}).get("name") for ind in raw_inds if ind.get("disease")]
                            
                            drugs.append({
                                "name": d_name,
                                "id": d_id,
                                "indications": indication_ids,
                                "indication_names": indication_names
                            })
                            seen_drugs.add(d_id)
                            if len(drugs) >= limit:
                                break
                return drugs
        except Exception as e:
            print(f"Error fetching drugs for target {symbol}: {e}")
        return []

if __name__ == "__main__":
    client = OpenTargetsClient()
    # Disease: Alzheimer (EFO_0000249)
    efo = client.get_efo_id("Alzheimer")
    print(f"EFO for Alzheimer: {efo}")
    if efo:
        print(f"Drugs for Alzheimer: {client.get_drugs_by_disease(efo, limit=5)}")
    # Target: EGFR
    print(f"Drugs for EGFR: {client.get_drugs_by_target('EGFR', limit=5)}")
