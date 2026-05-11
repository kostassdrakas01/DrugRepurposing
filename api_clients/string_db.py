import requests
from typing import List, Dict

class StringDBClient:
    """Client for fetching protein-protein interactions from String-DB."""
    
    URL = "https://string-db.org/api/json/network"
    
    def __init__(self):
        self.session = requests.Session()

    def get_interactors(self, symbol: str, limit: int = 5, required_score: int = 700) -> List[Dict]:
        """
        Fetches high-confidence interactors for a gene symbol.
        required_score: 400 (medium), 700 (high), 900 (highest)
        """
        params = {
            "identifiers": symbol,
            "species": 9606,  # Human
            "limit": limit,
            "required_score": required_score,
            "caller_identity": "KEGG-ID_Discovery_Engine"
        }
        
        try:
            r = self.session.get(self.URL, params=params, timeout=10)
            if r.status_code == 200:
                data = r.json()
                interactors = []
                seen = {symbol.upper()}
                
                for interaction in data:
                    # String-DB returns pairs (p1, p2)
                    p1 = interaction.get("preferredName_A", "").upper()
                    p2 = interaction.get("preferredName_B", "").upper()
                    score = interaction.get("score")
                    
                    # We want the 'other' protein in the pair
                    other = p2 if p1 == symbol.upper() else p1
                    
                    if other and other not in seen:
                        interactors.append({
                            "symbol": other,
                            "score": score
                        })
                        seen.add(other)
                        
                return interactors
        except Exception as e:
            print(f"Error fetching String-DB interactors for {symbol}: {e}")
        return []
