import requests
import pandas as pd
from typing import Dict, Optional, List

class GTExClient:
    BASE_URL = "https://gtexportal.org/api/v2"

    def __init__(self):
        self.session = requests.Session()
        self.gene_id_cache = {}

    def get_gencode_id(self, gene_symbol: str) -> Optional[str]:
        """Resolves a gene symbol to a Gencode ID."""
        if gene_symbol in self.gene_id_cache:
            return self.gene_id_cache[gene_symbol]

        try:
            params = {"geneId": gene_symbol}
            response = self.session.get(f"{self.BASE_URL}/reference/gene", params=params, timeout=15)
            if response.status_code == 200:
                data = response.json().get("data", [])
                if data:
                    gencode_id = data[0].get("gencodeId")
                    self.gene_id_cache[gene_symbol] = gencode_id
                    return gencode_id
        except Exception:
            pass
        return None

    def get_expression(self, gene_symbol: str, tissue: str) -> float:
        """Returns the median expression (TPM) for a gene in a specific tissue."""
        gencode_id = self.get_gencode_id(gene_symbol)
        if not gencode_id:
            return 0.0

        try:
            params = {
                "gencodeId": gencode_id,
                "tissueSiteDetailId": tissue,
                "datasetId": "gtex_v8"
            }
            response = self.session.get(f"{self.BASE_URL}/expression/medianGeneExpression", params=params, timeout=15)
            if response.status_code == 200:
                data = response.json().get("data", [])
                if data:
                    return data[0].get("median", 0.0)
        except Exception:
            pass
        return 0.0

    def get_tissues(self) -> List[str]:
        """Returns a list of available tissue IDs."""
        try:
            response = self.session.get(f"{self.BASE_URL}/dataset/tissueSiteDetail", timeout=15)
            if response.status_code == 200:
                return [t["tissueSiteDetailId"] for t in response.json().get("data", [])]
        except Exception:
            pass
        return []
