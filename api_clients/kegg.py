import requests
import io
from Bio.KEGG.KGML import KGML_parser
from typing import List, Dict, Optional

class KEGGClient:
    BASE_URL = "https://rest.kegg.jp"

    def __init__(self):
        self.session = requests.Session()
        self.cache = {
            "cpd:C00001": "H2O",
            "cpd:C00002": "ATP",
            "cpd:C00008": "ADP"
        }
        self.kgml_cache = {}
        self.gene_pathways_cache = {}

    def _get(self, endpoint: str, operation: str, argument: str) -> str:
        """Helper for KEGG REST API calls with timeouts."""
        url = f"{self.BASE_URL}/{operation}/{argument}"
        try:
            response = self.session.get(url, timeout=10)
            if response.status_code == 200:
                return response.text
        except Exception:
            pass
        return ""

    def get_pathways_for_gene(self, kegg_gene_id: str) -> List[str]:
        """Finds pathways associated with a KEGG Gene ID."""
        if kegg_gene_id in self.gene_pathways_cache:
            return self.gene_pathways_cache[kegg_gene_id]

        result = self._get("pathway", "link", f"pathway/{kegg_gene_id}")
        pathways = []
        for line in result.strip().split("\n"):
            if "\t" in line:
                pathway_id = line.split("\t")[1].replace("path:", "")
                pathways.append(pathway_id)

        self.gene_pathways_cache[kegg_gene_id] = pathways
        return pathways

    def get_pathway_category(self, pathway_id: str) -> str:
        """Maps pathway ID to a high-level KEGG category based on official numbering."""
        pid_num = "".join(filter(str.isdigit, pathway_id))
        if not pid_num: return "Unknown"
        
        # Ensure we have a 5-digit number
        if len(pid_num) < 5:
            pid_num = pid_num.zfill(5)
            
        val = int(pid_num)
        
        if val < 2000: return "Metabolism"
        if val < 4000: return "Genetic Information"
        if val < 5000:
            if 4720 <= val <= 4730: return "Nervous System"
            return "Signal Transduction"
        if val < 6000:
            if 5010 <= val <= 5022: return "Neurodegenerative Diseases"
            return "Human Diseases"
        return "Biological Systems"

    def get_entity_name(self, entity_id):
        """Resolves technical IDs with caching and batch support."""
        if not entity_id: return "Unknown"
        
        if isinstance(entity_id, list):
            to_fetch = [eid for eid in entity_id if eid not in self.cache]
            if to_fetch:
                # Batch 10 at a time
                for i in range(0, len(to_fetch), 10):
                    chunk = to_fetch[i:i+10]
                    raw = self._get("get", "get", "+".join(chunk))
                    records = raw.split("///")
                    for record in records:
                        e_found = None
                        for line in record.split("\n"):
                            if line.startswith("ENTRY"):
                                parts = line.split()
                                if len(parts) > 1:
                                    raw_id = parts[1]
                                    for c_id in chunk:
                                        if raw_id in c_id:
                                            e_found = c_id
                                            break
                            if line.startswith("NAME") and e_found:
                                name = line.replace("NAME", "").strip().split(";")[0]
                                self.cache[e_found] = name
                                break
            return [self.cache.get(eid, eid) for eid in entity_id]

        if entity_id in self.cache: return self.cache[entity_id]
        
        if ":" not in entity_id: return entity_id
        result = self._get("get", "get", entity_id)
        name = entity_id
        for line in result.split("\n"):
            if line.startswith("NAME"):
                name = line.replace("NAME", "").strip().split(";")[0]
                break
        self.cache[entity_id] = name
        return name

    def get_pathway_info(self, pathway_id: str) -> Dict:
        """Fetches metadata for a pathway."""
        result = self._get("get", "get", pathway_id)
        info = {"compounds": [], "genes": []}
        for line in result.split("\n"):
            if line.startswith("NAME"):
                info["name"] = line.replace("NAME", "").strip()
            if line.startswith("DESCRIPTION"):
                info["description"] = line.replace("DESCRIPTION", "").strip()
        return info

    def get_kgml_pathway(self, pathway_id: str):
        """Fetches and parses KGML for a pathway (with caching)."""
        if pathway_id in self.kgml_cache: return self.kgml_cache[pathway_id]
        
        url = f"{self.BASE_URL}/get/{pathway_id}/kgml"
        try:
            response = self.session.get(url, timeout=10)
            if response.status_code == 200:
                pathway = KGML_parser.read(response.text)
                self.kgml_cache[pathway_id] = pathway
                return pathway
        except Exception:
            pass
        return None

    def get_map_url(self, pathway_id: str, highlighted_genes: List[str]) -> str:
        pid = pathway_id.replace("path:", "")
        base = f"https://www.kegg.jp/kegg-bin/show_pathway?{pid}"
        if highlighted_genes:
            return f"{base}+{'+'.join(highlighted_genes)}"
        return base
