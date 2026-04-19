from Bio import Entrez

class NCBIClient:
    def __init__(self, email="your.email@example.com"):
        Entrez.email = email

    def get_pathway_summary(self, pathway_name):
        """
        Attempts to find a summary for a pathway using NCBI Entrez.
        Search in 'pmc' or 'mesh' for general descriptions.
        """
        try:
            # Search PubMed Central for a general overview
            search_query = f"{pathway_name} pathway summary[Title/Abstract]"
            handle = Entrez.esearch(db="pmc", term=search_query, retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if record["IdList"]:
                pmcid = record["IdList"][0]
                fetch_handle = Entrez.efetch(db="pmc", id=pmcid, rettype="xml", retmode="text")
                # This returns full text which might be too much.
                # A better way is to search in Entrez Gene if the pathway is specific.
                # But for pathways, searching BioSystems or similar might be better.
                # NCBI BioSystems is deprecated.
                fetch_handle.close()
                return f"NCBI Summary available at: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmcid}"
            
            return "No specific NCBI summary found. Refer to KEGG description."
        except Exception as e:
            return f"Error fetching NCBI summary: {e}"

if __name__ == "__main__":
    client = NCBIClient()
    print(client.get_pathway_summary("Regulation of lipolysis in adipocytes"))
