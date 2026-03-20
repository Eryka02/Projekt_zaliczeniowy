import ssl
from Bio import Entrez

ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "erykasworczuk993@gmail.com"


def fetch_sequence_from_ncbi(query):

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        return None

    seq_id = record["IdList"][0]

    handle = Entrez.efetch(
        db="nucleotide",
        id=seq_id,
        rettype="fasta",
        retmode="text"
    )

    fasta_data = handle.read()
    handle.close()

    lines = fasta_data.splitlines()

    sequence = "".join(
        [l.strip() for l in lines if not l.startswith(">")]
    )

    return {
        "query": query,
        "id": seq_id,
        "sequence": sequence
    }