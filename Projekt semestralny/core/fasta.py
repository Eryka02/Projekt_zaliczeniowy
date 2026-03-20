import os


def parse_fasta_file(path):

    with open(path) as f:
        lines = f.readlines()

    name = os.path.basename(path)

    sequence = "".join(
        [l.strip() for l in lines if not l.startswith(">")]
    )

    return name, sequence