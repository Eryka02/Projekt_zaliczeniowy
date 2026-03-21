from core.analysis import count_motif

class AnalysisService:

    def __init__(self, sequences):
        self.sequences = sequences

    def count_all(self, motifs):
        results = []

        for motif in motifs:
            row = {
                "motif": motif,
                "counts": []
            }

            for seq in self.sequences:
                count = count_motif(seq["sequence"], motif)
                row["counts"].append(count)

            results.append(row)

        return results