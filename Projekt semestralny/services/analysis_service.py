from core.analysis import get_positions_single

class AnalysisService:
    def __init__(self, sequences):
        self.sequences = sequences

    def analyze_all(self, motifs):
        result = {}

        for seq in self.sequences:
            seq_name = seq["name"]
            seq_text = seq["sequence"].upper()
            result[seq_name] = {}

            for motif in motifs:
                motif = motif.upper()
                positions = get_positions_single(seq_text, motif).tolist()
                result[seq_name][motif] = {
                    "count": len(positions),
                    "positions": positions
                }

        return result