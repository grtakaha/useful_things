def translate(transcript):
    """
    Takes in a nucleotide (DNA - T, not U) transcript and returns its protein translation.

        Parameters:
            transcript (str): Transcript in DNA nucleotides (T, not U).

        Returns:
            protein (str): Protein translation of transcript.
    """

    # Standard codon table adapted from Wikipedia.
    # Saved in https://github.com/grtakaha/useful_things.
    # Stop codons represented by "-".
    codons = {
        "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
        "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "TAT":"Y", "TAC":"Y", "TAA":"-", "TAG":"-",
        "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "TGT":"C", "TGC":"C", "TGA":"-", "TGG":"W",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
        "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
        }    
    
    protein = ""
    caps_transcript = transcript.upper()
    i = 0
    while i < len(transcript):
        codon = caps_transcript[i:i+3]
        # Won't output anything for the end of the sequence if it's not a multiple of 3.
        if len(codon) == 3:
            aa = codons.get(codon)
            # If this is not a recognized codon, then use "_" as a placeholder.
            if not aa:
                aa = "_"
            protein += aa
        i += 3

    return protein