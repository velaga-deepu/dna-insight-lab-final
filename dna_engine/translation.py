# DNA to RNA conversion

def dna_to_rna(sequence):
    return sequence.replace("T", "U")


# RNA codon table

CODON_TABLE = {
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L",

    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y",
    "UAA": "STOP", "UAG": "STOP",

    "UGU": "C", "UGC": "C",
    "UGA": "STOP", "UGG": "W",

    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",

    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M",

    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",

    "AGU": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


# RNA to protein conversion

def rna_to_protein(rna):

    protein = []

    for i in range(0, len(rna) - 2, 3):

        codon = rna[i:i + 3]

        amino_acid = CODON_TABLE.get(codon)

        if amino_acid == "STOP":
            break

        if amino_acid:
            protein.append(amino_acid)

    return "".join(protein)


# Amino acid full names

AMINO_ACID_NAMES = {
    "F": "Phenylalanine",
    "L": "Leucine",
    "S": "Serine",
    "Y": "Tyrosine",
    "C": "Cysteine",
    "W": "Tryptophan",
    "P": "Proline",
    "H": "Histidine",
    "Q": "Glutamine",
    "R": "Arginine",
    "I": "Isoleucine",
    "M": "Methionine",
    "T": "Threonine",
    "N": "Asparagine",
    "K": "Lysine",
    "V": "Valine",
    "A": "Alanine",
    "D": "Aspartic Acid",
    "E": "Glutamic Acid",
    "G": "Glycine",
    "STOP": "Stop Codon"
}


# Codon analysis

def analyze_codons(rna):

    codons = []

    for i in range(0, len(rna) - 2, 3):

        codon = rna[i:i + 3]

        amino_acid = CODON_TABLE.get(codon)

        if amino_acid:
            codons.append({
                "codon": codon,
                "amino_acid": amino_acid,
                "amino_acid_name": AMINO_ACID_NAMES[amino_acid]
            })

    return codons

# ORF detection

def find_orfs(sequence):

    stop_codons = ["TAA", "TAG", "TGA"]
    sequence = sequence.upper()
    orfs = []

    # Forward strand
    for frame in range(3):

        for i in range(frame, len(sequence) - 2, 3):

            codon = sequence[i:i + 3]

            if codon == "ATG":

                for j in range(i + 3, len(sequence) - 2, 3):

                    stop_codon = sequence[j:j + 3]

                    if stop_codon in stop_codons:

                        orfs.append({
                            "strand": "+",
                            "frame": frame + 1,
                            "start": i + 1,
                            "end": j + 3,
                            "sequence": sequence[i:j + 3],
                            "length": len(sequence[i:j + 3])
                        })

                        break

    # Reverse complement
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    reverse_sequence = "".join(
        complement[base] for base in sequence[::-1]
    )

    # Reverse strand
    for frame in range(3):

        for i in range(frame, len(reverse_sequence) - 2, 3):

            codon = reverse_sequence[i:i + 3]

            if codon == "ATG":

                for j in range(i + 3, len(reverse_sequence) - 2, 3):

                    stop_codon = reverse_sequence[j:j + 3]

                    if stop_codon in stop_codons:

                        orfs.append({
                            "strand": "-",
                            "frame": frame + 1,
                            "start": i + 1,
                            "end": j + 3,
                            "sequence": reverse_sequence[i:j + 3],
                            "length": len(reverse_sequence[i:j + 3])
                        })

                        break

    return orfs