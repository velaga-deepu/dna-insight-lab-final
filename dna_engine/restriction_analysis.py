RESTRICTION_ENZYMES = {
    "EcoRI": {
        "site": "GAATTC"
    },
    "BamHI": {
        "site": "GGATCC"
    },
    "HindIII": {
        "site": "AAGCTT"
    },
    "HaeIII": {
        "site": "GGCC"
    }
}


def find_restriction_sites(sequence, enzyme_name):

    sequence = sequence.strip().upper()

    if not sequence:
        return {
            "error": True,
            "error_message": "Please enter a DNA sequence."
        }

    for base in sequence:
        if base not in "ATGC":
            return {
                "error": True,
                "error_message":
                    "Invalid DNA sequence. Use only A, T, G and C."
            }

    if enzyme_name not in RESTRICTION_ENZYMES:
        return {
            "error": True,
            "error_message":
                "Please select a valid restriction enzyme."
        }

    recognition_site = RESTRICTION_ENZYMES[enzyme_name]["site"]

    positions = []
    cut_offsets = {
        "EcoRI": 1,
        "BamHI": 1,
        "HindIII": 1,
        "HaeIII": 2
    }

    cut_offset = cut_offsets[enzyme_name]

    for i in range(
        len(sequence) - len(recognition_site) + 1
    ):

        if sequence[
            i:i + len(recognition_site)
        ] == recognition_site:

            positions.append(i + 1)

    # Calculate fragment sizes using actual enzyme cut positions
    cut_positions = [
    position + cut_offset - 1 for position in positions
]
    fragment_sizes = []

    if cut_positions:
        fragment_sizes.append(cut_positions[0])

        for i in range(1, len(cut_positions)):
            fragment_sizes.append(
                cut_positions[i] - cut_positions[i - 1]
            )

        fragment_sizes.append(
            len(sequence) - cut_positions[-1]
        )
    else:
        fragment_sizes.append(len(sequence))

    return {
        "error": False,
        "sequence": sequence,
        "enzyme": enzyme_name,
        "recognition_site": recognition_site,
        "positions": positions,
        "site_count": len(positions),
        "fragment_sizes": fragment_sizes
    }