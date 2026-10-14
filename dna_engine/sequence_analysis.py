def clean_sequence(sequence):
    return sequence.upper()


def analyze_sequence(sequence):

    sequence = clean_sequence(sequence)

    # Empty input
    if not sequence:
        return {
            "error": True,
            "error_message": "Please enter a DNA sequence."
        }

    # Check for spaces
    if " " in sequence:
        return {
            "error": True,
            "error_message": "DNA sequence can't contain spaces."
        }

    # Check for numbers
    if any(character.isdigit() for character in sequence):
        return {
            "error": True,
            "error_message": "DNA sequence can't contain numbers."
        }

    # Check for special characters
    if any(not character.isalpha() for character in sequence):
        return {
            "error": True,
            "error_message": "DNA sequence can't contain special characters."
        }

    # Check for letters other than A, T, G and C
    if any(character not in "ATGC" for character in sequence):
        return {
            "error": True,
            "error_message": "DNA sequence can only contain A, T, G and C."
        }

    # ---------------------------------------------------------
    # Basic DNA calculations
    # ---------------------------------------------------------

    length = len(sequence)

    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")

    # GC and AT percentages
    gc_content = ((g_count + c_count) / length) * 100
    at_content = ((a_count + t_count) / length) * 100

    # ---------------------------------------------------------
    # Return analysis result
    # ---------------------------------------------------------

    return {
        "error": False,

        # Original sequence
        "sequence": sequence,

        # Length
        "length": length,

        # Individual nucleotide counts
        "A": a_count,
        "T": t_count,
        "G": g_count,
        "C": c_count,

        # GC / AT percentages
        "GC": gc_content,
        "AT": at_content,

        # Dictionary used by the Virtual DNA Lab template
        "counts": {
            "A": a_count,
            "T": t_count,
            "G": g_count,
            "C": c_count
        },

        # Names used by the Virtual DNA Lab template
        "gc_content": round(gc_content, 2),
        "at_content": round(at_content, 2),

        # Compatibility aliases
        "gc_percentage": round(gc_content, 2),
        "at_percentage": round(at_content, 2)
    }


def reverse_complement(sequence):

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    reverse_sequence = sequence[::-1]

    return "".join(
        complement[base]
        for base in reverse_sequence
    )