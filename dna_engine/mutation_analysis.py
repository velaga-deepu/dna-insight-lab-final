from dna_engine.translation import CODON_TABLE, AMINO_ACID_NAMES


def compare_sequences(reference, sample):

    reference = reference.upper()
    sample = sample.upper()

    differences = []

    if not reference or not sample:
        return {
            "error": True,
            "error_message": "Both DNA sequences are required."
        }

    valid_bases = "ATGC"

    if any(base not in valid_bases for base in reference):
        return {
            "error": True,
            "error_message":
                "Reference sequence can only contain A, T, G and C."
        }

    if any(base not in valid_bases for base in sample):
        return {
            "error": True,
            "error_message":
                "Sample sequence can only contain A, T, G and C."
        }

    if len(reference) != len(sample):
        return {
            "error": True,
            "error_message":
                "Reference and sample sequences must have the same length."
        }

    for i in range(len(reference)):

        if reference[i] != sample[i]:

            differences.append({
                "position": i + 1,
                "reference": reference[i],
                "sample": sample[i],
                "type": "Substitution",
                "change":
                    reference[i] + " → " + sample[i]
            })

    return {
        "error": False,
        "reference": reference,
        "sample": sample,
        "total_mutations": len(differences),
        "differences": differences
    }


def analyze_mutation_effect(reference, sample):

    result = compare_sequences(reference, sample)

    if result["error"]:
        return result

    mutation_details = []

    for mutation in result["differences"]:

        position = mutation["position"]

        codon_start = ((position - 1) // 3) * 3

        reference_codon = reference[
            codon_start:codon_start + 3
        ]

        sample_codon = sample[
            codon_start:codon_start + 3
        ]

        # Check whether both sequences contain a complete codon
        if len(reference_codon) == 3 and len(sample_codon) == 3:

            reference_rna_codon = reference_codon.replace(
                "T", "U"
            )

            sample_rna_codon = sample_codon.replace(
                "T", "U"
            )

            reference_amino_acid = CODON_TABLE.get(
                reference_rna_codon
            )

            sample_amino_acid = CODON_TABLE.get(
                sample_rna_codon
            )

            reference_amino_acid_name = AMINO_ACID_NAMES.get(
                reference_amino_acid,
                "Unknown"
            )

            sample_amino_acid_name = AMINO_ACID_NAMES.get(
                sample_amino_acid,
                "Unknown"
            )

            if reference_amino_acid == sample_amino_acid:

                effect = "Synonymous change"

            else:

                effect = "Amino-acid change"

            mutation["reference_codon"] = reference_codon

            mutation["sample_codon"] = sample_codon

            mutation["reference_amino_acid"] = (
                reference_amino_acid
            )

            mutation["sample_amino_acid"] = (
                sample_amino_acid
            )

            mutation["reference_amino_acid_name"] = (
                reference_amino_acid_name
            )

            mutation["sample_amino_acid_name"] = (
                sample_amino_acid_name
            )

            mutation["effect"] = effect

        else:

            mutation["reference_codon"] = (
                "⚠️ Codon effect cannot be determined"
            )

            mutation["sample_codon"] = (
                "⚠️ Codon effect cannot be determined"
            )

            mutation["reference_amino_acid"] = None

            mutation["sample_amino_acid"] = None

            mutation["reference_amino_acid_name"] = (
                "⚠️ Codon effect cannot be determined"
            )

            mutation["sample_amino_acid_name"] = (
                "⚠️ Codon effect cannot be determined"
            )

            mutation["effect"] = (
                "⚠️ Codon effect cannot be determined"
            )

        mutation_details.append(mutation)

    result["differences"] = mutation_details

    return result