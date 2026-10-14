# =========================================================
# DNA GLOBAL SEQUENCE ALIGNMENT
# Needleman-Wunsch Algorithm
# =========================================================


def needleman_wunsch(
    sequence1,
    sequence2,
    match_score=2,
    mismatch_score=-1,
    gap_penalty=-2
):
    """
    Perform global alignment using the Needleman-Wunsch algorithm.

    Parameters:
        sequence1       First DNA sequence
        sequence2       Second DNA sequence
        match_score     Score for matching bases
        mismatch_score  Score for different bases
        gap_penalty     Penalty for insertion/deletion

    Returns:
        alignment information including:
        - aligned sequences
        - alignment score
        - matches
        - mismatches
        - gaps
        - differences
        - similarity
    """

    sequence1 = sequence1.upper()
    sequence2 = sequence2.upper()

    rows = len(sequence1)
    columns = len(sequence2)

    # -----------------------------------------------------
    # Create scoring matrix
    # -----------------------------------------------------

    matrix = [
        [0 for _ in range(columns + 1)]
        for _ in range(rows + 1)
    ]

    # -----------------------------------------------------
    # Initialize first column
    # -----------------------------------------------------

    for i in range(1, rows + 1):

        matrix[i][0] = (
            i * gap_penalty
        )

    # -----------------------------------------------------
    # Initialize first row
    # -----------------------------------------------------

    for j in range(1, columns + 1):

        matrix[0][j] = (
            j * gap_penalty
        )

    # -----------------------------------------------------
    # Fill scoring matrix
    # -----------------------------------------------------

    for i in range(1, rows + 1):

        for j in range(1, columns + 1):

            if sequence1[i - 1] == sequence2[j - 1]:

                diagonal_score = (
                    matrix[i - 1][j - 1]
                    + match_score
                )

            else:

                diagonal_score = (
                    matrix[i - 1][j - 1]
                    + mismatch_score
                )

            up_score = (
                matrix[i - 1][j]
                + gap_penalty
            )

            left_score = (
                matrix[i][j - 1]
                + gap_penalty
            )

            matrix[i][j] = max(
                diagonal_score,
                up_score,
                left_score
            )

    # -----------------------------------------------------
    # Traceback
    # -----------------------------------------------------

    aligned1 = []
    aligned2 = []

    i = rows
    j = columns

    while i > 0 or j > 0:

        # Diagonal move
        if i > 0 and j > 0:

            if sequence1[i - 1] == sequence2[j - 1]:

                diagonal_score = (
                    matrix[i - 1][j - 1]
                    + match_score
                )

            else:

                diagonal_score = (
                    matrix[i - 1][j - 1]
                    + mismatch_score
                )

            if matrix[i][j] == diagonal_score:

                aligned1.append(
                    sequence1[i - 1]
                )

                aligned2.append(
                    sequence2[j - 1]
                )

                i -= 1
                j -= 1

                continue

        # Gap in sequence 2
        if i > 0:

            up_score = (
                matrix[i - 1][j]
                + gap_penalty
            )

            if matrix[i][j] == up_score:

                aligned1.append(
                    sequence1[i - 1]
                )

                aligned2.append("-")

                i -= 1

                continue

        # Gap in sequence 1
        if j > 0:

            left_score = (
                matrix[i][j - 1]
                + gap_penalty
            )

            if matrix[i][j] == left_score:

                aligned1.append("-")

                aligned2.append(
                    sequence2[j - 1]
                )

                j -= 1

                continue

    # Reverse traceback results
    aligned1.reverse()
    aligned2.reverse()

    aligned1 = "".join(aligned1)
    aligned2 = "".join(aligned2)

    # -----------------------------------------------------
    # Analyze alignment
    # -----------------------------------------------------

    matches = 0
    mismatches = 0
    gaps = 0

    differences = []

    for position in range(
        len(aligned1)
    ):

        base1 = aligned1[position]
        base2 = aligned2[position]

        alignment_position = position + 1

        # Match
        if (
            base1 != "-"
            and base2 != "-"
            and base1 == base2
        ):

            matches += 1

        # Gap
        elif base1 == "-" or base2 == "-":

            gaps += 1

            differences.append({

                "position":
                    alignment_position,

                "sequence1":
                    base1,

                "sequence2":
                    base2,

                "type":
                    "Insertion/Deletion"

            })

        # Mismatch
        else:

            mismatches += 1

            differences.append({

                "position":
                    alignment_position,

                "sequence1":
                    base1,

                "sequence2":
                    base2,

                "type":
                    "Base mismatch"

            })

    alignment_length = len(aligned1)

    # -----------------------------------------------------
    # Similarity
    # -----------------------------------------------------

    if alignment_length > 0:

        similarity = (
            matches /
            alignment_length
        ) * 100

    else:

        similarity = 0

    # -----------------------------------------------------
    # Status
    # -----------------------------------------------------

    if similarity == 100:

        status = "Exact Match"

    elif similarity >= 90:

        status = "Very High Similarity"

    elif similarity >= 70:

        status = "High Similarity"

    elif similarity >= 50:

        status = "Moderate Similarity"

    else:

        status = "Low Similarity"

    # -----------------------------------------------------
    # Alignment visualization
    # -----------------------------------------------------

    alignment_line = []

    for position in range(
        alignment_length
    ):

        if (
            aligned1[position] != "-"
            and aligned2[position] != "-"
            and aligned1[position]
            == aligned2[position]
        ):

            alignment_line.append("|")

        elif (
            aligned1[position] == "-"
            or aligned2[position] == "-"
        ):

            alignment_line.append(" ")

        else:

            alignment_line.append(".")

    alignment_line = "".join(
        alignment_line
    )

    return {

        "aligned_sequence1":
            aligned1,

        "aligned_sequence2":
            aligned2,

        "alignment_line":
            alignment_line,

        "alignment_score":
            matrix[rows][columns],

        "alignment_length":
            alignment_length,

        "matches":
            matches,

        "mismatches":
            mismatches,

        "gaps":
            gaps,

        "similarity":
            similarity,

        "status":
            status,

        "differences":
            differences

    }