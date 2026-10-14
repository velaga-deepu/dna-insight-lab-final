from flask import Flask, render_template, request, redirect, url_for, jsonify
import json
import sqlite3
from datetime import datetime

from dna_engine.sequence_analysis import (
    analyze_sequence,
    reverse_complement
)

from dna_engine.restriction_analysis import (
    find_restriction_sites
)

from dna_engine.translation import (
    dna_to_rna,
    rna_to_protein,
    analyze_codons,
    find_orfs
)

from dna_engine.mutation_analysis import (
    compare_sequences,
    analyze_mutation_effect
)

from dna_engine.alignment import (
    needleman_wunsch
)


# =========================================================
# FLASK APP
# =========================================================

app = Flask(__name__)

DATABASE = "dna_reports.db"
MAX_DNA_LENGTH = 10000


def sequence_length_error(sequence):
    if len(sequence) > MAX_DNA_LENGTH:
        return {
            "error": True,
            "error_title": "⚠️ Sequence too long",
            "error_message": (
                f"Your sequence contains {len(sequence):,} bp.\n"
                f"The maximum supported length is {MAX_DNA_LENGTH:,} bp.\n\n"
                "Please provide a shorter sequence."
            )
        }

    return None


# =========================================================
# DATABASE SETUP
# =========================================================

def init_database():

    connection = sqlite3.connect(DATABASE)

    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS saved_analyses (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            module TEXT NOT NULL,
            title TEXT NOT NULL,
            details TEXT NOT NULL,
            created_at TEXT NOT NULL
        )
        """
    )

    connection.commit()
    connection.close()


init_database()


# =========================================================
# DASHBOARD
# =========================================================

@app.route("/dashboard")
def dashboard():

    return render_template("dashboard.html")


# =========================================================
# HOME / DNA ANALYZER
# =========================================================

@app.route("/", methods=["GET", "POST"])
def home():

    result = None

    if request.method == "POST":

        sequence = request.form.get(
            "sequence",
            ""
        ).strip().upper()

        length_error = sequence_length_error(sequence)

        if length_error:
            return render_template(
                "index.html",
                error=length_error["error_message"],
                error_title=length_error["error_title"]
            )

        if not sequence:

            result = {
                "error": True,
                "error_message":
                    "Invalid input! DNA sequence is required."
            }

        elif any(base not in "ATGC" for base in sequence):

            result = {
                "error": True,
                "error_message":
                    "Invalid input! DNA sequence can contain only A, T, G and C."
            }

        else:

            result = analyze_sequence(sequence)

            result["reverse_complement"] = reverse_complement(
                sequence
            )

    return render_template(
        "index.html",
        result=result
    )
 

# =========================================================
# FASTA UPLOAD
# =========================================================

@app.route("/fasta/upload", methods=["POST"])
def fasta_upload():

    file = request.files.get("fasta_file")

    if not file or file.filename == "":

        return render_template(
            "index.html",
            result={
                "error": True,
                "error_message":
                    "Invalid input! Please select a FASTA file."
            }
        )

    try:

        content = file.read().decode("utf-8")

    except UnicodeDecodeError:

        return render_template(
            "index.html",
            result={
                "error": True,
                "error_message":
                    "Invalid input! The FASTA file could not be read."
            }
        )

    lines = content.splitlines()

    sequence_parts = []

    for line in lines:

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            continue

        sequence_parts.append(line)

    sequence = "".join(sequence_parts).upper()

    if not sequence:

        return render_template(
            "index.html",
            result={
                "error": True,
                "error_message":
                    "Invalid input! The FASTA file does not contain a DNA sequence."
            }
        )

    if any(base not in "ATGC" for base in sequence):

        return render_template(
            "index.html",
            result={
                "error": True,
                "error_message":
                    "Invalid input! The FASTA sequence can contain only A, T, G and C."
            }
        )

    result = analyze_sequence(sequence)

    result["reverse_complement"] = reverse_complement(
        sequence
    )

    return render_template(
        "index.html",
        result=result
    )


# =========================================================
# MUTATION ANALYSIS
# =========================================================

@app.route("/mutation", methods=["GET", "POST"])
def mutation():

    mutation_result = None

    if request.method == "POST":

        reference = request.form.get("reference","" ).strip().upper()

        sample = request.form.get("sample","").strip().upper()

        reference_length_error = sequence_length_error(reference)
        sample_length_error = sequence_length_error(sample)

        if reference_length_error:
         return render_template(
        "mutation.html",
        mutation_result=reference_length_error
    )

        if sample_length_error:
         return render_template(
        "mutation.html",
        mutation_result=sample_length_error
    )

        if not reference or not sample:

            mutation_result = {
                "error": True,
                "error_message":
                    "Invalid input! Both reference and sample DNA sequences are required."
            }

        elif any(base not in "ATGC" for base in reference):

            mutation_result = {
                "error": True,
                "error_message":
                    "Invalid input! Reference DNA can contain only A, T, G and C."
            }

        elif any(base not in "ATGC" for base in sample):

            mutation_result = {
                "error": True,
                "error_message":
                    "Invalid input! Sample DNA can contain only A, T, G and C."
            }

        else:

            mutation_result = analyze_mutation_effect(
                reference,
                sample
            )

    return render_template(
        "mutation.html",
        mutation_result=mutation_result
    )


# =========================================================
# DNA MATCHING
# NEEDLEMAN-WUNSCH GLOBAL ALIGNMENT
# =========================================================

@app.route("/matching", methods=["GET", "POST"])
def matching():

    matching_result = None

    if request.method == "POST":

        sequence1 = request.form.get(
            "sequence1",
            ""
        ).strip().upper()

        sequence2 = request.form.get(
            "sequence2",
            ""
        ).strip().upper()

        sequence1_length_error = sequence_length_error(sequence1)
        sequence2_length_error = sequence_length_error(sequence2)

        if sequence1_length_error:
            matching_result = sequence1_length_error

        elif sequence2_length_error:
            matching_result = sequence2_length_error

        elif not sequence1 or not sequence2:

            matching_result = {
                "error": True,
                "error_message":
                    "Invalid input! Both DNA sequences are required."
            }

        elif any(base not in "ATGC" for base in sequence1):

            matching_result = {
                "error": True,
                "error_message":
                    "Invalid input! Sequence 1 can contain only A, T, G and C."
            }

        elif any(base not in "ATGC" for base in sequence2):

            matching_result = {
                "error": True,
                "error_message":
                    "Invalid input! Sequence 2 can contain only A, T, G and C."
            }

        else:

            alignment = needleman_wunsch(
                sequence1,
                sequence2,
                match_score=2,
                mismatch_score=-1,
                gap_penalty=-2
            )

            matching_result = {

                "error": False,

                "sequence1": sequence1,
                "sequence2": sequence2,

                "length1": len(sequence1),
                "length2": len(sequence2),

                "aligned_sequence1":
                    alignment["aligned_sequence1"],

                "aligned_sequence2":
                    alignment["aligned_sequence2"],

                "alignment_line":
                    alignment["alignment_line"],

                "alignment_score":
                    alignment["alignment_score"],

                "alignment_length":
                    alignment["alignment_length"],

                "matching_count":
                    alignment["matches"],

                "difference_count":
                    alignment["mismatches"] +
                    alignment["gaps"],

                "matches":
                    alignment["matches"],

                "mismatches":
                    alignment["mismatches"],

                "gaps":
                    alignment["gaps"],

                "similarity":
                    alignment["similarity"],

                "status":
                    alignment["status"],

                "differences":
                    alignment["differences"]
            }

    return render_template(
        "matching.html",
        matching_result=matching_result
    )


# =========================================================
# VIRTUAL DNA LAB
# =========================================================
@app.route("/lab")
def lab():

    return render_template(
        "lab.html",
        mutation_result=None
    )

# =========================================================
# LAB - EXPERIMENT 2
# DNA SEQUENCE ANALYSIS
# =========================================================

@app.route("/lab/analyze", methods=["POST"])
def lab_analyze():

    sequence = request.form.get("sequence","").strip().upper()

    length_error = sequence_length_error(sequence)

    if length_error:
     return render_template(
        "lab.html",
        result=length_error
    )
    if not sequence:

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence is required."
        }

    elif any(base not in "ATGC" for base in sequence):

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence can contain only A, T, G and C."
        }

    else:

        result = analyze_sequence(sequence)

        result["reverse_complement"] = reverse_complement(
            sequence
        )
    return render_template(
    "lab.html",
    lab_result=result,
    mutation_result=None
)


# =========================================================
# LAB - EXPERIMENT 3
# DNA → RNA → PROTEIN
# =========================================================

@app.route("/lab/translation", methods=["POST"])
def lab_translation():

    sequence = request.form.get("sequence","").strip().upper()
    length_error = sequence_length_error(sequence)

    if length_error:
     return render_template(
        "lab.html",
        translation_result=length_error
    )
    if not sequence:

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence is required."
        }

    elif any(base not in "ATGC" for base in sequence):

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence can contain only A, T, G and C."
        }

    else:

        rna = dna_to_rna(sequence)

        protein = rna_to_protein(
            rna
        )

        codons = analyze_codons(
            rna
        )

        orfs = find_orfs(
            sequence
        )

        result = {

            "error": False,

            "dna": sequence,

            "rna": rna,

            "protein": protein,

            "codons": codons,

            "orfs": orfs
        }

    return render_template(
        "lab.html",
        translation_result=result
    )


# =========================================================
# LAB - EXPERIMENT 6
# INTERACTIVE TRANSLATION LAB
#
# IMPORTANT:
# This route is an additional alias for Experiment 6.
# Therefore BOTH URLs work:
#
# /lab/translation
# /lab/experiment6
# =========================================================

@app.route("/lab/experiment6", methods=["POST"])
def lab_experiment6():

    sequence = request.form.get("sequence","").strip().upper()
    length_error = sequence_length_error(sequence)

    if length_error:
     return render_template(
        "lab.html",
        translation_result=length_error
    )

    if not sequence:

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence is required."
        }

    elif any(base not in "ATGC" for base in sequence):

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence can contain only A, T, G and C."
        }

    else:

        rna = dna_to_rna(
            sequence
        )

        protein = rna_to_protein(
            rna
        )

        codons = analyze_codons(
            rna
        )

        orfs = find_orfs(
            sequence
        )

        result = {

            "error": False,

            "dna": sequence,

            "rna": rna,

            "protein": protein,

            "codons": codons,

            "orfs": orfs
        }

    return render_template(
        "lab.html",
        experiment6_result=result
    )


# =========================================================
# LAB - EXPERIMENT 7
# COMPLEMENT / REVERSE COMPLEMENT
# =========================================================

@app.route("/lab/complement", methods=["POST"])
def lab_complement():

    sequence = request.form.get("sequence","").strip().upper()
    length_error = sequence_length_error(sequence)

    if length_error:
     return render_template(
        "lab.html",
        complement_result=length_error
    )
    if not sequence:

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence is required."
        }

    elif any(base not in "ATGC" for base in sequence):

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence can contain only A, T, G and C."
        }

    else:

        complement = ""

        mapping = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G"
        }

        for base in sequence:

            complement += mapping[base]

        result = {

            "error": False,

            "sequence": sequence,

            "complement": complement,

            "reverse_complement":
                complement[::-1]
        }

    return render_template(
        "lab.html",
        complement_result=result
    )


# =========================================================
# LAB - EXPERIMENT 5
# RESTRICTION SITE ANALYSIS
# =========================================================

@app.route("/lab/restriction", methods=["POST"])
def lab_restriction():

    sequence = request.form.get("sequence","").strip().upper()
    length_error = sequence_length_error(sequence)

    if length_error:
     return render_template(
        "lab.html",
        restriction_result=length_error
    )
    enzyme = request.form.get(
        "enzyme",
        ""
    )

    if not sequence:

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence is required."
        }

    elif any(base not in "ATGC" for base in sequence):

        result = {
            "error": True,
            "error_message":
                "Invalid input! DNA sequence can contain only A, T, G and C."
        }

    else:

        result = find_restriction_sites(
            sequence,
            enzyme
        )

    return render_template(
        "lab.html",
        restriction_result=result
    )


# =========================================================
# LAB - EXPERIMENT 4
# VARIANT / MUTATION ANALYSIS
# =========================================================

@app.route("/lab/variant", methods=["POST"])
def lab_variant():

    reference = request.form.get("reference","").strip().upper()
    sample = request.form.get("sample","").strip().upper()
    reference_length_error = sequence_length_error(reference)

    if reference_length_error:
        return render_template(
            "lab.html",
            variant_result=reference_length_error
        )

    sample_length_error = sequence_length_error(sample)

    if sample_length_error:
        return render_template(
            "lab.html",
            variant_result=sample_length_error
        )

    if not reference or not sample:

        result = {
            "error": True,
            "error_message":
                "Invalid input! Both reference and sample DNA sequences are required."
        }

    elif any(base not in "ATGC" for base in reference):

        result = {
            "error": True,
            "error_message":
                "Invalid input! Reference DNA can contain only A, T, G and C."
        }

    elif any(base not in "ATGC" for base in sample):

        result = {
            "error": True,
            "error_message":
                "Invalid input! Sample DNA can contain only A, T, G and C."
        }

    else:

        result = analyze_mutation_effect(
            reference,
            sample
        )

    return render_template(
        "lab.html",
        variant_result=result
    )

@app.route("/lab/mutation", methods=["POST"])
def lab_mutation():
    sequence = request.form.get("mutationSequence", "").strip().upper()
    mutation_type = request.form.get("mutationType", "").strip().lower()
    position_text = request.form.get("mutationPosition", "").strip()
    new_base = request.form.get("mutationBase", "").strip().upper()

    length_error = sequence_length_error(sequence)

    if length_error:
        return render_template(
            "lab.html",
            mutation_result={
                "error": True,
                "title": length_error["error_title"],
                "message": length_error["error_message"]
            }
        )
    # Validate DNA sequence
    if not sequence:
        return render_template(
            "lab.html",
            mutation_result={
                "error": True,
                "message": "Invalid input! DNA sequence cannot be empty."
            }
        )

    if not all(base in "ATGC" for base in sequence):
        return render_template(
            "lab.html",
            mutation_result={
                "error": True,
                "message": "Invalid input! DNA sequence can contain only A, T, G and C."
            }
        )

    if mutation_type not in ["substitution", "insertion", "deletion"]:
        return render_template(
            "lab.html",
            mutation_result={
                "error": True,
                "message": "Invalid input! Please select a mutation type."
            }
        )

    # Validate position
    try:
        position = int(position_text)
    except ValueError:
        return render_template(
            "lab.html",
            mutation_result={
                "error": True,
                "message": "Invalid input! Position must be a number."
            }
        )

    if position < 1 or position > len(sequence):
        return render_template(
            "lab.html",
            mutation_result={
                "error": True,
                "message": f"Invalid input! Position must be between 1 and {len(sequence)}."
            }
        )

    # SUBSTITUTION
    if mutation_type == "substitution":

        if new_base not in "ATGC":
            return render_template(
                "lab.html",
                mutation_result={
                    "error": True,
                    "message": "Invalid input! Replacement base must be A, T, G or C."
                }
            )

        old_base = sequence[position - 1]

        mutated_sequence = (
            sequence[:position - 1]
            + new_base
            + sequence[position:]
        )

        result = {
            "error": False,
            "mutation_type": "Substitution",
            "original_sequence": sequence,
            "mutated_sequence": mutated_sequence,
            "position": position,
            "original_base": old_base,
            "new_base": new_base,
            "description": f"{old_base} → {new_base}"
        }

    # INSERTION
    elif mutation_type == "insertion":

        if new_base not in "ATGC":
            return render_template(
                "lab.html",
                mutation_result={
                    "error": True,
                    "message": "Invalid input! Inserted base must be A, T, G or C."
                }
            )

        mutated_sequence = (
            sequence[:position - 1]
            + new_base
            + sequence[position - 1:]
        )

        result = {
            "error": False,
            "mutation_type": "Insertion",
            "original_sequence": sequence,
            "mutated_sequence": mutated_sequence,
            "position": position,
            "original_base": "-",
            "new_base": new_base,
            "description": f"Inserted {new_base}"
        }

    # DELETION
    else:

        deleted_base = sequence[position - 1]

        mutated_sequence = (
            sequence[:position - 1]
            + sequence[position:]
        )

        result = {
            "error": False,
            "mutation_type": "Deletion",
            "original_sequence": sequence,
            "mutated_sequence": mutated_sequence,
            "position": position,
            "original_base": deleted_base,
            "new_base": "-",
            "description": f"Deleted {deleted_base}"
        }

    return render_template(
        "lab.html",
        mutation_result=result
    )


# =========================================================
# GENETICS DATA
# =========================================================

def get_conditions():

    return [

        {
            "id": "sickle_cell",
            "name": "Sickle Cell Anemia",
            "gene": "HBB",
            "description":
                "An inherited blood disorder associated with changes in the HBB gene."
        },

        {
            "id": "cystic_fibrosis",
            "name": "Cystic Fibrosis",
            "gene": "CFTR",
            "description":
                "An inherited condition associated with changes in the CFTR gene."
        },

        {
            "id": "beta_thalassemia",
            "name": "Beta-Thalassemia",
            "gene": "HBB",
            "description":
                "An inherited blood condition associated with variants in the HBB gene."
        },

        {
            "id": "pku",
            "name": "Phenylketonuria (PKU)",
            "gene": "PAH",
            "description":
                "An inherited metabolic condition associated with changes in the PAH gene."
        },

        {
            "id": "huntington",
            "name": "Huntington's Disease",
            "gene": "HTT",
            "description":
                "An inherited neurological condition associated with CAG repeat expansion in the HTT gene."
        },

        {
            "id": "duchenne",
            "name": "Duchenne Muscular Dystrophy",
            "gene": "DMD",
            "description":
                "An inherited muscle condition associated with changes in the DMD gene."
        },

        {
            "id": "hypertension",
            "name": "High Blood Pressure",
            "gene": "Multiple genes",
            "description":
                "A multifactorial condition influenced by genetic, lifestyle and environmental factors."
        },

        {
            "id": "type2_diabetes",
            "name": "Type 2 Diabetes",
            "gene": "Multiple genes",
            "description":
                "A multifactorial condition influenced by genetic, lifestyle and environmental factors."
        }

    ]


# =========================================================
# GENETICS
# =========================================================

@app.route("/genetics")
def genetics():

    conditions = get_conditions()

    return render_template(
        "genetics.html",
        conditions=conditions
    )


# =========================================================
# GENETICS ANALYSIS
# =========================================================

@app.route("/genetics/analyze", methods=["POST"])
def genetics_analyze():

    condition_id = request.form.get("condition_id", "").strip()

    # Backward compatibility with the older form field name.
    if not condition_id:
        condition_id = request.form.get("condition", "").strip()

    reference = request.form.get("reference", "").strip().upper()
    sample = request.form.get("sample", "").strip().upper()

    reference_length_error = sequence_length_error(reference)
    sample_length_error = sequence_length_error(sample)

    conditions = get_conditions()
    condition = next(
        (item for item in conditions if item["id"] == condition_id),
        None
    )

    if condition is None:
        analysis_result = {
            "error": True,
            "error_message": "Invalid input! Please select a genetic condition."
        }

    elif reference_length_error:
        analysis_result = reference_length_error

    elif sample_length_error:
        analysis_result = sample_length_error

    elif condition_id in {"hypertension", "type2_diabetes"}:
        analysis_result = {
            "error": False,
            "condition_id": condition_id,
            "condition": condition["name"],
            "special_type": "multifactorial",
            "message": "Educational genetics result",
            "explanation": condition["description"]
        }

    elif not reference or not sample:
        analysis_result = {
            "error": True,
            "condition_id": condition_id,
            "error_message": "Invalid input! Both reference and sample DNA sequences are required."
        }

    elif any(base not in "ATGC" for base in reference):
        analysis_result = {
            "error": True,
            "condition_id": condition_id,
            "error_message": "Invalid input! Reference DNA can contain only A, T, G and C."
        }

    elif any(base not in "ATGC" for base in sample):
        analysis_result = {
            "error": True,
            "condition_id": condition_id,
            "error_message": "Invalid input! Sample DNA can contain only A, T, G and C."
        }

    elif condition_id == "huntington":
        reference_repeats = reference.count("CAG")
        sample_repeats = sample.count("CAG")
        difference = sample_repeats - reference_repeats
        analysis_result = {
            "error": False,
            "condition_id": condition_id,
            "condition": condition["name"],
            "special_type": "repeat",
            "reference_repeats": reference_repeats,
            "sample_repeats": sample_repeats,
            "repeat_difference": difference,
            "result": "CAG repeat comparison completed for educational exploration."
        }

    else:
        differences = []
        max_length = max(len(reference), len(sample))

        for index in range(max_length):
            ref_base = reference[index] if index < len(reference) else "-"
            sample_base = sample[index] if index < len(sample) else "-"

            if ref_base != sample_base:
                if ref_base == "-":
                    mutation_type = "Insertion"
                elif sample_base == "-":
                    mutation_type = "Deletion"
                else:
                    mutation_type = "Substitution"

                differences.append({
                    "position": index + 1,
                    "reference": ref_base,
                    "sample": sample_base,
                    "type": mutation_type,
                    "change": f"{ref_base} → {sample_base}"
                })

        analysis_result = {
            "error": False,
            "condition_id": condition_id,
            "condition": condition["name"],
            "total_mutations": len(differences),
            "differences": differences
        }

    # Genetics analysis is consumed by the page with fetch(), so the page
    # itself is never reloaded and unsaved inputs remain in the browser.
    analysis_result["gene"] = condition.get("gene", "") if condition else ""
    return jsonify(analysis_result)


# =========================================================
# SAVE ANALYSIS REPORT
# =========================================================

@app.route("/report/save", methods=["POST"])
def save_analysis_report():

    module = request.form.get(
        "module",
        ""
    ).strip()

    title = request.form.get(
        "title",
        ""
    ).strip()

    details_text = request.form.get(
        "details",
        ""
    )

    if not module or not title or not details_text:

        return (
            "Invalid report data. Please try again.",
            400
        )

    if len(details_text) > 500000:

        return (
            "Report data is too large.",
            413
        )

    try:

        details = json.loads(
            details_text
        )

        if not isinstance(details, dict):

            return (
                "Invalid report format.",
                400
            )

    except (
        json.JSONDecodeError,
        TypeError
    ):

        return (
            "Could not read the analysis results.",
            400
        )

    created_at = datetime.now().strftime(
        "%d %B %Y, %I:%M %p"
    )

    connection = sqlite3.connect(
        DATABASE
    )

    try:

        connection.execute(
            """
            INSERT INTO saved_analyses
                (module, title, details, created_at)
            VALUES (?, ?, ?, ?)
            """,
            (
                module[:100],
                title[:200],
                json.dumps(details),
                created_at
            )
        )

        connection.commit()

    finally:

        connection.close()

    return redirect(
        url_for("saved_reports")
    )


# =========================================================
# SAVED REPORTS LIST
# =========================================================

@app.route("/report/saved")
def saved_reports():

    connection = sqlite3.connect(
        DATABASE
    )

    connection.row_factory = sqlite3.Row

    try:

        rows = connection.execute(
            """
            SELECT id, module, title, created_at
            FROM saved_analyses
            ORDER BY id DESC
            """
        ).fetchall()

        saved_analyses = [
            dict(row)
            for row in rows
        ]

    finally:

        connection.close()

    return render_template(
        "saved_reports.html",
        saved_analyses=saved_analyses
    )


# =========================================================
# VIEW ONE SAVED ANALYSIS
# =========================================================
@app.route("/report/view/<int:analysis_id>")
def view_saved_analysis(analysis_id):

    connection = sqlite3.connect(DATABASE)
    connection.row_factory = sqlite3.Row

    try:
        row = connection.execute(
            """
            SELECT id, module, title, details, created_at
            FROM saved_analyses
            WHERE id = ?
            """,
            (analysis_id,)
        ).fetchone()
    finally:
        connection.close()

    if row is None:
        return ("Saved analysis not found.", 404)

    try:
        details = json.loads(row["details"])
    except (json.JSONDecodeError, TypeError):
        return ("Saved analysis data could not be read.", 500)

    analysis = {
        "id": row["id"],
        "module": row["module"],
        "title": row["title"],
        "details": details,
        "created_at": row["created_at"]
    }

    module = row["module"].strip().lower()

    if module == "dna matching":
        return render_template(
            "saved_matching.html",
            analysis=analysis
        )

    elif module == "mutation analysis":
        return render_template(
            "saved_mutation.html",
            analysis=analysis
        )

    else:
        return render_template(
            "combined_report.html",
            analyses=[analysis],
            generated_at=row["created_at"]
        )
# =========================================================
# COMBINED REPORT
# =========================================================

@app.route("/report/generate", methods=["POST"])
def generate_combined_report():

    selected_ids = request.form.getlist(
        "analysis_ids"
    )

    if not selected_ids:

        return (
            "Please select at least one analysis.",
            400
        )

    try:

        selected_ids = [
            int(value)
            for value in selected_ids
        ]

    except ValueError:

        return (
            "Invalid analysis selection.",
            400
        )

    placeholders = ",".join(
        "?"
        for _ in selected_ids
    )

    connection = sqlite3.connect(
        DATABASE
    )

    connection.row_factory = sqlite3.Row

    try:

        rows = connection.execute(
            f"""
            SELECT id, module, title, details, created_at
            FROM saved_analyses
            WHERE id IN ({placeholders})
            ORDER BY id DESC
            """,
            selected_ids
        ).fetchall()

    finally:

        connection.close()

    analyses = []

    for row in rows:

        try:

            details = json.loads(
                row["details"]
            )

        except (
            json.JSONDecodeError,
            TypeError
        ):

            details = {}

        analyses.append(
            {
                "id": row["id"],
                "module": row["module"],
                "title": row["title"],
                "details": details,
                "created_at": row["created_at"]
            }
        )

    return render_template(
        "combined_report.html",
        analyses=analyses,
        generated_at=datetime.now().strftime(
            "%d %B %Y, %I:%M %p"
        )
    )


# =========================================================
# RUN APPLICATION
# =========================================================
@app.route("/helix")
def helix():
    return render_template("helix.html")
if __name__ == "__main__":

    app.run(
        host="127.0.0.1",
        port=8000,
        debug=True
    )