# DNA Insight Lab

A Flask web application that works as a virtual DNA laboratory. Analyse DNA
sequences, simulate mutations, align and compare sequences, explore genetics
concepts, and save your results as reports.

**3rd place, healthcare track, Electrolython 2026 (CMRIT).**

> For educational purposes only. This tool does not provide medical or
> clinical diagnosis.

## Features
- **DNA analyzer:** analyse a sequence and get its reverse complement, with
  input validation (A, T, G, C only) and a 10,000 bp length limit with clear
  error messages
- **FASTA upload:** analyse sequences straight from a FASTA file
- **Mutation analysis:** compare a reference and a sample sequence and see the
  effect of each variant
- **DNA matching:** Needleman-Wunsch global alignment (match +2, mismatch -1,
  gap -2) reporting similarity, matches, mismatches and gaps
- **Virtual DNA Lab:** hands-on experiments covering sequence analysis,
  DNA -> RNA -> protein translation with codon analysis and ORF finding,
  variant analysis, restriction site analysis, complement and reverse
  complement, and an interactive mutation simulator (substitution, insertion,
  deletion)
- **Genetics explorer:** educational modules for 8 conditions, including
  sickle cell anemia, cystic fibrosis, beta-thalassemia, PKU, Huntington's
  disease (CAG repeat comparison) and Duchenne muscular dystrophy, plus
  multifactorial conditions such as hypertension and type 2 diabetes. Results
  load without refreshing the page, so entered data is never lost
- **Saved reports:** save analyses to a SQLite database, reopen them later,
  and generate a combined report from any selection

## Built With
- Python, Flask
- SQLite
- HTML, CSS, JavaScript
- A custom `dna_engine` package for the bioinformatics logic

## Project Structure
```
dna-insight-lab-final/
│
├── app.py
├── requirements.txt
├── dna_engine/
│   ├── sequence_analysis.py
│   ├── restriction_analysis.py
│   ├── translation.py
│   ├── mutation_analysis.py
│   └── alignment.py
├── templates/
└── static/
```

## Installation & Running
```
git clone https://github.com/velaga-deepu/dna-insight-lab-final.git
cd dna-insight-lab-final
pip install -r requirements.txt
python app.py
```
Then open `http://127.0.0.1:8000` in your browser. The SQLite database is
created automatically on first run.

## Purpose
Built for Electrolython 2026 to make core genetics concepts such as
mutations, translation and sequence comparison interactive and easy to
explore, then improved using the judges' feedback.

## 🏆 Recognition
**3rd place, Healthcare track, Electrolython 2026 (CMRIT)**

Originally built for the hackathon. This repo is the final version, improved
using the judges' feedback.
