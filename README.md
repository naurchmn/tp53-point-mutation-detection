# Deteksi Mutasi Titik (Point Mutation) pada Gen TP53 Kanker

## Deskripsi

Proyek ini mengimplementasikan **Smith-Waterman Algorithm** untuk mendeteksi mutasi titik (point mutation) pada gen **TP53** — gen penekan tumor yang paling sering bermutasi pada kanker manusia.

## Jenis Mutasi Titik yang Dideteksi

- **Substitusi (Missense/Nonsense)**: Pergantian satu basa nukleotida
- **Insersi**: Penyisipan satu basa nukleotida
- **Delesi**: Penghapusan satu basa nukleotida

## Struktur Proyek

```
tp53-point-mutation-detection/
├── data/                  # Data sekuens FASTA
│   ├── tp53_reference.fasta   # Sekuens referensi TP53 (wildtype)
│   └── tp53_mutant.fasta      # Sekuens mutan untuk analisis
├── notebooks/
│   └── tp53_mutation_detection.ipynb  # Notebook utama
├── src/
│   ├── __init__.py
│   ├── smith_waterman.py      # Implementasi algoritma Smith-Waterman
│   └── mutation_detector.py   # Logika deteksi mutasi
├── pyproject.toml
└── README.md
```

## Cara Menjalankan

```bash
# Install dependencies
uv sync

# Jalankan Jupyter Notebook
uv run jupyter notebook
```

## Tech Stack

- **Python 3.12**
- **BioPython** — parsing FASTA dan utilitas bioinformatika
- **NumPy** — operasi matriks scoring
- **Matplotlib** — visualisasi alignment dan heatmap
- **Jupyter** — notebook interaktif

## Referensi

- Smith, T.F. & Waterman, M.S. (1981). Identification of common molecular subsequences.
- NCBI TP53 Gene: https://www.ncbi.nlm.nih.gov/gene/7157
