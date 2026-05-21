"""
Implementasi Smith-Waterman Algorithm untuk Local Sequence Alignment.

Algoritma ini digunakan untuk menemukan alignment lokal optimal antara
dua sekuens DNA, yang berguna untuk mendeteksi mutasi titik pada gen TP53.

Referensi:
    Smith, T.F. & Waterman, M.S. (1981).
    Identification of common molecular subsequences.
    Journal of Molecular Biology, 147(1), 195-197.
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class AlignmentResult:
    """Hasil alignment dari Smith-Waterman Algorithm."""
    aligned_seq1: str        # Sekuens 1 yang sudah di-align
    aligned_seq2: str        # Sekuens 2 yang sudah di-align
    score: float             # Skor alignment optimal
    score_matrix: np.ndarray # Matriks scoring H
    start_pos: tuple         # Posisi awal alignment (i, j) pada sekuens asli
    end_pos: tuple           # Posisi akhir alignment (i, j) pada sekuens asli
    identity: float          # Persentase kecocokan


def smith_waterman(
    seq1: str,
    seq2: str,
    match_score: int = 2,
    mismatch_penalty: int = -1,
    gap_penalty: int = -1,
) -> AlignmentResult:
    """
    Menjalankan Smith-Waterman local alignment antara dua sekuens.

    Parameters
    ----------
    seq1 : str
        Sekuens referensi (misalnya TP53 wildtype).
    seq2 : str
        Sekuens query (misalnya sampel pasien/mutan).
    match_score : int
        Skor untuk match (default: 2).
    mismatch_penalty : int
        Penalti untuk mismatch (default: -1).
    gap_penalty : int
        Penalti untuk gap/indel (default: -1).

    Returns
    -------
    AlignmentResult
        Objek berisi hasil alignment lengkap.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    m = len(seq1)
    n = len(seq2)

    # ============================================================
    # LANGKAH 1: Inisialisasi matriks scoring H (ukuran (m+1) x (n+1))
    # Baris pertama dan kolom pertama diisi 0
    # ============================================================
    H = np.zeros((m + 1, n + 1), dtype=int)

    # Matriks traceback: 0=stop, 1=diagonal, 2=atas, 3=kiri
    traceback = np.zeros((m + 1, n + 1), dtype=int)

    # ============================================================
    # LANGKAH 2: Pengisian matriks scoring (Matrix Filling)
    # H(i,j) = max(0, H(i-1,j-1)+s(a_i,b_j), H(i-1,j)+gap, H(i,j-1)+gap)
    # ============================================================
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Hitung skor untuk match/mismatch (diagonal)
            if seq1[i - 1] == seq2[j - 1]:
                diag = H[i - 1, j - 1] + match_score
            else:
                diag = H[i - 1, j - 1] + mismatch_penalty

            # Hitung skor untuk gap di seq2 (dari atas)
            up = H[i - 1, j] + gap_penalty

            # Hitung skor untuk gap di seq1 (dari kiri)
            left = H[i, j - 1] + gap_penalty

            # Ambil nilai maksimum (minimal 0 karena local alignment)
            H[i, j] = max(0, diag, up, left)

            # Simpan arah traceback
            if H[i, j] == 0:
                traceback[i, j] = 0  # Stop
            elif H[i, j] == diag:
                traceback[i, j] = 1  # Diagonal (match/mismatch)
            elif H[i, j] == up:
                traceback[i, j] = 2  # Atas (gap di seq2)
            else:
                traceback[i, j] = 3  # Kiri (gap di seq1)

            # Track posisi skor maksimum
            if H[i, j] > max_score:
                max_score = H[i, j]
                max_pos = (i, j)

    # ============================================================
    # LANGKAH 3: Traceback — rekonstruksi alignment optimal
    # Mulai dari posisi skor tertinggi, ikuti traceback sampai 0
    # ============================================================
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_pos
    end_pos = (i - 1, j - 1)  # Konversi ke 0-indexed

    while traceback[i, j] != 0 and i > 0 and j > 0:
        if traceback[i, j] == 1:  # Diagonal
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i, j] == 2:  # Atas (gap di seq2)
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:  # Kiri (gap di seq1)
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    start_pos = (i, j)  # Posisi awal alignment (0-indexed)

    # Balik urutan karena traceback berjalan mundur
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    # Hitung identity (persentase kecocokan)
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    identity = (matches / len(aligned_seq1) * 100) if aligned_seq1 else 0.0

    return AlignmentResult(
        aligned_seq1=aligned_seq1,
        aligned_seq2=aligned_seq2,
        score=max_score,
        score_matrix=H,
        start_pos=start_pos,
        end_pos=end_pos,
        identity=identity,
    )
