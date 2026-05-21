"""
Modul Deteksi Mutasi Titik (Point Mutation) pada Gen TP53.

Mendeteksi tiga jenis mutasi titik:
1. Substitusi — pergantian satu basa (missense/nonsense/silent)
2. Insersi  — penyisipan satu basa
3. Delesi   — penghapusan satu basa
"""

from dataclasses import dataclass
from enum import Enum
from .smith_waterman import AlignmentResult


class MutationType(Enum):
    """Jenis mutasi titik."""
    SUBSTITUTION = "Substitusi"
    INSERTION = "Insersi"
    DELETION = "Delesi"


@dataclass
class Mutation:
    """Representasi satu mutasi yang terdeteksi."""
    mutation_type: MutationType
    position: int           # Posisi pada sekuens referensi (1-indexed)
    ref_base: str           # Basa pada referensi
    mut_base: str           # Basa pada mutan
    context: str            # Konteks sekuens di sekitar mutasi

    def __str__(self) -> str:
        if self.mutation_type == MutationType.SUBSTITUTION:
            return (
                f"[{self.mutation_type.value}] Posisi {self.position}: "
                f"{self.ref_base} → {self.mut_base} | Konteks: ...{self.context}..."
            )
        elif self.mutation_type == MutationType.INSERTION:
            return (
                f"[{self.mutation_type.value}] Posisi {self.position}: "
                f"Insersi '{self.mut_base}' | Konteks: ...{self.context}..."
            )
        else:  # DELETION
            return (
                f"[{self.mutation_type.value}] Posisi {self.position}: "
                f"Delesi '{self.ref_base}' | Konteks: ...{self.context}..."
            )


def detect_mutations(
    alignment: AlignmentResult,
    context_size: int = 5,
) -> list[Mutation]:
    """
    Mendeteksi mutasi titik dari hasil alignment Smith-Waterman.

    Parameters
    ----------
    alignment : AlignmentResult
        Hasil alignment dari fungsi smith_waterman().
    context_size : int
        Jumlah basa di kiri-kanan mutasi untuk konteks (default: 5).

    Returns
    -------
    list[Mutation]
        Daftar mutasi yang terdeteksi.
    """
    mutations = []
    ref_aligned = alignment.aligned_seq1
    query_aligned = alignment.aligned_seq2

    # Posisi pada sekuens referensi (untuk pelaporan posisi asli)
    ref_pos = alignment.start_pos[0]  # 0-indexed start

    for i in range(len(ref_aligned)):
        ref_base = ref_aligned[i]
        query_base = query_aligned[i]

        if ref_base != '-':
            ref_pos += 1  # Increment posisi referensi

        if ref_base == query_base:
            continue  # Tidak ada mutasi

        # Ambil konteks sekuens
        ctx_start = max(0, i - context_size)
        ctx_end = min(len(ref_aligned), i + context_size + 1)
        context = ref_aligned[ctx_start:ctx_end]

        if ref_base == '-':
            # Gap di referensi = Insersi pada query
            mutations.append(Mutation(
                mutation_type=MutationType.INSERTION,
                position=ref_pos,
                ref_base='-',
                mut_base=query_base,
                context=context,
            ))
        elif query_base == '-':
            # Gap di query = Delesi pada query
            mutations.append(Mutation(
                mutation_type=MutationType.DELETION,
                position=ref_pos,
                ref_base=ref_base,
                mut_base='-',
                context=context,
            ))
        else:
            # Substitusi (basa berbeda)
            mutations.append(Mutation(
                mutation_type=MutationType.SUBSTITUTION,
                position=ref_pos,
                ref_base=ref_base,
                mut_base=query_base,
                context=context,
            ))

    return mutations


def summarize_mutations(mutations: list[Mutation]) -> dict:
    """
    Membuat ringkasan statistik dari mutasi yang terdeteksi.

    Returns
    -------
    dict
        Ringkasan berisi jumlah tiap jenis mutasi dan total.
    """
    summary = {
        "total": len(mutations),
        "substitusi": sum(1 for m in mutations if m.mutation_type == MutationType.SUBSTITUTION),
        "insersi": sum(1 for m in mutations if m.mutation_type == MutationType.INSERTION),
        "delesi": sum(1 for m in mutations if m.mutation_type == MutationType.DELETION),
    }
    return summary
