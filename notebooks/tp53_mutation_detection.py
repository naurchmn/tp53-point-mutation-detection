import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from Bio import SeqIO

# Tambahkan path root project agar bisa import modul src
sys.path.insert(0, os.path.abspath('..'))

from src.smith_waterman import smith_waterman, AlignmentResult
from src.mutation_detector import detect_mutations, summarize_mutations, Mutation, MutationType

print("✅ Semua library berhasil diimpor!")

# Load sekuens referensi
ref_record = next(SeqIO.parse("../data/tp53_reference.fasta", "fasta"))
ref_seq = str(ref_record.seq)

print(f"📋 Referensi: {ref_record.id}")
print(f"   Panjang : {len(ref_seq)} bp")
print(f"   Preview : {ref_seq[:80]}...")
print()

# Load sekuens mutan
mutant_records = list(SeqIO.parse("../data/tp53_mutant.fasta", "fasta"))

print(f"🧬 Jumlah sampel mutan: {len(mutant_records)}")
for rec in mutant_records:
    print(f"   - {rec.id} ({len(rec.seq)} bp)")

# Demonstrasi dengan sekuens pendek
demo_ref   = "ACGTACGT"
demo_query = "ACGAACGT"  # Substitusi T→A di posisi 4

demo_result = smith_waterman(demo_ref, demo_query)

print("=" * 50)
print("DEMO: Smith-Waterman Local Alignment")
print("=" * 50)
print(f"Referensi : {demo_ref}")
print(f"Query     : {demo_query}")
print(f"\nSkor alignment  : {demo_result.score}")
print(f"Identity        : {demo_result.identity:.1f}%")
print(f"\nHasil Alignment:")
print(f"  Ref  : {demo_result.aligned_seq1}")
print(f"         {''.join('|' if a == b else 'X' for a, b in zip(demo_result.aligned_seq1, demo_result.aligned_seq2))}")
print(f"  Query: {demo_result.aligned_seq2}")

# Visualisasi matriks scoring untuk demo
fig, ax = plt.subplots(figsize=(10, 8))

H = demo_result.score_matrix
im = ax.imshow(H, cmap='YlOrRd', aspect='auto')

# Label sumbu
ax.set_xticks(range(len(demo_query) + 1))
ax.set_yticks(range(len(demo_ref) + 1))
ax.set_xticklabels(['-'] + list(demo_query))
ax.set_yticklabels(['-'] + list(demo_ref))
ax.set_xlabel('Query Sequence', fontsize=12)
ax.set_ylabel('Reference Sequence', fontsize=12)
ax.set_title('Smith-Waterman Scoring Matrix (Demo)', fontsize=14, fontweight='bold')

# Tampilkan nilai di setiap sel
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        ax.text(j, i, str(H[i, j]), ha='center', va='center',
                fontsize=11, fontweight='bold',
                color='white' if H[i, j] > H.max() * 0.6 else 'black')

plt.colorbar(im, label='Alignment Score')
plt.tight_layout()
plt.show()

# Jalankan alignment dan deteksi mutasi untuk setiap sampel
results = {}

for record in mutant_records:
    mut_seq = str(record.seq)
    
    # Smith-Waterman alignment
    alignment = smith_waterman(ref_seq, mut_seq)
    
    # Deteksi mutasi
    mutations = detect_mutations(alignment)
    summary = summarize_mutations(mutations)
    
    results[record.id] = {
        'alignment': alignment,
        'mutations': mutations,
        'summary': summary,
    }
    
    print("=" * 60)
    print(f"🧬 Sampel: {record.id}")
    print(f"   Deskripsi: {record.description}")
    print(f"   Skor Alignment : {alignment.score}")
    print(f"   Identity       : {alignment.identity:.1f}%")
    print(f"   Total Mutasi   : {summary['total']}")
    print(f"     - Substitusi : {summary['substitusi']}")
    print(f"     - Insersi    : {summary['insersi']}")
    print(f"     - Delesi     : {summary['delesi']}")
    print()
    
    if mutations:
        print("   Mutasi yang terdeteksi:")
        for mut in mutations:
            print(f"     {mut}")
    print()

# Bar chart perbandingan skor alignment
sample_names = [name.replace('TP53_Mutant_Sample_', '') for name in results.keys()]
scores = [r['alignment'].score for r in results.values()]
identities = [r['alignment'].identity for r in results.values()]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Skor alignment
colors = plt.cm.RdYlGn(np.array(identities) / 100)
bars1 = ax1.bar(sample_names, scores, color=colors, edgecolor='black', linewidth=0.8)
ax1.set_title('Skor Alignment Smith-Waterman', fontsize=13, fontweight='bold')
ax1.set_ylabel('Skor')
ax1.set_xlabel('Sampel')
ax1.tick_params(axis='x', rotation=45)
for bar, score in zip(bars1, scores):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
             str(score), ha='center', va='bottom', fontweight='bold')

# Identity percentage
bars2 = ax2.bar(sample_names, identities, color=colors, edgecolor='black', linewidth=0.8)
ax2.set_title('Persentase Identity', fontsize=13, fontweight='bold')
ax2.set_ylabel('Identity (%)')
ax2.set_xlabel('Sampel')
ax2.set_ylim(0, 105)
ax2.tick_params(axis='x', rotation=45)
for bar, ident in zip(bars2, identities):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f'{ident:.1f}%', ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
plt.show()

# Stacked bar chart jenis mutasi per sampel
fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(sample_names))
width = 0.5

subs = [r['summary']['substitusi'] for r in results.values()]
ins = [r['summary']['insersi'] for r in results.values()]
dels = [r['summary']['delesi'] for r in results.values()]

ax.bar(x, subs, width, label='Substitusi', color='#e74c3c', edgecolor='black')
ax.bar(x, ins, width, bottom=subs, label='Insersi', color='#3498db', edgecolor='black')
ax.bar(x, dels, width, bottom=[s+i for s, i in zip(subs, ins)], label='Delesi', color='#2ecc71', edgecolor='black')

ax.set_xticks(x)
ax.set_xticklabels(sample_names, rotation=45, ha='right')
ax.set_ylabel('Jumlah Mutasi')
ax.set_title('Distribusi Jenis Mutasi Titik per Sampel', fontsize=13, fontweight='bold')
ax.legend()
ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

plt.tight_layout()
plt.show()

# Heatmap matriks scoring untuk sampel pertama (cropped region)
first_sample = list(results.keys())[0]
H_full = results[first_sample]['alignment'].score_matrix

# Crop area sekitar alignment (terlalu besar untuk visualisasi penuh)
start_i = max(0, results[first_sample]['alignment'].start_pos[0] - 10)
end_i = min(H_full.shape[0], results[first_sample]['alignment'].end_pos[0] + 10)
start_j = max(0, results[first_sample]['alignment'].start_pos[1] - 10)
end_j = min(H_full.shape[1], results[first_sample]['alignment'].end_pos[1] + 10)

H_crop = H_full[start_i:end_i, start_j:end_j]

fig, ax = plt.subplots(figsize=(12, 10))
im = ax.imshow(H_crop, cmap='inferno', aspect='auto', interpolation='nearest')
ax.set_title(f'Scoring Matrix Heatmap — {first_sample}\n(Cropped region sekitar alignment)',
             fontsize=13, fontweight='bold')
ax.set_xlabel('Query Position')
ax.set_ylabel('Reference Position')
plt.colorbar(im, label='Alignment Score')
plt.tight_layout()
plt.show()

print(f"Ukuran matriks penuh: {H_full.shape}")
print(f"Region ditampilkan: [{start_i}:{end_i}, {start_j}:{end_j}]")
print(f"Skor maksimum: {H_full.max()}")

def visualize_alignment(alignment: AlignmentResult, sample_name: str, window: int = 80):
    """
    Menampilkan alignment dalam format yang mudah dibaca,
    mirip output BLAST.
    """
    ref_al = alignment.aligned_seq1
    query_al = alignment.aligned_seq2
    
    print(f"\n{'='*70}")
    print(f"Alignment Detail: {sample_name}")
    print(f"Score: {alignment.score} | Identity: {alignment.identity:.1f}%")
    print(f"{'='*70}")
    
    for start in range(0, len(ref_al), window):
        end = min(start + window, len(ref_al))
        ref_chunk = ref_al[start:end]
        query_chunk = query_al[start:end]
        
        # Baris pencocokan
        match_line = ''.join(
            '|' if r == q else ('.' if r == '-' or q == '-' else 'X')
            for r, q in zip(ref_chunk, query_chunk)
        )
        
        print(f"\nRef   {start+1:>5d}  {ref_chunk}  {end}")
        print(f"             {match_line}")
        print(f"Query {start+1:>5d}  {query_chunk}  {end}")

# Tampilkan alignment untuk setiap sampel
for sample_name, result in results.items():
    visualize_alignment(result['alignment'], sample_name)

print("\n" + "=" * 70)
print("RINGKASAN ANALISIS MUTASI TITIK GEN TP53")
print("=" * 70)
print(f"\n📊 Jumlah sampel dianalisis: {len(results)}")
print(f"📏 Panjang sekuens referensi: {len(ref_seq)} bp")
print(f"⚙️  Algoritma: Smith-Waterman (Local Alignment)")
print(f"   Parameter: Match=+2, Mismatch=-1, Gap=-1")

print("\n" + "-" * 70)
print(f"{'Sampel':<30} {'Skor':>8} {'Identity':>10} {'Mutasi':>8}")
print("-" * 70)

for name, result in results.items():
    short_name = name.replace('TP53_Mutant_Sample_', '')
    al = result['alignment']
    s = result['summary']
    print(f"{short_name:<30} {al.score:>8} {al.identity:>9.1f}% {s['total']:>8}")

print("-" * 70)

total_mutations = sum(r['summary']['total'] for r in results.values())
total_subs = sum(r['summary']['substitusi'] for r in results.values())
total_ins = sum(r['summary']['insersi'] for r in results.values())
total_dels = sum(r['summary']['delesi'] for r in results.values())

print(f"\n📈 Total mutasi terdeteksi: {total_mutations}")
print(f"   Substitusi : {total_subs}")
print(f"   Insersi    : {total_ins}")
print(f"   Delesi     : {total_dels}")
print("\n" + "=" * 70)