import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from src.smith_waterman import smith_waterman
from src.mutation_detector import detect_mutations, summarize_mutations

def main():
    parser = argparse.ArgumentParser(description="Deteksi Mutasi Titik pada Gen TP53")
    parser.add_argument("--ref", type=str, required=True, help="Path ke file FASTA referensi")
    parser.add_argument("--query", type=str, required=True, help="Path ke file FASTA sampel/query")
    parser.add_argument("--context", type=int, default=5, help="Ukuran konteks basa di sekitar mutasi")
    
    args = parser.parse_args()
    
    if not Path(args.ref).exists():
        print(f"Error: File referensi {args.ref} tidak ditemukan.")
        sys.exit(1)
        
    if not Path(args.query).exists():
        print(f"Error: File query {args.query} tidak ditemukan.")
        sys.exit(1)
        
    # Load sequences
    ref_record = next(SeqIO.parse(args.ref, "fasta"))
    query_record = next(SeqIO.parse(args.query, "fasta"))
    
    ref_seq = str(ref_record.seq)
    query_seq = str(query_record.seq)
    
    print(f"Memproses alignment: {query_record.id} terhadap {ref_record.id}...")
    
    # Run algorithm
    alignment = smith_waterman(ref_seq, query_seq)
    mutations = detect_mutations(alignment, context_size=args.context)
    summary = summarize_mutations(mutations)
    
    # Output results
    print("\n" + "="*50)
    print("HASIL ALIGNMENT & DETEKSI MUTASI")
    print("="*50)
    print(f"Skor Alignment : {alignment.score}")
    print(f"Identity       : {alignment.identity:.2f}%")
    print(f"Total Mutasi   : {summary['total']}")
    print(f"  - Substitusi : {summary['substitusi']}")
    print(f"  - Insersi    : {summary['insersi']}")
    print(f"  - Delesi     : {summary['delesi']}")
    
    if mutations:
        print("\nDetail Mutasi:")
        for i, mut in enumerate(mutations, 1):
            print(f"  {i}. {mut}")

if __name__ == "__main__":
    main()
