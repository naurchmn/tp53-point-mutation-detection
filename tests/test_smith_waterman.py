import pytest
from src.smith_waterman import smith_waterman
from src.mutation_detector import detect_mutations, MutationType

def test_smith_waterman_exact_match():
    seq1 = "ATGC"
    seq2 = "ATGC"
    result = smith_waterman(seq1, seq2)
    
    assert result.score == 8  # 4 matches * 2
    assert result.aligned_seq1 == "ATGC"
    assert result.aligned_seq2 == "ATGC"
    assert result.identity == 100.0

def test_smith_waterman_substitution():
    seq1 = "ATGC"
    seq2 = "ATTC"  # G substituted with T
    result = smith_waterman(seq1, seq2)
    
    assert result.aligned_seq1 == "ATGC"
    assert result.aligned_seq2 == "ATTC"
    
    mutations = detect_mutations(result)
    assert len(mutations) == 1
    assert mutations[0].mutation_type == MutationType.SUBSTITUTION
    assert mutations[0].ref_base == "G"
    assert mutations[0].mut_base == "T"
    assert mutations[0].position == 3

def test_smith_waterman_insertion():
    seq1 = "ATGC"
    seq2 = "ATAGC" # A inserted after T
    result = smith_waterman(seq1, seq2)
    
    assert result.aligned_seq1 == "AT-GC"
    assert result.aligned_seq2 == "ATAGC"
    
    mutations = detect_mutations(result)
    assert len(mutations) == 1
    assert mutations[0].mutation_type == MutationType.INSERTION
    assert mutations[0].mut_base == "A"
    assert mutations[0].position == 2 # Position in reference sequence where it happened

def test_smith_waterman_deletion():
    seq1 = "ATGC"
    seq2 = "ATC" # G deleted
    result = smith_waterman(seq1, seq2)
    
    assert result.aligned_seq1 == "ATGC"
    assert result.aligned_seq2 == "AT-C"
    
    mutations = detect_mutations(result)
    assert len(mutations) == 1
    assert mutations[0].mutation_type == MutationType.DELETION
    assert mutations[0].ref_base == "G"
    assert mutations[0].position == 3
