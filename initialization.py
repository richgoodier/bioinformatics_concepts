from typing import List, Tuple
import numpy as np
import os
import csv

def create_reference_genome(reference_length: int) -> np.ndarray:
    """
    Creates an nd.array of nucleotide base n long as a reference genome

    Args:
        reference_length (int): size of reference genome

    Returns:
        nd.array: reference genome where each index is one nucleotide
    """
    return np.random.choice(["A", "T", "G", "C"], reference_length)

def index_reference_genome(reference_genome: np.ndarray, kmer_length: int = 3) -> dict:
    """
    Takes a reference genome and indexes it, returning the index

    Args:
        referenc_genome (nd.array): reference genome
        kmer_length (int): kmer length to index

    Returns:
        dict: index of genome by kmer
    """
    reference_index = {}
    reference_length = len(reference_genome)
    
    for i in range(reference_length - (kmer_length - 1)):
        kmer = "".join(list(reference_genome[i: i + kmer_length]))
        if kmer in reference_index:
            reference_index[kmer].append(i)
        else:
            reference_index[kmer] = [i]
    
    return reference_index

