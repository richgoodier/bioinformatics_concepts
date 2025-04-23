import numpy as np
import matplotlib.pyplot as plt

import random

def generate_non_overlapping_sites(
    genome_length: int, binding_site_length: int, count: int = 4
) -> tuple[list[int], list[int]]:
    """
    Generate two groups of non-overlapping binding site start indices.

    Args:
        genome_length (int): Total length of the genome.
        binding_site_length (int): Length of each binding site.
        count (int): Total number of non-overlapping sites (must be even).

    Returns:
        tuple[list[int], list[int]]: Two lists of binding site start positions.
    """
    assert count % 2 == 0, "Count must be even to divide into two equal groups."
    max_start = genome_length - binding_site_length
    used = set()
    sites = []

    while len(sites) < count:
        candidate = random.randint(0, max_start)
        if all(abs(candidate - s) >= binding_site_length for s in used):
            sites.append(candidate)
            used.add(candidate)

    mid = count // 2
    return sites[:mid], sites[mid:]




def create_reference_genome_chip(
    genome_length: int,
    binding_site_length: int,
    ideal_site_locations: list[int],
    good_site_locations: list[int]
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate a random reference genome and insert ideal and good binding sites.

    Args:
        genome_length (int): Total length of the reference genome.
        binding_site_length (int): Length of each binding site motif.
        ideal_site_locations (list[int]): Start positions to insert ideal (perfect match) binding sites.
        good_site_locations (list[int]): Start positions to insert good (imperfect match) binding sites.

    Returns:
        tuple[np.ndarray, np.ndarray]: 
            - The generated reference genome as a NumPy array of nucleotides.
            - The ideal binding site as a NumPy array for comparison.
    """

    
    nucleotides = ["A", "T", "G", "C"]
    ideal_binding_site = np.random.choice(nucleotides, binding_site_length)
    mutations = np.random.choice(nucleotides, 3)
    good_binding_site = np.copy(ideal_binding_site)
    for i in range(3):
        base_location = np.random.randint(0, binding_site_length)
        good_binding_site[base_location] = mutations[i]
    
    reference_genome = np.random.choice(["A", "T", "G", "C"], genome_length)
    for i in good_site_locations:
        reference_genome[i:i + binding_site_length] = good_binding_site
    for i in ideal_site_locations:
        reference_genome[i:i + binding_site_length] = ideal_binding_site
    
    return reference_genome, ideal_binding_site


def create_reads_chip(
    reference_genome: np.ndarray,
    binding_site: np.ndarray,
    binding_site_length: int,
    antibody_specificity: int = 4
) -> np.ndarray:
    """
    Simulate ChIP-seq reads based on binding site similarity and antibody specificity.

    Args:
        reference_genome (np.ndarray): Array of nucleotides representing the genome.
        binding_site (np.ndarray): The ideal binding site to detect.
        binding_site_length (int): Length of the binding site to scan for.
        antibody_specificity (int): Integer from 1 (poor) to 4 (best), indicating how selective the antibody is.
                                    Higher values produce sharper peaks centered on perfect binding sites.
                                    Typical dropdown menu labels: [1 → "poor", 2 → "good", 3 → "better", 4 → "best"]

    Returns:
        np.ndarray: Read coverage map with jittered alignment counts per position.
    """

    genome_length = len(reference_genome)

    # Map affinities
    affinities = []
    for i in range(genome_length - binding_site_length):
        affinity = (reference_genome[i:i+binding_site_length] == binding_site).sum()
        affinities.append(affinity)

    raw_affinities = np.array(affinities)
    normalized_affinities = raw_affinities / binding_site_length
    soft_affinities = np.clip(normalized_affinities, a_min=0.05, a_max=None)
    
    # Create reads
    num_reads = 100
    reads_per_position = []
    affinity_freq = soft_affinities ** antibody_specificity
    for i in range(genome_length - binding_site_length):
        affinity = np.array([affinity_freq[i] for _ in range(num_reads)])
        reads = (np.random.random(num_reads) < affinity).sum()
        reads_per_position.append(reads)

    reads_per_position = np.array(reads_per_position)

    # Jitter reads
    k = 9  # spread, aka number of buckets
    mu = (k - 1) / 2
    ends = int(mu)

    sigma = k / 6  # adjust spread; ~99.7% of values within [0, k-1]
    len_norm_rpp = len(reads_per_position) + 2 * ends
    jittered_rpp = np.zeros(len_norm_rpp, dtype=int)
    
    for i in range(len(reads_per_position)):
        samples = np.random.normal(loc=mu, scale=sigma, size=reads_per_position[i])
        buckets, _ = np.histogram(samples, bins=np.arange(k + 1))
        jittered_rpp[i:i+k] += buckets

    return jittered_rpp


def plot_read_map(read_map: np.ndarray) -> plt.Figure:
    """
    Plot a bar graph showing the number of reads aligned to each position in the genome.

    Args:
        read_map (np.ndarray): Array representing the number of reads aligned at each genomic position.

    Returns:
        matplotlib.figure.Figure: Bar plot of read alignment across the genome.
    """

    x = [i for i in range(len(read_map))]
    y_max = int(read_map.max() * 1.25)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(x, read_map, width=1, align='edge')
    ax.set_xlim(0, len(read_map) + 2)
    ax.set_ylim(0, y_max)
    ax.set_title("Positive Reads Aligned to Reference Genome (after Gaussian Spread)")
    ax.set_xlabel("Base Index")
    ax.set_ylabel("Reads")
    
    return fig