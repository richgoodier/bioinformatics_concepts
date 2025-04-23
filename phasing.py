import numpy as np
import matplotlib.pyplot as plt


def generate_sequence(length: int = 50) -> np.ndarray:
    """
    Generate a random nucleotide sequence of the given length.

    Args:
        length (int): Length of the sequence.

    Returns:
        np.ndarray: Array of nucleotide characters.
    """
    nucleotides = ["A", "T", "G", "C"]
    return np.random.choice(nucleotides, length)

def simulate_one_read(sequence: np.ndarray, error_rate: float = 0.1) -> np.ndarray:
    """
    Simulate a sequencing read with phasing and prephasing errors.

    Args:
        sequence (np.ndarray): The original reference sequence (1D array of bases).
        error_rate (float): Probability of a phasing error per cycle (default: 0.1).

    Returns:
        np.ndarray: Simulated read with potential insertions ('-') and base shifts.
    """
    read = []
    sequence_index = 0
    read_distance = len(sequence)
    for i in range(read_distance):
        if np.random.random() < error_rate:
            if np.random.random() < 0.5:
                read.append('-')
                #print('lag')
            else:
                if sequence_index + 1 < len(sequence):
                    sequence_index +=1
                read.append(sequence[sequence_index])
                #print('jump')
        else:
            read.append(sequence[sequence_index])
            if sequence_index + 1 < len(sequence):
                sequence_index += 1

    return np.array(read)

def calculate_read_values(reads: np.ndarray) -> np.ndarray:
    """
    Count the frequency of each base (G, A, T, C) at every read cycle.

    Args:
        reads (np.ndarray): 2D array of shape (num_reads, read_length).

    Returns:
        np.ndarray: Array of shape (read_length, 4) with base counts ordered as [G, A, T, C].
    """
    read_len = reads.shape[1] # (100, 45) 100 reads of len 45
    read_values = np.zeros((read_len, 4), dtype=int)
    for i in range(read_len):
        counts = {
            base: np.count_nonzero(reads[:, i] == base)
            for base in ["G", "A", "T", "C"]
        }
        read_values[i] = [counts["G"], counts["A"], counts["T"], counts["C"]]
    
    return read_values

def generate_consensus_sequence(read_values: np.ndarray) -> np.ndarray:
    """
    Generate a consensus sequence by selecting the most common base at each cycle.

    Args:
        read_values (np.ndarray): Array of base counts per cycle (shape: read_len × 4).

    Returns:
        np.ndarray: Consensus sequence (1D array of base characters).
    """
    base_labels = np.array(["G", "A", "T", "C"])

    # Get the index of the max base per cycle
    max_indices = np.argmax(read_values, axis=1)  # read_values shape: (read_len, counts)

    # Convert indices to base characters
    consensus_sequence = base_labels[max_indices]
    
    return consensus_sequence

def record_misreads(sequence: np.ndarray, consensus_sequence: np.ndarray) -> tuple[np.ndarray, int | None]:
    """
    Track cumulative misreads and the index of the first mismatch.

    Args:
        sequence (np.ndarray): True reference sequence.
        consensus_sequence (np.ndarray): Inferred sequence from base signals.

    Returns:
        tuple[np.ndarray, int]: 
            - Array of cumulative misreads per cycle.
            - Index of first misread (or None if all correct).
    """
    accurate_read = (sequence == consensus_sequence)
    accumulated_misreads = np.zeros(len(sequence), dtype=int)
    num_misreads = 0
    first_misread_index = None
    for i, read in enumerate(accurate_read):
        if not read:
            num_misreads += 1
            if num_misreads == 1:
                first_misread_index = i
        accumulated_misreads[i] = num_misreads

    return accumulated_misreads, first_misread_index

def plot_Illumina_read(
    sequence: np.ndarray,
    read_values: np.ndarray,
    consensus_sequence: np.ndarray,
    accumulated_misreads: np.ndarray,
    error_rate: float
) -> plt.Figure:
    """
    Plot base signal intensities and cumulative misreads over sequencing cycles.

    Displays the consensus sequence along the bottom, color-coded by accuracy.

    Args:
        sequence (np.ndarray): True sequence.
        read_values (np.ndarray): Base signal counts per cycle (shape: read_len × 4).
        consensus_sequence (np.ndarray): Most frequent base at each cycle.
        accumulated_misreads (np.ndarray): Cumulative mismatches at each position.
        error_rate (float): Error rate used in the simulation.

    Returns:
        matplotlib.figure.Figure: The generated plot figure.
    """
    G = read_values[:, 0]
    A = read_values[:, 1]
    T = read_values[:, 2]
    C = read_values[:, 3]

    fig, ax = plt.subplots(figsize=(15, 4))
    x = range(len(sequence))

    ax.plot(x, G, color='orange', label='G')
    ax.plot(x, A, color='green', label='A')
    ax.plot(x, T, color='red', label='T')
    ax.plot(x, C, color='blue', label='C')
    ax.plot(x, accumulated_misreads, color='black', label='accumulated misreads', linestyle=':')

    ax.set_xlabel("Read Cycle")
    ax.set_ylabel("Fluorescence Intensity")
    ax.set_title(f"Base Composition Over Read Cycles (error rate = {error_rate})", pad=30)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=5)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_xlim((-1, len(sequence)))
    ax.set_ylim(bottom=-6, top=105)
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)

    # === Add color-coded consensus bases at y = -5 ===
    base_colors = {'G': 'orange', 'A': 'green', 'T': 'red', 'C': 'blue'}
    #for i, base in enumerate(consensus_sequence):
    #    plt.text(i, -2, base, color=base_colors[base], ha='center', va='center', fontsize=5)

    for i, base in enumerate(consensus_sequence):
        color = 'black' if base == sequence[i] else base_colors[base]
        ax.text(i, -3, base, color=color, ha='center', va='center', fontsize=7)
    
    return fig


def simulate_phasing_errors(actual: np.ndarray, error_rate: float) -> np.ndarray:
    """
    Simulate phasing errors by introducing increasing mismatches over cycles.

    Args:
        actual (np.ndarray): The true sequence.
        error_rate (float): Base error rate; mismatch probability increases with cycle.

    Returns:
        np.ndarray: Simulated read with phasing errors.
    """
    read = actual.copy()
    length = len(actual)
    for i in range(length):
        error_chance = 1 - np.exp(-error_rate * i)  # grows with position
        if np.random.rand() < error_chance:
            read[i] = np.random.choice([b for b in ["A", "T", "G", "C"] if b != actual[i]])
    return read


def plot_phasing_degradation(actual: np.ndarray, read: np.ndarray, error_rate: float) -> plt.Figure:
    """
    Create a plot showing cumulative mismatch accumulation due to phasing errors.

    Args:
        actual (np.ndarray): The true sequence.
        read (np.ndarray): The simulated read.
        error_rate (float): The phasing error rate used.

    Returns:
        plt.Figure: A matplotlib figure object.
    """
    mismatches = np.array([actual[i] != read[i] for i in range(len(actual))], dtype=int)
    cumulative = np.cumsum(mismatches)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(range(1, len(actual) + 1), cumulative, color='crimson', marker='o')
    ax.set_title(f"Phasing Error Accumulation (error rate = {error_rate})")
    ax.set_xlabel("Read Cycle")
    ax.set_ylabel("Cumulative Mismatches")
    ax.grid(True, linestyle="--", alpha=0.5)
    return fig
