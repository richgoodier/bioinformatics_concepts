import numpy as np
import matplotlib.pyplot as plt

def create_reads(reference_genome: np.ndarray, read_length: int, num_reads: int) -> np.ndarray:
    reference_length = len(reference_genome)
    reads = []
    for _ in range(num_reads):
        start = np.random.randint(0, reference_length - read_length + 1)
        read = reference_genome[start:start+read_length]
        reads.append(read)

    return np.array(reads)

def align_reads(reference_genome: np.ndarray, reference_index: dict[str, list[int]], kmer_length: int, reads:np.ndarray) -> list[int]:
    reference_length = len(reference_genome)
    num_reads, read_length = reads.shape
    read_starts = np.zeros(num_reads, dtype='int')

    for r, read in enumerate(reads):
        read_kmer = "".join(reads[r][0:kmer_length])
        possible_starts = []
        for i in reference_index[read_kmer]:
            if i + read_length > reference_length:
                continue
            possible_starts.append(sum(reference_genome[i:i+read_length] == read))
        best_start = np.argmax(possible_starts)
        read_starts[r] = reference_index[read_kmer][best_start]
    
    return read_starts

def create_scaffold(reference_length: np.ndarray, reads: np.ndarray, read_starts: list[int]) -> np.ndarray:
    num_reads, read_length = reads.shape
    scaffold = np.array(["N" for _ in range(reference_length)])
    for i in range(num_reads):
        scaffold[read_starts[i]: read_starts[i] + read_length] = reads[i]
    
    return scaffold

def calculate_coverage(read_length: int, num_reads: int, reference_length: int):
    return read_length * num_reads / reference_length

def count_unread_bases(reference_genome: np.ndarray, scaffold: np.ndarray) -> int:
    reference_length = len(reference_genome)
    return reference_length - sum(scaffold == reference_genome)

def color_sequence(sequence: str) -> str:
    """Format sequence so that all characters are blue except 'N', which remains black."""
    return "".join(
        f'<span style="color: black;">{char}</span>' if char == "N" 
        else f'<span style="color: blue;">{char}</span>' 
        for char in sequence
    )

def plot_reads(read_length: int, read_starts: np.ndarray, scaffold: np.ndarray, reference_length: int):
    """_summary_

    Args:
        read_length (int): _description_
        read_starts (np.ndarray): _description_
        reference_length (int): _description_

    Returns:
        _type_: _description_
    """
    num_reads = len(read_starts)
    
    fig, ax = plt.subplots(figsize=(10, 4))

    # Plot each read
    for i, start in enumerate(read_starts):
        ax.plot([start, start + read_length], [num_reads - i, num_reads - i], color='blue', linewidth=0.5)

    # Plot scaffold
    scaffold_colors = np.where(scaffold == "N", "gray", "blue")
    for i in range(reference_length - 1):
        ax.plot([i, i + 1], [-0.5, -0.5], color=scaffold_colors[i], linewidth=4)

    # Formatting the plot
    ax.set_xlim(0, reference_length)
    ax.set_ylim(-1, num_reads + 1)
    if num_reads >= 10:
        yticks = list(np.linspace(1, num_reads, min(11, num_reads), dtype=int))
    else:
        yticks = list(range(1, num_reads + 1))
    ax.set_yticks(yticks)

    ax.set_xlabel("Reference Genome Position (Gray are unread bases)")
    ax.set_ylabel("Read")
    ax.set_title("Aligned Reads to Reference Genome")
    ax.grid(True, linestyle="--", alpha=0.5)
    
    return fig

if __name__ == '__main__':
    main()