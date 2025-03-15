import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter


def simulate_sanger(seq_len: int, dd_ratio: float, num_reactions: int) -> np.ndarray:
    fragment_counts = np.zeros(seq_len, dtype=int)
    
    # Generate termination positions directly using a geometric distribution
    termination_sites = np.random.geometric(p=dd_ratio, size=num_reactions)
    
    # Filter out terminations beyond sequence length
    termination_sites = termination_sites[termination_sites <= seq_len]

    # Count occurrences of each termination site
    unique, counts = np.unique(termination_sites, return_counts=True)
    fragment_counts[unique - 1] = counts  # Adjust indexing since positions start at 1

    return fragment_counts, termination_sites

def plot_fragment_counts(fragment_counts, dd_ratio, seq_len):
    fig, ax = plt.subplots(figsize=(10, 5))
    
    ax.bar(range(1, seq_len + 1), fragment_counts, color="blue", alpha=0.7, edgecolor="black")
    ax.set_xlabel("Fragment Length (Termination Position)")
    ax.set_ylabel("Count of Fragments")
    ax.set_title(f"Sanger Sequencing Fragment Distribution (ddNTP/dNTP = {dd_ratio})")
    
    return fig

def plot_gel_electrophoresis(fragment_counts, seq_len, dd_ratio):
    """
    Simulate gel electrophoresis for Sanger sequencing results.

    Parameters:
        fragment_counts (np.ndarray): Number of occurrences of each fragment length.
        seq_len (int): Maximum sequence length.
        dd_ratio (float): ddNTP/dNTP ratio.

    Returns:
        fig, ax: Matplotlib figure and axis.
    """
    # Create gel matrix
    gel_height = 600  # Increased height for better separation
    gel_width = 300   
    gel = np.zeros((gel_height, gel_width))  

    # Log-based migration but adjusted to avoid excessive compression at small sizes
    fragment_sizes = np.arange(1, seq_len + 1)
    log_migration = np.log(seq_len + 1 - fragment_sizes) ** 1.2  # Exaggerate differences in small fragment spacing
    log_migration = (log_migration - np.min(log_migration)) / (np.max(log_migration) - np.min(log_migration))  
    log_migration = (log_migration * 0.75 + 0.15) * gel_height  

    # Normalize fragment intensity to adjust brightness
    fragment_intensity = fragment_counts / np.max(fragment_counts) if np.max(fragment_counts) > 0 else fragment_counts

    # Generate bands with adjusted thickness
    lane_x = gel_width // 2  
    lane_width = 50  # Slightly wider bands for better visibility

    for i, (size, intensity) in enumerate(zip(fragment_sizes, fragment_intensity)):
        if intensity > 0:
            y_pos = int(gel_height - log_migration[i])  
            band_thickness = max(1, int(3 * intensity))  # Thicker bands for more intense signals
            gel[y_pos - band_thickness:y_pos + band_thickness, lane_x - lane_width//2 : lane_x + lane_width//2] += intensity  

    # Apply Gaussian blur for realistic band diffusion
    gel = gaussian_filter(gel, sigma=2)  # Reduced blur for clarity

    # Improved colormap for better contrast
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["black", "deepskyblue", "white"])  

    # Plot gel
    fig, ax = plt.subplots(figsize=(2, 5))
    ax.imshow(gel, cmap=cmap, aspect="auto", extent=[0, gel_width, 0, gel_height])

    # Formatting for realism
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"Sanger Gel\n(ddNTP/dNTP = {dd_ratio})", fontsize=10, color="white", pad=10)
    ax.set_frame_on(False)  

    # Dark background mimicking UV transilluminator
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    return fig
