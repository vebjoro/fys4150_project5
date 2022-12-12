import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

# Fontsize settings
label_fontsize = 16
ticks_fontsize = 16


def plot_prob(P, filename):
    """
    Plot the probability deviation for a given system.
    """

    P = np.array(P).T[0] # Convert to numpy array
    n = len(P) # Number of iterations

    # Plot
    fig = plt.figure(figsize=(6, 4.5))
    plt.plot(range(n), P, color="#8a1629", alpha=0.8, linewidth=2.0)
    plt.xlabel('Number of iterations', fontsize=label_fontsize)
    plt.ylabel('Total probability deviation', fontsize=label_fontsize)
    plt.xticks(fontsize=ticks_fontsize)
    plt.yticks(fontsize=ticks_fontsize)
    plt.grid()
    ax = plt.gca()
    ax.set_facecolor("#e6e6e6")
    plt.tight_layout()

    # Save figure
    plt.savefig(f"./figs/{filename}.pdf")

# Load data for system without double slit barrier
P0 = pa.mat()
P0.load("./data/system1_P_0.bin")

# Load data for system with double slit barrier
P1 = pa.mat()
P1.load("./data/system1_P_1.bin")

# Plot
plot_prob(P0, "system1_P_0")
plot_prob(P1, "system1_P_1")