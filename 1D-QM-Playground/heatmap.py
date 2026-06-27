import os
import numpy as np
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

x0 = 10  # initial position of the wavepacket

if __name__ == "__main__":
    # plot settings
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 12,
        'figure.dpi': 1000,
        'savefig.dpi': 1000
    })

    timeDir = os.getcwd() + "/cpp/build/timesteps"
    pointfiles = sorted(
        (os.path.join(timeDir, p) for p in os.listdir(timeDir) if p.startswith("Timestep_")),
        key=lambda s: int(os.path.basename(s).split("_")[-1])
    )

    data = []
    for path in tqdm(pointfiles):
        # tqdm.write(f"processing {path}")
        step = int(os.path.basename(path).split("_")[-1])
        with open(path, "r") as f:
            points = f.readlines()
            
        # Extract data
        x = [p.split()[0] for p in points]
        re = [float(p.split()[1]) for p in points]
        im = [float(p.split()[2]) for p in points]

        # Compute probability density
        d = [r**2+i**2 for r, i in zip(re, im)]

        data.append(d)

    # Create heatmap
    data = np.array(data)
    data = data[::-1]     
    plt.tight_layout()
    sns.heatmap(data, cmap="cool")

    # Set x-ticks from 0 to 100 with step 10
    plt.xlabel("r")
    physical = np.arange(0, 101, 10)   # desired x-values
    indices = physical * 3             # convert to array indices
    plt.xticks(ticks=indices, labels=physical)
    
    # Set y-ticks from 0 to 1000 with step 100, with 0 at the bottom
    N = data.shape[0]
    tick_labels = np.arange(0, 1001, 100)              # 0, 100, ..., 1000
    tick_positions = np.linspace(N-1, 0, len(tick_labels))
    plt.yticks(ticks=tick_positions, labels=tick_labels)

    # Set labels and title
    plt.ylabel("Time")
    plt.title(fr"$|\psi(x,t)|^2$ for $x_0 = {x0}$", pad=20)
    plt.savefig(f"{timeDir}/heatmap_x{x0}.png", bbox_inches="tight")