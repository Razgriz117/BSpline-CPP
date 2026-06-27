import os
import cv2
import glob
import numpy as np
from tqdm import tqdm
import matplotlib
import matplotlib.pyplot as plt

# for parallelization
matplotlib.use("Agg")
from concurrent.futures import ProcessPoolExecutor, as_completed

x0 = 10  # initial position of the wavepacket

# Plot a single time step file
def plotOneTimestep(path, color=None, label=None):
    with open(path, "r") as f:
        lines = f.readlines()
        x = [float(l.split()[0]) for l in lines]
        real = [float(l.split()[1]) for l in lines]
        im = [float(l.split()[2]) for l in lines]

        # Compute probability density
        y = [r**2 + i**2 for r, i in zip(real, im)]

        plt.ylim(0, 0.13)
        plt.xlim(0, 60)
        plt.xlabel("r")
        plt.plot(x, y, color=color, label=label)

        plt.title(fr"$|\psi(x,t)|^2$ for $x_0 = {x0}$", pad=20)

# Plot a single curve from a single file
def plot(path, color=None, label=None):
    with open(path, "r") as f:
        lines = f.readlines()
        x = [float(l.split()[0]) for l in lines]
        y = [float(l.split()[1]) for l in lines]
        
        plt.plot(x, y, color=color, label=label)

# Process a single time step file: plot and save the figure
def process_timestep(path, output_dir):
    step = path.split("_")[-1]
    plt.clf()
    plotOneTimestep(path, label=step)
    out_path = os.path.join(output_dir, f"wave_{step}.png")
    plt.savefig(out_path)
    return out_path

if __name__ == "__main__":
    timeDir = os.getcwd() + "/cpp/build/timesteps"

    # Plot all curves described in every single time step file in timeDir and animate
    
    # MAKE FRAMES
    pointfiles = [os.path.join(timeDir, p) for p in os.listdir(timeDir) if p.startswith("Timestep_")]

    results = [] # Collects all paths returned by process_timestep()
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_timestep, path, timeDir) for path in pointfiles]
        for f in tqdm(as_completed(futures), total=len(futures)):
            results.append(f.result())

    # ANIMATE FRAMES
    frame_paths = sorted(results)

    # Get frame size and initialize video writer
    frame = cv2.imread(frame_paths[0])
    height, width, layers = frame.shape
    out = cv2.VideoWriter(f'{timeDir}/time_evolution_x{x0}.mp4', cv2.VideoWriter_fourcc(*'mp4v'), 45, (width, height))

    # Write frames
    for path in frame_paths:
        frame = cv2.imread(path)
        out.write(frame)

    out.release()

    # Plot a single curve from a single file in timeDir
    # filename = f"{timeDir}/B_G_reconstructed"
    # plt.clf()
    # plot(filename)
    # plt.savefig(f"{filename}_plotted.png")

