import os
import glob
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from natsort import natsorted

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Set global domain length for all cases
L = 1.0

def make_animation(snapshot_dir, output_filename, title):
    files = glob.glob(os.path.join(snapshot_dir, "wave_timestep_*.txt"))
    
    if not files:
        print(f"No snapshot files found in {snapshot_dir}")
        return
    files = natsorted(files)

    # Load sample to get spatial grid
    y_sample = np.loadtxt(files[0])
    Nx = len(y_sample)
    x = np.linspace(-L/2, L/2, Nx)
    wave_frames = [np.loadtxt(f) for f in files]

    fig, ax = plt.subplots()
    line, = ax.plot(x, wave_frames[0])
    ax.set_xlim(-L/2, L/2)
    ax.set_ylim(-1.1 * np.max(np.abs(wave_frames)), 1.1 * np.max(np.abs(wave_frames)))
    ax.set_xlabel("x")
    ax.set_ylabel("Amplitude")
    ax.set_title(title)
    ax.grid(True)

    writer = PillowWriter(fps=20)
    print(f"Rendering frames for {output_filename}...")
    start_time = time.time()

    with writer.saving(fig, output_filename, dpi=100):
        total = len(wave_frames)
        for i, frame in enumerate(wave_frames):
            line.set_ydata(frame)
            writer.grab_frame()
            if (i + 1) % (total // 10) == 0:
                print(f"{int((i + 1) / total * 100)}% complete (frame rendering)")

    print("Rendering complete. Finalizing and saving the animation...")
    print(f"Successfully saved: {output_filename}")
    print(f"Total time elapsed: {time.time() - start_time:.2f} seconds\n")
    plt.close(fig)

# Generate animations
make_animation("Finite_Diff_Snapshots", "wave_animation_fd.gif", "Wave Evolution\n(Finite Difference Method)")
make_animation("SineFourier_Snapshots", "wave_animation_fourier.gif", "Wave Evolution\n(Fourier Series Method)")
time.sleep(3)
