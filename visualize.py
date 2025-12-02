#!/usr/bin/env python3
"""
Visualization tool for grid-based FMM results
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys


def plot_particles(csv_file, output_file=None):
    """Plot particle positions from CSV file"""
    # Read data
    data = np.loadtxt(csv_file, delimiter=',', skiprows=1)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    xs = data[:, 0]
    ys = data[:, 1]

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 10))

    ax.scatter(xs, ys, s=10, alpha=0.6, c='blue', edgecolors='black', linewidths=0.5)

    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)
    ax.set_title(f'Particle Distribution: {csv_file}', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150)
        print(f"Saved plot to {output_file}")
    else:
        plt.show()

    plt.close()


def animate_trajectory():
    """Create animation from frame files"""
    import matplotlib.animation as animation

    # Find all frame files
    frame_files = sorted(Path('.').glob('frame_*.csv'))

    if not frame_files:
        print("No frame files found!")
        return

    print(f"Found {len(frame_files)} frames")

    # Read first frame to get limits
    data0 = np.loadtxt(frame_files[0], delimiter=',', skiprows=1)
    xmin, xmax = data0[:, 0].min(), data0[:, 0].max()
    ymin, ymax = data0[:, 1].min(), data0[:, 1].max()
    margin = 0.1 * max(xmax - xmin, ymax - ymin)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(xmin - margin, xmax + margin)
    ax.set_ylim(ymin - margin, ymax + margin)
    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    scatter = ax.scatter([], [], s=10, alpha=0.6, c='blue',
                        edgecolors='black', linewidths=0.5)
    title = ax.text(0.5, 1.02, '', transform=ax.transAxes,
                   ha='center', fontsize=16)

    def init():
        scatter.set_offsets(np.empty((0, 2)))
        return scatter, title

    def update(frame_idx):
        data = np.loadtxt(frame_files[frame_idx], delimiter=',', skiprows=1)
        scatter.set_offsets(data)
        title.set_text(f'Frame {frame_idx} / {len(frame_files)-1}')
        return scatter, title

    ani = animation.FuncAnimation(fig, update, init_func=init,
                                 frames=len(frame_files),
                                 interval=200, blit=True)

    # Save animation
    ani.save('trajectory.gif', writer='pillow', fps=5, dpi=100)
    print("Saved animation to trajectory.gif")

    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Visualize FMM results')
    parser.add_argument('--file', type=str, help='CSV file to plot')
    parser.add_argument('--animate', action='store_true',
                       help='Create animation from frame_*.csv files')
    parser.add_argument('--output', type=str, help='Output file for plot')

    args = parser.parse_args()

    if args.animate:
        animate_trajectory()
    elif args.file:
        plot_particles(args.file, args.output)
    else:
        # Default: plot first frame if exists
        if Path('frame_00000.csv').exists():
            plot_particles('frame_00000.csv', 'particles.png')
        else:
            print("Usage:")
            print("  python visualize.py --file frame_00000.csv")
            print("  python visualize.py --animate")


if __name__ == '__main__':
    main()
