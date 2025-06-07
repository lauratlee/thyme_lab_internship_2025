# plots the % similarities of a given csv alignment summary file and saves it as an image

import sys, re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# Load CSV
csv_file = sys.argv[1]

if "super_outputs" in csv_file:
    folder = "super_outputs"
else if "align_outputs" in csv_file:
    folder = "align_outputs"
else:
    folder = "unknown"

match = re.search(r"_([0-9.]+)\.csv$", filename)
if match:
    threshold = float(match.group(1))  # or keep as string with match.group(1)
else:
    print("Pattern not found.")
    threshold = 0.0

    
df = pd.read_csv(csv_file)

if "% similarity" not in df.columns:
    raise ValueError("Column '% similarity' not found in the CSV file.")

similarities = df["% similarity"]

# Create histogram data
counts, bins = np.histogram(similarities, bins=20, range=(0, 100))
bin_centers = 0.5 * (bins[1:] + bins[:-1])

# Create a colormap from blue (low) to red (high)
norm = mcolors.Normalize(vmin=0, vmax=100)
cmap = plt.cm.viridis

# Plot bars with color gradient
plt.figure(figsize=(10, 5))
for i in range(len(counts)):
    plt.bar(bin_centers[i], counts[i], width=(bins[1] - bins[0]),
            color=cmap(norm(bin_centers[i])))

plt.xlabel('% Similarity')
plt.ylabel('Frequency')
plt.title('Distribution of % Similarity')
plt.grid(True)

# Save and open
output_path = os.join("..", folder, f"plot_{threshold}.png")
plt.savefig(output_path)
print(f"Saved gradient histogram to {output_path}")
