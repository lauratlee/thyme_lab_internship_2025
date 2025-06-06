# plots the % similarities of a given csv alignment summary file and saves it as an image

import sys
import pandas as pd
import matplotlib.pyplot as plt

csv_file = sys.argv[1]

# Read the CSV file
df = pd.read_csv(csv_file)

# Check if the column exists
if "% similarity" not in df.columns:
    raise ValueError("Column '% similarity' not found in the CSV file.")

# Extract the values from the column
similarities = df["% similarity"]

# Plot the values
plt.figure(figsize=(10, 5))
plt.plot(similarities, marker='o', linestyle='-', color='blue')
plt.title("Percent Similarity")
plt.xlabel("Index")
plt.ylabel("% Similarity")
plt.grid(True)

# Save the plot as an image file
output_path = "similarity_plot.png"
plt.savefig(output_path)
print(f"Plot saved as '{output_path}'")
