# plots the % similarities of a given csv alignment summary file and saves it as an image
# run from directory that has the summaries

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys, re

file_path = sys.argv[1]

match = re.search(r'(\d+\.\d+)', file_path)
if match:
    threshold = match.group(1)

# Load the CSV file
df = pd.read_csv(file_path)

# Clean up column names (remove extra spaces)
df.columns = df.columns.str.strip()

# Number of data points
n_points = df['% similarity'].dropna().shape[0]

# Create the distribution plot
plt.figure(figsize=(10, 6))
bins = list(range(60, 105, 5)
sns.histplot(df['% similarity'], bins=bins, color='skyblue', kde=False)
plt.title('Distribution of % Similarity')
plt.xlabel('% Similarity')
plt.ylabel('Frequency')
plt.grid(True)

# Add number of data points text at top right
plt.text(
    x=0.95, y=1.05, 
    s=f'N = {n_points}', 
    ha='right', va='top', 
    transform=plt.gca().transAxes,
    fontsize=12,
    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7)
)

plt.tight_layout()

# Save the plot to a file
plt.savefig(f"similarity_distribution_{threshold}.png", dpi=300)  # High-quality PNG
plt.close()  # Close the figure to free up memory
print(f"saved to similarity_distribution_{threshold}.png")

