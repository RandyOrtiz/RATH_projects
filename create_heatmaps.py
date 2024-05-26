import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def create_heatmap(file_path, output_image):
    # Load the data
    data = pd.read_csv(file_path, sep='\t')

    # Set the index to the first column (Name)
    data.set_index(data.columns[0], inplace=True)

    # Calculate dynamic figure size
    num_rows, num_cols = data.shape
    fig_width = num_cols  # Make the width dynamic based on the number of columns
    fig_height = (num_rows * 0.5) / 2  # Half the height

    # Generate the heatmap with dynamic figure size and smaller font sizes
    plt.figure(figsize=(fig_width, fig_height))
    sns.set(font_scale=0.5)  # Set the font scale to be smaller
    heatmap = sns.heatmap(data, annot=False, cmap='viridis', linewidths=.5, cbar_kws={"shrink": .75})

    # Trim x-axis labels to the first 10 characters
    x_labels = [label[:10] for label in data.columns]
    heatmap.set_xticklabels(x_labels, fontsize=6, rotation=45, ha='right')
    heatmap.set_yticklabels(heatmap.get_yticklabels(), fontsize=6)

    # Remove axis labels
    heatmap.set_xlabel('')
    heatmap.set_ylabel('')

    plt.tight_layout()

    # Save the heatmap to a file
    plt.savefig(output_image, dpi=300, bbox_inches='tight')  # Ensure layout fits
    plt.close()

if __name__ == "__main__":
    file_paths = sys.argv[1:]
    for file_path in file_paths:
        output_image = file_path.replace('.txt', '.png')
        create_heatmap(file_path, output_image)
