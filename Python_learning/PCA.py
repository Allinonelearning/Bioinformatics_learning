import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Load the data
file_path = './Data/PCA_RAWDATA.csv'  # Update the file path if necessary
df = pd.read_csv(file_path)

# Select the relevant columns for PCA
X = df[['L..sativa', 'L..serriola', 'L..saligna', 'L..virosa', 'other']]

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)




# Define categories and their respective colors and markers
categories = ['L..sativa', 'L..serriola', 'L..saligna', 'L..virosa', 'other']
colors = ['red', 'blue', 'green', 'orange', 'purple']
markers = ['o', '^', 's', 'D', '*']

# Create a PCA plot
plt.figure(figsize=(8, 6))

# Plot points with different shapes and colors based on categories
for i, category in enumerate(categories):
    plt.scatter(X_pca[df[category] > 0, 0], X_pca[df[category] > 0, 1], 
                color=colors[i], marker=markers[i], label=category)

# Labeling axes based on explained variance ratio
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% of variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% of variance)')
plt.title('PCA of Lactuca Species Data')
plt.grid(True)
plt.legend()
plt.savefig('pca_plot_by_category1.png')
plt.show()