# PCA以国家为维度
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Data input (replace with your dataset)
file_path = './Data/PCA_RAWDATA.csv'  # 替换为你的CSV文件路径
df = pd.read_csv(file_path)

# 将DataFrame转换为字典（行索引为键，列为值）
data_dict = df.to_dict(orient='list') 
data_dict.update({"Country": data_dict.pop("Unnamed: 0")})
# print(data_dict)

# Convert the data dictionary into a pandas DataFrame
df = pd.DataFrame(data_dict)

# Isolate the PCA input data (dropping the 'Country' column)
X = df.drop(columns='Country')
# print(X)

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

# Create a PCA plot
plt.figure(figsize=(8,6))
plt.scatter(X_pca[:, 0], X_pca[:, 1], c='blue')

# Annotating points with country names
for i, country in enumerate(df['Country']):
    plt.annotate(country, (X_pca[i, 0], X_pca[i, 1]), fontsize=9)

# Labeling axes based on explained variance ratio
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% of variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% of variance)')
plt.title('PCA of Country Data')
plt.grid(True)
plt.savefig('pca_plot.png')
plt.show()

