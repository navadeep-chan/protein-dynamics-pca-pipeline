import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

INPUT_FILE = "PCA_projection.csv"
OUTPUT_CSV = "PCA_clusters.csv"
OUTPUT_PLOT = "PCA_plot.png" 

#  LOAD THE PCA DATA----------------------------------------------------------

if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(f"Error: '{INPUT_FILE}' not found in the working directory.")

print(f"Loading data from {INPUT_FILE}...")
df = pd.read_csv(INPUT_FILE)

if not {'PC1', 'PC2'}.issubset(df.columns):
    raise ValueError("The input CSV must contain 'PC1' and 'PC2' columns.")


#  PERFORM K-MEANS CLUSTERING (4 CLASSES)-------------------------------

print("Performing K-Means clustering (k=4)...")

pca_features = df[['PC1', 'PC2']]

kmeans = KMeans(n_clusters=4, random_state=40)
cluster_labels = kmeans.fit_predict(pca_features)

df['Cluster'] = cluster_labels

print(f"Saving clustered data to {OUTPUT_CSV}...")
df.to_csv(OUTPUT_CSV, index=False)


# CREATE AND SAVE THE 2D PLOT----------------------------------------------------------
print("Generating PCA scatter plot...")

plt.figure(figsize=(10, 8))

colors = ['red', 'blue', 'yellow', 'green']

for cluster_id in range(4):
    cluster_data = df[df['Cluster'] == cluster_id]
    
    plt.scatter(
        cluster_data['PC1'], 
        cluster_data['PC2'], 
        c=colors[cluster_id], 
        label=f'Class {cluster_id}', 
        edgecolor='black', 
        s=60,              
        alpha=0.8          
    )

centroids = kmeans.cluster_centers_
plt.scatter(
    centroids[:, 0], 
    centroids[:, 1], 
    c='black', 
    marker='X', 
    s=200, 
    label='Centroids'
)

plt.title('PCA of PDB Structures with K-Means Clustering (4 Classes)', fontsize=14)
plt.xlabel('Principal Component 1', fontsize=12)
plt.ylabel('Principal Component 2', fontsize=12)
plt.legend(title="Clusters")
plt.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300)
print(f"Plot saved successfully as {OUTPUT_PLOT}")

plt.show()

print("\nClustering and plotting completed successfully.")