from rdkit import Chem
from rdkit.Chem import SDWriter, AllChem
from rdkit.ML.Cluster import Butina
import numpy as np
import os

# Load molecules from SDF
sdf_file = r"D:\Research Dr. Kireev\Clustering project\Fragment clustering\CIB1 protein-FRASE hits.sdf"
supplier = Chem.SDMolSupplier(sdf_file, removeHs=True)
output_dir = r"D:\Research Dr. Kireev\Clustering project\Fragment clustering"
os.makedirs(output_dir, exist_ok=True)
output_sdf = os.path.join(output_dir, "aligned_clustered_fragments.sdf")

# Function to compute geometric center of a molecule
def compute_centroid(mol):
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    centroid = np.mean([conf.GetAtomPosition(i) for i in range(num_atoms)], axis=0)
    return centroid

# Compute centroids
centroids = []
valid_mols = []
for mol in supplier:
    if mol is not None:
        centroid = compute_centroid(mol)
        centroids.append(centroid)
        valid_mols.append(mol)

# Compute Euclidean distance matrix
n_mols = len(centroids)
distance_matrix = []
for i in range(n_mols):
    for j in range(i):
        dist = np.linalg.norm(centroids[i] - centroids[j])  # Euclidean distance
        distance_matrix.append(dist)

# Set a clustering threshold (adjust as needed)
threshold = 15.0  # Distance in Angstroms

# Apply Butina clustering
clusters = Butina.ClusterData(distance_matrix, n_mols, threshold, isDistData=True)

# Compute mean centroid for each cluster
cluster_centroids = {}
for cluster_idx, cluster in enumerate(clusters):
    cluster_points = np.array([centroids[i] for i in cluster])
    mean_centroid = np.mean(cluster_points, axis=0)
    cluster_centroids[cluster_idx] = mean_centroid

# Function to align a molecule to a new centroid
def align_to_centroid(mol, new_centroid):
    conf = mol.GetConformer()
    mol_centroid = compute_centroid(mol)
    translation_vector = new_centroid - mol_centroid
    
    # Apply translation to each atom
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        new_pos = pos + translation_vector
        conf.SetAtomPosition(i, new_pos)


writer = SDWriter(output_sdf)

# Assign clusters, align molecules, and write to SDF
for cluster_idx, cluster in enumerate(clusters):
    mean_centroid = cluster_centroids[cluster_idx]
    
    for mol_idx in cluster:
        mol = valid_mols[mol_idx]
        align_to_centroid(mol, mean_centroid)  # Align fragment to cluster centroid
        mol.SetProp("Cluster_ID", str(cluster_idx + 1))  # Assign cluster number
        writer.write(mol)

writer.close()
print(f"Aligned clustered fragments saved to: {output_sdf}")