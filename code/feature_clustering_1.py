# code to generate clusters with singletons

import numpy as np
from scipy.spatial.distance import pdist, squareform

# Function to perform Butina clustering based on a threshold distance
def butina_clustering(coordinates, threshold):
    coordinates = np.array(coordinates)
    
    # Calculate the pairwise distance matrix (condensed format) for all points
    dist_matrix = pdist(coordinates)
    
    # Convert the condensed distance matrix to a square form (full 2D matrix)
    dist_matrix = squareform(dist_matrix)
    
    # Initialize the clusters list and start with the first point as the first cluster
    clusters = []
    clusters.append([0])  # First cluster containing the first point (index 0)
    
    # Iterate through the remaining points to assign them to clusters
    for i in range(1, len(coordinates)):
        is_in_cluster = False
        # Check each existing cluster for the minimum distance from point i
        for cluster in clusters:
            min_distance = np.min(dist_matrix[i, cluster])  # Get the minimum distance to the cluster
            if min_distance < threshold:  # If it's within the threshold, add the point to the cluster
                cluster.append(i)
                is_in_cluster = True
                break
        
        # If the point wasn't added to any existing cluster, create a new cluster
        if not is_in_cluster:
            clusters.append([i])

    return clusters

# Function to extract 3D coordinates from an SDF file
def extract_coordinates(sdf_file_path):
    coordinates = []  # List to store coordinates
    with open(sdf_file_path, 'r') as file:
        lines = file.readlines()  # Read the file lines

        # Process each line to extract coordinates
        for line in lines:
            if len(line.split()) >= 4 and line.split()[3].isalpha() and line.split()[3] != "M":
                # Extract x, y, z coordinates from the first 3 columns
                x, y, z = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
                coordinates.append((x, y, z))  # Add the coordinates to the list
                
    return coordinates  # Return the list of coordinates

# Function to calculate centroids of each cluster
def calculate_centroids(clusters, coordinates):
    centroids = []  # List to store the centroid of each cluster
    
    # Iterate over each cluster to calculate its centroid
    for cluster in clusters:
        cluster_coords = np.array([coordinates[i] for i in cluster])  # Extract the coordinates of the cluster points
        centroid = np.mean(cluster_coords, axis=0)  # Calculate the mean (centroid) of the cluster coordinates
        centroids.append(centroid)  # Add the centroid to the list
    
    return centroids  # Return the list of centroids

# Function to generate a Mol block format for a centroid
def generate_mol_block(fragment_name, feature, element):
    # Extract the x, y, z coordinates from the feature dictionary
    x, y, z = feature["coordinates"]

    # Format the coordinates to 4 decimal places
    x_str = f"{x:.4f}"
    y_str = f"{y:.4f}"
    z_str = f"{z:.4f}"

    # Determine appropriate spacing for the coordinates based on their length
    def get_spacing(coord_str):
        if len(coord_str.split('.')[0]) == 1:
            return "    "
        elif len(coord_str.split('.')[0]) == 2:
            return "   "
        else:
            return "  "
    
    # Apply appropriate spacing for x, y, and z coordinates
    spacing_x = get_spacing(x_str)
    spacing_y = get_spacing(y_str)
    spacing_z = get_spacing(z_str)
    
    # Construct the Mol block format for the centroid
    mol_block = f"""{fragment_name}
     RDKit          3D

  1  0  0  0  0  0  0  0  0  0999 V2000
{spacing_x}{x:.4f}{spacing_y}{y_str}{spacing_z}{z:.4f} {element}   0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <PDBID_FRAGID>  (1) 
{fragment_name}-{feature["type"]}

>  <feature_type>  (1) 
{feature["type"]}

$$$$
"""
    return mol_block  # Return the generated Mol block

# Function to write the generated centroids to an SDF file
def write_centroids_to_sdf(centroids, output_file_path):
    with open(output_file_path, 'w') as sdf_file:
        # Iterate through the centroids and write each one as a Mol block
        for idx, centroid in enumerate(centroids):
            fragment_name = f"Fragment_{idx + 1}"  # Create a unique fragment name
            feature = {
                "coordinates": centroid,  # Store the centroid coordinates
                "type": "centroid_aromatic",  # Set the type as 'centroid_family_type'
                "family": "cluster_aromatic"  # Set the family as 'cluster_family_type'
            }
            element = "N"  # Specify the element according to "family to element"
            '''family_to_element = {
            "Donor": "P",
            "Acceptor": "Ne",
            "Aromatic": "N",
            "Hydrophobe": "Cl",
            "LumpedHydrophobe": "Cl",
            "PosIonizable": "O",
            "NegIonizable": "S",
            }'''
            mol_block = generate_mol_block(fragment_name, feature, element)  # Generate the Mol block
            sdf_file.write(mol_block)  # Write the Mol block to the output file

# Path to an input SDF file
sdf_file_path = r"D:\Research Dr. Kireev\Clustering project\output_families\Aromatic_features.sdf"
output_sdf_file_path = r"D:\Research Dr. Kireev\Clustering project\output_families\centroids_Aromatic_features.sdf"

# Extract the coordinates from the SDF file
coordinates = extract_coordinates(sdf_file_path)

# Set a threshold distance for clustering (e.g., 3.0 Ångströms)
threshold = 3.0  # You can adjust this value based on your needs

# Call the butina_clustering function to generate clusters based on the threshold
clusters = butina_clustering(coordinates, threshold)

# Calculate centroids for the generated clusters
centroids = calculate_centroids(clusters, coordinates)

# Write the centroids to the SDF file
write_centroids_to_sdf(centroids, output_sdf_file_path)

# Print a success message with the output file path
print(f"Centroids written to {output_sdf_file_path}")