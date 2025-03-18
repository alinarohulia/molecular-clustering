# This code generates 6 files for 6 feature types:
# Hydrogen Bond Donor, Hydrogen Bond Acceptor, Aromatic Rings, Hydrophobes + Lumped Hydrophobes, Positive Ionizable, Negative Ionizable

import os
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

# Paths for input SDF file
input_sd_file = "D:/Research Dr. Kireev/Clustering project/CIB1 protein-FRASE hits.sdf"
output_dir = "D:/Research Dr. Kireev/Clustering project/output_families/"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Assign chemical elements to feature families (different elements to visualize atoms in different colors in Maestro Schrodinger)
family_to_element = {
    "Donor": "P",
    "Acceptor": "Ne",
    "Aromatic": "N",
    "Hydrophobe": "Cl",
    "LumpedHydrophobe": "Cl",
    "PosIonizable": "O",
    "NegIonizable": "S",
}

# Build feature factory
fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

def generate_mol_block(fragment_name, feature, element):
    x, y, z = feature["coordinates"]

    x_str = f"{x:.4f}"
    y_str = f"{y:.4f}"
    z_str = f"{z:.4f}"

    # coordinates in mol blocks shoud be separated with spaces according to amount of digits from the beginning of the string
    # this code accounts for single and double digits, negative and positive numbers for x, y, and z coordinates
    if len(x_str.split('.')[0]) == 1:
        spacing_x = "    "
    elif len(x_str.split('.')[0]) == 2:
        spacing_x = "   "
    else:
        spacing_x = "  "

    if len(y_str.split('.')[0]) == 1:
        spacing_y = "    "
    elif len(y_str.split('.')[0]) == 2:
        spacing_y = "   "
    else:
        spacing_y = "  "

    if len(z_str.split('.')[0]) == 1:
        spacing_z = "    "
    elif len(z_str.split('.')[0]) == 2:
        spacing_z = "   "
    else:
        spacing_z = "  "

    mol_block = f"""{fragment_name}
     RDKit          3D

  1  0  0  0  0  0  0  0  0  0999 V2000
{spacing_x}{x:.4f}{spacing_y}{y_str}{spacing_z}{z:.4f} {element}   0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <PDBID_FRAGID>  (1) 
{fragment_name}-{feature["type"]}

>  <feature_family>  (1) 
{feature["family"]}

>  <feature_type>  (1) 
{feature["type"]}

$$$$
"""
    return mol_block


# Write output .SDF file
with Chem.SDMolSupplier(input_sd_file) as supplier:
    feature_files = {}

    for mol in supplier:
        if mol is None:
            continue

        fragment_name = mol.GetProp("_Name")
        feats = factory.GetFeaturesForMol(mol)

        features = {}

        for feat in feats:
            family = feat.GetFamily()
            coords = feat.GetPos()

            if family == "ZnBinder":  # Skip ZnBinder family
                continue

            # Combine "Hydrophobe" and "LumpedHydrophobe" under "Hydrophobe"
            if family in ["Hydrophobe", "LumpedHydrophobe"]:
                family = "Hydrophobe"

            if family not in features:
                features[family] = []

            features[family].append({
                "id": feat.GetId(),
                "type": feat.GetType(),
                "family": feat.GetFamily(),
                "coordinates": (coords.x, coords.y, coords.z),
            })

        # Write each family to its own file
        for family, features in features.items():
            output_file = os.path.join(output_dir, f"{family}_features.sdf")

            # Open file in append mode so we can write multiple fragments
            if family not in feature_files:
                feature_files[family] = open(output_file, "w")

            for feature in features:
                element = family_to_element.get(family, "C") # Assign element "C" as a default value for family
                mol_block = generate_mol_block(fragment_name, feature, element)
                feature_files[family].write(mol_block)

for file in feature_files.values():
    file.close()


print("\nFeature files have been successfully generated in the following directory:")
print(output_dir)
print("\nGenerated files:")
for family in feature_files.keys():
    print(f"- {os.path.join(output_dir, f'{family}_features.sdf')}")