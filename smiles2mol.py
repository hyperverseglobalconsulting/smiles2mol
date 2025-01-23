from rdkit import Chem
from rdkit.Chem import Draw

def draw_molecule_from_smiles(smiles, output_file="molecule.png"):
    """
    Draws a molecule from its SMILES string and saves it as an image.

    Parameters:
        smiles (str): The SMILES string of the molecule.
        output_file (str): The name of the output image file.
    """
    try:
        # Convert SMILES string to a molecule
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError("Invalid SMILES string")

        # Draw the molecule and save it to a file
        img = Draw.MolToImage(molecule)
        img.save(output_file)
        print(f"Molecule image saved to {output_file}")
    except Exception as e:
        print(f"Error: {e}")

# Example SMILES strings
smiles_input = "CCO"  # Ethanol
output_image = "ethanol.png"

# Draw and save the molecule
draw_molecule_from_smiles(smiles_input, output_image)

