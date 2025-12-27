from typing import Dict, Any
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Lipinski

class MoleculeTools:
    @staticmethod
    def compute_descriptors(smiles: str) -> Dict[str, Any]:
        """
        Calculate the necessary physicochemical parameters.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"valid": False, "smiles": smiles}

        try:
            return {
                "valid": True,
                "smiles": Chem.MolToSmiles(mol, canonical=True), # Canonicalize
                "mw": Descriptors.MolWt(mol),
                "logp": Descriptors.MolLogP(mol),
                "hbd": Lipinski.NumHDonors(mol),
                "hba": Lipinski.NumHAcceptors(mol),
                "tpsa": Descriptors.TPSA(mol),
                "rot_bonds": Descriptors.NumRotatableBonds(mol),
                "qed": QED.qed(mol)
            }
        except:
            return {"valid": False, "smiles": smiles}