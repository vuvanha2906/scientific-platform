from typing import List
import random
from rdkit import Chem


class GeneratorAgent:
    @staticmethod
    def _mutate_smiles(smiles: str) -> str:
        """
        Perform simple chemical mutation (Deterministic logic).
        Use simple carbon addition/subtraction techniques for demonstration.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return smiles

        try:
            # Switch to RWMol for editing.
            rw_mol = Chem.RWMol(mol)
            num_atoms = rw_mol.GetNumAtoms()

            # Tactics:
            # 1. Add a methyl group (-CH3) to a random atom if possible.
            # 2. Or return to the original if too small.

            if num_atoms > 0:
                idx = random.randint(0, num_atoms - 1)
                atom = rw_mol.GetAtomWithIdx(idx)

                # Only add carbon or nitrogen with an equal valence.
                if atom.GetSymbol() in ['C', 'N'] and atom.GetExplicitValence() < atom.GetTotalValence():
                    new_idx = rw_mol.AddAtom(Chem.Atom(6))  # Add Carbon
                    rw_mol.AddBond(idx, new_idx, Chem.BondType.SINGLE)

            new_smiles = Chem.MolToSmiles(rw_mol)
            return new_smiles
        except:
            return smiles  # Safe fallback

    @staticmethod
    def mutate_batch(seeds: List[str], n_variants: int) -> List[str]:
        """Create n variants for each seed."""
        generated = []
        for seed in seeds:
            # Add the original seed to preserve the good results.
            generated.append(seed)
            for _ in range(n_variants):
                new_smi = GeneratorAgent._mutate_smiles(seed)
                if new_smi:
                    generated.append(new_smi)
        return list(set(generated))  # Deduplicate
