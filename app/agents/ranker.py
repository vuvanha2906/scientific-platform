from typing import List, Dict, Any


class RankerAgent:
    @staticmethod
    def rank_and_select(molecules: List[Dict[str, Any]], top_k: int) -> List[Dict[str, Any]]:
        """
        Sort the list of elements based on 'final_score' in descending order.
        Select the top_k elements.
        """
        # Sort in descending order by score
        sorted_mols = sorted(molecules, key=lambda x: x.get("final_score", 0), reverse=True)

        # Cut off the top
        return sorted_mols[:top_k]