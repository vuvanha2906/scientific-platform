from typing import Dict, Any


class ScreeningFilter:
    @staticmethod
    def evaluate(props: Dict[str, Any], max_violations_allowed: int) -> Dict[str, Any]:
        """
        Apply the Lipinski Rules and calculate the score.
        Score = QED - 0.1 * violations
        """
        violations = 0
        violation_details = []

        # Rule 1: MW <= 500
        if props["mw"] > 500:
            violations += 1
            violation_details.append("MW > 500")

        # Rule 2: LogP <= 5
        if props["logp"] > 5:
            violations += 1
            violation_details.append("LogP > 5")

        # Rule 3: HBD <= 5
        if props["hbd"] > 5:
            violations += 1
            violation_details.append("HBD > 5")

        # Rule 4: HBA <= 10
        if props["hba"] > 10:
            violations += 1
            violation_details.append("HBA > 10")

        # Rule 5: TPSA <= 140
        if props["tpsa"] > 140:
            violations += 1
            violation_details.append("TPSA > 140")

        # Scoring Logic
        # Score = QED - 0.1 * violations
        base_qed = props.get("qed", 0)
        penalty = 0.1 * violations
        final_score = base_qed - penalty

        is_accepted = violations <= max_violations_allowed

        return {
            "violations_count": violations,
            "violation_details": violation_details,
            "final_score": final_score,
            "is_accepted": is_accepted
        }