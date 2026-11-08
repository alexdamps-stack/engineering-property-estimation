"""Check the approved A2 budget registry against immutable reference budgets."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "validation/tolerances.json"


def main():
    registry = json.loads(PATH.read_text())
    assert registry["status"] == "approved for the listed fixed reference sets"
    for row in registry["reference_sets"]:
        fixture = json.loads((ROOT / row["fixture"]).read_text())
        assert fixture[row["budget_key"]] == row["budgets"], row["fixture"]
    import hashlib
    process = json.loads((ROOT / "validation/fixtures/process-v1.json").read_text())
    for name, digest in process["input_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    assert {r["capability"] for r in registry["deferred"]} == {"critical-point/retrograde envelope continuation"}
    print("PASS: approved budgets exactly match every registered reference fixture, including bounded A13; critical/retrograde extensions deferred")


if __name__ == "__main__":
    main()
