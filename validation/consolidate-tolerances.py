"""Check the approved A2 budget registry against immutable reference budgets."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "validation/tolerances.json"


def main():
    registry = json.loads(PATH.read_text())
    assert registry["status"] == "approved for the listed fixed reference sets"
    registered = [row["fixture"] for row in registry["reference_sets"]]
    assert len(registered) == len(set(registered)), "duplicate registry rows"
    # Every committed reference fixture must be registered with its budgets.
    committed = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / "validation/fixtures").glob("*.json"))
    assert sorted(registered) == committed, ("unregistered or missing fixtures", sorted(set(committed) ^ set(registered)))
    for row in registry["reference_sets"]:
        fixture = json.loads((ROOT / row["fixture"]).read_text())
        if row["budget_key"] is None:
            # The fixture declares no budget; the registered one lives in the
            # named Rust test, which must exist and still name the fixture.
            source, function = row["budget_source"].split("::")
            text = (ROOT / source).read_text()
            assert f"fn {function}(" in text, row["budget_source"]
            assert Path(row["fixture"]).name in text, (row["fixture"], source)
            assert "tolerances" not in fixture and "tolerance" not in fixture, row["fixture"]
            continue
        assert fixture[row["budget_key"]] == row["budgets"], row["fixture"]
    import hashlib
    process = json.loads((ROOT / "validation/fixtures/process-v1.json").read_text())
    for name, digest in process["input_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    assert {r["capability"] for r in registry["deferred"]} == {"critical-point/retrograde envelope continuation"}
    print(f"PASS: approved budgets exactly match all {len(registered)} committed reference fixtures, including bounded A13; critical/retrograde extensions deferred")


if __name__ == "__main__":
    main()
