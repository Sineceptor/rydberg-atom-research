"""Recompute a selected study without replacing delivered results."""
import argparse
from pathlib import Path
from run_support import run_managed


def calculate(output_root, section="all"):
    summaries = {}
    if section in ("all", "foundations"):
        import run_august_experiments
        import phase_diagram_bridge
        summaries["august"] = run_august_experiments.run(output_root)
        summaries["bridge"] = phase_diagram_bridge.run(output_root)
    if section in ("all", "chains"):
        import run_extension
        import parity_diagnostics
        summaries["chains"] = run_extension.run(output_root)
        summaries["parity"] = parity_diagnostics.run(output_root)
    if section in ("all", "optional"):
        import kagome_pilot
        import kagome_shape_control
        summaries["kagome"] = kagome_pilot.run(output_root)
        summaries["kagome_shapes"] = kagome_shape_control.run(output_root)
    import circuit_controls
    summaries["circuits"] = circuit_controls.run(output_root)
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--section", choices=["all", "foundations", "chains", "optional"], default="all")
    args = parser.parse_args()
    run_managed(lambda destination: calculate(destination, args.section), args.section)


if __name__ == "__main__":
    main()
