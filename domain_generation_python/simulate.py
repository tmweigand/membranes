"""simulate.py

Unified entry point for the full membrane simulation pipeline.
Select a stage via the --stage argument:

    python simulate.py --stage polymerize
    python simulate.py --stage equilibrate
    python simulate.py --stage hydrate

Universal parameters (shared across stages) are defined once in UNIVERSAL_PARAMS.
Stage-specific parameters are in their respective PARAMS dicts below.
"""

import argparse

from membranes.domain_generation.polymerization import PolymerizationSimulation
from membranes.domain_generation.equilibration import EquilibrationSimulation
from membranes.domain_generation.hydration import HydrationSimulation

# ──────────────────────────────────────────────────────────────────────
# Universal parameters (shared across all stages)
# ──────────────────────────────────────────────────────────────────────

UNIVERSAL_PARAMS = dict(
    in_dir="rv",  # Force field / input directory key ("rv" or other)
    multiple=0.1,  # System size scaling factor
)

# ──────────────────────────────────────────────────────────────────────
# Stage-specific parameters
# ──────────────────────────────────────────────────────────────────────

POLYMERIZE_PARAMS = dict(
    **UNIVERSAL_PARAMS,
    xlink=0.80,  # Target cross-linking degree (0–1)
    temperature=300,  # Simulation temperature (K)
    bond_frequency=50,  # Steps between bond/react attempts
    bond_distance=(0.0, 5.0),  # (min, max) bond formation distance (Å)
    stabilization=0.03,  # Bond/react stabilization factor
    max_cycles=375,  # Maximum NVT/NPT cycles per polymerization stage
)

EQUILIBRATE_PARAMS = dict(
    **UNIVERSAL_PARAMS,
    input_data="logs/term_final.lmps",  # Output from polymerization
    rho_target=1.24,  # Experimental PA density (g/cm³)
    nsteps=500,  # Base MD steps per stage
    initial_seed=58447419,  # Velocity seed for first NVT run
)

HYDRATE_PARAMS = dict(
    **UNIVERSAL_PARAMS,
    mult=0.1,  # System size scaling factor (used by HydrationSimulation)
    rand=1,  # Random seed multiplier
    input_data="logs/equil_polymer.lmps",  # Output from equilibration
    feed_pressure_atm=0.5,  # Feed-side applied pressure (atm)
    perm_pressure_atm=0.5,  # Permeate-side applied pressure (atm)
    hydration_steps=3000,  # Steps for initial hydration run
    production_steps=5000,  # Steps for production run
)

# ──────────────────────────────────────────────────────────────────────
# Stage runners
# ──────────────────────────────────────────────────────────────────────


def run_polymerize():
    sim = PolymerizationSimulation(**POLYMERIZE_PARAMS)

    print("\n=== Stage 1: Packing molecules ===")
    sim.pack_molecules()

    print("\n=== Stage 2: Minimization ===")
    sim.minimize()

    print("\n=== Stage 3: PA polymerization ===")
    max_bonds = sim.polymerize_pa()

    print("\n=== Stage 4: Adding hydroxide ===")
    sim.add_hydroxide(max_bonds)

    print("\n=== Stage 5: OH polymerization ===")
    sim.polymerize_oh()

    print("\n=== Stage 6: Cleanup ===")
    sim.cleanup()

    print("\nPolymerization complete.")


def run_equilibrate():
    sim = EquilibrationSimulation(**EQUILIBRATE_PARAMS)

    print("\n=== Stage 1: Load and prepare ===")
    sim.load()

    print("\n=== Stage 2: Density equilibration loop ===")
    sim.equilibrate()

    print("\n=== Stage 3: Final equilibration ===")
    sim.finalize()

    print("\nEquilibration complete.")


def run_hydrate():
    sim = HydrationSimulation(**HYDRATE_PARAMS)

    print("\n=== Stage 1: Load and prepare ===")
    sim.load()

    print("\n=== Stage 2: Unwrap membrane ===")
    sim.unwrap()

    print("\n=== Stage 3: Add water ===")
    sim.add_water()

    print("\n=== Stage 4: Add pistons ===")
    sim.add_pistons()

    print("\n=== Stage 5: Run hydration ===")
    sim.run_hydration()

    print("\nHydration simulation complete.")


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────

STAGES = {
    "polymerize": run_polymerize,
    "equilibrate": run_equilibrate,
    "hydrate": run_hydrate,
}

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run a membrane simulation stage.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(f"  {s}" for s in STAGES),
    )
    parser.add_argument(
        "--stage",
        choices=STAGES,
        required=True,
        metavar="STAGE",
        help="Simulation stage to run: %(choices)s",
    )
    args = parser.parse_args()
    STAGES[args.stage]()
