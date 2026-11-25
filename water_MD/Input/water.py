#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

# Simulation settings
temperature = 300
equil_steps = 250000
production_steps = 1000000
dump_frequency = 1000

input_dir = os.path.dirname(os.path.abspath(__file__))

forcefields = {
    "riley": {
        "pair_style": "lj/charmm/coul/long",
        "mass_o": 16.000000,
        "mass_h": 1.008000,
        "epsilon_o": 0.155865,
        "sigma_o": 3.178,
        "epsilon_h": 0.0,
        "sigma_h": 0.0,
        "bond_k": 450.0,
        "bond_r0": 1.0118,
        "angle_k": 55.0,
        "angle_theta": 108.15,
    },
    "tip3p_orig": {
        "pair_style": "lj/cut/coul/long",
        "mass_o": 15.9994,
        "mass_h": 1.00800,
        "epsilon_o": 0.1521,
        "sigma_o": 3.1507,
        "epsilon_h": 0.0,
        "sigma_h": 1.0,
        "bond_k": 450.0,
        "bond_r0": 0.9572,
        "angle_k": 55.0,
        "angle_theta": 104.52,
    },
}

parser = argparse.ArgumentParser()
parser.add_argument(
    "--forcefield",
    required=True,
    choices=("riley", "tip3p_orig"),
)
parser.add_argument("--pressure", required=True, type=float)
parser.add_argument("--output", required=True)
args = parser.parse_args()

# Select forcefield parameters
parameters = forcefields[args.forcefield]

# Prepare output directory
output_dir = os.path.abspath(args.output)
os.makedirs(output_dir, exist_ok=True)

data_file = os.path.join(output_dir, "water.data")
in_file = os.path.join(input_dir, "in.water")
data_generator = os.path.join(input_dir, "make_water_data.py")

# Convert pressure from bar to atm,
pressure_atm = args.pressure * 0.986923

# Generate the initial water data file
subprocess.run(
    [
        sys.executable,
        data_generator,
        "--forcefield",
        args.forcefield,
        "--output",
        data_file,
    ],
    check=True,
)

# LAMMPS variables
variables = {
    "forcefield": args.forcefield,
    "output": output_dir,
    "data": data_file,
    "pair_style": parameters["pair_style"],
    "mass_o": parameters["mass_o"],
    "mass_h": parameters["mass_h"],
    "epsilon_o": parameters["epsilon_o"],
    "sigma_o": parameters["sigma_o"],
    "epsilon_h": parameters["epsilon_h"],
    "sigma_h": parameters["sigma_h"],
    "bond_k": parameters["bond_k"],
    "bond_r0": parameters["bond_r0"],
    "angle_k": parameters["angle_k"],
    "angle_theta": parameters["angle_theta"],
    "temperature": temperature,
    "pressure_atm": pressure_atm,
    "equil_steps": equil_steps,
    "production_steps": production_steps,
    "dump_frequency": dump_frequency,
}

command = ["srun", "lmp", "-in", in_file]

for key, value in variables.items():
    command.extend(["-var", key, str(value)])

subprocess.run(command, check=True, cwd=output_dir)
