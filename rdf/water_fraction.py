"""membrane_domain.py"""

import glob
import numpy as np
from mpi4py import MPI
import pmmoto
import rdf_helpers
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc_size = comm.Get_size()

logger = pmmoto.logger

# This maps from from the lammps input file id and charge
# to a unique id. Atom_map.txt relates ids to atom types
atom_id_charge_map = {
    (1, 0.6797): 1,
    (1, 0.743425): 2,
    (3, -0.23): 3,
    (3, -0.1956): 4,
    (3, -0.1565): 5,
    (3, 0.014): 6,
    (3, 0.1716): 7,
    (4, -0.587509): 8,
    (5, 0.10745): 9,
    (5, 0.131): 10,
    (5, 0.1816): 11,
    (7, -0.4621): 12,
    (7, -0.398375): 13,
    (8, 0.23105): 14,
    (12, -0.5351): 15,
    (14, 0.4315): 16,
}


def initialize_domain(voxels):
    """
    Initialize the membrane domain
    """

    if proc_size == 8:
        subdomains = (2, 2, 2)
    elif proc_size == 27:
        subdomains = (3, 3, 3)
    elif proc_size == 64:
        subdomains = (4, 4, 4)
    elif proc_size == 125:
        subdomains = (5, 5, 5)
    elif proc_size == 216:
        subdomains = (6, 6, 6)
    else:
        print("oops")
        subdomains = (0, 0, 0)

    # Full domain with reservoirs
    box = [
        [0.0, 176],
        [0.0, 176],
        # [-35, 65],
        # [-287, 237],
        [-100, 100],
    ]

    sd = pmmoto.initialize(
        voxels=voxels,
        box=box,
        rank=rank,
        subdomains=subdomains,
        boundary_types=(
            (pmmoto.BoundaryType.PERIODIC, pmmoto.BoundaryType.PERIODIC),
            (pmmoto.BoundaryType.PERIODIC, pmmoto.BoundaryType.PERIODIC),
            (pmmoto.BoundaryType.END, pmmoto.BoundaryType.END),
        ),
        verlet_domains=(20, 20, 20),
        inlet=((0, 0), (0, 0), (1, 0)),
        outlet=((0, 0), (0, 0), (0, 1)),
    )

    return sd


def generate_bounded_rdf():
    """
    Generate radii
    """
    atom_folder = "rdf/bridges_results/bins_extended/"
    atom_map, rdf = pmmoto.io.data_read.read_binned_distances_rdf(atom_folder)

    bounded_rdf = {}
    for _id, _rdf in rdf.items():
        bounded_rdf[_id] = pmmoto.domain_generation.rdf.BoundedRDF.from_rdf(
            _rdf, 1.0e-3
        )

    return bounded_rdf


def determine_radii(bounded_rdf, pmf_value):
    """
    Collect the radii given a pmf cutoff
    """
    radii = {}
    for atom_id, _rdf in bounded_rdf.items():
        radii[atom_id] = _rdf.interpolate_radius_from_pmf(pmf_value)

    return radii


def determine_uff_radii():
    """
    Collect the radii given a pmf cutoff
    """
    atom_folder = "rdf/bridges_results/bins_extended/"
    atom_map, _ = pmmoto.io.data_read.read_binned_distances_rdf(atom_folder)
    radii = {}
    for atom_id, atom_data in atom_map.items():
        radii[atom_id] = (
            list(pmmoto.particles.uff_radius(atom_names=atom_data["element"]).values())[
                0
            ]
            + 1.4
        )

    return radii


def generate_membrane_domain(pmf_value, subdomain, membrane_file):
    """
    Test for generating a radial distribution function from LAMMPS data
    """
    bounded_rdf = generate_bounded_rdf()

    radii = determine_radii(bounded_rdf, pmf_value)

    pm = pmmoto.domain_generation.gen_pm_atom_file(
        subdomain=subdomain,
        lammps_file=membrane_file,
        atom_radii=radii,
        type_map=atom_id_charge_map,
        kd=False,
        add_periodic=True,
    )

    pm_morph = pmmoto.filters.morphological_operators.dilate(subdomain, pm.img, 1.4)

    porosity = pmmoto.analysis.average.average_image_along_axis(
        subdomain, pm_morph, (0, 1)
    )

    cc, _ = pmmoto.filters.connected_components.connect_components(
        img=pm.img, subdomain=subdomain
    )

    connections = pmmoto.filters.connected_components.inlet_outlet_connections(
        subdomain=subdomain, labeled_img=cc
    )

    if not connections:
        return porosity, np.zeros_like(porosity), None

    connected = np.where(cc == connections[0], 1, 0)

    conected_morph = pmmoto.filters.morphological_operators.dilate(
        subdomain, connected, 1.4
    )

    connected_porosity = pmmoto.analysis.average.average_image_along_axis(
        subdomain, conected_morph, (0, 1)
    )

    return porosity, connected_porosity, conected_morph


def generate_water_domain(water_radii, subdomain, water_file):
    """
    Test for generating a radial distribution function from LAMMPS data
    """
    pm = pmmoto.domain_generation.gen_pm_atom_file(
        subdomain=subdomain,
        lammps_file=water_file,
        atom_radii=water_radii,
        kd=False,
        add_periodic=True,
    )

    porosity = pmmoto.analysis.average.average_image_along_axis(
        subdomain, pm.img, (0, 1)
    )

    # pmmoto.io.output.save_img_data_parallel(
    #     "data_out/water_domain",
    #     subdomain,
    #     pm.img,  # additional_img={"edt": _edt}
    # )

    return porosity


def count_water_locations(subdomain, water_file, img=None):
    """
    Count the water atoms
    """
    water_positions, water_atom_type, _, _ = pmmoto.io.data_read.read_lammps_atoms(
        water_file
    )

    water_radii = {15: 1.4, 16: 0.6}

    water = pmmoto.particles.initialize_atoms(
        subdomain,
        water_positions,
        water_radii,
        water_atom_type,
        by_type=True,
        trim_within=True,
    )

    # hydrogen = water.return_list(16)
    oxygen = water.return_list(15)

    num_bins = int(subdomain.domain.voxels[2] / 10)
    start = subdomain.domain.box[2][0]
    end = subdomain.domain.box[2][1]

    water_bins = pmmoto.analysis.bins.Bin(start, end, num_bins, name="oxygen")

    coords = oxygen.return_coordinates()

    if img is not None:
        zeros = np.zeros((coords.shape[0], 1))
        coords = np.hstack((coords, zeros))

        for n, atom in enumerate(coords):
            index = subdomain.get_img_index(atom[0:3])
            if index and img[index] == 1:
                coords[n, 3] = 1

        coords = coords[coords[:, 3] == 1]

    pmmoto.analysis.bins.count_locations(
        coordinates=coords,
        dimension=2,
        bin=water_bins,
    )

    return water_bins


def save_water_locations(subdomain, water_file, img=None):
    """
    Save the locations of water molecules
    """
    water_positions, water_atom_type, _, _ = pmmoto.io.data_read.read_lammps_atoms(
        water_file
    )

    water_radii = {15: 1.4, 16: 0.6}

    water = pmmoto.particles.initialize_atoms(
        subdomain,
        water_positions,
        water_radii,
        water_atom_type,
        by_type=True,
        trim_within=True,
    )

    hydrogen = water.return_list(16)
    oxygen = water.return_list(15)

    save_based_on_img("data_out/water_oxygen", subdomain, img, oxygen)
    save_based_on_img("data_out/water_hydrogen", subdomain, img, hydrogen)


def save_based_on_img(file_name, subdomain, img, particle):

    particle_array = particle.return_np_array()

    particle_out = np.zeros_like(particle_array)
    particle_out[:, 0:3] = particle_array[:, 0:3]

    if img is not None:
        for n, atom in enumerate(particle_array):
            index = subdomain.get_img_index(atom[0:3])
            if index and img[index] == 1:
                particle_out[n, 3] = 1

    particle_out = particle_out[particle_out[:, 3] == 1]

    pmmoto.io.output.save_particle_data(file_name, subdomain, particle_out)


if __name__ == "__main__":

    bridges = False

    # Grab Files
    if bridges:
        voxels_in = (3520, 3520, 4000)
        membrane_files, _ = rdf_helpers.get_bridges_files()
    else:
        voxels_in = (1500, 1500, 1500)
        membrane_files = glob.glob("data/membrane_data/membranedata.100030000")
        water_files = glob.glob("data/water_data/pressuredata.100030000")

    sd = initialize_domain(voxels_in)
    _water_file = water_files[0]
    _membrane_file = membrane_files[0]

    upper_pmf_data = 17.315

    pmf = upper_pmf_data  # 1.40  # 3.61  # upper_pmf_data

    if sd.rank == 0:
        fig, ax1 = plt.subplots(figsize=(10, 5))

    porosity, connected_porosity, connect_img = generate_membrane_domain(
        pmf, sd, _membrane_file
    )

    if sd.rank == 0:
        coords = sd.domain.get_coords(
            sd.domain.box, sd.domain.voxels, sd.domain.resolution
        )

        print(
            "Porosity entire domain",
            porosity.shape,
            np.mean(porosity),
            np.mean(connected_porosity),
        )
        data = np.column_stack((coords[2], porosity, connected_porosity))
        np.savetxt(
            f"porosity_output_{pmf}.txt",
            data,
            fmt="%.6f",
            delimiter="\t",
            header="A\tB\tC",
            comments="",
        )

    if sd.rank == 0:
        ax1.plot(coords[2], porosity, "k", label="Total Pore Space")
        ax1.plot(
            coords[2],
            connected_porosity,
            "b",
            label="Connected Pore Space",
        )

    water_radii = {15: 1.4, 16: 0.65}
    # water_fraction = generate_water_domain(water_radii, sd, _water_file)
    water_bins = count_water_locations(sd, _water_file, img=None)
    if np.mean(connected_porosity) > 1.0e-6:
        connected_water_bins = count_water_locations(sd, _water_file, img=connect_img)

    if sd.rank == 0:
        # ax1.plot(sd.domain.coordinates[2], 1 - water_fraction, label="WATER")
        # plt.legend()
        bin_volume = 176 * 176 * water_bins.width
        volume_water = 29.7  # Angstroms cubed
        ax1.set_ylim(bottom=0.0)

        ax1.set_xlabel("z-coordinate (Å)", fontsize=16)
        ax1.set_ylabel("Volume Fraction", fontsize=16)

        ax1.plot(
            water_bins.centers,
            water_bins.values * volume_water / bin_volume,
            "k",
            linestyle="dashed",
            label="Total Water",
        )

        ax1.plot(
            connected_water_bins.centers,
            connected_water_bins.values * volume_water / bin_volume,
            "b",
            linestyle="dashed",
            label="Connected Water",
        )

        ax1.grid(True)

        ax1.tick_params(axis="x", labelsize=16)
        ax1.tick_params(axis="y", labelsize=16)

        ax1.axvline(x=-35, color="gray", linestyle="dotted", linewidth=1.5)
        ax1.axvline(x=60, color="gray", linestyle="dotted", linewidth=1.5)

        plt.legend(loc="lower left", fontsize=12)
        plt.savefig(
            "data_out/water_plotsss.pdf",
            dpi=300,
            bbox_inches="tight",
        )

    if sd.rank == 0:
        bin_volume = 176 * 176 * water_bins.width
        volume_water = 29.7  # Angstroms cubed
        if np.mean(connected_porosity) > 1.0e-6:
            print(
                "Water domain average",
                water_bins.values.shape,
                np.mean(water_bins.values * volume_water / bin_volume),
                np.mean(connected_water_bins.values * volume_water / bin_volume),
            )
        else:
            print(
                "Water domain average",
                np.mean(water_bins.values * volume_water / bin_volume),
                0.0,
            )
        if np.mean(connected_porosity) > 1.0e-6:
            data = np.column_stack(
                (
                    water_bins.centers,
                    water_bins.values * volume_water / bin_volume,
                    connected_water_bins.values * volume_water / bin_volume,
                )
            )
            np.savetxt(
                f"water_output_{pmf}.txt",
                data,
                fmt="%.6f",
                delimiter="\t",
                header="A\tB\tC",
                comments="",
            )
