import numpy as np
import pmmoto


def generate_bounded_rdf():
    """
    Generate radii
    """
    atom_folder = "rdf/bridges_results/bins_extended/"
    atom_map, rdf = pmmoto.io.data_read.read_binned_distances_rdf(atom_folder)

    bounded_rdf = {}
    for _id, _rdf in rdf.items():
        bounded_rdf[_id] = pmmoto.domain_generation.rdf.Bounded_RDF.from_rdf(
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
        print(atom_id, bounded_rdf[atom_id].name, radii[atom_id])

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


if __name__ == "__main__":
    bounded_rdf = generate_bounded_rdf()
    uff_radii = determine_uff_radii()
    # pmf_radii = determine_radii(bounded_rdf, 1.40)
    # pmf_radii2 = determine_radii(bounded_rdf, 3.61)
    pmf_radii2 = determine_radii(bounded_rdf, 17.315)

    # print(pmf_radii)

    # sum = 0
    # sum2 = 0
    # count = 0
    # for key in uff_radii.keys() & pmf_radii.keys() & pmf_radii2.keys():
    #     u_r = uff_radii[key]
    #     p_r = pmf_radii[key]
    #     p_r2 = pmf_radii2[key]
    #     print((key, u_r, p_r, p_r2))
    #     sum += abs(u_r - p_r)
    #     sum2 += abs(p_r - p_r2)
    #     count += 1

    # print(f"Average difference from equilibrium: {sum/count}")
    # print(f"Difference between G values : {sum2/count}")
