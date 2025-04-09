import cProfile
import pstats
import os
import sys
import generate_rdf


def main():
    # Create profiling output directory if it doesn't exist
    profile_dir = os.path.join("profiling", "output")
    os.makedirs(profile_dir, exist_ok=True)

    # Setup profiler
    profiler = cProfile.Profile()

    # Profile the main function
    profiler.enable()
    generate_rdf.generate_rdf(bridges=False)
    profiler.disable()

    # Save stats to file
    stats_path = os.path.join(profile_dir, "rdf_stats.prof")
    profiler.dump_stats(stats_path)

    # Print some basic stats to console
    stats = pstats.Stats(profiler)
    stats.sort_stats("cumulative")
    stats.print_stats(20)

    print(f"\nProfile data saved to: {stats_path}")
    print("To visualize with snakeviz, run:")
    print(f"snakeviz {stats_path}")


if __name__ == "__main__":
    main()
