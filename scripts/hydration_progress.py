import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
#Process output files for hydration to check progress towards desired density
initial_file_path = 'hydration_outputs/hydration_output.out'
restart_file_path = 'hydration_outputs/hydration_restart_output.out'
key_line = '   Step          Temp          Press         v_press     v_actualdensity   c_totalvol  '

timesteps = []
densities = []
OPEN_INITIAL = False
if OPEN_INITIAL:
     with open(initial_file_path, 'r') as file:
        start_track = False
        skip_counter = 0
        
        for line in file:
            if key_line in line:
                start_track = True
                skip_counter = -1
                continue

            if start_track:
                skip_counter += 1
                if skip_counter % 5 == 0:
                    #Starts tracking every 5 lines after the header line is found
                    values = line.split()
                    if len(values) >= 5:
                        timesteps.append(values[0])
                        densities.append(values[4])
             
OPEN_RESTART = True
if OPEN_RESTART:
    with open(restart_file_path, 'r') as file:
        start_track = False
        skip_counter = 0
        
        for line in file:
            if key_line in line:
                start_track = True
                skip_counter = -1
                continue

            if start_track:
                skip_counter += 1
                if skip_counter % 4 == 0:
                    #Starts tracking every 4 lines after the header line is found
                    values = line.split()
                    if len(values) >= 5:
                        timesteps.append(values[0])
                        densities.append(values[4])

print("Timesteps: ", timesteps)
print("Density: ", densities)
timesteps = np.array(timesteps, dtype=np.float64)  # Convert to float or int if appropriate
densities = np.array(densities, dtype=np.float64)
#PLOT RESULTS
PLOT = True
if PLOT:
    plt.figure(figsize=(10,6))
    plt.plot(timesteps, densities, linestyle='-', linewidth=2,color='dodgerblue')
    plt.title('Water Density over Time During Hydration')
    plt.xlabel('Timestep')
    plt.ylabel('Density (g/cm^3)')
    plt.xlim(min(timesteps), max(timesteps))
    plt.ylim(min(densities), max(densities))
    plt.locator_params(axis='x', nbins=10)  # Reduce number of x-axis ticks
    plt.locator_params(axis='y', nbins=10)
    plt.tight_layout
    plt.savefig('hydration_outputs/density_plot.png')
                
