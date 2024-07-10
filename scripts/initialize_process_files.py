import glob
import os
import sys
import time

lammps_in_files = glob.glob('./bridges_data/*')
num_files = len(lammps_in_files)

directory_status_file = 'director_process_status.txt'
os.path.isfile(directory_status_file)

# if os.path.isfile(directory_status_file):
#     print("File alreay exists!")
#     sys.exit()

out_file = open(directory_status_file, 'w+', encoding="utf-8")

out_file.write(f'{num_files} \n')
for file in lammps_in_files:
    out_string = f'{file.split("/")[-1]},{0},{time.time()} \n'
    out_file.write(out_string)
    
out_file.close()