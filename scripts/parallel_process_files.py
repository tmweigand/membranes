import glob
import os
import sys
import time
from mpi4py import MPI

class CurrentFile:
    def __init__(self,file):
        self.file = file
        self.read_time = time.time()
        self.write_time = os.path.getmtime(file)


    def select_file(self,status,attempts = 10):
        """
        Select file unti no error
        """
        result = None
        for n in range(0,attempts):
            try:
                # connect
                result = self._select_file(status)
            except:
                pass
        # other code that uses result but is not involved in getting it

    def _select_file(self,status):
        """
        Determine the next file of specfied status.
            0: To be processed
            1: Processing 
            2: Completed
        """
        process_files = open(self.file, 'r', encoding="utf-8")
        line = process_files.readline()
        num_files = int(line)
        for n in range(0,num_files):

            line = process_files.readline()
            split = line.split(",")
            
            # next_file = split[0]
            _status = int(split[1])

            mod_time = float(split[2])
            # assert( mod_time > self.read_time)
            
            if status == _status:
                found = True
                print("Success")
                output = n
                break
        else:
            print("Errror")
            output = None

        process_files.close()   
        return output

def replace_line(file_name, line_num, next_file, read_time):
    
    lines = open(file_name, 'r', encoding="utf-8").readlines()
    split = lines[line_num].split(",")

    out_value = f'{next_file},{1},{time.time()} \n'
    lines[line_num] = out_value

    last_time = os.path.getmtime(file_name)
    
    if last_time != read_time:
        print("Error File Changed. Try again.")
    
    out = open(file_name, 'w', encoding="utf-8")
    out.writelines(lines)
    out.close()




def my_function():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    lammps_in_files = glob.glob('./bridges_data/*')
    directory_status_file = 'director_process_status.txt'

    if rank == 0:
        f = CurrentFile(directory_status_file)
        f.select_file(0)
        # next_file,line_number,read_time = select_file(directory_status_file)

        # replace_line(directory_status_file, line_number, next_file, read_time)




        
if __name__ == "__main__":
    my_function()
    MPI.Finalize()