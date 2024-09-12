import glob
import sys
import random
from time import sleep
from mpi4py import MPI
import parallel_process_files

def my_function(id):

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    directory = '../data/test_data/*'
    print(directory)
    directory_status_file = 'process/directory_process_status.txt'

    process_completed_file = "process/process_" + id + "_.txt"
    out_process = open(process_completed_file, 'w+', encoding="utf-8")

    if rank == 0:
        f = parallel_process_files.CurrentFile(directory_status_file,id)
        while (f.process_file is not None) or (f.init is True):
            f.select_file(0)
            out_process.write(parallel_process_files.data_string(f.process_file,1,id))
            
            # Perform Job
            sleep(random.random()*2)

            success = f.close_file()
            out_process.write(parallel_process_files.data_string(f.process_file,2,id))


    out_process.close()


if __name__ == "__main__":
    my_function(sys.argv[-1])
    MPI.Finalize()
