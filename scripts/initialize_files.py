import sys
import parallel_process_files

def my_function():

    directory = '../data/test_data/*'
    directory_status_file = './process/director_process_status.txt'
    fd = parallel_process_files.FileDirectory(directory,directory_status_file)
    fd.initialize_list()
    
if __name__ == "__main__":
    my_function()
