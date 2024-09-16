import sys
import parallel_process_files

def my_function(name=None):
    directory = '../data/test_data/*'
    if name is None:
        directory_status_file = './process/directory_process_status.txt'
    else:
        directory_status_file = f'./process/directory_process_status.{name}.txt'
    fd = parallel_process_files.FileDirectory(directory,directory_status_file)
    fd.initialize_list()
    
if __name__ == "__main__":
    if len(sys.argv) > 1:
        name = sys.argv[1]
    else:
        name = None
    my_function(name)
