import sys
import parallel_process_files

def my_function():

    directory = './data/sample_bridges_data/*'
    directory_status_file = './process/directory_process_status.txt'
    df = parallel_process_files.FileDirectory(directory,directory_status_file)
    df.check_progress()

if __name__ == "__main__":
    my_function()
