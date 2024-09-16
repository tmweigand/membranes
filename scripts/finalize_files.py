import sys
import parallel_process_files

def my_function(flag=None):

    directory = './data/sample_bridges_data/*'
    if flag is None:
        directory_status_file = './process/directory_process_status.txt'
    else:
        directory_status_file = f'./process/directory_process_status.{flag}.txt'
    df = parallel_process_files.FileDirectory(directory,directory_status_file)
    df.check_progress(flag)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        flag = sys.argv[1]
    else:
        flag = None
    my_function(flag)
