import sys
import parallel_process_files

def my_function(job_name):

    directory = './data/sample_bridges_data/*'
    if job_name is None:
        print('Error: Must pass through job_name that has already been initialized')
        directory_status_file = './process/directory_process_status.txt'
    else:
        directory_status_file = f'./{job_name}/process/directory_process_status.txt'
        
    df = parallel_process_files.FileDirectory(directory,directory_status_file)
    df.check_progress(job_name)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        job_name = sys.argv[1]
    else:
        job_name = None
    my_function(job_name)
