import sys
import parallel_process_files
import os

def my_function(job_name):
    directory = '../data/test_data/*'
    if job_name is None:
        print('Error: Must pass through job batch name for initialize files')
        # directory_status_file = './process/directory_process_status.txt'
    else:
        #Creates job folder with given job_name
        job_folder = f'./{job_name}/process'
        os.makedirs(job_folder, exist_ok='True')

        #Creates directory process status file in created job_folder
        directory_status_file = os.path.join(job_folder, 'directory_process_status.txt')

    fd = parallel_process_files.FileDirectory(directory,directory_status_file)
    fd.initialize_list()
    
if __name__ == "__main__":
    if len(sys.argv) > 1:
        job_name = sys.argv[1]
    else:
        job_name = None
    my_function(job_name)
