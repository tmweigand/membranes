import glob
import os
import sys
import time

def data_string(file,status,job_id):
    """
    Output for data in files
    
    File - Status - Time last moidified - Job ID
    
    """
    out_string = (f'{file.split("/")[-1]} \t'
                    f'{status} \t' 
                    f'{time.time()} \t'
                    f'{job_id}\t'
                    f'\n'
                )

    return out_string


def match_file(main_data,match):
    """
    Get the line number of file in main database
    """
    line_out = None
    for line_number, line in enumerate(main_data):
        if match in line:
            line_out = line_number
                
    return line_out

class FileDirectory:
    """
    TO make processing multiple files on a remote server 
    this create a list of jobs to perform. We first generate the datafiles
    to be processed. 
    """
    def __init__(self,directory,out_file):
        self.directory = directory
        self.directory_file = out_file
        self.files = glob.glob(self.directory)
        self.num_files = len(self.files)

    def initialize_list(self):
        """
        Initialize the list of files to process
        """
        os.path.isfile(self.directory_file)
        # if os.path.isfile(self.directory_file):
        #     print("File already exists!")
        #     sys.exit()

        out_file = open(self.directory_file, 'w+', encoding="utf-8")

        out_file.write(f'{self.num_files} \n')
        for file in self.files:
            out_file.write(data_string(file,0,0))
    
        out_file.close()

    def check_progress(self):
        """
        Ensure that all files have been processed
        """
        files = glob.glob('./process/*')
        files.remove(self.directory_file)
        main_counts = self.process_main(self.directory_file)
        main_data = self.gen_stats(self.directory_file,main = True)

        job_data = {}
        job_counts = {}
        for file in files:
            job_data.update(self.gen_stats(file,main = False))
            job_counts.update(self.process_job(file))

        for job,counts in job_counts.items():
            for status,status_count in counts.items():
                if status == 2:
                    if status_count != main_counts[job][status]:
                        print(f'Error! Main and Job {job} Counts for Status {status} differ!')
                        print(f'Main Count: {main_counts[job][status]:}')
                        print(f'Job {job}  Count: {status_count}')

        for job,data in job_data.items():
            for file,status in zip(data['files'],data['status']):
                if status == 2:
                    line_number = match_file(main_data[job]['files'],file)
                    main_status = main_data[job]['status'][line_number]
                    if main_status != status:
                        print(f'Error! Status is different for {file}')

    def process_main(self,file):
        """
        Check the main file to determine status
        """
        stats = self.gen_stats(file,main = True)
        counts = {}
        for id,job in stats.items():
            counts[id] = {0:0,1:0,2:0}
            for n in job['status']:
                counts[id][n] += 1

        return counts

    def process_job(self,file):
        """
        Coolect the information from the project files
        """
        stats = self.gen_stats(file,main = False)
        counts = {}
        for id,job in stats.items():
            counts[id] = {0:0,1:0,2:0}
            for n in job['status']:
                counts[id][n] += 1
        
        return counts

        

    def gen_stats(self,file,main = False):
        """
        Read in the main and ID files
        """

        stats = {}
        
        read_file = open(file, 'r', encoding="utf-8")
        
        if main:
            read_file.readline()
        
        lines = read_file.readlines()
        for line in lines:
            split = line.split('\t')

            file_name = split[0]
            status = int(split[1])
            id = split[3]

            if id not in stats:
                stats[id] = {'files':(),
                             'status':()
                            }
            
            stats[id]['files'] += (file_name,)
            stats[id]['status'] += (status,)

            
        return stats





class CurrentFile:
    def __init__(self,file,id):
        self.file = file
        self.id = id
        self.init = True
        self.write_time = os.path.getmtime(file)
        self.read_time = None
        self.num_files = None
        self.process_file = None
        self.line_number = None


    def select_file(self,status,attempts = 10):
        """
        Try to select file 
        """
        self.line_number = None
        for _ in range(0,attempts):
            
            self._select_file(status)
            
            if self.line_number is not None:
                
                assert(self.line_number <= self.num_files)
                file_found,file = self.replace_line(change_status = status + 1)
                
                if file_found:
                    self.process_file = file
                    break
        else:
            print("Attempting to read from file and failed....exiting")
            sys.exit()

    def read_file(self):
        """
        Open the file as read
        """
        self.read_time = time.time()
        read_file = open(self.file, 'r', encoding="utf-8")

        if self.init:
            self.get_num_files(read_file)
            self.init = False

        return read_file

    def get_num_files(self,read_file):
        """
        Get the total number of files listed
        """
        line = read_file.readline()
        try:
            self.num_files = int(line)
        except ValueError:
            print(f'First line {line} is not an int!')

    def _select_file(self,status):
        """
        Determine the next file of specfied status.
            0: To be processed
            1: Processing 
            2: Completed
        """
        read_file = self.read_file()
        line = read_file.readline()

        # Find file of given status
        for line_number in range(1,self.num_files+1):
            assert(line_number <= self.num_files)
            line = read_file.readline()
            split = line.split('\t')
            
            if len(split) != 5:
                continue
            else:
                _status = int(split[1])
                if status == _status:
                    mod_time = float(split[2])
                    assert( mod_time < self.read_time)
                    self.line_number = line_number
                    break

        read_file.close()

    def replace_line(self,change_status):
        """
        Replace the identified line in file
        """
        read_file = self.read_file()
        lines = read_file.readlines()
        read_file.close()

        write_complete = False
        process_file = None

        if self.line_number <= len(lines):
            split = lines[self.line_number].split('\t')
            process_file = split[0]
            old_status = int(split[1])

            if old_status != change_status:

                lines[self.line_number] = data_string(process_file,change_status,self.id)
                last_time = os.path.getmtime(self.file)
                
                if last_time <= self.read_time:
                    out = open(self.file, 'w', encoding="utf-8")
                    out.writelines(lines)
                    out.close()
                    write_complete = True
        
        return write_complete,process_file

    def close_file(self,status = 2,attempts = 20):
        """
        Change the status to be completed
        """
        write = False
        for _ in range(0,attempts):
            write,_ = self.replace_line(change_status = status)
            if write:
                break
        return write

