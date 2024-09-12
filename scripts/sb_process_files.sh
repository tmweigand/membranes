#!/bin/bash

#checks if initialize_files.py has already been ran once before
if [ -f "../process/director_process_status.txt" ]; then
    echo "Initialization already complete."
else
    python initialize_files.py

    if [ $? -ne 0 ]; then
        echo "Error Initializing"
        exit 1
    fi
fi

# Run process files with slurm job id as identifier argument 
python process_files.py JOB.%j

if [ $? -eq 0 ]; then
    echo "Processing Complete"
else
    echo "Error Processing"
    exit 1
fi