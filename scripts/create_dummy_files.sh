#!/bin/bash

# Number of files to create
n=100 # Change this to however many files you want to create

# Loop to create files
mkdir ../data/test_data
for i in $(seq 1 $n)
do
  cd ../data/test_data/
  touch "file_$i.txt"
done