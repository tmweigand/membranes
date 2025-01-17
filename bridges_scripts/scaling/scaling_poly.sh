#!/bin/bash

# Array of values
values=(1 2 4 8 16 32)

# Loop through the values
for val in "${values[@]}"; do
    echo "Processing value: $val"
    ./pack_poly_1_of_2.sh "$val"
    # Add your commands here that use $val
done
