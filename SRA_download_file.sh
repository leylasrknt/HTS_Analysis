#!/bin/bash

# Path to the directory where you want to store the downloaded data
experimentDir="/path/to/your/data_SRA"

# Replace 'YOUR_SRA_ID' with the specific SRA experiment ID you want to download
experimentID="YOUR_SRA_ID"

# Path to the prefetch tool (SRA Toolkit)
prefetch="/path/to/sratoolkit.3.0.6-ubuntu64/bin/prefetch"

# Download the data from SRA using prefetch
echo "Downloading data from SRA..."
$prefetch $experimentID --output-directory $experimentDir

# Check if the download was successful
if [ $? -eq 0 ]; then
    echo "Download completed successfully."
else
    echo "Error: Download failed. Please check your SRA experiment ID and network connection."
fi
