#!/bin/bash
# Preparing images for candidates selection

# Trim the images so they are the same size for source finding.
echo "Trimming images."
cd ../
python3 trim.py

# make fake sources and determine the completeness.
python3 inject_false_sources.py

# Run source finding
cd shell_scripts/
sh find_sources.sh		# Sources from DECAM.
sh find_mock_sources.sh		# Sources from DECAM + Injected images.

# Remove all sources that are not within the DECAM area (remove edge cases).
echo "Only selecting source that are within the DECam area."
cd ../
python3 regionfy_catalog.py 
#python3 mask_sex_catalog.py

echo "Images are prepared for candidate search."
