# Preparing images for candidates selection

# Unpack the .fits.fz files.
sh unzip_fz_files.sh

# Mask the very negative values from the images.
echo "Removing negative pixels."
cd ../
python3 mask_negative_pixels.py

# Swarp the images from multi-extension to one mosaic.
cd shell_scripts/
sh swarp_files.sh

# Trim the images so they are the same size for source finding.
echo "Trimming images."
cd ../
python3 trim.py

# Run source finding
cd shell_scripts/
sh find_sources.sh

# Remove all sources that are not within the DECAM area (remove edge cases)
echo "Only selecting source that are within the DECam area."
cd ../
python3 regionfy_catalog.py

echo "Images are prepared for candidate search."
