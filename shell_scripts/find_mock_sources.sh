#!/bin/bash

# Script to copy the same config file used for the narrow band
# and then update it so that the setting are the same for running on
# the injected narrow band image.

# Change into the directory.
cd ../../correct_stacks/N964/

# Remove the the sex files if they exist.

rm default_n964_false_135.sex
rm default_n964_false.sex

# First copy the original .sex file with a new name.
cp default_n964.sex default_n964_false.sex
cp default_n964.sex default_n964_false_135.sex

# Update the name of the output catalog.
sed -i 's/n964.cat/n964_false.cat/' default_n964_false.sex

# Update the config file for 1.35"
# apertures = 5.2 pixels (pixel scale = 0.26"/pix)

sed -i 's/n964.cat/n964_false_135.cat/g' default_n964_false_135.sex
sed -i 's/7.4/5.2/g' default_n964_false_135.sex

# Run sextractor on both apertures.
echo "Running SEXTRACTOR"
sex n964.injected.fits -c default_n964_false.sex
sex n964.injected.fits -c default_n964_false_135.sex
