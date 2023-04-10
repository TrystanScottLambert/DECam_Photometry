#!/bin/bash

# Script to copy the same config file used for the narrow band
# and then update it so that the setting are the same for running on
# the injected narrow band image.

DATA_PATH='../../correct_stacks/N964/'
IMAGE='n964.injected.fits'
NORMAL_CAT='default_n964.sex'
NEW_CAT='default_n964_false.sex'

# First copy the original .sex file with a new name.
cp $DATA_PATH$NORMAL_CAT $DATA_PATH$NEW_CAT

# Update the name of the output catalog.
cd $DATA_PATH
sed -i 's/n964.cat/n964_false.cat/' $NEW_CAT

# run sextractor on the image
sex $IMAGE -c $NEW_CAT

# Update the config file for 1.35"
# apertures = 5.2 pixels (pixel scale = 0.26"/pix)

sed -i 's/7.4/5.2/' $NEW_CAT
sed -i 's/n964_false.cat/n964_false_135.cat/' $NEW_CAT

# Run sextractor for 1.35" apertures.
sex $IMAGE -c $NEW_CAT

# Remove the clutter
rm $NEW_CAT

