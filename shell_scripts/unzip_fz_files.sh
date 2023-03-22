# Script to unzip the .fz files.

# First, remove the current files (which may have been edidted).
echo "Removing current c4d.fits files"
rm ../../correct_stacks/i/c4d*.fits
rm ../../correct_stacks/z/c4d*.fits
rm ../../correct_stacks/N964/c4d*.fits

# Unpack all of the .fz files
echo "Unpacking all fits.fz files."
funpack ../../correct_stacks/i/*.fits.fz
funpack ../../correct_stacks/z/*.fits.fz
funpack ../../correct_stacks/N964/*.fits.fz
