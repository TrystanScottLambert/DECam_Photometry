# Script which will swarp all of the c4d.fits files into one mosiac.

echo "Swarping images using the default.swarp file in each directory."
echo "i-band:"
cd ../../correct_stacks/i/
swarp c4d_211021_003940_osj_i_vik1.fits

echo "z-band:"
cd ../z/
swarp c4d_210831_053503_osj_z_vik1.fits

echo "N964-band:"
cd ../N964/
swarp c4d_210831_050404_osj_N964_vik1.fits
