#!/bin/bash

# Script which will swarp all of the c4d.fits files into one mosiac.

echo "Swarping images using the default.swarp file in each directory."
echo "i-band:"
cd ../../correct_stacks/i/
weight_name='c4d*w*.fits'
swarp c4d_*j_i_*.fits -weight_image $weight_name

echo "z-band:"
cd ../z/
weight_name='c4d*w*.fits'
swarp c4d_*j_z_*.fits -weight_image $weight_name

echo "N964-band:"
cd ../N964/
weight_name='c4d*w*.fits'
swarp c4d_*j_N964_*.fits  -weight_image $weight_name
