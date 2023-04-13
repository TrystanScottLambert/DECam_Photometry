#!/bin/bash
# Script which runs sextractor on the images

echo "Running sextractor on i-band."
cd ../../correct_stacks/N964/
sex n964.fits,i.fits -c default_i.sex

echo "Running sextractor on z-band."
sex n964.fits,z.fits -c default_z.sex

echo "Running sextractor on N964-band."
sex n964.fits -c default_n964.sex

echo "Running sextractor on N964-band with 1.35 apertures."
sex n964.fits -c default_n964_135.sex
