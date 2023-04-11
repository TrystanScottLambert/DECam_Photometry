#!/bin/bash
#test script

cd ../../correct_stacks/N964/
sex n964.fits -c default_n964.sex
sex n964.injected.fits -c default_n964_false.sex
