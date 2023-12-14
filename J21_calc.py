"""Calculating the J21 isotropic UV intensity at the Lyman limit"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.units.quantity import Quantity



def f_lyman_limit(ab_mag: float) -> Quantity:
	"""Calculates the flux density at the lyman limit given
	the extinction corrected AB magnitude at rest 1450."""

	f0 = 3631*u.Jy
	ratio = (912/1450)
	return f0*(10**(-0.4*ab_mag))*(ratio**0.99)

'''def f_lyman_limit(ab_mag: float) -> Quantity:
	"""Working out the lyman limit using Fan+ 2001"""
	return (10**((ab_mag + 48.6)/(-2.5)))*(u.erg /(u.cm**2)/u.s/u.Hz)'''

def f_q(ab_mag: float, redshift: float, radius_of_sphere: Quantity) -> Quantity:
	"""Calculates the Fq flux density using the ab mag at
	rest 1450 within a sphere of given radius"""

	cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
	lum_distance = cosmo.luminosity_distance(redshift)
	radius_mpc = radius_of_sphere.to(u.Mpc)
	f_lya = f_lyman_limit(ab_mag)
	numerator = (lum_distance**2)*f_lya
	denominator = (radius_mpc**2)
	janskies = numerator/denominator 
	return janskies.to(u.erg /(u.cm**2)/u.s/u.Hz)

def j_21(fq: Quantity) -> Quantity:
	val = (fq/(4*np.pi*u.sr))
	return val*1e21

if __name__ == '__main__':
	# Checking calulation is consistent
	#ab = 19.93
	#z = 4.874
	#flux_density = f_q(ab, z, 770*u.kpc)
	#print(flux_density)
	#print(j_21(flux_density))
	#print('-------------------------------------')
	ab = 21.17
	z = 6.9
	print('Stuff at 2Mpc')
	flux_density = f_q(ab, z, 2*u.Mpc)
	print(flux_density)
	print(j_21(flux_density))
	print('-------------------------------------')
	print('Stuff at 5Mpc')
	flux_density = f_q(ab, z, 5*u.Mpc)
	print(flux_density)
	print(j_21(flux_density))

	

