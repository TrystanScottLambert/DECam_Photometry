"""
Chiara plot but looking at the densities of the galaxies. 
"""

from typing import Literal
from dataclasses import dataclass
from astropy.units.quantity import Quantity
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import pylab as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from plotting import start_plot, end_plot


@dataclass
class Observation:
    """All the data related to observations of QSO"""
    qso_name: str
    qso_redshift: float
    survey_area: Quantity
    overdensity_factor: float
    survey_type: Literal['LAE', 'LBG']

    @property
    def comoving_area(self) -> Quantity:
        """The on sky area of the survey in co-moving units"""
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        arc_kpc_squared = cosmo.arcsec_per_kpc_comoving(self.qso_redshift)**2
        area_cmpc = self.survey_area/arc_kpc_squared
        return area_cmpc.to(u.Mpc**2)

# Inputs
stiavelli_2005 = Observation('SDSS_J1030+0524', 6.28, 11.3*(u.arcmin**2), 2, 'LBG')
husband_2013_2130 = Observation('SDSS_J2130+0026', 4.951, 50*(u.arcmin**2), 2, 'LBG')
morselli_2014_1048 = Observation('SDSS_J1048+4637', 6.2, 0.14*(u.deg**2), 2.09, 'LBG')
champagne_2023_2348 = Observation('SDSS_2348-3054', 6.9, 7.29*(u.arcmin**2), 2.3, 'LBG')
morselli_2014_1148 = Observation('SDSS J1148+5251', 6.41, 0.14*(u.deg**2), 2.32, 'LBG')
morselli_2014_1411 = Observation('SDSS_J1411+1217', 5.95, 0.14 * (u.deg**2), 2.79, 'LBG')
champagne_2023_2054 = Observation('SDSS_J2054-0005', 6.04, 7.29*(u.arcmin**2), 2.9, 'LBG')
balmaverde_2017 = Observation('SDSS J1030+0524', 6.28, 0.17 * (u.deg**2), 3.4, 'LBG')
morselli_2014_1030 = Observation('SDSS_J1030+0524', 6.28, 0.14 * (u.deg**2), 3.72, 'LBG')
husband_2013_1204 = Observation('SDSS_J1204-0021', 5.086, 49*(u.arcmin**2), 4, 'LBG')
husband_2013_0338 = Observation('SDSS_J0338+0021', 5.027, 49*(u.arcmin**2), 6, 'LBG')
zheng_2006 = Observation('SDSS_J0836+0054', 5.8, 5*(u.arcmin**2), 6, 'LBG')
utsumi_2010 = Observation('CFHQS_J2329-0301', 6.43, 0.219*(u.deg**2), 7, 'LBG')
champagne_2023_3051 = Observation('VIK_J030516.92-315056.0', 6.6, 7.9*(u.arcmin**2), 9.8, 'LBG')
lambert_2023 = Observation('VIK_2348-3054', 6.9, 2.87*(u.deg**2), 11, 'LBG')
kim_2009_1030 = Observation('SDSS J1030+0524', 6.28, 11.3 * (u.arcmin**2), 1.75, 'LBG')
kim_2009_1630 = Observation('SDSS J1630+4012', 6.05, 11.3 * (u.arcmin**2), 1.375, 'LBG')


goto_2017 = Observation('CFHQS_J2329-0301', 6.43, 0.219 * (u.deg**2), 1, 'LBG')
muzzucchelli_2017 = Observation('PSO_J2329-1601', 6.4, 37 * (u.arcmin**2), 1, 'LAE')
willot_2005_1030 = Observation('SDSS_J1030+0524', 6.28, 30 *(u.arcmin**2), 1, 'LBG')
willot_2005_1048 = Observation('SDSS_J1048+4637', 6.23, 30 *(u.arcmin**2), 1, 'LBG')
willot_2005_1148 = Observation('SDSS_J1148+5251', 6.42, 30 *(u.arcmin**2), 1, 'LBG')
banados_2013 = Observation('ULAS_J0203+0012', 5.72, 46.23 * (u.arcmin**2), 1, 'LAE')
simpson_2014 = Observation('ULAS_J1120+0641', 7.08, 13.125 * (u.arcmin**2), 1, 'LBG')
ota_2018 = Observation('VIK_J0305-3150', 6.61, 0.2*(u.deg**2),1, 'LAE')
kim_2009_1048 = Observation('SDSS_J1048+4637', 6.23, 11.3 * (u.arcmin**2), 1, 'LBG')
kim_2009_1148 = Observation('SDSS_J1148+5251', 6.4, 11.3*(u.arcmin**2), 0.375, 'LBG')
kim_2009_1306 = Observation('SDSS J1306+0356', 5.99, 11.3*(u.arcmin**2), 0.125, 'LBG')
kikuta_2017_0807_lae = Observation('SDSS J0807+1328', 4.885, 0.25*(u.deg**2), 1, 'LAE')
kikuta_2017_0807_lbg = Observation('SDSS J0807+1328', 4.885, 0.25*(u.deg**2), 1, 'LBG')
kikuta_2017_1113_lae = Observation('SDSS J1113+0253', 4.87, 0.25*(u.deg**2), 1, 'LAE')
kikuta_2017_1113_lbg = Observation('SDSS J1113+0253', 4.87, 0.25*(u.deg**2), 1, 'LBG')



surveys = [
    stiavelli_2005,
    husband_2013_0338,
    husband_2013_1204,
    husband_2013_2130,
    morselli_2014_1030,
    morselli_2014_1048,
    morselli_2014_1148,
    morselli_2014_1411,
    champagne_2023_2054,
    champagne_2023_2348,
    champagne_2023_3051,
    balmaverde_2017,
    zheng_2006,
    utsumi_2010,
]

no_surveys = [
    goto_2017,
    muzzucchelli_2017,
    willot_2005_1030,
    willot_2005_1048,
    willot_2005_1148,
    banados_2013,
    simpson_2014,
    ota_2018,
    kikuta_2017_0807_lbg,
    kikuta_2017_0807_lae,
    kikuta_2017_1113_lae,
    kikuta_2017_1113_lbg,
]

yes_redshifts = [obs.qso_redshift for obs in surveys]
no_redshifts = [obs.qso_redshift for obs in no_surveys]
all_redshifts = yes_redshifts + no_redshifts

norm = Normalize(vmin=min(all_redshifts), vmax=max(all_redshifts))
cmap = plt.cm.get_cmap('inferno_r')
sm = ScalarMappable(norm=norm, cmap=cmap)


start_plot(r'Search Area [cMpc$^{2}$]', 'Overdensity factor')
for survey in surveys:
    color = sm.to_rgba(survey.qso_redshift)
    if survey.survey_type == 'LBG':
        SHAPE = 's'
    else:
        SHAPE = 'd'
    plt.scatter(survey.comoving_area, survey.overdensity_factor, s=20, c=color, zorder=99, marker=SHAPE)


for survey in no_surveys:
    color = sm.to_rgba(survey.qso_redshift)
    if survey.survey_type == 'LBG':
        SHAPE = 's'
    else:
        SHAPE = 'd'
    plt.scatter(survey.comoving_area, survey.overdensity_factor, s=30, facecolor='None', edgecolors=color, zorder=99, marker=SHAPE)

plt.axhline(1, ls=':', color='k', lw=1, alpha=0.4)
plt.errorbar(lambert_2023.comoving_area, lambert_2023.overdensity_factor, ms=10, c='m', marker='*', yerr=2, capsize=2, ecolor='k')
#plt.errorbar(lambert_2023.comoving_area, 10, ms=10, c='m', marker='*', yerr=7, capsize=2, ecolor='k')
plt.axvspan(10**(2.5), 10**(3.5), color='r', alpha=0.2, label='Protocluster area')
plt.legend(frameon=False)
plt.xscale('log')
cbar = plt.colorbar(sm, label='QSO Redshift')
end_plot('plots/chiara_lambert.png')
