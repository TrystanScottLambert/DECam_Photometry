"""
Zero points of the CDFS observations
"""

from zero_points import ZeroPoint, ZeroPoints

zero_points_cdfs = ZeroPoints(
    i_band = ZeroPoint(31.34957577696574, 1.22),
    z_band = ZeroPoint(31.529631523941692, 1.14),
    n964_band = ZeroPoint(29.0775366272342, 1.18)
)
