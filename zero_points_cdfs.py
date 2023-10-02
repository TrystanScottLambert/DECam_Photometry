"""
Zero points of the CDFS observations
"""

from zero_points import ZeroPoint, ZeroPoints

zero_points_cdfs = ZeroPoints(
    i_band = ZeroPoint(31.328244356052533, 1.22),
    z_band = ZeroPoint(31.434538659814653, 1.14),
    n964_band = ZeroPoint(29.064151501608553, 1.18)
)
