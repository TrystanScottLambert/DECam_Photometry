"""
Zero points of the CDFS observations
"""

from zero_points import ZeroPoint, ZeroPoints

zero_points_cdfs = ZeroPoints(
    i_band = ZeroPoint(31.332476348946987, 1.22),
    z_band = ZeroPoint(31.438599090248633, 1.14),
    n964_band = ZeroPoint(29.066102291892573, 1.18)
)
