"""
Zero points of the CDFS observations
"""

from zero_points import ZeroPoint, ZeroPoints

zero_points_cdfs = ZeroPoints(
    i_band = ZeroPoint(31.65954985394665, 1.22),
    z_band = ZeroPoint(31.744812537112118, 1.14),
    n964_band = ZeroPoint(29.128647505175802, 1.18)
)
