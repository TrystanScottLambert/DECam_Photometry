"""
Zero points of the CDFS observations
"""

from zero_points import ZeroPoint, ZeroPoints

zero_points_cdfs = ZeroPoints(
    i_band = ZeroPoint(31.32494051892022, 1.22),
    z_band = ZeroPoint(31.42510074230392, 1.14),
    n964_band = ZeroPoint(29.06573971755607, 1.18)
)
