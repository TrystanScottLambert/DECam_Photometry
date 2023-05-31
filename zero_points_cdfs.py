"""
Zero points of the CDFS observations
"""

from zero_points import ZeroPoint, ZeroPoints

zero_points_cdfs = ZeroPoints(
    i_band = ZeroPoint(31.658628725529308, 1.22),
    z_band = ZeroPoint( 31.751747297788803, 1.14),
    n964_band = ZeroPoint(29.128647505175802, 1.18)
)
