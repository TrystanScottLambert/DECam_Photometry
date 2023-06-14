"""
Zero points of the CDFS observations
"""

from zero_points import ZeroPoint, ZeroPoints

zero_points_cdfs = ZeroPoints(
    i_band = ZeroPoint(31.29834347060596, 1.22),
    z_band = ZeroPoint(31.43450767279138, 1.14),
    n964_band = ZeroPoint(29.08108350751592, 1.18)
)
