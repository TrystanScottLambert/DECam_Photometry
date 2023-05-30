"""
Module to store the zero points taking into account the k constant.
"""

from dataclasses import dataclass
from k_constant import calculate_k_constant_mag

class ZeroPoint:
    """Base class for zero points"""
    def __init__(self, prime_value, seeing):
        self.prime_value = prime_value
        self.seeing = seeing

    def mag_correct(self, aperture_radius:float) -> float:
        """Determines the mag with the k constant correction."""
        return self.prime_value + calculate_k_constant_mag(aperture_radius, self.seeing)

@dataclass
class ZeroPoints:
    """Storing all the information for the zero points."""
    i_band: ZeroPoint
    z_band: ZeroPoint
    n964_band: ZeroPoint

zero_points = ZeroPoints(
    
    i_band = ZeroPoint(30.92423685413585, 1.17),
    z_band = ZeroPoint(30.56990938472794, 1.23),
    n964_band = ZeroPoint(29.058777664420884, 1.47)
    )
