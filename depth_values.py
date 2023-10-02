""" Data classes for the depth parameters"""

from dataclasses import dataclass

@dataclass
class FilterDepthInfo:
    """
    Storing the sigma values for the depth of a particular band
    and the exponential parameters which were used in fitting that
    data.
    """
    sigma_1: float
    sigma_2: float
    sigma_3: float
    sigma_4: float
    sigma_5: float
    exponential_fit_params: list[float, float]

@dataclass
class DepthInfo:
    """
    All the depth info for a particular survey.
    """
    i_band: FilterDepthInfo
    z_band: FilterDepthInfo
    n_band: FilterDepthInfo

OUR_DEPTH = DepthInfo(
    i_band = FilterDepthInfo(27.148597376286858, 26.39565702944552, 25.955215161263357, 25.642716682604185, 25.40032403118172, [71477524534.81918, -0.9205871135313299]),
    z_band = FilterDepthInfo(27.00302357678343, 26.250174093690596, 25.809785377393975, 25.497324610597758, 25.254961210768744, [62700533749.60914, -0.9206982220567864]),
    n_band = FilterDepthInfo(26.41050514582186, 25.659170228114288, 25.21966747577294, 24.90783531040671, 24.66595949172676, [38162619550.192345, -0.9225541955042275]),
)

CDFS_DEPTH = DepthInfo(
    i_band = FilterDepthInfo(28.742923861348057, 27.983866787968037, 27.539846864133583, 27.224809714588023, 26.980447917044017, [250603829396.6726, -0.9131687258685562]),
    z_band = FilterDepthInfo(28.46646348454595, 27.71066199480761, 27.268546465321496, 26.954860505069266, 26.71154677136477, [217759094293.51932, -0.9171021623679402]),
    n_band = FilterDepthInfo(26.552147283420556, 25.786783049106962, 25.339073672640346, 25.02141881479336, 24.77502656494586, [27757748648.752213 -0.9056435478482744]),
)
