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
    i_band = FilterDepthInfo(27.14799344881014, 26.39522705157127, 25.95488693738356, 25.642460654332393, 25.40012400217406, [71851556429.43224, -0.9207998432214697]),
    z_band = FilterDepthInfo(27.002434075864002, 26.249755521486453, 25.80946679207857, 25.497076967108896, 25.254768594035557, [63021315517.44405, -0.9209073070152256]),
    n_band = FilterDepthInfo(26.40512447904803, 25.65384528287756, 25.2143751255459, 24.902566086707086, 24.66070820635542, [38042322869.9511, -0.9226226203163284]),
)

CDFS_DEPTH = DepthInfo(
    i_band = FilterDepthInfo(28.741277354413807, 27.983494488925896, 27.540219928926444, 27.225711623437984, 26.981760029213174, [261517625433.30716, -0.9147042142654559]),
    z_band = FilterDepthInfo(28.4641835053049, 27.70983988603929, 27.26857715611063, 26.955496266773682, 26.712651862533065, [228548624459.7423, -0.9188745856096189]),
    n_band = FilterDepthInfo(26.543819446284726, 25.780819519558644, 25.334493174370895, 25.017819592832556, 24.772188480022432, [29679701475.91775, -0.908449865171198]),
)
