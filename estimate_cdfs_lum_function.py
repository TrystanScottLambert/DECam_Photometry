"""
Script to estimate lyman alpha lum function of cdfs
by assuming an erf function. 
"""

def completeness(narrowbands):
    values = []
    for narrowband in narrowbands:
        if narrowband < 24.5:
            values.append(0.63)
        else:
            values.append(0)
    return values

EFFECTIVE_VOLUME = 1.29e6 # co moving Mpc