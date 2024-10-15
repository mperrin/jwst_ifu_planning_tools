from jwst_gtvt.jwst_tvt import Ephemeris
from jwst_gtvt.display_results import get_visibility_windows
from astropy.time import Time
from astroquery.simbad import Simbad
import numpy as np
import warnings
from erfa import ErfaWarning

warnings.simplefilter("ignore", category=ErfaWarning)

def get_v3PA_range(target_name, ra=None, dec=None, start=None, end=None, cycle=4):
    """

    :param target_name: SIMBAD name of target - if applicable
    :param ra: (optional) ra of target in HH:MM:SS format - if no SIMBAD name
    :param dec: (optional) dec of target in HH:MM:SS format - if no SIMBAD name
    :param start: (optional) start date of period of interest in YYYY-MM-DDTHH:MM:SS format eg. 2024-06-01T00:00:00
    :param end: (optional) end date of period of interest in YYYY-MM-DDTHH:MM:SS format eg. 2024-06-01T00:00:00
    :param cycle: JWST cycle
    :return: start, end, ra, dec, v3pa range
    """
    if (ra is None) or (dec is None):
        try:
            target = Simbad.query_object(target_name)
            ra = target["RA"].data[0].replace(" ", ":")
            dec = target["DEC"].data[0].replace(" ", ":")
            # print(f"Coordinates for {target_name} are: {ra}, {dec}")
        except:
            return None
    else:
        ra = str(ra)
        dec = str(dec)


    if (start is None) or (end is None):
        year = 2021 + cycle
        start = Time(f"{year}-06-01T00:00:00")
        end = Time(f"{year+1}-06-01T00:00:00")

    eph = Ephemeris(start_date=start, end_date=end)
    eph.get_fixed_target_positions(ra, dec)
    dataframe = eph.dataframe
    dataframe['times'] = Time(dataframe['MJD'], format='mjd').datetime

    df = dataframe.loc[dataframe['in_FOR'] == True]
    # These indices allow us to get the visible regions regions
    window_indices = get_visibility_windows(df.index.tolist())
    if len(window_indices) != 2:
        return None

    min_PA_data, max_PA_data = [], []
    v3pa_range = []
    for starti, endi in window_indices:
        data_to_plot = df.loc[starti:endi]
        min_PA_data = np.round(np.min(data_to_plot['V3PA_min_pa_angle']), 1)
        max_PA_data = np.round(np.max(data_to_plot['V3PA_max_pa_angle']), 1)
        v3pa_range.append([min_PA_data, max_PA_data])
    return start, end, ra, dec, v3pa_range