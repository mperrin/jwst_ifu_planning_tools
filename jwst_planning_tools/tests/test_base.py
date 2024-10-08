from jwst_planning_tools.tools import mrs_planning_tool
from jwst_planning_tools.util import get_v3PA_range
import numpy as np
import webbpsf
import os
from miricoord.mrs import mrs_tools as mt

# test getting the coordinates and V3 PA range

def test_get_V3PA_from_target():
    target_name = "HR 8799"
    _, _, _, _, V3PA_range = get_v3PA_range(target_name=target_name, ra=None, dec=None, start=None, end=None,
                                            cycle=3)
    # print(~np.isnan(np.array(V3PA_range)).any())
    assert ~np.isnan(np.array(V3PA_range)).any()

def test_webbpsf_MRS():
    miri = webbpsf.MIRI()
    miri.mode = "IFU"
    miri.options['ifu_broadening'] = "empirical"
    miri.include_si_wfe = False
    miri.band = "1A"
    miri.options['source_offset_x'] = 0.
    miri.options['source_offset_y'] = 0.
    psf = miri.calc_psf(monochromatic=5.3*1e-6, oversample=7, fov_arcsec=6, add_distortion=True)
    print(miri.pixelscale, miri.band, miri.aperturename, miri._get_aperture_rotation(miri.aperturename),
          miri._get_pixelscale_from_apername(miri.aperturename))

def test_miricoord_transform():
    print(mt.v2v3_to_xyideal(-503.378, -318.999))


if __name__ == "__main__":
    test_get_V3PA_from_target()
    test_webbpsf_MRS()
    test_miricoord_transform()