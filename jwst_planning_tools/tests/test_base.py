from jwst_planning_tools.tools import mrs_planning_tool
from jwst_planning_tools.util import get_v3PA_range
import numpy as np

# test getting the coordinates and V3 PA range

def test_get_V3PA_from_target():
    target_name = "HR 8799"
    _, _, _, _, V3PA_range = get_v3PA_range(target_name=target_name, ra=None, dec=None, start=None, end=None,
                                            cycle=3)

    assert np.isnan(V3PA_range).any()


if __name__ == "__main__":
    test_get_V3PA_from_target()