from jwst_planning_tools.tools import mrs_planning_tool, simulate_geometry
"""
Tutorial to run the mrs planning tool.

author : Polychronis Patapis
updated : 2024-10-08

useage:
    python tutorial.py

Input V3PA parameters for the instrument of choice:
    NIRCam : 0.0°
    MIRI : 4.84°
    NIRISS : 0.56°
    NIRSpec : 138.5°
Parameters for mrs_planning_tool:
    planets : list((separation ["], parallactic angle [degree], contrast)) (whereistheplanet.com)
    target_name : str, SIMBAD resolvable name
    band : str, MRS sub channel to use
    primary : int, 0 will center the scene on the host, 1-N will be the first to Nth planet
    vscale_im : float, image scaling for color scale
    vscale_res : float, residual image scaling for color scale
    jwst_cycle : int, when will the observation occur
"""
if __name__ == "__main__":
    mrs_planning_tool(planets=[(1.72, 72.8, 4e-4), (0.96, 345, 8e-4), (0.7, 244, 8e-4)], 
                      target_name="HR 8799",
                      band="2A",
                      primary=0, 
                      vscale_im=0.001,
                      vscale_res=2e-4,
                      sign="ext", 
                      jwst_cycle=4)
