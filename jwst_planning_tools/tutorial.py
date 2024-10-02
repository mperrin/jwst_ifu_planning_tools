from jwst_planning_tools.tools import mrs_planning_tool, simulate_geometry

if __name__=="__main__":
    # specify V3 PA on your own
    # simulate_geometry(planets=[(1.7, 72, 1e-4), (0.96, 344, 1e-4), (0.7, 243, 1e-4)], v3pa=80, band="1A", offset=(0., -0.),
    #                   which="4pt", sign="ext", system_name="HR8799-", primary=0, webbpsf_plot=True)

    # automatically fetch V3PA range
    mrs_planning_tool(planets=[(1.72, 72.8, 4e-4), (0.96, 345, 8e-4), (0.7, 244, 8e-4)], target_name="HR 8799",
                      band="2A",
                      primary=0, vscale_im=0.001,
                      vscale_res=2e-4,
                      sign="ext", jwst_cycle=3)