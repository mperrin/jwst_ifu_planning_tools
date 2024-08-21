import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import pandas as pd
from astropy.io import fits
import astropy.units as u

#### MRS ###

def pa_to_MRS_pa(pa, v3pa, band):
    """
    Transforms the position angle on Sky to position angle on the MRS FoV (counterclockwise of beta)
    Args:
        pa: PA in degrees with respect to North
        v3pa: PA in degrees of JWST V3 axis during time of observation (V3PA)
        band: MRS band (1A-4C) or (1SHORT - 4LONG)

    Returns: PA in degrees with respect to MRS beta

    """
    mrs_rot = {"1": 8.3, "2": 8.2, "3": 7.6, "4": 8.4}  # from Patapis+23
    north = -(v3pa + 180 + mrs_rot[band[0]])  # the 180 comes from the orientation of alpha/beta with respect to V3
    return north + pa

def albe_to_v2v3(alpha, beta):
    theta = 188 * np.pi / 180
    v2 = alpha * np.cos(theta) + beta * np.sin(theta)
    v3 = beta * np.cos(theta) - alpha * np.sin(theta)
    return v2, v3

def dither_pattern(ax, center, sign, ch="1", which="4pt", color="C0"):
    w, h = fov_MRS[ch]

    d1 = Rectangle(xy=(dither_off[sign][0][0] - w * 0.5 + center[0], dither_off[sign][0][1] - h * 0.5 + center[1]),
                   width=w, height=h, fill=False, color=color, label=sign)
    d2 = Rectangle(xy=(dither_off[sign][1][0] - w * 0.5 + center[0], dither_off[sign][1][1] - h * 0.5 + center[1]),
                   width=w, height=h, fill=False, color=color)
    d3 = Rectangle(xy=(dither_off[sign][2][0] - w * 0.5 + center[0], dither_off[sign][2][1] - h * 0.5 + center[1]),
                   width=w, height=h, fill=False, color=color)
    d4 = Rectangle(xy=(dither_off[sign][3][0] - w * 0.5 + center[0], dither_off[sign][3][1] - h * 0.5 + center[1]),
                   width=w, height=h, fill=False, color=color)
    ax.add_patch(d1)
    ax.add_patch(d2)
    if which == "4pt":
        ax.add_patch(d3)
        ax.add_patch(d4)
    return


def albe_v2v3_off(al, be):  # offset requirement in V2/V3
    theta = 11 * np.pi / 180
    return al * np.cos(theta) + be * np.sin(theta), -al * np.sin(theta) + be * np.cos(theta)

def create_dithered_scene(planet_angles, planet_sep, v3_pa, planet_names, dither_sign, center=None, which="4pt",
                          band="1A"):
    fig, ax = plt.subplots()
    if center is None:
        center = [0, 0]
    else:
        ax.scatter(center[0], center[1], c="red", marker="+", label="offset")
    ax.scatter(0, 0, marker="+")
    ax.set_xlim([-4 + center[0], 4 + center[0]])
    ax.set_ylim(-4 + center[0], 4 + center[0])
    ax.set_xlabel("MRS alpha [arcsec]")
    ax.set_ylabel("MRS beta [arcsec]")

    # add planets
    for i, r in enumerate(planet_sep):
        v3PA = np.linspace(v3_pa[0], v3_pa[1], endpoint=True)
        MRSPAS = pa_to_MRS_pa(pa=planet_angles[i], v3pa=v3PA, band=band)
        beta, alpha = planet_sep[i] * np.cos(-MRSPAS * np.pi / 180), planet_sep[i] * np.sin(-MRSPAS * np.pi / 180)
        ax.scatter(alpha, beta, label=planet_names[i])

        dither_pattern(ax=ax, center=[np.median(alpha), np.median(beta)], sign=dither_sign, which=which, ch=band[0])

    #     plt.plot((0, 1.5*np.cos(PA.mean()+NE*np.pi/180)), (0, 1.5*np.sin(PA.mean()+NE*np.pi/180)), "red", alpha=0.5, label="N")
    #     plt.plot((0, 1.5*np.cos(PA.mean()+(90-NE)*np.pi/180)), (0, 1.5*np.sin(PA.mean()+(90-NE)*np.pi/180)), "green", alpha=0.5, label="E")

    ax.legend()
    return


### LRS ###


def pa_to_LRS_pa(pa, v3pa):
    """
        Transforms the position angle on Sky to position angle on the LRS FoV
        Args:
            pa: PA in degrees with respect to North
            v3pa: PA in degrees of JWST V3 axis during time of observation (V3PA)
        Returns: PA in degrees with respect to LRS x, y coordinates

        """
    north = -(v3pa + 4.8)
    return north + pa


def v2v3_to_XidlYidl(v2, v3):
    return -v2, v3






