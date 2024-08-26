import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.ndimage import shift
from matplotlib.collections import PatchCollection
import pandas as pd
from astropy.io import fits
import astropy.units as u

#### MRS ###

fov_MRS = {"1": (3.3,3.7), "2":(4.0 , 4.8), "3":(5.6, 6.2)}
dither_off = {"neg": [(-1.078, 0.528), (0.98, -0.44), (1.078, -0.528), (-0.980, 0.44)],
              "pos": [(1.078, 0.528), (-0.98, -0.44), (-1.078, -0.528), (0.980, 0.44)],
              "ext":[(0.4*0.18, 0.4*0.18), (-0.4*0.18, 0.4*0.18), (0.4*0.18,-0.4*0.18), (-0.4*0.18,-0.4*0.18)],
             "ext-all":[(-0.09740,0.46396), (0.18792, -0.51062), (0.18845, 0.46396), (-0.09800, -0.51062)]}

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


def simulate_geometry(planets, v3pa, band, which, sign, offset=None, system_name=None, primary=None,
                      webbpsf_plot=False):
    """

    :param planets: list of tuples with (separation, pa, flux (optional)) for each planet
    :param v3pa:
    :param band:
    :param which:
    :param sign:
    :param offset:
    :param system_name:
    :param primary:
    :param webbpsf:
    :return:
    """
    nplanets = np.shape(planets)[0]
    coordinates = np.zeros((nplanets+1, 2))
    contrasts = np.ones(nplanets+1)
    for i, p in enumerate(planets):
        separation, pa, contrast = p
        mrspa = pa_to_MRS_pa(pa, v3pa, band)
        beta, alpha = separation * np.cos(-mrspa * np.pi / 180), separation * np.sin(-mrspa * np.pi / 180)
        coordinates[i+1, 0] = alpha
        coordinates[i+1, 1] = beta
        contrasts[i+1] = contrast

    if primary is not None:
        try:
            coordinates[:, 0] -= coordinates[primary, 0]
            coordinates[:, 1] -= coordinates[primary, 1]
            contrasts[:] /= contrasts[primary]
        except IndexError:
            print("Invalid index for planet. Setting star as primary")
    else:
        try:
            contrasts[:] /= contrasts[1]
        except IndexError:
            print("No planet in list.")

    if offset is not None:
        coordinates[:, 0] -= offset[0]
        coordinates[:, 1] -= offset[1]

    if system_name is None:
        system_name = "Tattoine"

    names = [system_name+n for n in ["A", "b", "c", "d", "e"][:nplanets + 1]]

    plt.figure()

    plt.scatter(coordinates[0, 0], coordinates[0, 1], s=100, marker="o", color="gold", label=names[0])
    for i in np.arange(1, nplanets+1):
        plt.scatter(coordinates[i, 0], coordinates[i, 1], marker="o", color=f"C{i}", label=names[i])

    plt.gca().set_aspect('equal')

    # plt.gca().arrow(2.3, -0.5, northx, northy, head_width=0.1, head_length=0.1, fc='k', ec='k')
    # plt.gca().arrow(2.3, -0.5, eastx, easty, head_width=0.1, head_length=0.1, fc='k', ec='k')
    # plt.text(1.6, -0.1, "N")
    # plt.text(2, -1, "E")
    dither_pattern(plt.gca(), [0., 0.], sign=sign, ch=band[0], which=which)

    # plt.xlim([-3, 4])
    # plt.ylim([-2, 6])
    plt.legend()
    plt.suptitle(f"V3 PA: {v3pa}")
    plt.xlabel("MRS alpha [arcsec]")
    plt.ylabel("MRS beta [arcsec]")
    plt.show()

    if webbpsf_plot:
        simfov = 2*fov_MRS[band[0]][0]
        miri = webbpsf.MIRI()
        miri.image_mask = "MIRI-IFU_3"
        miri.filter = "DSHORT"
        band = "3A"
        xoff = coordinates[0, 0]
        yoff = coordinates[0, 1]
        miri.options['source_offset_x'] = xoff
        miri.options['source_offset_y'] = yoff
        # produce PSF for each object given fluxes in Jy
        starpsf = miri.calc_psf(monochromatic=miri.wavelength, outfile=None, add_distortion=False,
                                fov_arcsec=simfov)

        miri.options['source_offset_x'] = xoff + np.random.normal(0, 0.1*miri._IFU_pixelscale[f"Ch{band[0]}"][1])
        miri.options['source_offset_y'] = yoff + np.random.normal(0, 0.1*miri._IFU_pixelscale[f"Ch{band[0]}"][1])
        # produce PSF for each object given fluxes in Jy
        refpsf = miri.calc_psf(monochromatic=miri.wavelength, outfile=None, add_distortion=False,
                                fov_arcsec=simfov)
        print("Placing star at ", miri.options['source_offset_x'], miri.options['source_offset_y'])
        xoff = coordinates[1, 0]
        yoff = coordinates[1, 1]
        miri.options['source_offset_x'] = xoff
        miri.options['source_offset_y'] = yoff
        print("Placing planet at ", coordinates[1, :])
        planetpsf = miri.calc_psf(monochromatic=miri.wavelength, outfile=None, add_distortion=False,
                                fov_arcsec=simfov)
        # ref_psf = shift(starpsf[3].data, shift=np.random.normal(0, 0.1, size=2))
        plt.figure()

        plt.imshow(starpsf[3].data*10000+planetpsf[3].data, origin="lower", cmap="plasma", vmax=1,
                   extent=[-simfov/2, simfov/2, -simfov/2, simfov/2])
        plt.scatter(coordinates[1, 0], coordinates[1, 1], marker="+", color=f"red")
        dither_pattern(plt.gca(), [0., 0.], sign=sign, ch=band[0], which=which, color="white")
        plt.gca().set_aspect('equal')
        plt.suptitle(f"V3 PA: {v3pa} - Science")
        plt.xlabel("MRS alpha [arcsec]")
        plt.ylabel("MRS beta [arcsec]")
        plt.show()

        plt.figure()

        plt.imshow(starpsf[3].data * contrasts[0] + planetpsf[3].data - refpsf[3].data*contrasts[0], origin="lower",
                   cmap="RdBu_r", vmin=-0.1, vmax=0.1,
                   extent=[-simfov / 2, simfov / 2, -simfov / 2, simfov / 2])
        plt.scatter(coordinates[1, 0], coordinates[1, 1], marker="+", color=f"red")
        dither_pattern(plt.gca(), [0., 0.], sign=sign, ch=band[0], which=which, color="white")
        plt.gca().set_aspect('equal')
        plt.suptitle(f"V3 PA: {v3pa} - Residuals")
        plt.xlabel("MRS alpha [arcsec]")
        plt.ylabel("MRS beta [arcsec]")
        plt.show()

if __name__=="__main__":

    import sys

    sys.path.append("/Users/polychronispatapis/Documents/PhD/MIRI/MRS_PSF/paper_analysis/")
    from optimise_cube_webbpsf import *


    simulate_geometry(planets=[(4.1, 37.4, 6*1e-3)], v3pa=60, band="3A",
                      which="4pt", sign="neg", offset=None, system_name="Tatooine-", primary=1, webbpsf_plot=True)