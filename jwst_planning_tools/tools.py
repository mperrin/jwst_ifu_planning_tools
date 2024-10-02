import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.ndimage import shift
import webbpsf
from matplotlib.collections import PatchCollection
import pandas as pd
from astropy.io import fits
import astropy.units as u
from photutils.aperture import CircularAperture
from jwst_planning_tools.util import get_v3PA_range
from miricoord.mrs import mrs_tools as mt
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
    ## elisabeth maths
    north = -(v3pa + 180 - mrs_rot[band[0]])  # the 180 comes from the orientation of alpha/beta with respect to V3
    return north + pa

def albe_to_v2v3(alpha, beta):
    theta = 188 * np.pi / 180
    v2 = alpha * np.cos(theta) + beta * np.sin(theta)
    v3 = beta * np.cos(theta) - alpha * np.sin(theta)
    return v2, v3

def dither_pattern(ax, center, sign, ch="1", which="4pt", color="C0"):
    w, h = fov_MRS[ch]

    d1 = Rectangle(xy=(dither_off[sign][0][0] - w * 0.5 + center[0], dither_off[sign][0][1] - h * 0.5 + center[1]),
                   width=w, height=h, fill=False, color=color, label=f"{which} {sign}")
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
                      webbpsf_plot=False, vscale_im=None, vscale_res=None):
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

    if vscale_im is None:
        vscale_im = 0.1
    if vscale_res is None:
        vscale_res = 0.1

    if webbpsf_plot:
        simfov = 2*fov_MRS[band[0]][0]
        miri = webbpsf.MIRI()
        miri.image_mask = "MIRI-IFU_3"
        miri.filter = "DSHORT"
        miri.band = band
        xoff = coordinates[0, 0]
        yoff = coordinates[0, 1]
        miri.options['source_offset_x'] = xoff
        miri.options['source_offset_y'] = yoff
        # produce PSF for each object given fluxes in Jy
        starpsf = miri.calc_psf(monochromatic=miri.wavelength, outfile=None, add_distortion=False,
                                fov_arcsec=simfov)
        system = starpsf[3].data*contrasts[0]  # initialise system with star
        system += np.random.normal(0, 0.1 * system)
        # miri.options['source_offset_x'] = xoff + np.random.normal(0, 0.1*miri._IFU_pixelscale[f"Ch{band[0]}"][1])
        # miri.options['source_offset_y'] = yoff + np.random.normal(0, 0.1*miri._IFU_pixelscale[f"Ch{band[0]}"][1])
        # produce PSF for each object given fluxes in Jy
        # refpsf = miri.calc_psf(monochromatic=miri.wavelength, outfile=None, add_distortion=False,
        #                         fov_arcsec=simfov)
        print("Placing star at ", miri.options['source_offset_x'], miri.options['source_offset_y'])
        for j in range(nplanets):
            xoff = coordinates[j+1, 0]
            yoff = coordinates[j+1, 1]
            miri.options['source_offset_x'] = xoff
            miri.options['source_offset_y'] = yoff
            print("Placing planet at ", coordinates[j+1, :])
            planetpsf = miri.calc_psf(monochromatic=miri.wavelength, outfile=None, add_distortion=False,
                                    fov_arcsec=simfov)
            system += planetpsf[3].data*contrasts[j+1]


        refpsf = shift(starpsf[3].data*contrasts[0], shift=np.random.normal(0, 0.01, size=2))


        fig, ax = plt.subplots(1, 3, figsize=(12, 4), sharey=True, squeeze=True)

        ax[0].scatter(coordinates[0, 0], coordinates[0, 1], s=100, marker="o", color="gold", label=names[0])
        for i in np.arange(1, nplanets + 1):
            ax[0].scatter(coordinates[i, 0], coordinates[i, 1], marker="o", color=f"C{i}", label=names[i])

        # plt.gca().set_aspect('equal')

        # plt.gca().arrow(2.3, -0.5, northx, northy, head_width=0.1, head_length=0.1, fc='k', ec='k')
        # plt.gca().arrow(2.3, -0.5, eastx, easty, head_width=0.1, head_length=0.1, fc='k', ec='k')
        # plt.text(1.6, -0.1, "N")
        # plt.text(2, -1, "E")
        dither_pattern(ax[0], [0., 0.], sign=sign, ch=band[0], which=which)
        ax[0].set_xlim([-simfov/2, simfov/2])
        ax[0].set_ylim([-simfov/2, simfov/2])
        ax[0].set_aspect('equal')
        ax[0].legend()
        ax[0].set_title(f"Geometry")
        ax[0].set_xlabel("MRS alpha [arcsec]")
        ax[0].set_ylabel("MRS beta [arcsec]")


        ax[1].imshow(system, origin="lower", cmap="plasma", vmax=vscale_im,
                   extent=[-simfov/2, simfov/2, -simfov/2, simfov/2])
        for i in np.arange(1, nplanets + 1):
            ap = CircularAperture((coordinates[i, 0], coordinates[i, 1]), r=0.3)
            ap.plot(ax[1], color="red")
            # ax[0].scatter(coordinates[i, 0], coordinates[i, 1], marker="o", color=f"re", label=names[i])
        dither_pattern(ax[1], [0., 0.], sign=sign, ch=band[0], which=which, color="white")
        ax[1].set_aspect('equal')
        ax[1].set_title(f"Science")
        ax[1].set_xlabel("MRS alpha [arcsec]")


        ax[2].imshow(system - refpsf, origin="lower",
                   cmap="RdBu_r", vmin=-vscale_res, vmax=vscale_res,
                   extent=[-simfov / 2, simfov / 2, -simfov / 2, simfov / 2])
        for i in np.arange(1, nplanets + 1):
            ap = CircularAperture((coordinates[i, 0], coordinates[i, 1]), r=0.3)
            ap.plot(ax[2], color="red")
        dither_pattern(ax[2], [0., 0.], sign=sign, ch=band[0], which=which, color="black")
        ax[2].set_aspect('equal')
        ax[2].set_title(f"Residuals")
        ax[2].set_xlabel("MRS alpha [arcsec]")

        fig.suptitle(f"V3 PA: {v3pa}, wavelength: {miri.wavelength:.1f}")
        plt.show()
    else:
        fig = plt.figure()
        plt.scatter(coordinates[0, 0], coordinates[0, 1], s=100, marker="o", color="gold", label=names[0])
        for i in np.arange(1, nplanets + 1):
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
    return fig


def mrs_planning_tool(planets, target_name, band="1A", primary=0, which="4pt", sign="pos",
                      ra=None, dec=None, start_time=None, end_time=None, jwst_cycle=4, webbpsf_plot=True,
                      vscale_im=None, vscale_res=None):
    # get system V3PA range in given time period of JWST cycle
    _, _, _, _, V3PA_range = get_v3PA_range(target_name=target_name, ra=ra, dec=dec, start=start_time, end=end_time,
                                            cycle=jwst_cycle)
    v3pa_msg = ""
    for i in range(len(V3PA_range)):
        v3pa_msg += f"{V3PA_range[0][i]}-{V3PA_range[1][i]}\n"

    while True:
        # Output the available V3 PA ranges
        print("The available V3 PA ranges for the target are:")
        print(v3pa_msg)

        # Take input from the user
        try:
            v3_pa = float(input("Enter V3 PA (float): "))
            offsetx = float(input("Enter offset x (tuple): "))
            offsety = float(input("Enter offset y (tuple): "))
        except ValueError:
            print("Invalid input. Please enter valid float values.")
            continue  # Loop again if input is invalid
        offset = (offsetx, offsety)
        v2off, v3off = mt.abtov2v3(-offsetx, -offsety, channel=band)
        offsetxyidl = np.round(mt.v2v3_to_xyideal(v2off, v3off), 2)
        # Perform actions with the inputs (e.g., calculate or process)
        print(f"You entered V3 PA: {v3_pa} and offset: {offset}")
        print("Simulating geometry")
        fig = simulate_geometry(planets=planets, v3pa=v3_pa, band=band, offset=offset,
                          which=which, sign=sign, system_name=target_name + "-", primary=primary,
                          webbpsf_plot=webbpsf_plot, vscale_im=vscale_im, vscale_res=vscale_res)
        # Ask if the user wants to continue
        saveplot = input("Do you want to save the plot? (y/n): ").lower()
        if saveplot == 'y':
            savepath = input("provide path: ").lower()
            fig.savefig(savepath+f"_{target_name}_{band}_V3_{v3_pa}_offsetXYidl_{offsetxyidl}.png", dpi=300)
        # Ask if the user wants to continue
        should_continue = input("Do you want to continue? (y/n): ").lower()
        if should_continue != 'y':
            print("Exiting...")
            break




