# ***********************************************************************
# Author           : Willem A. Ligtendag, De GISFabriek
# Created          : 07-14-2019
#
# Last Modified By : Willem A. Ligtendag, De GISFabriek
# Last Modified On : 07-15-2019
# ***********************************************************************
# Python PORT from C version of RDNAPTRANS
# ***********************************************************************
import constants
import helpers
import grdfile
import resources


def etrs2rd(phi_etrs, lambda_etrs, h_etrs):
    """
    etrs2rd.

    convert ETRS89 coordinates to RD coordinates.

    Parameters
    ----------
    phi_etrs : double
        Latitude
    lambda_etrs : double
        Longitude
    h_etrs : double
        Ellipsoidal height

    Returns
    -------
    tuple of double, double, double, int
        (x in RD, y in RD, z in NAP, errorcode)

    """
    x_amersfoort_bessel, y_amersfoort_bessel, z_amersfoort_bessel = helpers.geographic2cartesian(
        constants.PHI_AMERSFOORT_BESSEL, constants.LAMBDA_AMERSFOORT_BESSEL, constants.H_AMERSFOORT_BESSEL,
        constants.A_BESSEL, constants.INV_F_BESSEL)
    x_amersfoort_etrs = x_amersfoort_bessel + constants.TX_BESSEL_ETRS
    y_amersfoort_etrs = y_amersfoort_bessel + constants.TY_BESSEL_ETRS
    z_amersfoort_etrs = z_amersfoort_bessel + constants.TZ_BESSEL_ETRS

    x_etrs, y_etrs, z_etrs = helpers.geographic2cartesian(phi_etrs, lambda_etrs, h_etrs, constants.A_ETRS,
                                                          constants.INV_F_ETRS)
    x_bessel, y_bessel, z_bessel = helpers.sim_trans(x_etrs, y_etrs, z_etrs, constants.TX_ETRS_BESSEL,
                                                     constants.TY_ETRS_BESSEL, constants.TZ_ETRS_BESSEL,
                                                     constants.ALPHA_ETRS_BESSEL, constants.BETA_ETRS_BESSEL,
                                                     constants.GAMMA_ETRS_BESSEL, constants.DELTA_ETRS_BESSEL,
                                                     x_amersfoort_etrs, y_amersfoort_etrs, z_amersfoort_etrs)
    phi_bessel, lambda_bessel, h_bessel = helpers.cartesian2geographic(x_bessel, y_bessel, z_bessel,
                                                                       constants.A_BESSEL, constants.INV_F_BESSEL)
    x_pseudo_rd, y_pseudo_rd = helpers.rd_projection(phi_bessel, lambda_bessel)
    x_rd, y_rd, error = helpers.rd_correction(x_pseudo_rd, y_pseudo_rd)
    return x_rd, y_rd, h_bessel, error


def rd2etrs(x_rd, y_rd, nap):
    """
    rd2etrs.

    convert RD coordinates to ETRS89 coordinates.

    Parameters
    ----------
    x_rd : double
        RD x coordinate
    y_rd : double
        RD y coordinate
    nap : double
        Height in NAP

    Returns
    -------
    tuple of double, double, double, int
        (Latitude, Longitude, Ellipsoidal height, errorcode)

    """
    x_amersfoort_bessel, y_amersfoort_bessel, z_amersfoort_bessel = helpers.geographic2cartesian(
        constants.PHI_AMERSFOORT_BESSEL, constants.LAMBDA_AMERSFOORT_BESSEL, constants.H_AMERSFOORT_BESSEL,
        constants.A_BESSEL, constants.INV_F_BESSEL)
    h_bessel = nap + constants.MEAN_GEOID_HEIGHT_BESSEL
    x_pseudo_rd, y_pseudo_rd, error = helpers.inv_rd_correction(x_rd, y_rd)
    phi_bessel, lambda_bessel = helpers.inv_rd_projection(x_pseudo_rd, y_pseudo_rd)
    x_bessel, y_bessel, z_bessel = helpers.geographic2cartesian(phi_bessel, lambda_bessel, h_bessel, constants.A_BESSEL,
                                                                constants.INV_F_BESSEL)
    x_etrs, y_etrs, z_etrs = helpers.sim_trans(x_bessel, y_bessel, z_bessel,
                                               constants.TX_BESSEL_ETRS, constants.TY_BESSEL_ETRS,
                                               constants.TZ_BESSEL_ETRS,
                                               constants.ALPHA_BESSEL_ETRS, constants.BETA_BESSEL_ETRS,
                                               constants.GAMMA_BESSEL_ETRS,
                                               constants.DELTA_BESSEL_ETRS,
                                               x_amersfoort_bessel, y_amersfoort_bessel, z_amersfoort_bessel)
    phi_etrs, lambda_etrs, h_etrs = helpers.cartesian2geographic(x_etrs, y_etrs, z_etrs, constants.A_ETRS,
                                                                 constants.INV_F_ETRS)
    return phi_etrs, lambda_etrs, h_etrs, error


def etrs2nap(phi, lmbd, h):
    """
    etrs2nap.

    convert ellipsoidal ETRS89 height to NAP height.

    Parameters
    ----------
    phi : double
        Latitude
    lmbd : double
        Longitude
    h : double
        Ellipsoidal height

    Returns
    -------
    tuple of double, int
        (NAP height, errorcode)

    """
    n, err, error_message = grdfile.grid_interpolation(lmbd, phi, resources.gridFileGeoid, "geoid")
    nap = h - n + 0.0088
    return nap, err


def nap2etrs(phi, lmbd, nap):
    """
    nap2etrs.

    convert NAP height to ellipsoidal ETRS89 height.

    Parameters
    ----------
    phi : double
        Latitude
    lmbd : double
        Longitude
    NAP: double
        Height in NAP

    Returns
    -------
    tuple of double, int
        (Ellipsoidal ETRS89 height, errorcode)

    """
    n, err, error_message = grdfile.grid_interpolation(lmbd, phi, resources.gridFileGeoid, "geoid")
    h = nap + n - 0.0088
    return h, err


def etrs2rdnap(phi, lmbd, hg):
    """
    etrs2rdnap.

    convert ETRS89 coordinates to RD and NAP coordinates.

    Parameters
    ----------
    phi : double
        Latitude
    lmbd : double
        Longitude
    hg : double
        Ellipsoidal height

    Returns
    -------
    tuple of double, double, double
        (x in RD, y in RD, z in NAP)

    """
    x_rd, y_rd, h_bessel, error = etrs2rd(phi, lmbd, hg)
    h_geoid, err = etrs2nap(phi, lmbd, hg)
    if err == 3:
        nap = h_bessel
    else:
        nap = h_geoid
    return x_rd, y_rd, nap


def rdnap2etrs(x_rd, y_rd, nap):
    """
    rdnap2etrs.

    convert RD and NAP coordinates to ETRS89 coordinates.

    Parameters
    ----------
    x_rd : double
        RD x coordinate
    y_rd : double
        RD y coordinate
    nap : double
        Height in NAP

    Returns
    -------
    tuple of double, double, double
        (Latitude, Longitude, Ellipsoidal height)

    """
    phi, lmbd, h_etrs_sim, error = rd2etrs(x_rd, y_rd, nap)
    h_etrs_geoid, error = nap2etrs(phi, lmbd, nap)
    if error == 3:
        h = h_etrs_sim
    else:
        h = h_etrs_geoid
    return phi, lmbd, h
