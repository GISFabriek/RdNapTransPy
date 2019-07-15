# ***********************************************************************
# Author           : Willem A. Ligtendag, De GISFabriek
# Created          : 07-14-2019
#
# Last Modified By : Willem A. Ligtendag, De GISFabriek
# Last Modified On : 07-15-2019
# ***********************************************************************
# Python PORT from C version of RDNAPTRANS
# ***********************************************************************
import math
import constants
import grdfile
import resources


def deg_sin(alpha):
    return math.sin(alpha / 180.0 * math.pi)


def deg_cos(alpha):
    return math.cos(alpha / 180.0 * math.pi)


def deg_tan(alpha):
    return math.tan(alpha / 180.0 * math.pi)


def deg_asin(a):
    return math.asin(a) * 180.0 / math.pi


def deg_atan(a):
    return math.atan(a) * 180.0 / math.pi


def atanh(a):
    return 0.5 * math.log((1.0 + a) / (1.0 - a))


def deg_min_sec2decimal(deg, mins, sec):
    dec_deg = (deg + mins / 60.0 + sec / 3600.0)
    return dec_deg


def decimal2deg_min_sec(dec_deg):
    deg = int(dec_deg)
    mins = int((dec_deg - deg) * 60.0)
    sec = ((dec_deg - deg) * 60.0 - mins) * 60.0
    return deg, mins, sec


def geographic2cartesian(phi, lmbd, h, a, inv_f):
    f = 1.0 / inv_f
    ee = f * (2.0 - f)
    n = a / math.sqrt(1.0 - ee * pow(deg_sin(phi), 2))
    x = (n + h) * deg_cos(phi) * deg_cos(lmbd)
    y = (n + h) * deg_cos(phi) * deg_sin(lmbd)
    z = (n * (1.0 - ee) + h) * deg_sin(phi)
    return x, y, z


def cartesian2geographic(x, y, z, a, inv_f):
    f = 1.0 / inv_f
    ee = f * (2.0 - f)
    rho = math.sqrt(x * x + y * y)
    n = 0
    phi = 0
    diff = 90
    while diff > constants.DEG_PRECISION:
        previous = phi
        n = a / math.sqrt(1.0 - ee * math.pow(deg_sin(phi), 2))
        phi = deg_atan(z / rho + n * ee * deg_sin(phi) / rho)
        diff = math.fabs(phi - previous)

    lmbd = deg_atan(y / x)
    h = rho * deg_cos(phi) + z * deg_sin(phi) - n * (1.0 - ee * pow(deg_sin(phi), 2))
    return phi, lmbd, h


def sim_trans(x_in, y_in, z_in, tx, ty, tz, alpha, beta, gamma, delta, xa, ya, za):
    a = math.cos(gamma) * math.cos(beta)
    b = math.cos(gamma) * math.sin(beta) * math.sin(alpha) + math.sin(gamma) * math.cos(alpha)
    c = -math.cos(gamma) * math.sin(beta) * math.cos(alpha) + math.sin(gamma) * math.sin(alpha)
    d = -math.sin(gamma) * math.cos(beta)
    e = -math.sin(gamma) * math.sin(beta) * math.sin(alpha) + math.cos(gamma) * math.cos(alpha)
    f = math.sin(gamma) * math.sin(beta) * math.cos(alpha) + math.cos(gamma) * math.sin(alpha)
    g = math.sin(beta)
    h = -math.cos(beta) * math.sin(alpha)
    i = math.cos(beta) * math.cos(alpha)

    x = x_in - xa
    y = y_in - ya
    z = z_in - za

    x_out = (1.0 + delta) * (a * x + b * y + c * z) + tx + xa
    y_out = (1.0 + delta) * (d * x + e * y + f * z) + ty + ya
    z_out = (1.0 + delta) * (g * x + h * y + i * z) + tz + za
    return x_out, y_out, z_out


def rd_projection(phi, lmbd):
    f = 1 / constants.INV_F_BESSEL
    ee = f * (2 - f)
    e = math.sqrt(ee)
    eea = ee / (1.0 - ee)

    phi_amersfoort_sphere = deg_atan(deg_tan(constants.PHI_AMERSFOORT_BESSEL) / math.sqrt(
        1 + eea * math.pow(deg_cos(constants.PHI_AMERSFOORT_BESSEL), 2)))
    lambda_amersfoort_sphere = constants.LAMBDA_AMERSFOORT_BESSEL

    r1 = constants.A_BESSEL * (1 - ee) / math.pow(math.sqrt(1 - ee *
                                                            math.pow(deg_sin(constants.PHI_AMERSFOORT_BESSEL), 2)), 3)
    r2 = constants.A_BESSEL / math.sqrt(1.0 - ee * pow(deg_sin(constants.PHI_AMERSFOORT_BESSEL), 2))
    r_sphere = math.sqrt(r1 * r2)

    n = math.sqrt(1 + eea * pow(deg_cos(constants.PHI_AMERSFOORT_BESSEL), 4))
    q_amersfoort = math.atanh(deg_sin(constants.PHI_AMERSFOORT_BESSEL)) - e * atanh(
        e * deg_sin(constants.PHI_AMERSFOORT_BESSEL))
    w_amersfoort = math.log(deg_tan(45 + 0.5 * phi_amersfoort_sphere))
    m = w_amersfoort - n * q_amersfoort

    q = atanh(deg_sin(phi)) - e * atanh(e * deg_sin(phi))
    w = n * q + m
    phi_sphere = 2 * deg_atan(math.exp(w)) - 90
    delta_lambda_sphere = n * (lmbd - lambda_amersfoort_sphere)
    sin_half_psi_squared = pow(deg_sin(0.5 * (phi_sphere - phi_amersfoort_sphere)), 2) + \
                           pow(deg_sin(0.5 * delta_lambda_sphere), 2) * deg_cos(phi_sphere) * \
                           deg_cos(phi_amersfoort_sphere)
    sin_half_psi = math.sqrt(sin_half_psi_squared)
    cos_half_psi = math.sqrt(1 - sin_half_psi_squared)
    tan_half_psi = sin_half_psi / cos_half_psi
    sin_psi = 2 * sin_half_psi * cos_half_psi
    cos_psi = 1 - 2 * sin_half_psi_squared
    sin_alpha = deg_sin(delta_lambda_sphere) * (deg_cos(phi_sphere) / sin_psi)
    cos_alpha = (deg_sin(phi_sphere) - deg_sin(phi_amersfoort_sphere) * cos_psi) / (deg_cos(phi_amersfoort_sphere) *
                                                                                    sin_psi)
    r = 2 * constants.SCALE_RD * r_sphere * tan_half_psi
    x_rd = r * sin_alpha + constants.X_AMERSFOORT_RD
    y_rd = r * cos_alpha + constants.Y_AMERSFOORT_RD
    return x_rd, y_rd


def inv_rd_projection(x_rd, y_rd):
    f = 1 / constants.INV_F_BESSEL
    ee = f * (2 - f)
    e = math.sqrt(ee)
    eea = ee / (1.0 - ee)

    phi_amersfoort_sphere = deg_atan(deg_tan(constants.PHI_AMERSFOORT_BESSEL) /
                                     math.sqrt(1 + eea * math.pow(deg_cos(constants.PHI_AMERSFOORT_BESSEL), 2)))

    r1 = constants.A_BESSEL * (1 - ee) / pow(math.sqrt(1 - ee * pow(deg_sin(constants.PHI_AMERSFOORT_BESSEL), 2)), 3)
    r2 = constants.A_BESSEL / math.sqrt(1.0 - ee * pow(deg_sin(constants.PHI_AMERSFOORT_BESSEL), 2))
    r_sphere = math.sqrt(r1 * r2)

    n = math.sqrt(1 + eea * math.pow(deg_cos(constants.PHI_AMERSFOORT_BESSEL), 4))
    q_amersfoort = atanh(deg_sin(constants.PHI_AMERSFOORT_BESSEL)) - e * atanh(e *
                                                                               deg_sin(constants.PHI_AMERSFOORT_BESSEL))
    w_amersfoort = math.log(deg_tan(45 + 0.5 * phi_amersfoort_sphere))
    m = w_amersfoort - n * q_amersfoort

    r = math.sqrt(pow(x_rd - constants.X_AMERSFOORT_RD, 2) + pow(y_rd - constants.Y_AMERSFOORT_RD, 2))
    sin_alpha = (x_rd - constants.X_AMERSFOORT_RD) / r
    if r < constants.PRECISION:
        sin_alpha = 0
    cos_alpha = (y_rd - constants.Y_AMERSFOORT_RD) / r
    if r < constants.PRECISION:
        cos_alpha = 1
    psi = 2 * deg_atan(r / (2 * constants.SCALE_RD * r_sphere))
    phi_sphere = deg_asin(cos_alpha * deg_cos(phi_amersfoort_sphere) * deg_sin(psi) +
                          deg_sin(phi_amersfoort_sphere) * deg_cos(psi))
    delta_lambda_sphere = deg_asin((sin_alpha * deg_sin(psi)) / deg_cos(phi_sphere))
    lmbd = delta_lambda_sphere / n + constants.LAMBDA_AMERSFOORT_BESSEL
    w = atanh(deg_sin(phi_sphere))
    q = (w - m) / n
    phi = 0
    diff = 90
    while diff > constants.DEG_PRECISION:
        previous = phi
        phi = 2 * deg_atan(math.exp(q + 0.5 * e * math.log((1 + e * deg_sin(phi)) / (1 - e * deg_sin(phi))))) - 90
        diff = math.fabs(phi - previous)
    return phi, lmbd


def rd_correction(x_ps_rd, y_ps_rd):
    dx, xerror, error_message = grdfile.grid_interpolation(x_ps_rd, y_ps_rd, resources.gridFileDx, "x")
    dy, yerror, error_message = grdfile.grid_interpolation(x_ps_rd, y_ps_rd, resources.gridFileDy, "y")
    x_rd = x_ps_rd - dx
    y_rd = y_ps_rd - dy
    error = ""
    if xerror != 0:
        error =xerror
    else:
        error = yerror
    return x_rd, y_rd, error


def inv_rd_correction(x_rd, y_rd):
    x_ps_rd = x_rd
    y_ps_rd = y_rd
    dx, xerror, error_message = grdfile.grid_interpolation(x_ps_rd, y_ps_rd, resources.gridFileDx, "x")
    dy, yerror, error_message = grdfile.grid_interpolation(x_ps_rd, y_ps_rd, resources.gridFileDy, "y")
    x_pseudo_rd = x_rd + dx
    y_pseudo_rd = y_rd + dy
    error = ""
    if xerror != 0:
        error = xerror
    else:
        error = yerror
    return x_pseudo_rd, y_pseudo_rd, error
