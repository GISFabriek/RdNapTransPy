# ***********************************************************************
# Author           : Willem A. Ligtendag, De GISFabriek
# Created          : 07-14-2019
#
# Last Modified By : Willem A. Ligtendag, De GISFabriek
# Last Modified On : 07-15-2019
# ***********************************************************************
# Python PORT from C version of RDNAPTRANS
# ***********************************************************************
import binascii
import struct
import math
import array
import constants


def grid_interpolation(x, y, gridfile, gridfileAction):
    decoded = binascii.a2b_base64(gridfile)
    error_message = ""
    size_x, size_y, min_x, max_x, min_y, max_y, min_value, max_value, error = read_grd_file_header(decoded)
    if error != 0:
        return 0, error, "reading grid file header failed"
    step_size_x = (max_x - min_x) / (size_x - 1)
    step_size_y = (max_y - min_y) / (size_y - 1)
    if (x <= (min_x + step_size_x) or x >= (max_x - step_size_x) or
            y <= (min_y + step_size_y) or y >= (max_y - step_size_y)):
        if gridfileAction == "x":
            error = 1
        if gridfileAction == "y":
            error = 2
        if gridfileAction == "geoid":
            error = 3
        return 0, error, "Outside bounding box"
    ddx = (x - min_x) / step_size_x - math.floor((x - min_x) / step_size_x)
    ddy = 1 - ((y - min_y) / step_size_y - math.floor((y - min_y) / step_size_y))
    array.array('i')
    record_number = array.array('i', (0 for i in range(0, 16)))
    record_number[5] = int((x - min_x) / step_size_x + math.floor((y - min_y) / step_size_y) * size_x)
    record_number[0] = record_number[5] - size_x - 1
    record_number[1] = record_number[5] - size_x
    record_number[2] = record_number[5] - size_x + 1
    record_number[3] = record_number[5] - size_x + 2
    record_number[4] = record_number[5] - 1
    record_number[6] = record_number[5] + 1
    record_number[7] = record_number[5] + 2
    record_number[8] = record_number[5] + size_x - 1
    record_number[9] = record_number[5] + size_x
    record_number[10] = record_number[5] + size_x + 1
    record_number[11] = record_number[5] + size_x + 2
    record_number[12] = record_number[5] + 2 * size_x - 1
    record_number[13] = record_number[5] + 2 * size_x
    record_number[14] = record_number[5] + 2 * size_x + 1
    record_number[15] = record_number[5] + 2 * size_x + 2
    array.array('f')
    record_value = array.array('f', (0 for i in range(0, 16)))
    for i in range(16):
        record_value[i] = read_grd_file_body(decoded, record_number[i])
        if record_value[i] > max_value + constants.PRECISION or record_value[i] < min_value - constants.PRECISION:
            if gridfileAction == "x":
                error = 1
            if gridfileAction == "y":
                error = 2
            if gridfileAction == "geoid":
                error = 3
            return 0, error, "Outside validity area"
    array.array('d')
    f = array.array('d', (0 for i in range(0, 4)))
    array.array('d')
    g = array.array('d', (0 for i in range(0, 4)))
    array.array('d')
    gfac = array.array('d', (0 for i in range(0, 16)))
    f[0] = -0.5 * ddx + ddx * ddx - 0.5 * ddx * ddx * ddx
    f[1] = 1.0 - 2.5 * ddx * ddx + 1.5 * ddx * ddx * ddx
    f[2] = 0.5 * ddx + 2.0 * ddx * ddx - 1.5 * ddx * ddx * ddx
    f[3] = -0.5 * ddx * ddx + 0.5 * ddx * ddx * ddx
    g[0] = -0.5 * ddy + ddy * ddy - 0.5 * ddy * ddy * ddy
    g[1] = 1.0 - 2.5 * ddy * ddy + 1.5 * ddy * ddy * ddy
    g[2] = 0.5 * ddy + 2.0 * ddy * ddy - 1.5 * ddy * ddy * ddy
    g[3] = -0.5 * ddy * ddy + 0.5 * ddy * ddy * ddy

    gfac[12] = f[0] * g[0]
    gfac[8] = f[0] * g[1]
    gfac[4] = f[0] * g[2]
    gfac[0] = f[0] * g[3]
    gfac[13] = f[1] * g[0]
    gfac[9] = f[1] * g[1]
    gfac[5] = f[1] * g[2]
    gfac[1] = f[1] * g[3]
    gfac[14] = f[2] * g[0]
    gfac[10] = f[2] * g[1]
    gfac[6] = f[2] * g[2]
    gfac[2] = f[2] * g[3]
    gfac[15] = f[3] * g[0]
    gfac[11] = f[3] * g[1]
    gfac[7] = f[3] * g[2]
    gfac[3] = f[3] * g[3]

    value = 0.0
    for i in range(16):
        value += gfac[i] * record_value[i]
    return value, error, error_message


def read_grd_file_header(decoded_string):
    id_string = decoded_string[0:4].decode("utf-8")
    if id_string != "DSBB":
        return 0, 0, 0, 0, 0, 0, 0, 0, -1
    temp = decoded_string[4: 6]
    size_x = int.from_bytes(temp, byteorder='little', signed=True)
    temp = decoded_string[6:8]
    size_y = int.from_bytes(temp, byteorder='little', signed=True)
    min_x = get_double(decoded_string, 8)
    max_x = get_double(decoded_string, 16)
    min_y = get_double(decoded_string, 24)
    max_y = get_double(decoded_string, 32)
    min_value = get_double(decoded_string, 40)
    max_value = get_double(decoded_string, 48)
    return size_x, size_y, min_x, max_x, min_y, max_y, min_value, max_value, 0


def get_slice(lst, rng):
    for i in rng:
        yield lst[i]


def get_double(decoded_string, start):
    end = start + 8
    temp = decoded_string[start:end]
    unpacked = struct.unpack('d', temp)
    result = unpacked[0]
    return result


def get_float(decoded_string, start):
    end = start + 4
    temp = decoded_string[start:end]
    unpacked = struct.unpack('f', temp)
    result = unpacked[0]
    return result


def read_grd_file_body(decoded_string, record_number):
    record_length = 4
    header_length = 56

    start = header_length + record_number * record_length
    result = get_float(decoded_string, start)
    return result
