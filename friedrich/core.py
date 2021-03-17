import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import interp2d
from shapely.geometry import Polygon
import rasterio
import rioxarray


def distance_between_two_points(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def IDW(window, coords_x, coords_y, center, power=4.5):

    window = window.T
    weight_list = []

    for i in pd.unique(coords_x):
        for j in pd.unique(coords_y):
            d = distance_between_two_points(center, (i, j))
            weight_list.append(1 / d ** power)

    # flatten window array and remove nans
    window = window.flatten()
    nan_indeces = np.argwhere(np.isnan(window))
    window = np.delete(window, nan_indeces)
    weight_list = np.delete(weight_list, nan_indeces)

    # return IDW
    return np.dot(window, weight_list) / np.nansum(weight_list)


def nearest_neighbor(window, coords_x, coords_y, center):

    weight_list = []

    for i in pd.unique(coords_x):
        for j in pd.unique(coords_y):
            d = distance_between_two_points(center, (i, j))
            weight_list.append(d)

    min_distance_value_index = np.argmin(weight_list)

    window = window.flatten()

    return window[min_distance_value_index]


def scipy_methods(window, coords_x, coords_y, center, kind="linear"):

    if window.shape[0] != window.shape[1]:
        return np.nan
    else:
        fun = interp2d(
            x=pd.unique(coords_x), y=pd.unique(coords_y), z=window, kind=kind
        )
        return fun(center[0], center[1])[0]


def find_optimal_IDW_power(
    positions, data, dx=3, dy=3, spacing=1, powers=np.arange(0, 20, 0.1)
):
    diffs = []

    for power in powers:
        window_IDW = []
        scipy_cubic = []

        for p in positions:
            window = data[p[0] - dx : p[0] + dx, p[1] - dy : p[1] + dy]

            if np.isnan(window).any():
                continue
            else:
                center = dx - spacing / 2, dy - spacing / 2

                coords = []
                for x in np.arange(0, window.shape[0], 1):
                    for y in np.arange(0, window.shape[1], 1):
                        coords.append((x, y))

                coords_x = np.array(coords)[:, 0]
                coords_y = np.array(coords)[:, 1]

                IDW_return = IDW(window, coords_x, coords_y, center, power=power)
                scipy_cubic_return = scipy_methods(
                    window, coords_x, coords_y, center, kind="cubic"
                )

                window_IDW.append(IDW_return)
                scipy_cubic.append(scipy_cubic_return)

        diff = np.array(window_IDW) - np.array(scipy_cubic)

        diffs.append(diff)

    idx = np.nanargmin(np.abs(diffs))
    print("Optimal IDW power:", powers[idx])
    return powers[idx]


def find_optimal_georaster_IDW_power(
    raster_file_name, positions, offset=3, spacing=1, powers=np.arange(0, 20, 0.1)
):

    source = rasterio.open(raster_file_name)
    spacing = source.res[0]
    dx = offset * spacing
    dy = offset * spacing

    diffs = []

    for power in powers:
        window_IDW = []
        scipy_cubic = []

        for p in positions:
            UL = (p[0] - dx, p[1] + dy)
            UR = (p[0] + dx, p[1] + dy)
            LR = (p[0] + dx, p[1] - dy)
            LL = (p[0] - dx, p[1] - dy)

            window_polygon = Polygon([UL, UR, LR, LL, UL])
            window, transform = rasterio.mask.mask(
                source, [window_polygon], crop=True, all_touched=True
            )
            window = window.squeeze()
            window = replace_and_fill_nodata_value(window, source.nodata, np.nan)
            if np.isnan(window).any():
                continue
            else:
                coords = []
                for x in np.arange(0, window.shape[0], 1):
                    for y in np.arange(0, window.shape[1], 1):
                        coords.append(transform * (x, y))
                coords = np.array(coords)
                coords_x = pd.unique(coords[:, 0]) + spacing / 2
                coords_y = pd.unique(coords[:, 1]) - spacing / 2

                IDW_return = IDW(window, coords_x, coords_y, p, power=power)
                scipy_cubic_return = scipy_methods(
                    window, coords_x, coords_y, p, kind="cubic"
                )

                window_IDW.append(IDW_return)
                scipy_cubic.append(scipy_cubic_return)
                break

        diff = np.array(window_IDW) - np.array(scipy_cubic)

        diffs.append(diff)

    idx = np.nanargmin(np.abs(diffs))
    print("Optimal IDW power:", powers[idx])
    return powers[idx]


def distance_along_transect(list_of_coordinate_tuples):

    cumulative_distance_along_transect = 0
    distances = [0]
    c = 0
    for coordinate in list_of_coordinate_tuples:
        if c < len(list_of_coordinate_tuples) - 1:
            x1 = coordinate[0]
            y1 = coordinate[1]
            x2 = list_of_coordinate_tuples[c + 1][0]
            y2 = list_of_coordinate_tuples[c + 1][1]

            distance = distance_between_two_points((x1, y1), (x2, y2))
            cumulative_distance_along_transect += distance

            distances.append(cumulative_distance_along_transect)

            c += 1

    return distances


def replace_and_fill_nodata_value(array, nodata_value, fill_value):
    if np.isnan(nodata_value):
        masked_array = np.nan_to_num(array, nan=fill_value)
    else:
        mask = array == nodata_value
        masked_array = np.ma.masked_array(array, mask=mask)
        masked_array = np.ma.filled(masked_array, fill_value=fill_value)

    return masked_array


def interpolate(raster_file_name, positions, offset=3, kind="IDW"):

    source = rasterio.open(raster_file_name)
    spacing = source.res[0]
    dx = offset * spacing
    dy = offset * spacing

    results = []

    if kind == "IDW":
        IDW_power = find_optimal_georaster_IDW_power(
            raster_file_name, positions, offset=offset, powers=np.arange(0, 20, 0.1)
        )
    for p in positions:
        UL = (p[0] - dx, p[1] + dy)
        UR = (p[0] + dx, p[1] + dy)
        LR = (p[0] + dx, p[1] - dy)
        LL = (p[0] - dx, p[1] - dy)

        window_polygon = Polygon([UL, UR, LR, LL, UL])
        window, transform = rasterio.mask.mask(
            source, [window_polygon], crop=True, all_touched=True
        )
        window = window.squeeze()
        window = replace_and_fill_nodata_value(window, source.nodata, np.nan)
        if np.isnan(window).all():
            results.append(np.nan)

        else:
            coords = []
            for x in np.arange(0, window.shape[0], 1):
                for y in np.arange(0, window.shape[1], 1):
                    coords.append(transform * (x, y))
            coords = np.array(coords)
            coords_x = pd.unique(coords[:, 0]) + spacing / 2
            coords_y = pd.unique(coords[:, 1]) - spacing / 2

            if kind == "IDW":
                results.append(IDW(window, coords_x, coords_y, p, power=IDW_power))
            elif kind == "mean":
                results.append(np.nanmean(window))
            elif kind == "rio_NN":
                results.append(
                    list(source.sample([p,]))[
                        0
                    ][0]
                )
            elif kind == "NN":
                results.append(nearest_neighbor(window, coords_x, coords_y, p))
            elif kind == "linear":
                results.append(
                    scipy_methods(window, coords_x, coords_y, p, kind="linear")
                )
            elif kind == "cubic":
                results.append(
                    scipy_methods(window, coords_x, coords_y, p, kind="cubic")
                )
            elif kind == "quintic":
                results.append(
                    scipy_methods(window, coords_x, coords_y, p, kind="quintic")
                )

    if kind == "rio_NN":
        results = replace_and_fill_nodata_value(
            np.array(results), source.nodata, np.nan
        )

    return np.array(results)
