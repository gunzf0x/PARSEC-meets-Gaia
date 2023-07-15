from astropy.table import Table, Column
import argparse
from dataclasses import dataclass, field
from pathlib import Path
import sys
import copy
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import re


# Conversion factors from PARSEC for passbands 
A_GBP_sub_Av: float  = 1.08337
A_GRP_sub_Av: float  = 0.63439


@dataclass(kw_only=True)
class isochrone_data:
    name: str
    object_type: str
    log_age: float
    MH: float
    extinction: float
    distance: float
    MSTO_color: float
    MSTO_magnitude: float
    shifted_color: float
    shifted_magnitude: float
    first_n_elem_isochrone: int
    last_n_elem_isochrone: int


@dataclass(kw_only=True)
class contour_area:
    colors: list[float] = field(default_factory=list)
    magnitudes: list[float] = field(default_factory=list)


@dataclass(kw_only=True)
class CMD_point:
    color: float
    magnitude: float


def parse_flags() -> argparse.Namespace:
    """
    Parse flags from user
    """
    # Create the parser
    parser = argparse.ArgumentParser(description='Plot Zero Age Main-Sequence (ZAMS) and an isochrone to GaiaDR3 data.')
    # Add flags/options
    parser.add_argument('-n', '--object-name', type=str, help='Object name to extract data from', required=True)
    parser.add_argument('-i', '--isochrone-filename', type=str, help='Filename containing isochrone that fits the current object', required=True)
    parser.add_argument('-g', '--gaia-data', type=str, help='Filename containing original Gaia DR3 data to extract Blue Straggler Stars from', required=True)
    parser.add_argument('-z', '--zams-isochrone', type=str, help='Filename containing Zero-Age Main Sequence data', required=True)
    parser.add_argument('-d', '--data-filename', type=str, help = 'Filename containing the data for the isolated isochrone', required=True)
    parser.add_argument('--add-mag-MSTO', type=float, help = 'Quantity to add to the identified MS-turnoff magnitude', required=True)
    parser.add_argument('--subs-mag-MSTO', type=float, help='Quantity to substract to the identified MSTO magnitude', required=True)
    parser.add_argument('--add-color-MSTO', type=float, help='Quantity to substract to MSTO color')
    parser.add_argument('--first-n-elements-isochrone', type=int, default=35, help='Select the first N elements of the isochrone')
    parser.add_argument('--last-n-elements-isochrone', type=int, default=350, help='Select the last N elements of the isochrone')
    parser.add_argument('--first-n-elements-ZAMS', type=int, default=3, help='Select the first N elements for the ZAMS isochrone')
    parser.add_argument('--last-n-elements-ZAMS', type=int, default=12, help='Select the last N elements for the ZAMS isochrone')
    parser.add_argument('--skip-confirmation', action='store_true', help="Do not ask if you want to extract BSS, just does it")
    parser.add_argument('--shift-n-times', type=float, default=2.0, help="Shift ZAMS N times more in color than the shifted isochrone")
    parser.add_argument('--show-only-contour', action='store_true', help="Only show lines that will set the region to select Blue Stragglers")
    parser.add_argument('--text-x-coord', type=float, help="X (color) coordinate to plot text in CMD")
    parser.add_argument('--text-y-coord', type=float, help="Y (color) coordinate to plot text in CMD")
    parser.add_argument('--save-data', action='store_true', help="Sava Stragglers Stars and the region created to select them")
    parser.add_argument('--save-figures', action='store_true', help="Save plots generated within the directory containing 'Gaia Data'")
    parser.add_argument('--figure-format', type=str, default="pdf", help="Format to save plots generated")
    # Parse the command-line arguments
    args = parser.parse_args()
    return args


def check_args(args: argparse.Namespace, data_file: isochrone_data)->None:
    if args.add_color_MSTO is None:
        setattr(args, 'add_color_MSTO', data_file.shifted_color)
    return


def ask_to_continue(text_to_ask: str, max_attempts=10) -> bool | None:
    attempts = 0
    while attempts < max_attempts:
        user_input = input(f"{text_to_ask} (Y/N): ")

        if re.match(r'^\s*y', user_input, re.IGNORECASE):
            return True
        elif re.match(r'^\s*n', user_input, re.IGNORECASE):
            return False
        elif not user_input.strip():
            return True
        else:
            print("Invalid input. Please enter 'Y' or 'N'.")
        attempts += 1
    print(f"Maximum number of attempts reached ({max_attempts!r}). Returning system error.")
    sys.exit(1)


# Pass model data (PARSEC) to empiric data (Gaia)
def pass_derredened_color_to_apparent_color(derredened_color, a_1: float, a_2: float, extinction: float):
    return derredened_color + ((a_1 - a_2)*extinction)


def get_apparent_magnitude_from_abs_magnitude_and_dist(distance: float, abs_magnitude: float, extinction: float)->float:
    return (5*(np.log10(distance)-1)) + extinction + abs_magnitude


def get_data_from_datafile(args: argparse.Namespace)-> isochrone_data | None:
    path_object_datafile = Path(args.data_filename)
    if not path_object_datafile.exists():
        print(f"[!] {str(path_object_datafile)!r} file does not exists. Run previous steps to create it and retry")
        print("    Exiting...")
        sys.exit(1)
    with path_object_datafile.open('r') as f:
        lines = f.readlines()
    for line in lines:
        if line.split()[0].lower() == args.object_name.lower():
            return isochrone_data(name=line.split()[0], object_type=line.split()[3],
                                  log_age=float(line.split()[4]), MH=float(line.split()[5]),
                                  extinction=float(line.split()[6]), distance=float(line.split()[7]),
                                  MSTO_color=float(line.split()[8]), MSTO_magnitude=float(line.split()[9]),
                                  shifted_color=float(line.split()[10]), shifted_magnitude=float(line.split()[11]),
                                  first_n_elem_isochrone=int(line.split()[12]),
                                  last_n_elem_isochrone=int(line.split()[13]))
    print(f"[!] Data for {args.object_name!r} could not be found in {str(path_object_datafile)!r}")
    print("Exiting...")
    sys.exit(1)


def get_data_from_files(args: argparse.Namespace) -> (Table, Table, Table):
    gaia_data = Table.read(args.gaia_data, format="ascii.ecsv")
    ZAMS_data = Table.read(args.zams_isochrone, format="ascii")
    isochrone_fit_data = Table.read(args.isochrone_filename, format="ascii")
    return gaia_data, ZAMS_data, isochrone_fit_data


def fit_simple_line(x1: float, x2: float, y1: float, y2: float):
    try:
        slope = (y2-y1)/(x2-x1)
    except ZeroDivisionError:
        return None, None
    intercept = y1 - (slope * x1)
    return slope, intercept


def get_intersection_between_isochrone_and_line_ZAMS(args: argparse.Namespace, original_data_isochrone: Table, 
                                                     data_file: isochrone_data, cut_value: float) -> CMD_point | None:
    # Create a copy, so we avoid to modify the original data by accident
    data_isochrone: Table = copy.deepcopy(original_data_isochrone)
    data_isochrone = data_isochrone[args.first_n_elements_ZAMS:-args.last_n_elements_ZAMS]
    # Compare colors until they reach the value we are looking for in "color" axis
    for index, data in enumerate(data_isochrone):
        apparent_magnitude_data = get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, data['G_BPmag'], data_file.extinction)
        if index == 0:
            continue
        try:
            prev_apparent_magnitude = get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, data_isochrone[index-1]['G_BPmag'], data_file.extinction) # y1
            current_apparent_magnitude = get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, data_isochrone[index]['G_BPmag'], data_file.extinction) # y2
            if prev_apparent_magnitude <= cut_value <= current_apparent_magnitude:
                useful_index = index
                break
        except IndexError:
            print(f"[!] Error when trying to get index value (index '{index}' and '{index-1}' for ZAMS isochrone)")
            print("     Maybe try to change '--first-n-elements-ZAMS' and '--last-n-elements-ZAMS' for ZAMS isochrone? (see '-h')")
            sys.exit(1)
        if index == len(data_isochrone)-1:
            print("[!] No value found")
    prev_color = data_isochrone[useful_index-1]['G_BPmag'] - data_isochrone[useful_index-1]['G_RPmag'] 
    prev_apparent_color = (pass_derredened_color_to_apparent_color(prev_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction)) - data_file.shifted_color * args.shift_n_times # x1
    current_color = data_isochrone[useful_index]['G_BPmag'] - data_isochrone[useful_index]['G_RPmag']
    current_apparent_color = pass_derredened_color_to_apparent_color(current_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction) - data_file.shifted_color * args.shift_n_times # x2
    slope, intercept = fit_simple_line(prev_apparent_color, current_apparent_color, prev_apparent_magnitude, current_apparent_magnitude)
    if slope is None and intercept is None:
        print("[!] Slope and intercept are 'None'")
        sys.exit(1)
    return CMD_point(color=(cut_value-intercept)/slope, magnitude=cut_value)


def get_intersection_between_isochrone_and_line_fit(args: argparse.Namespace, original_data_isochrone: Table, 
                                                    data_file: isochrone_data, cut_value: float) -> CMD_point | None:
    # Create a copy, so we avoid to modify the original data by accident
    data_isochrone: Table = copy.deepcopy(original_data_isochrone)
    data_isochrone = data_isochrone[args.first_n_elements_isochrone:-args.last_n_elements_isochrone]
    isochrone_color_derredened = data_isochrone['G_BPmag'] - data_isochrone['G_RPmag']
    isochrone_color = np.asarray([pass_derredened_color_to_apparent_color(og_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction) for og_color in isochrone_color_derredened]) - data_file.shifted_color
    isochrone_app_mag = np.asarray([get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, og_mag, data_file.extinction) for og_mag in data_isochrone['G_BPmag']]) - data_file.shifted_magnitude
    if len(isochrone_color) != len(data_isochrone):
        print("[!] Warning! List does not have the same length in 'get_intersection' for 'fit' isochrone")
        sys.exit(1)
    for index in range(0, len(data_isochrone)):
        if index == 0:
            continue
        try: 
            prev_apparent_magnitude = isochrone_app_mag[index-1]
            current_apparent_magnitude = isochrone_app_mag[index]
        except IndexError:
            print(f"[!] Error when trying to get index value (index '{index}' and '{index-1}' for 'fit' isochrone)")
            print("    Maybe try to change '--first-n-elements-isochrone' and '--last-n-elements-isochrone' for fitting isochrone? (see '-h')")
            sys.exit(1)
        if (prev_apparent_magnitude >= cut_value) and (cut_value >= current_apparent_magnitude):
            useful_index = index
            break
        if index == len(data_isochrone)-1:
            print("[!] No value found")
    prev_apparent_color = isochrone_color[useful_index-1]
    current_apparent_color = isochrone_color[useful_index]    
    slope, intercept = fit_simple_line(prev_apparent_color, current_apparent_color, prev_apparent_magnitude, current_apparent_magnitude)
    if slope is None and intercept is None:
        print("[!] Slope and intercept are 'None'")
        sys.exit(1)
    return CMD_point(color=(cut_value-intercept)/slope, magnitude=cut_value)


def get_intersection_between_isochrone_and_line_fit_right_mid(args:argparse.Namespace, original_data_isochrone: Table, 
                                                              data_file: isochrone_data, cut_value: float,
                                                              MSTO_point: CMD_point) -> CMD_point | None:
    """
    Now we want to get the intersection between the magnitude (Y-axis) and the fitting isochrone
    """
    data_isochrone: Table = copy.deepcopy(original_data_isochrone)
    data_isochrone = data_isochrone[args.first_n_elements_isochrone:-args.last_n_elements_isochrone]
    isochrone_color_derredened = data_isochrone['G_BPmag'] - data_isochrone['G_RPmag']
    isochrone_color = np.asarray([pass_derredened_color_to_apparent_color(og_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction) for og_color in isochrone_color_derredened]) - data_file.shifted_color
    isochrone_app_mag = np.asarray([get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, og_mag, data_file.extinction) for og_mag in data_isochrone['G_BPmag']]) - data_file.shifted_magnitude
    mask_data = isochrone_app_mag <= MSTO_point.magnitude + args.add_mag_MSTO
    isochrone_app_mag = isochrone_app_mag[mask_data]
    isochrone_color = isochrone_color[mask_data]
    if len(isochrone_color) != len(isochrone_app_mag):
        print("[!] Warning! List does not have the same length in 'get_intersection' for 'fit' isochrone (mid-right vertex)")
        print(f"    Size magnitude list: {len(isochrone_app_mag)}; size colors list: {len(isochrone_color)}")
        sys.exit(1)
    for index in range(0, len(isochrone_color)):
        if index == 0:
            continue
        try: 
            prev_color = isochrone_color[index-1]
            current_color = isochrone_color[index]
        except IndexError:
            print(f"[!] Error when trying to get index value (index '{index}' and '{index-1}' for 'fit' isochrone)")
            print("    Maybe try to change '--first-n-elements-isochrone' and '--last-n-elements-isochrone' for fitting isochrone? (see '-h')")
            sys.exit(1)
        if (prev_color <= cut_value) and (cut_value <= current_color):
            useful_index = index
            break
        if index == len(isochrone_color)-1:
            print("[!] No value found")
    prev_apparent_magnitude = isochrone_app_mag[useful_index-1]
    current_apparent_magnitude = isochrone_app_mag[useful_index]
    slope, intercept = fit_simple_line(prev_color, current_color, prev_apparent_magnitude, current_apparent_magnitude)
    if slope is None and intercept is None:
        print("[!] Slope and intercept are 'None'")
        sys.exit(1)
    return CMD_point(color=cut_value, magnitude=slope*cut_value + intercept)


def get_data_inside_a_range(args: argparse.Namespace, original_data: Table, isochrone_type: str,
                            MSTO_point: CMD_point, data_file: isochrone_data):
    data_table = copy.deepcopy(original_data)
    if isochrone_type.lower() == "zams":
        data_table = data_table[args.first_n_elements_ZAMS:-args.last_n_elements_ZAMS]
        ZAMS_derredened_color = data_table['G_BPmag'] - data_table['G_RPmag']
        apparent_color_isochrone = np.asarray([pass_derredened_color_to_apparent_color(og_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction) for og_color in ZAMS_derredened_color]) - data_file.shifted_color * args.shift_n_times
        apparent_magnitude_isochrone = np.asarray([get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, og_mag, data_file.extinction) for og_mag in data_table['G_BPmag']])
        top_value = MSTO_point.magnitude - args.subs_mag_MSTO
        bottom_value = MSTO_point.magnitude + args.add_mag_MSTO
        # Filter by top value (delete every star brighter than the vertical top line)
        mask_data = top_value < apparent_magnitude_isochrone
        apparent_magnitude_isochrone = apparent_magnitude_isochrone[mask_data]
        apparent_color_isochrone = apparent_color_isochrone[mask_data]
        # Filter by bottom value (delete every object fainter than the vertical bottom line)
        mask_data = apparent_magnitude_isochrone < bottom_value
        apparent_magnitude_isochrone = apparent_magnitude_isochrone[mask_data]
        apparent_color_isochrone = apparent_color_isochrone[mask_data]
    elif isochrone_type.lower() == "fit":
        data_table = data_table[args.first_n_elements_isochrone:-args.last_n_elements_isochrone]
        fit_isochrone_og_color = data_table['G_BPmag'] - data_table['G_RPmag']
        apparent_color_isochrone = np.asarray([pass_derredened_color_to_apparent_color(og_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction) for og_color in fit_isochrone_og_color]) - data_file.shifted_color
        apparent_magnitude_isochrone = np.asarray([get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, og_mag, data_file.extinction) for og_mag in data_table['G_BPmag']]) - data_file.shifted_magnitude
        bottom_value = MSTO_point.magnitude + args.add_mag_MSTO
        right_side_value = MSTO_point.color + data_file.shifted_color
        # Filter by all the values above the MSTO point
        mask_data = apparent_magnitude_isochrone < bottom_value
        apparent_magnitude_isochrone = apparent_magnitude_isochrone[mask_data]
        apparent_color_isochrone = apparent_color_isochrone[mask_data]
        # Filter by everything with a color lower than the max allowed right side
        mask_data = apparent_color_isochrone < right_side_value
        apparent_magnitude_isochrone = apparent_magnitude_isochrone[mask_data]
        apparent_color_isochrone = apparent_color_isochrone[mask_data]
    else:
        print(f"[!] You have provided an invalid data type in 'get_data_inside_range' function ({isochrone_type!r})")
        sys.exit(1)
    return apparent_color_isochrone, apparent_magnitude_isochrone


def set_selection_contour_area(args: argparse.Namespace, data_file: isochrone_data,
                               ZAMS_isochrone: Table, fit_isochrone: Table) -> contour_area:
    contour = contour_area()
    # Get the MSTO-point from in apparent magnitude (PARSEC is in absolute magnitude) with a non-derredened color
    MSTO_point = CMD_point(color=pass_derredened_color_to_apparent_color(data_file.MSTO_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction), 
                           magnitude=get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, data_file.MSTO_magnitude, data_file.extinction))
    # 1) Get top right vertex
    top_right_vertex = CMD_point(color=MSTO_point.color + args.add_color_MSTO, magnitude=MSTO_point.magnitude - args.subs_mag_MSTO)
    contour.colors.append(top_right_vertex.color)
    contour.magnitudes.append(top_right_vertex.magnitude)
    top_left_vertex = get_intersection_between_isochrone_and_line_ZAMS(args, ZAMS_isochrone, data_file, MSTO_point.magnitude - args.subs_mag_MSTO)
    # 2) Get top left vertex
    contour.colors.append(top_left_vertex.color)
    contour.magnitudes.append(top_left_vertex.magnitude)
    # 3) Get the ZAMS that is at the left side, between the vertical lines
    left_ZAMS_color, left_ZAMS_magnitude = get_data_inside_a_range(args, ZAMS_isochrone, "zams", MSTO_point, data_file)
    for c, m in zip(left_ZAMS_color, left_ZAMS_magnitude):
        contour.colors.append(c)
        contour.magnitudes.append(m)
    # 4) Get bottom left vertex
    bottom_left_vertex = get_intersection_between_isochrone_and_line_ZAMS(args, ZAMS_isochrone, data_file, MSTO_point.magnitude + args.add_mag_MSTO)
    contour.colors.append(bottom_left_vertex.color)
    contour.magnitudes.append(bottom_left_vertex.magnitude)
    get_intersection_between_isochrone_and_line_fit(args, fit_isochrone, data_file, MSTO_point.magnitude + args.add_mag_MSTO)
    # 5) Get bottom right vertex
    bottom_right_vertex = get_intersection_between_isochrone_and_line_fit(args, fit_isochrone, data_file, MSTO_point.magnitude + args.add_mag_MSTO)
    contour.colors.append(bottom_right_vertex.color)
    contour.magnitudes.append(bottom_right_vertex.magnitude)
    # 6) Get the shifted isochrone upper than the MSTO and colder than the right-side in the CMD
    right_fit_color, right_fit_magnitude = get_data_inside_a_range(args, fit_isochrone, "fit", MSTO_point, data_file)
    for c, m in zip(right_fit_color, right_fit_magnitude):
        contour.colors.append(c)
        contour.magnitudes.append(m)
    # 7) Get the mid-right point. Where the fitting isochrone touches the color-limit right side
    mid_right_vertex = get_intersection_between_isochrone_and_line_fit_right_mid(args, fit_isochrone, data_file, MSTO_point.color + args.add_color_MSTO , MSTO_point)
    contour.colors.append(mid_right_vertex.color)
    contour.magnitudes.append(mid_right_vertex.magnitude)
    # Finally, enclose the area repeating the first point
    contour.colors.append(contour.colors[0])
    contour.magnitudes.append(contour.magnitudes[0])
    return contour


def plot_before_extracting(args: argparse.Namespace, gaia_data: Table, ZAMS_isochrone: Table,
                           fit_isochrone: Table, data_file: isochrone_data, add_BSS_region: bool = False,
                           plot_BSS: bool = False, text_number: int = 1)->None | Table:
    text_size = 30
    fig, ax = plt.subplots()
    plt.scatter(gaia_data['bp_rp'], gaia_data['phot_bp_mean_mag'], color="black", s= 2.)
    # Pass theorical data (PARSEC) to empirical data (Gaia) for ZAMS
    ZAMS_derredened_color = ZAMS_isochrone['G_BPmag'] - ZAMS_isochrone['G_RPmag']
    ZAMS_apparent_color = [pass_derredened_color_to_apparent_color(og_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction) for og_color in ZAMS_derredened_color]
    ZAMS_apparent_magnitude = [get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, og_mag, data_file.extinction) for og_mag in ZAMS_isochrone['G_BPmag']]
    # Pass theorical data (PARSEC) to empirical data (Gaia) for the fitting isochrone
    fit_isochrone_derredened_color = fit_isochrone['G_BPmag'] - fit_isochrone['G_RPmag']
    fit_apparent_color = [pass_derredened_color_to_apparent_color(og_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction) for og_color in fit_isochrone_derredened_color]
    fit_apparent_magnitude = [get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, og_mag, data_file.extinction) for og_mag in fit_isochrone['G_BPmag']]
    # Cut the data so they can be inside a nice plot
    fixed_ZAMS_apparent_color = ZAMS_apparent_color[args.first_n_elements_ZAMS:-args.last_n_elements_ZAMS]
    fixed_ZAMS_apparent_magnitude = ZAMS_apparent_magnitude[args.first_n_elements_ZAMS:-args.last_n_elements_ZAMS]
    fixed_fitting_isochrone_color = fit_apparent_color[args.first_n_elements_isochrone:-args.last_n_elements_isochrone]
    fixed_fitting_isochrone_magnitude = fit_apparent_magnitude[args.first_n_elements_isochrone:-args.last_n_elements_isochrone]
    # Plot shifted ZAMS and isochrone
    shifted_isochrone_color = np.asarray(fixed_fitting_isochrone_color) - data_file.shifted_color
    shifted_isochrone_mag = np.asarray(fixed_fitting_isochrone_magnitude) - data_file.shifted_magnitude
    shifted_ZAMS_color = np.asarray(fixed_ZAMS_apparent_color) - data_file.shifted_color * args.shift_n_times
    shifted_ZAMS_mag = np.asarray(fixed_ZAMS_apparent_magnitude)
    if not args.show_only_contour:
        plt.plot(fixed_ZAMS_apparent_color, fixed_ZAMS_apparent_magnitude, linestyle="-", color="lightsalmon")
        plt.plot(fixed_fitting_isochrone_color, fixed_fitting_isochrone_magnitude, linestyle="-", color="orangered")
    plt.plot(shifted_ZAMS_color, shifted_ZAMS_mag, linestyle = "--", color = "orange")
    plt.plot(shifted_isochrone_color, shifted_isochrone_mag, linestyle = "--", color = "red")
    plt.axhline(y=get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, data_file.MSTO_magnitude, data_file.extinction)+args.add_mag_MSTO, 
               color='grey', linestyle='--')
    plt.axhline(y=get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, data_file.MSTO_magnitude, data_file.extinction)-args.subs_mag_MSTO, 
               color='grey', linestyle='--')
    plt.axvline(x=pass_derredened_color_to_apparent_color(data_file.MSTO_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction)+args.add_color_MSTO,
                linestyle='--', color = 'grey')
    plt.plot(pass_derredened_color_to_apparent_color(data_file.MSTO_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction),
             get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, data_file.MSTO_magnitude, data_file.extinction),
             marker="X", markersize=12, color="magenta",  markeredgecolor="black")
    if args.text_x_coord is not None and args.text_y_coord is not None and text_number == 1:
        text = r'${\rm NGC2141}$'
        text +='\n'
        text+= r'${\log ({\rm Age}/{\rm yr})} = $' + str(round(data_file.log_age, 2))
        text+= '\n'
        text += r"$[{\rm M}/{\rm H}] = $"+str(round(data_file.MH, 2))
        text += '\n'
        text += r'${\rm A}_{\rm V} = $' + str(round(data_file.extinction, 3))
        plt.text(args.text_x_coord, args.text_y_coord, text, size=text_size, color = 'gray')
        # 1.7 19.4
    if add_BSS_region:

        contour = set_selection_contour_area(args, data_file, ZAMS_isochrone, fit_isochrone)
        if text_number == 2:
            plt.plot(contour.colors, contour.magnitudes, marker='o', color = "dodgerblue", linestyle="none")

            if args.text_x_coord is not None and args.text_y_coord is not None:
                text = f"Contour points: {len(contour.colors)}"
                plt.text(args.text_x_coord, args.text_y_coord, text, size=text_size-2, color="gray")
            
        if plot_BSS:
            plt.plot(contour.colors, contour.magnitudes, linestyle="-", color="royalblue")
            plt.fill(contour.colors, contour.magnitudes, color='dodgerblue', alpha=0.3)
            BSS = get_Blue_Stragglers(gaia_data, contour)
            if args.text_x_coord is not None and args.text_y_coord is not None and text_number==3:
                text = f"Stars inside\narea: {len(BSS)}"
                plt.text(args.text_x_coord, args.text_y_coord, text, size=text_size-2, color = "gray")

            color_bss = BSS['phot_bp_mean_mag'] - BSS['phot_rp_mean_mag']
            magnitude_BSS = BSS['phot_bp_mean_mag']
            MSTO_point = CMD_point(color=pass_derredened_color_to_apparent_color(data_file.MSTO_color, A_GBP_sub_Av, A_GRP_sub_Av, data_file.extinction), 
                           magnitude=get_apparent_magnitude_from_abs_magnitude_and_dist(data_file.distance, data_file.MSTO_magnitude, data_file.extinction))
            stragglers = get_Straggler_type(args, BSS, MSTO_point)
            if text_number == 3:
                color_stragglers_for_plot = "blue"
            else:
                color_stragglers_for_plot = get_colors_to_plot_for_stragglers(stragglers)
            color_stragglers = stragglers['phot_bp_mean_mag'] - stragglers['phot_rp_mean_mag']
            magnitude_stragglers = stragglers['phot_bp_mean_mag']
            plt.scatter(color_stragglers, magnitude_stragglers, color=color_stragglers_for_plot)
            n_bss, n_yss, n_rss = get_number_of_stragglers_type(stragglers)
            if args.text_x_coord is not None and args.text_y_coord is not None and text_number==4:
                text = f"BSS: {n_bss}\nYSS: {n_yss}\nRSS: {n_rss}"
                plt.text(args.text_x_coord+0.1, args.text_y_coord, text, size = text_size+1, color = "gray")
            
    plt.gca().invert_yaxis()
    plt.xlabel(r"${\rm G}_{\rm BP} - {\rm G}_{\rm RP}$", size=33)
    plt.ylabel(r"${\rm G}_{\rm BP}$")
    # We change the fontsize of minor ticks label 
    ax = plt.gca()
    # Increase the size of tick labels
    ax.tick_params(axis='both', which='major', labelsize=27)
    ax.tick_params(axis='both', which='minor', labelsize=27)
    if args.save_figures:
        data_path = Path(args.gaia_data)
        path_save_data = data_path.parent / f"{text_number}_{args.object_name.lower().replace(' ', '_')}_BSS_progress.{args.figure_format}"
        plt.savefig(path_save_data, format=args.figure_format)
        print(f"[+] Figure saved at {str(data_path.parent)!r}")
    plt.show()
    plt.close()
    if add_BSS_region and plot_BSS:
        return contour, stragglers
    else:
        return None, None


def get_number_of_stragglers_type(stragglers_data: Table) -> (float, float, float):
    bss, yss, rss = 0, 0, 0
    for star in stragglers_data:
        if star['straggler_type'] == "blue":
            bss +=1
        elif star['straggler_type'] == "yellow":
            yss += 1
        elif star['straggler_type'] == "red":
            rss +=1
    return bss, yss, rss


def select_straggler_type(args: argparse.Namespace, data, MSTO_point: CMD_point):
    """
    Selects the straggler type (Blue, Yellow or Red)
    """
    color = data['phot_bp_mean_mag'] - data['phot_rp_mean_mag']
    magnitude = data['phot_bp_mean_mag']
    if (color > MSTO_point.color) and (magnitude < (MSTO_point.magnitude-(args.subs_mag_MSTO/2.0))):
        return "red"
    elif (color > MSTO_point.color) and (magnitude >= (MSTO_point.magnitude-(args.subs_mag_MSTO/2.0))):
        return "yellow"
    else:
        return "blue"


def get_colors_to_plot_for_stragglers(data: Table, if_color_blue="blue", if_color_red="crimson", if_color_yellow="gold")-> list[str] | None:
    colors_to_plot: list[str] = []
    for data_it in data:
        if data_it['straggler_type'] == "blue":
            colors_to_plot.append(if_color_blue)
        if data_it['straggler_type'] == "yellow":
            colors_to_plot.append(if_color_yellow)
        if data_it['straggler_type'] == "red":
            colors_to_plot.append(if_color_red)
    if len(colors_to_plot) != len(data):
        print(f"[!] Error. Color list to plot and data does not have the same length ({len(colors_to_plot)!r} and {len(data)!r}, respectively)")
        sys.exit(1)
    return colors_to_plot


def get_Straggler_type(args: argparse.Namespace, original_data: Table, MSTO_point: CMD_point,
                       new_column_name: str ="straggler_type")-> Table | None:
    data_table = copy.deepcopy(original_data)
    data_to_add: list[str] = []
    for data in data_table:
        data_to_add.append(select_straggler_type(args, data, MSTO_point))
    if len(data_to_add) != len(data_table):
        print("[!] Warning! Straggler type list and Possible BSS list does not have the same length")
        sys.exit(1)
    new_col_with_straggler_type = Column(data=data_to_add, name=new_column_name, dtype="S7")
    data_table.add_column(new_col_with_straggler_type)
    return data_table


def get_Blue_Stragglers(original_gaia_data: Table, contour: contour_area)->Table:
    # Create a copy since the data will be modified, so we can compare later with original data
    gaia_data = copy.deepcopy(original_gaia_data)
    # Get the vertices coordinates for the Blue Straggler region
    list_of_points_contour: list[Point] = []
    for c_contour, m_contour in zip(contour.colors, contour.magnitudes):
        list_of_points_contour.append(Point(c_contour, m_contour))
    list_of_points_contour = np.asarray(list_of_points_contour)
    area_bss = Polygon(list_of_points_contour)
    # Get the coordinates for the data in Color Magnitud Diagram
    list_of_points_data: list[Point] = []
    for data in gaia_data:
        list_of_points_data.append(Point(data['phot_bp_mean_mag'] - data['phot_rp_mean_mag'], data['phot_bp_mean_mag']))
    mask_bss = np.asarray(area_bss.contains(list_of_points_data))
    BSS_stars = gaia_data[mask_bss]
    return BSS_stars


def get_contour_data_content(contour: contour_area) -> str:
    data_to_write = "# N_point color G_BP_apparent\n"
    counter=0
    for c, m in zip(contour.colors, contour.magnitudes):
        counter+=1
        data_to_write += f"{counter} {c} {m}\n"
    return data_to_write


def save_contour(args: argparse.Namespace, contour: contour_area)->None:
    isochrone_data_path = Path(args.gaia_data)
    path_to_save_data = isochrone_data_path.parent / f"{args.object_name.lower().replace(' ','_')}_contour_candidate_BSS.dat"
    if not path_to_save_data.exists():
        data_to_write = get_contour_data_content(contour)
        with path_to_save_data.open('w') as f:
            f.write(data_to_write)
        print(f"[+] Contour data saved at {str(path_to_save_data)!r}")
    else:
        print(f"[!] Warning. File {str(path_to_save_data)!r} already exists")
        want_to_continue = ask_to_continue("Do you want to replace the file?")
        if want_to_continue:
            data_to_write = get_contour_data_content(contour)
            with path_to_save_data.open('w') as f:
                f.write(data_to_write)
            print(f"[+] Contour data saved at {str(path_to_save_data)!r}")
        else:
            print("[!] Contour data was not saved")
    return


def save_stragglers_stars(args: argparse.Namespace, data_stragglers: Table)->None:
    isochrone_data_path = Path(args.gaia_data)
    path_to_save_data = isochrone_data_path.parent / f"{args.object_name.lower().replace(' ','_')}_possible_BSS.dat"
    if not path_to_save_data.exists():
        data_stragglers.write(path_to_save_data, format="ascii.ecsv")
        print(f"[+] Stragglers data saved at {str(path_to_save_data)!r}")
    else:
        print(f"[!] Warning. File {str(path_to_save_data)!r} already exists")
        want_to_continue = ask_to_continue("Do you want to replace the file?")
        if want_to_continue:
            data_stragglers.write(path_to_save_data, format="ascii.ecsv", overwrite=want_to_continue)
            print(f"[+] Straggler data saved at {str(path_to_save_data)!r}")
        else:
            print("[!] Straggler data was not saved")
    return


def save_stragglers_and_contour(args: argparse.Namespace, contour: contour_area,
                                stragglers_data: Table) -> None:
    if not args.save_data:
        return
    save_contour(args, contour)
    save_stragglers_stars(args, stragglers_data)
    return


def main()->None:
    # Get flags from user
    args = parse_flags()
    # Get data from file containing data
    data_inside_file = get_data_from_datafile(args)
    # Check some values have been provided. If not, use the ones provided into the file containing data
    check_args(args, data_inside_file)
    # Get the object data from Gaia DR3, ZAMS generated and isochrone fit
    gaia_data, zams_isochrone, fit_isochrone = get_data_from_files(args)
    _,_ = plot_before_extracting(args, gaia_data, zams_isochrone, fit_isochrone, data_inside_file, text_number=1)
    _,_ = plot_before_extracting(args, gaia_data, zams_isochrone, fit_isochrone, data_inside_file, add_BSS_region=True, text_number=2)
    _,_ = plot_before_extracting(args, gaia_data, zams_isochrone, fit_isochrone, data_inside_file, add_BSS_region=True, plot_BSS = True, text_number=3)
    contour, stragglers = plot_before_extracting(args, gaia_data, zams_isochrone, fit_isochrone, data_inside_file, add_BSS_region=True, plot_BSS=True, text_number=4)
    save_stragglers_and_contour(args, contour, stragglers)


if __name__ == "__main__":
    main()
