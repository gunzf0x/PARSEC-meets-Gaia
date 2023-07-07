from dataclasses import dataclass, field
from pathlib import Path
import matplotlib.pyplot as plt
import sys
import argparse
import requests
from astropy.table import Table
import numpy as np


# Conversion factors from PARSEC for passbands 
A_GBP_sub_Av: float  = 1.08337
A_GRP_sub_Av: float  = 0.63439

@dataclass(kw_only=True)
class isochrone_combinations:
    isochrone_id: int
    log_age: float
    MH: float
    lines_containing_data: list[str] = field(default_factory=list)


@dataclass(kw_only=True)
class data_to_plot:
    id: int
    log_age: float
    MH: float
    color: list[float] = field(default_factory=list)
    magnitude: list[float] = field(default_factory=list)


@dataclass(kw_only=True)
class CMD_coords:
    color: float
    magnitude: float


@dataclass(kw_only=True)
class data_source:
    study_url: str
    distance: float 
    A_v: float

def check_if_user_has_provided_arguments()->None:
     # Check if any command-line arguments were passed
    if len(sys.argv) == 1:
        print("[-] No arguments provided. Use -h or --help for more information.")
        sys.exit(1)
    return


def parse_flags():
    """
    Parse flags from user
    """
    # Create the parser
    parser = argparse.ArgumentParser(description='Plot isochrones from PARSEC to Gaia DR3 data.')

    # Add flags/options
    parser.add_argument('-i', '--isochrones-filename', type=str, help='Filename containing PARSEC isochrones.', required=True)
    parser.add_argument('-f', '--data-filename', type=str, help='Filename containing GAIA data')
    parser.add_argument('-d', '--distance', type = float, help='Distance to the cluster or object(in pc)')
    parser.add_argument('-av', '--extinction', type=float, help = "A_v extinction for cluster (in mag)")
    parser.add_argument('-dm', '--d-modulus', type=float, help="Distance modulus (mag)")
    parser.add_argument('-s', '--select-isochrone', type=int, help='Show a specific isochrone from the set of isochrones identified')
    parser.add_argument('-t', '--get-turnoff-candidate', type=int, help="Get a turnoff candidate")
    parser.add_argument('--first-n-elements', type=int, default=30, help='Select the first N elements of the isochrone')
    parser.add_argument('--last-n-elements', type=int, default=260, help='Select the last N elements of the isochrone')
    parser.add_argument('--save-isochrones', action='store_true', help="Save different PARSEC isochrones found in separated files")
    parser.add_argument('--print-isochrone-details', action='store_true', help="Print some details about the isochrones found in PARSEC file")
    parser.add_argument('--cluster-type', type=str, help="Cluster Type: {Globular Cluster, GC, OpenCluster, OC, Other}")
    parser.add_argument('--shift-color', type=float, help="When one isochrone is selected, create a second one shifting its color")
    parser.add_argument('--shift-mag', type=float, help="When one isochrone is selected, create a second one shifting its color")

    # Parse the command-line arguments
    args = parser.parse_args()
    return args


def check_arguments_provided(args)->None:
    if args.data_filename is not None:
        if args.distance is None and args.d_modulus is None:
            print("[-] You have provided data but you have not provided distance.")
            print("    Neither a distance ('--distance') or distance modulus ('--d-modulus')")
            print("    Distance is required to fit the data to PARSEC Isochrones")
            sys.exit(1)
        if args.extinction is None:
            print("[-] You have provided data but you have not provided extinction A_v ('--extinction')")
            print("    Extinction is required to fit data to PARSEC Isochrones")
            sys.exit(1)
    return


def read_online_dat_file(url: str) -> list[str] | None:
    response = requests.get(url)
    if response.status_code != 200:
        print(f"[-] Something happened when trying to request data to {url!r}: status code {response.status_code}")
        sys.exit(1)
    data = response.text
    lines = data.split('\n')
    return lines


def read_PARSEC_isochrone_file(args):
    file_path = args.isochrones_filename
    print("[+] Reading file containing isochrones...")
    lines = []
    counter_isochrones = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                if line.startswith('# Zini') or line.startswith("#isochrone terminated"):
                    lines.append(line.strip())
                    if line.startswith("# Zini"):
                        counter_isochrones += 1
                else:
                    continue
            else:
                lines.append(line.strip())
    print(f"[+] Total isochrones found in original file: {counter_isochrones}")
    return lines


def get_different_log_age_and_MH_combinations(lines: list[str], header_file : str) -> list[isochrone_combinations]: 
    id_isochrone: int = 1
    list_of_isochrones = []
    for number_line, line in enumerate(lines):
        if line.startswith("# Zini") or line.startswith("#isochrone terminated"):
            if id_isochrone != 1:
                list_of_isochrones.append(current_isochrone_combination)
            line_number_with_header = number_line
            continue
        if number_line == line_number_with_header + 1:
            log_age = float(line.split()[2])
            MH = float(line.split()[1])
            current_isochrone_combination = isochrone_combinations(isochrone_id = id_isochrone,
                                                                   log_age = log_age,
                                                                   MH = MH)
            current_isochrone_combination.lines_containing_data.append(header_file)
            id_isochrone += 1
        try:
            current_log_age = float(line.split()[2])
            current_MH = float(line.split()[1]) 
            if current_log_age == current_isochrone_combination.log_age and current_MH == current_isochrone_combination.MH:
                current_isochrone_combination.lines_containing_data.append(line)
        except ValueError:
            continue
    print(f"[+] Isochrones extracted: {len(list_of_isochrones)} extracted")
    return list_of_isochrones


def print_obtained_isochrones(list_of_isochrones: list[isochrone_combinations]) -> None:
    if len(list_of_isochrones) == 0:
        print("[-] No ischrones to print")
        return
    for isochrone in list_of_isochrones:
        print(f"    {isochrone.isochrone_id} - log10 age: {isochrone.log_age} - [M/H]: {isochrone.MH} - lines containing data: {len(isochrone.lines_containing_data)}")
    return


def correct_color(color_to_correct, a_1: float, a_2: float):
    return color_to_correct - (a_1 - a_2)


def get_abs_magnitude_given_a_distance(distance: float, apparent_magnitude, Av):
    return apparent_magnitude - (5 * (np.log10(distance) - 1)) - Av


def get_abs_magnitude_given_distance_modulus(distance_modulus: float, apparent_magnitude, Av):
    return apparent_magnitude - distance_modulus - Av


def correct_data(args, magnitude_to_correct, color_to_correct, 
                 A_GBP_sub_Av: float, A_GRP_sub_Av: float):
    if args.distance is not None:
        absolute_magnitude = get_abs_magnitude_given_a_distance(args.distance, magnitude_to_correct, args.extinction)
    elif args.d_modulus is not None:
        absolute_magnitude = get_abs_magnitude_given_distance_modulus(args.d_modulus, magnitude_to_correct, args.extinction)
    color_corrected = correct_color(color_to_correct, A_GBP_sub_Av, A_GRP_sub_Av)
    return color_corrected, absolute_magnitude


def save_isochrones(args, list_of_isochrones: list[isochrone_combinations],
                    directory_containing_isochrones_data: str = "isochrones"
                    ) -> list[data_to_plot]:
    filename_containing_data = args.isochrones_filename
    list_of_data_to_plot = []
    parent_directory = Path(filename_containing_data).parent
    isochrones_folder = parent_directory / directory_containing_isochrones_data
    if not isochrones_folder.exists():
        print(f"[+] '{directory_containing_isochrones_data}' not found. Creating it...")
        isochrones_folder.mkdir()
    for isochrone in list_of_isochrones:
        data_for_future_plot = data_to_plot(id=isochrone.isochrone_id, log_age=isochrone.log_age, MH=isochrone.MH)
        for line_containing_data in isochrone.lines_containing_data:
            if line_containing_data.startswith("#"):
                continue
            try:
                G_RP = float(line_containing_data.split()[-1])
                G_BP = float(line_containing_data.split()[-2])
                color = G_BP - G_RP
                data_for_future_plot.color.append(color)
                data_for_future_plot.magnitude.append(G_BP)
            except ValueError:
                print(f"[-] Error trying when trying to convert '{line_containing_data.split()[-1]}' , '{line_containing_data.split()[-2]}' to floats")
        list_of_data_to_plot.append(data_for_future_plot) 
        filename_to_save_data = isochrones_folder / f"{isochrone.isochrone_id}_log_age_{str(isochrone.log_age).replace('.','_')}_{str(isochrone.MH).replace('.','_')}.dat"
        if args.save_isochrones:
            with open(filename_to_save_data, 'w') as f:
                for data in isochrone.lines_containing_data:
                    f.write(data+"\n")
            print(f"[+] Data succesfully written into '{isochrones_folder}' directory")
    return list_of_data_to_plot


def plot_multiple_isochrones(args, list_of_data_to_plot: list[data_to_plot]) -> None:
    colors_to_plot = ['red', 'green', 'blue', 'orange', 'purple', 'yellow', 'cyan', 'magenta']  
    markers_list = ['o', '^', '*', 's', 'd']
    _, ax = plt.subplots()
    for index, data in enumerate(list_of_data_to_plot):
        color = data.color
        marker_selected = markers_list[index % len(markers_list)]
        magnitude = data.magnitude
        color_to_plot = colors_to_plot[index % len(colors_to_plot)]  # Select color cyclically
        ax.plot(color[args.first_n_elements:-args.last_n_elements], magnitude[args.first_n_elements:-args.last_n_elements], 
                color=color_to_plot, marker=marker_selected, markersize = 5, 
                label = f"{data.id}) log Age: {data.log_age}; [M/H]: {data.MH}", linestyle='None')
    legend = plt.legend()
    for text in legend.get_texts():
        text.set_fontsize('large')
    plt.gca().invert_yaxis()
    plt.show()
    plt.close()


def plot_one_isochrone(args, list_of_data_to_plot: list[data_to_plot], data_index: int,
                       turnoff_points: CMD_coords | None,
                       data_color=None, data_magnitude=None) -> None:
    color = list_of_data_to_plot[data_index].color[args.first_n_elements:-args.last_n_elements]
    magnitude = list_of_data_to_plot[data_index].magnitude[args.first_n_elements:-args.last_n_elements]
    plt.plot(color, magnitude, color="orange")
    if args.shift_color and args.shift_mag:
        plt.plot(np.asarray(color)-args.shift_color, np.asarray(magnitude)-args.shift_mag, linestyle="--", color="magenta")
    if turnoff_points is not None:
        plt.axhline(y = turnoff_points.magnitude, color = 'gray', linestyle = '--', alpha=0.5)
        plt.axvline(x = turnoff_points.color, color = 'gray', linestyle = '--', alpha=0.5)

        plt.plot(turnoff_points.color, turnoff_points.magnitude, 'rX', markersize=10)
    if data_color is not None and data_magnitude is not None:
        plt.scatter(data_color, data_magnitude)
    plt.gca().invert_yaxis()
    plt.show()
    plt.close()


def plot_isochrones(args, list_of_data_to_plot: list[data_to_plot]) ->None:
    # If the user wants to plot only one plot, plot it
    if isinstance(args.select_isochrone, int):
        which_data_to_plot = args.select_isochrone - 1
        print(f"    [*] Isochrone number: {which_data_to_plot+1}")
        print(f"    [*] log10 (Age/yr): {list_of_data_to_plot[which_data_to_plot].log_age:.3f}")
        print(f"    [*] [M/H]: {list_of_data_to_plot[which_data_to_plot].MH:.3f}")
        if which_data_to_plot < len(list_of_data_to_plot) and which_data_to_plot >= 0:
            turnoff_point = get_MSTO_Turnoff(args, list_of_data_to_plot, which_data_to_plot)
            if args.data_filename:
                plot_data = True
                gaia_data = Table.read(args.data_filename, format="ascii.ecsv")
                # Get the parameters from Gaia DR3 data
                G_data = gaia_data['phot_g_mean_mag']
                G_BP_data = gaia_data['phot_bp_mean_mag']
                G_RP_data = gaia_data['phot_rp_mean_mag']
                color_data = gaia_data['bp_rp']
                # Correct data by extinction and pass it from apparent to absolute magnitude
                corrected_color, absolute_magnitude_data = correct_data(args, G_BP_data, color_data, A_GBP_sub_Av, A_GRP_sub_Av) 
                plot_one_isochrone(args, list_of_data_to_plot, which_data_to_plot, turnoff_point, 
                                   data_color=corrected_color, data_magnitude=absolute_magnitude_data)
            else:
                plot_one_isochrone(args, list_of_data_to_plot, which_data_to_plot, turnoff_point)
        else:
            print(f"[-] You want to get isochrone number {which_data_to_plot+1}. However, only up to {len(list_of_data_to_plot)} are available")
            sys.exit(1)
    else:
        plot_multiple_isochrones(args, list_of_data_to_plot)


def get_MSTO_Turnoff(args, data_it: list[data_to_plot], data_index: int) -> CMD_coords | None:
    """
    Get the MSTO point based on the inclination 
    """
    color = data_it[data_index].color[args.first_n_elements:-args.last_n_elements]
    magnitude = data_it[data_index].magnitude[args.first_n_elements:-args.last_n_elements]
    possible_MSTO_points = []
    for i in range(0, len(color)-2):
        color1, color2 = color[i], color[i+1]
        mag1, mag2 = magnitude[i], magnitude[i+1]
        try:
            slope = (mag2-mag1) / (color2-color1)
        except ZeroDivisionError:
            continue
        if slope < 0:
            possible_MSTO_points.append(CMD_coords(color=color1, magnitude=mag1))
    if len(possible_MSTO_points) == 0:
        print("[-] 0 Turnoff point candidates found")
        return
    if args.get_turnoff_candidate is not None:
        if args.get_turnoff_candidate > len(possible_MSTO_points):
            print(f"[-] Cannot return Turnoff Candidate number {args.get_turnoff_candidate} since only {len(possible_MSTO_points)} are available")
            sys.exit(1)
        color_returned = possible_MSTO_points[args.get_turnoff_candidate-1].color 
        magnitude_returned = possible_MSTO_points[args.get_turnoff_candidate-1].magnitude
        print(f"\n[+] Turnoff Point returned: ({color_returned:.3f}, {magnitude_returned:.3f})")
        return possible_MSTO_points[args.get_turnoff_candidate-1]
    else:
        return


def main():
    # Check if the user has provided arguments. Otherwise print a "help" message
    check_if_user_has_provided_arguments()
    # Get flags from user
    args = parse_flags()
    # Check arguments provided by the user
    check_arguments_provided(args)
    # Get PARSEC data
    lines_file = read_PARSEC_isochrone_file(args)
    header_file = lines_file[0]
    # Get all the different isochrones contained in ".dat" file from PARSEC service
    list_of_isochrones = get_different_log_age_and_MH_combinations(lines_file, header_file)
    # Print some details about the isoschrones found if the user wants to
    if args.print_isochrone_details:
        print_obtained_isochrones(list_of_isochrones)
    # If desired, save the different isochrones found
    data_to_plot_var = save_isochrones(args, list_of_isochrones)
    # Plot the isochrone
    plot_isochrones(args, data_to_plot_var)


if __name__ == "__main__":
    main()
