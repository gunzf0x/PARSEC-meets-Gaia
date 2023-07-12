import argparse
from dataclasses import dataclass, field
import sys
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from pathlib import Path


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
class ZAMS:
    id: int
    log_age: float
    MH: float
    color: float
    magnitude: float


def parse_flags() -> argparse.Namespace:
    """
    Parse flags from user
    """
    # Create the parser
    parser = argparse.ArgumentParser(description='Plot Zero Age Main-Sequence (ZAMS) and an isochrone to GaiaDR3 data.')
    # Add flags/options
    parser.add_argument('-i', '--isochrones-filename', type=str, help='Filename containing PARSEC isochrones to create ZAMS', required=True)
    parser.add_argument('-f', '--data-filename', type=str, help='Filename containing original data to compare ZAMS')
    parser.add_argument('-c', '--compare-isochrone', type=str, help='Filename containing ONE isochrone to compare with the ZAMS created. This isochrone should have the same metallicity as the ones provided to create he ZAMS')
    parser.add_argument('-d', '--distance', type = float, help='Distance to the cluster or object(in pc)')
    parser.add_argument('-av', '--extinction', type=float, help = "A_v extinction for cluster (in mag)")
    parser.add_argument('-dm', '--d-modulus', type=float, help="Distance modulus (mag)")
    parser.add_argument('--first-n-elements', type=int, default=0, help='Select the first N elements of the isochrone')
    parser.add_argument('--last-n-elements', type=int, default=3, help='Select the last N elements of the isochrone')
    parser.add_argument('--save-isochrones', action='store_true', help="Save ZAMS and isochrone cut found in files")
    parser.add_argument('--object-name-to-save', type = str, help='Object name to save the data. For example, if ou provide "NGC104", output will be saved as "ngc104_ZAMS.dat"')
    parser.add_argument('--show-first-n-elem-isochrone', type=int, default=32, help="Select the first N elements to show for isochrone provided to compare with ZAMS")
    parser.add_argument('--show-last-n-elem-isochrone', type=int, default=350, help="Select the last N elements to show for isochrone provided to compare with ZAMS")

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


def read_PARSEC_isochrone_file(args)->list[str]:
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
    print(f"[+] Total isochrones to generate Zero-Age Main Sequence: {counter_isochrones}")
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
    if len(list_of_isochrones) == 0:
        print("[-] No isochrones could be detected")
        sys.exit(1)
    print(f"[+] Isochrones extracted: {len(list_of_isochrones)} extracted")
    return list_of_isochrones


def get_ZAMS_data_by_label(args, isochrones_list = list[isochrone_combinations],
                                 label_var: int=1)->list[ZAMS] | None:
    list_of_data_to_plot = []
    for index, isochrone in enumerate(isochrones_list):
        data_for_future_plot = data_to_plot(id=isochrone.isochrone_id, log_age=isochrone.log_age, MH=isochrone.MH)
        if index == 0:
            static_MH = isochrone.MH
        else:
            if isochrone.MH != static_MH:
                print("[-] Error. All isochrones provided into the isochrone file must contain the same metallicity [M/H]")
                print(f"    Current isochrone has a metallicity of '{isochrones_list[index].MH:.3f}' while the accepted one is '{isochrones_list[index-1].MH:.3f}'")
                sys.exit(1)
        for line_containing_data in isochrone.lines_containing_data:
            if line_containing_data.startswith("#"):
                continue
            try:
                G_RP = float(line_containing_data.split()[-1])
                G_BP = float(line_containing_data.split()[-2])
                label = int(line_containing_data.split()[9])
                color = G_BP - G_RP
                if label == label_var: # '1' label means Main Sequence stars; useful to set ZAMS
                    data_for_future_plot.color.append(color)
                    data_for_future_plot.magnitude.append(G_BP)
                else:
                    continue
            except ValueError:
                print(f"[-] Error trying when trying to convert '{line_containing_data.split()[-1]}' , '{line_containing_data.split()[-2]}' to floats")
        list_of_data_to_plot.append(data_for_future_plot)
    list_ZAMS=[]
    for index, data in enumerate(list_of_data_to_plot):
        list_ZAMS.append(ZAMS(id=index+1, log_age=data.log_age, MH=data.MH, color=data.color[0], magnitude=data.magnitude[0]))
    if len(list_ZAMS) == 0:
        print("[-] No points could be retrieved to build ZAMS")
        sys.exit(1)
    return list_ZAMS[args.first_n_elements:-args.last_n_elements]
    

def plot_isochrones(args, data_plot_list: list[ZAMS])->None:
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
        plt.scatter(corrected_color, absolute_magnitude_data, s=4.)
    color = [c.color for c in data_plot_list]
    magnitude = [m.magnitude for m in data_plot_list]
    if args.compare_isochrone is not None:
        data_isochrone_compare = Table.read(args.compare_isochrone, format="ascii")[args.show_first_n_elem_isochrone:-args.show_last_n_elem_isochrone]
        plt.plot(data_isochrone_compare['G_BPmag']-data_isochrone_compare['G_RPmag'], data_isochrone_compare['G_BPmag'],
                 linestyle = '-', color="orangered")
    plt.plot(color, magnitude, 'r-')
    plt.gca().invert_yaxis()
    plt.show()
    plt.close()


def save_isochrones_to_file(args: argparse.Namespace, data_lines: list[str])->None:
    if not args.save_isochrones:
        print("[!] Data not saved")
        return
    # Select lines we are interested in. For ZAMS, it is the first line where "label=1" and the age has changed wrt previous line
    all_data = []
    data_to_write = []
    list_of_ages = []
    list_of_metallicities = []
    for index, line in enumerate(data_lines):
        if line.startswith("# Zini") and index == 0:
            all_data.append(line+'\n')
        elif line.startswith("# Z ini") and index != 0:
            continue
        elif line.startswith("#isochrone terminated"):
            continue
        else:
            all_data.append(line+'\n')

    previous_MH, previous_log_age = None, None
    for index, data in enumerate(all_data):
        if data.startswith("# Zini") and index == 0:
            data_to_write.append(data)
            continue
        if data.startswith("# Zini") and index != 0:
            continue
        if data.startswith("#isochrone terminated"):
            continue
        MH = float(data.split()[1])
        log_age = float(data.split()[2])
        label = int(data.split()[9])
        if MH != previous_MH and index > 1:
            print(f"[!] Error. All metallicities must be the same")
            print(f"    Metallicity in line number {index+1} is {MH:.2f}, while previous value is {previous_MH:.2f}")
            sys.exit(1)
        if log_age != previous_log_age and label == 1:
            data_to_write.append(data)
            previous_log_age = log_age
        previous_MH = MH
    # Save file
    file_path = Path(args.isochrones_filename)
    parent_directory = file_path.parent
    if args.object_name_to_save is not None:
        object_str = args.object_name_to_save.lower()
    else:
        print("[!] You have not provided a Object name to save files.")
        object_str = str(input("    Please enter your object name here: ")).lower()
    output_file_path = parent_directory / f"{object_str}_ZAMS.dat"
    with open(output_file_path, 'w') as f:
        for write_data in data_to_write:
            f.write(write_data)
    print(f"[+] Data saved as {str(output_file_path)!r}")
    return


def main():
    # Get arguments from user
    args = parse_flags()
    # Check arguments from user (e.g., if the user passes Gaia-based data, then distance and extinction must be provided)
    check_arguments_provided(args)
    # Read the PARSEC file containing the data to extract the Zero-Age Main Sequence. They must have the same metallicity
    data_lines = read_PARSEC_isochrone_file(args)
    # Get the different isochrones contained into the PARSEC file provided
    list_of_isochrones = get_different_log_age_and_MH_combinations(data_lines, data_lines[0])
    # Get ZAMS
    data_ZAMS = get_ZAMS_data_by_label(args, list_of_isochrones)
    # Plot ZAMS
    plot_isochrones(args, data_ZAMS)
    # Save ZAMS
    save_isochrones_to_file(args, data_lines)
    

if __name__ == "__main__":
    main()
