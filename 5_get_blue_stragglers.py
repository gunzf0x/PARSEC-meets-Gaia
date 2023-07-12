from astropy.table import Table
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
import copy
import matplotlib.pyplot as plt


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
    parser.add_argument('--subs-mag-MSTO', type=float, help='Quantity to substract to the identified MS-turnoff magnitude', required=True)
    parser.add_argument('--first-n-elements-isochrone', type=int, help='Select the first N elements of the isochrone')
    parser.add_argument('--last-n-elements-isochrone', type=int, help='Select the last N elements of the isochrone')
    parser.add_argument('--first-n-elements-ZAMS', type=int , help='Select the first N elements for the ZAMS isochrone')
    parser.add_argument('--last-n-elements-ZAMS', type=int, help='Select the last N elements for the ZAMS isochrone')
    parser.add_argument('--skip-confirmation', action='store_true', help="Do not ask if you want to extract BSS, just does it")
    # Parse the command-line arguments
    args = parser.parse_args()
    return args


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
    gaia_data.pprint()
    ZAMS_data = Table.read(args.zams_isochrone, format="ascii")
    isochrone_fit_data = Table.read(args.isochrone_filename, format="ascii")
    return gaia_data, ZAMS_data, isochrone_fit_data
    


def plot_before_extracting(args: argparse.Namespace, gaia_data: Table, ZAMS_isochrone: Table,
                           fit_isochrone: Table, data_file: isochrone_data)->None:
    plt.scatter(gaia_data['bp_rp'], gaia_data['phot_bp_mean_mag'])
    plt.gca().invert_yaxis()
    plt.show()
    plt.close()
    


def get_Blue_Stragglers(args: argparse.Namespace, data_file: isochrone_data, original_gaia_data: Table):
    # Create a copy since the data will be modified, so we can comparel ater with original data
    gaia_data = copy.deepcopy(original_gaia_data)
    
    # Condition 1: Every star must brighter than the MSTO (or slightly fainter)
    mask_lower_MSTO = gaia_data['phot_bp_mean_mag'] < (isochrone_data.MSTO_magnitude + args.add_mag_MSTO)
    gaia_data = gaia_data[mask_lower_MSTO]

    # Condition 2: Every star must cannot be brighter than ~2 mags above the turnoff. This might change bewteen objects, but
    # 2 magnitudes is the usual value
    mask_upper_MSTO = gaia_data['phot_bp_mean_mag'] > (isochrone_data.MSTO_magnitude - args.subs_mag_MSTO)
    gaia_data = gaia_data[mask_upper_MSTO]



def main()->None:
    # Get flags from user
    args = parse_flags()
    # Get data from file containing data
    data_inside_file = get_data_from_datafile(args)
    # Get the object data from Gaia DR3, ZAMS generated and isochrone fit
    gaia_data, zams_isochrone, fit_isochrone = get_data_from_files(args)
    print(data_inside_file)
    plot_before_extracting(args, gaia_data, zams_isochrone, fit_isochrone, data_inside_file)


if __name__ == "__main__":
    main()
