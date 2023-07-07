import requests
from pwn import log
import sys
import re
import argparse

PARSEC_webpage_url: str = "http://stev.oapd.inaf.it"

def parse_flags():
    """
    Parse flags from user
    """
    # Create the parser
    parser = argparse.ArgumentParser(description='Request data to PARSEC isochrone online service.')

    # Add flags/options
    parser.add_argument('--min-log-age', type=float, help='Min log10(age/yr) to request (dex)', required=True)
    parser.add_argument('--max-log-age', type = float, help='Max log10(age/yr) to request (dex)', required=True)
    parser.add_argument('--step-log-age', type=float, help = "Steps in log10(age/yr) (dex) (use '0' fora single value)", required=True)
    parser.add_argument('--min-log-m', type=float, help="Min log(M/H) (dex)", required=True)
    parser.add_argument('--max-log-m', type=float, help="Max log(M/H) (dex)", required=True)
    parser.add_argument('--step-log-m', type=float, help="Steps in log(M/H) (dex) (use '0' for a single value)")
    parser.add_argument('-o', '--outfile', type=str, help="Output filename. If not provided it is saved using the name provided by PARSEC service")
    parser.add_argument('--no-save', action="store_true", help="Do not save requested data")
    # Parse the command-line arguments
    args = parser.parse_args()
    return args


def check_args_provided(args)->None:
    """
    Check parameters provided by the user
    """
    # Check if the values for step provided is bigger than the difference between the min and max value
    difference_log_age = abs(args.max_log_age - args.min_log_age)
    if args.step_log_age > difference_log_age:
        print(f"[-] The step in age ({args.step_log_age!r}) cannot be bigger than the difference bewteen them ('{difference_log_age:.2f}')")
        sys.exit(1)
    difference_log_m =  abs(args.max_log_m - args.min_log_m)
    if args.step_log_m > difference_log_m:
        print(f"[-] The step in metallicity ({args.step_log_m!r}) cannot be bigger than the difference bewteen them ('{difference_log_m:.2f}')")
        sys.exit(1)
    # Check if the min value is biggerthan the max value. If so, exchange them
    if args.min_log_age > args.max_log_age:
        temp = args.min_log_age
        setattr(args, 'min_log_age', args.max_log_age)
        setattr(args, 'max_log_age', temp)
    if args.min_log_m > args.max_log_m:
        temp = args.min_log_m
        setattr(args, 'min_log_m', args.max_log_m)
        setattr(args, 'max_log_m', temp)
    return


def get_PARSEC_isochrone_from_webpage(args)-> str | None:
    """
    Make request to PARSEC webpage
    """
    p = log.progress("PARSEC webpage data")
    PARSEC_cgi_bin = f"{PARSEC_webpage_url}/cgi-bin/cmd"
    p.status(f"Requesting data to {PARSEC_cgi_bin!r}")
    # Payload containing 
    payload = {
        "cmd_version": "3.7",
        "track_omegai": "0.00",
        "track_parsec": "parsec_CAF09_v1.2S",
        "track_colibri": "parsec_CAF09_v1.2S_S_LMC_08_web",
        "track_postagb": "no",
        "n_inTPC": "10",
        "eta_reimers": "0.2",
        "kind_interp": "1",
        "kind_postagb": "-1",
        "photsys_file": "YBC_tab_mag_odfnew/tab_mag_gaiaEDR3.dat",
        "photsys_version": "odfnew",
        "dust_sourceM": "nodustM",
        "dust_sourceC": "nodustC",
        "kind_mag": "2",
        "kind_dust": "0",
        "extinction_av": "0.0",
        "extinction_coeff": "constant",
        "extinction_curve": "cardelli",
        "kind_LPV": "3",
        "imf_file": "tab_imf/imf_chabrier_lognormal_salpeter.dat",
        "isoc_agelow": "1.0e9",
        "isoc_ageupp": "1.0e10",
        "isoc_dage": "0.0",
        "isoc_isagelog": "1",
        "isoc_lagelow": str(args.min_log_age), # "6.6"
        "isoc_lageupp": str(args.max_log_age), # "10.13"
        "isoc_dlage": str(args.step_log_age), # "1.5"
        "isoc_zlow": "0.0152",
        "isoc_zupp": "0.03",
        "isoc_dz": "0.0",
        "isoc_ismetlog": "1",
        "isoc_metlow": str(args.min_log_m), #"-2"
        "isoc_metupp": str(args.max_log_m), # "0.3"
        "isoc_dmet": str(args.step_log_m),
        "output_kind": "0",
        "output_evstage": "1",
        "lf_maginf": "-15",
        "lf_magsup": "20",
        "lf_deltamag": "0.5",
        "sim_mtot": "1.0e4",
        "submit_form": "Submit",
        ".cgifields": [
            "kind_LPV",
            "track_parsec",
            "output_kind",
            "output_gzip",
            "extinction_curve",
            "dust_sourceM",
            "track_omegai",
            "track_colibri",
            "track_postagb",
            "photsys_version",
            "isoc_isagelog",
            "dust_sourceC",
            "isoc_ismetlog",
            "extinction_coeff"
        ]
    }
    # Send the POST request to PARSEC webpage
    response = requests.post(PARSEC_cgi_bin, data=payload)
    # Check the response
    if response.status_code == 200:
        p.success("Data succesfully requested")
    else:
        p.failure(f"Could not retrieve data from PARSEC webpage (status code {response.status_code})")
        sys.exit(1)
    # Return the response
    return response.text


def get_filename_name_from_response(response_text: str) -> str | None:
    """
    Get the filename containing the output data generated
    """
    # Define the regular expression pattern
    pattern = r".*output\d+.*"
    # Find all lines matching the pattern
    matching_lines = re.findall(pattern, response_text, re.MULTILINE)
    if len(matching_lines) == 1:
        line_with_outfile =  matching_lines[0]
    elif len(matching_lines) > 1:
        print("[!] More than 1 output file found. Returning the first one found")
        line_with_outfile =  matching_lines[0]
    else:
        print("[-] No output file found")
        sys.exit(1)
    pattern_file = r"output\d+\.dat"
    # Find the first occurrence of the pattern
    match = re.search(pattern_file, line_with_outfile)
    # Extract the matched string
    if match:
        matched_file = match.group()
        return matched_file
    else:
        print("[-] Pattern not found for filename.")
        sys.exit(1)


def create_valid_filename_to_save(filename: str)->str:
    """
    Check if the outfile name ends with ".dat"
    """
    if not filename.endswith(".dat"):
        filename = f"{filename}.dat"
        return filename
    else:
        return filename


def save_file(args, original_filename: str)->None:
    """
    Save data requested to PARSEC service
    """
    if args.no_save:
        log.info("'--no-save' flag enabled. Data will not be saved")
        return
    p = log.progress("Saving data")
    if args.outfile is not None:
        filename: str = create_valid_filename_to_save(args.outfile)
    else:
        filename: str = original_filename
    data_url = f"{PARSEC_webpage_url}/tmp/{original_filename}"
    # Make a request to the data generated
    response = requests.get(data_url)
    # Check the response
    if response.status_code == 200:
        p.status("Data succesfully requested")
    else:
        p.failure(f"Could not retrieve data from PARSEC webpage (status code {response.status_code})")
        sys.exit(1)
    # Write the content into a file
    with open(filename, "wb") as f:
        f.write(response.content)
    p.success(f"Data saved into file {filename!r}")


# MAIN
def main()->None:
    # Get user parameters from flags provided
    args = parse_flags()
    # Check if the parameters provided are correct
    check_args_provided(args)
    # Get the response text for the PARSEC website request
    response_text = get_PARSEC_isochrone_from_webpage(args)
    # Get the filename containing data
    file_containing_data = get_filename_name_from_response(response_text)
    save_file(args, file_containing_data)


if __name__ == "__main__":
    main()
