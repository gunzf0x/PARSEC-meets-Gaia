import argparse
from dataclasses import dataclass
import requests


@dataclass(kw_only=True)
class data_params:
    name: str
    distance: float
    e_distance: float # error in distance
    log_age: float
    e_log_age: float
    metallicity: float
    e_metallicity: float
    Av: float
    e_Av: float
    source: str


def parse_flags():
    """
    Parse flags from user
    """
    # Create the parser
    parser = argparse.ArgumentParser(description='Plot isochrones from PARSEC to Gaia DR3 data.')

    # Add flags/options
    parser.add_argument('-n', '--name', type=str, help='Cluster name to search for.', required=True)
    

    # Parse the command-line arguments
    args = parser.parse_args()
    return args


def get_globular_cluster(args) -> None | data_params:
    """
    Checks if the object provided is in Harris catalogue
    """
    harris_url: str = 'https://raw.githubusercontent.com/bersavosh/GC_cat/master/gc_cat.txt'
    response = requests.get(harris_url)
    if response.status_code != 200:
        print(f"[-] Something happened when trying to requests Harris' data. Status code: {response.status_code}")
        return
    data = response.text
    for number_line, line in enumerate(data.split('\n')):
        if len(line) == 0 or number_line == 0:
            continue
        cluster_name = line.split()[0]
        if args.name.lower().replace("_","").replace(" ","") == cluster_name.lower():
            distance = float(line.split()[6])
            metallicity = float(line.split()[11])
            Av = 3.1 * float(line.split()[12])
            source = "Harris (1996, 2010)"
            return data_params(name=cluster_name, distance=distance, e_distance= 0.0,
                               log_age = 0., e_log_age = 0., metallicity=metallicity, 
                               e_metallicity= 0.0, Av=Av, e_Av=0.0, source=source)
    print("[-] Object not found as a Globular Cluster")
    return


def get_open_custer(args) -> None | data_params:
    dias_url: str = 'https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/504/356/table2.dat'
    response = requests.get(dias_url)
    data = response.text
    if response.status_code != 200:
        print(f"[-] Something happened when trying to requests Dias' data. Status code: {response.status_code}")
        return
    for line in data.split('\n'):
        cluster_name = line.split()[0]
        if args.name.lower().replace("_","").replace(" ","") == cluster_name.lower().replace("_", ""):
            distance = float(line.split()[1])
            e_distance = float(line.split()[2])
            log_age = float(line.split()[3])
            e_log_age = float(line.split()[4])
            metallicity = float(line.split()[5])
            e_metallicity = float(line.split()[6])
            Av = float(line.split()[7])
            e_Av = float(line.split()[8][:4])
            source = "Dias et al. (2021)"
            return data_params(name=cluster_name, distance=distance, e_distance= e_distance,
                               log_age=log_age, e_log_age=e_log_age,
                               metallicity=metallicity, e_metallicity=e_metallicity,
                               Av=Av, e_Av=e_Av, source=source)
    print("[-] Object not found as a Open Cluster")
    return



def get_online_data(args) -> None | data_params:
    globular_cluster = get_globular_cluster(args)
    if globular_cluster is not None:
        return globular_cluster
    open_cluster = get_open_custer(args)
    if open_cluster is not None:
        return open_cluster
    return None


def print_object_details(data)->None:
    if data is not None:
        if "Harris" in data.source:
            print(f"[+] Name: {data.name}")
            print(f"    Distance: {data.distance*1000} pc")
            print(f"    Metallicity: {data.metallicity} dex")
            print(f"    Av: {data.Av:.3f} mag")
            print(f"    Source: {data.source}")
            return
        if "Dias" in data.source:
            print(f"[+] Name: {data.name.replace('_', '')}")
            print(f"    Distance: ({data.distance} ± {data.e_distance}) pc")
            print(f"    Log10 age: ({data.log_age}  ± {data.e_log_age}) dex")
            print(f"    Metallicity: ({data.metallicity} ± {data.e_metallicity}) dex")
            print(f"    Av: ({data.Av:.3f} ± {data.e_Av}) mag")
            print(f"    Source: {data.source}")
            return
    print(f"[-] Object {data.name!r} could not be identified as Open Cluster or Globular Cluster")
    return


def main()->None:
    args = parse_flags()
    data = get_online_data(args)
    print_object_details(data) 
    return


if __name__ == "__main__":
    main()
