import os
import shutil
from chemical_network_builder import ChemicalNetworkBuilder

def prepare_basic_param_file(template_path, output_path, ramses_repo, snapnum, elements_dict, sed_dir):
    """
    Copies and customizes a simulation parameter file.

    Args:
        template_path (str): Path to the base parameter file
        output_path (str): Output folder for the current run
        ramses_repo (str): Path to the Ramses-RT outputs
        snapnum (int): Snapshot number
        elements_dict (dict): Element name to max ion stage (e.g., {"C": 4, "Si": 4})
    """

    # Start with krome parameter file:
    write_param_krome_file(output_path, elements_dict)
    
    # Read the template
    with open(template_path, 'r') as f:
        lines = f.readlines()

    str_elements = "_".join(f"{el}_{n}" for el, n in elements_dict.items())

    # Apply simple replacements
    new_lines = []
    for line in lines:
        if "repository =" in line:
            new_lines.append(f"  repository = {ramses_repo} \n")
        elif "snapnum" in line:
            new_lines.append(f"  snapnum = {snapnum} \n")
        elif "output_path" in line:
            new_lines.append(f"  output_path = {output_path} \n")
        elif "csn_file" in line:
            new_lines.append(f"  csn_file = csn_{str_elements}_{snapnum} \n")
        elif "sed_dir" in line:
            new_lines.append(f"  sed_dir = {sed_dir} \n")
        elif "n_elements" in line:
            new_lines.append(f"  n_elements = {len(elements_dict)} \n")
        elif "krome_parameter_file" in line:
            new_lines.append(f"  krome_parameter_file = {os.path.join(output_path, "params_krome.dat")} \n")
        else:
            new_lines.append(line)

    # Write to destination
    dest_file = os.path.join(output_path, "params_global.dat")
    with open(dest_file, 'w') as f:
        f.writelines(new_lines)

    print(f"Parameter file written to {dest_file}")


def write_param_krome_file(output_path, elements_dict):
    """
    Prepares the parameter file for Krome, with element names and max ionization. Assumes solar abundances, for now

    Args:
        output_path (str): Output folder for the current run
        elements_dict (dict): Element name to max ion stage (e.g., {"C": 4, "Si": 4})
    """

    builder = ChemicalNetworkBuilder()

    with open(os.path.join(output_path, "params_krome.dat"), 'w') as f:
        for element, max_ion in elements_dict.items():
            Z = builder.get_atomic_number(element)
            f.write(f"{Z} {max_ion} -1.00\n")
    
