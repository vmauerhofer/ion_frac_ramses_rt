import os
import shutil
import subprocess
from chemical_network_builder import ChemicalNetworkBuilder

def generate_krome_module(elements_dict, krome_path, output_folder_name=None):
    """
    Automates the generation of a KROME Fortran module using custom chemical network builder.
    
    Args:
        elements_dict (dict): Element name to max ion stage (e.g., {"C": 4, "Si": 4})
        krome_path (str): Absolute or relative path to the KROME directory
        output_folder_name (str): Name of the subfolder to store the resulting f90 files and the chemical network
    """

    # By default, we store results in krome/chem_net_ramses_rt
    if output_folder_name is None:
        output_folder_name = os.path.join(krome_path, "chem_net_ramses_rt")
    # Create a new folder for every unique elements_dict
    str_elements = "_".join(f"{el}_{n}" for el, n in elements_dict.items())
    output_folder_name = os.path.join(output_folder_name, str_elements)
    os.makedirs(output_folder_name, exist_ok=True)

    # Moving the modified krome code to krome/src
    for file in os.listdir('./krome_files_to_change'):
        shutil.copy(os.path.join('./krome_files_to_change', file), os.path.join(krome_path,'src'))

    # Generate chemical network
    builder = ChemicalNetworkBuilder()
    network_output_path = os.path.join(output_folder_name, f"network_{str_elements}")
    builder.build_krome_network(elements_dict, network_output_path)

    # Clean old KROME build
    build_path = os.path.join(krome_path, "build")
    for file in os.listdir(build_path):
        os.remove(os.path.join(build_path, file))
    print("Cleared KROME build directory")

    # Run KROME’s Python wrapper
    krome_script = os.path.join(krome_path, "krome")
    # Have to change working directory to run krome, then we come back
    original_dir = os.getcwd()
    os.chdir(krome_path)
    subprocess.run([
        "python3", krome_script, "-n", network_output_path,
        "-useN", "-compact", "-photoBins=1"], check=True, input=b"\n")  # Simulate "press any key" 
    os.chdir(original_dir)
    print("KROME build complete")

    # Move results to output_folder_name
    for file in os.listdir(build_path):
        shutil.copy(os.path.join(build_path, file), output_folder_name)
    print(f"Output files moved to {output_folder_name}")

    # Compile the combined krome fortran files and the local src_f90 files
    compile_fortran_code(makefile_dir="src_f90", krome_f90_dir=output_folder_name)


def compile_fortran_code(makefile_dir=".", krome_f90_dir=None):
    env = os.environ.copy()
    if krome_f90_dir:
        env["KROME_F90_DIR"] = krome_f90_dir

    try:
        subprocess.run(["make", "-C", makefile_dir], check=True, env=env)
        print("Fortran code compiled successfully.")
    except subprocess.CalledProcessError as e:
        print("Compilation failed. Please check the Makefile and paths.")
        raise e
