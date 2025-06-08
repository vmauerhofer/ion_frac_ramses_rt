import pandas as pd
import numpy as np
import os

class ChemicalNetworkBuilder:
    
    element_names = [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe'
    ]

    # === init ===
    def __init__(self, path=None):
        
        if path is None:
            try:
                base_dir = os.path.dirname(os.path.abspath(__file__))
            except NameError:
                raise RuntimeError("Cannot auto-detect data path: __file__ is not defined. Please provide the path explicitly.")
            self.path = os.path.join(base_dir, "reaction_rates")
        else:
            self.path = path

        # Load the Badnell Radiative Recombination (RR) data into a pandas DataFrame
        self.rr_df = pd.read_csv(
            f"{self.path}/Radiative_rec/clist",
            sep=r'\s+',
            skiprows=3,
            names=["Z", "N", "M", "W", "A", "B", "T0", "T1", "C", "T2"]
        )

        # Load the Badnell Dielectronic Recombination (DR) data into a pandas DataFrame
        # Column names
        base_cols = ["Z", "N", "M", "W"]
        C_cols = [f"C{i}" for i in range(1, 10)]
        E_cols = [f"E{i}" for i in range(1, 10)]

        # Read the two files
        df_c = pd.read_csv(
            f"{self.path}/Dielectronic_rec/clist_c",
            sep=r'\s+',
            skiprows=3,
            names=base_cols + C_cols
        )
        df_e = pd.read_csv(
            f"{self.path}/Dielectronic_rec/clist_E",
            sep=r'\s+',
            skiprows=3,
            names=base_cols + E_cols
        )

        # Merge the dataframes on Z, N, M, W
        self.dr_df = pd.merge(df_c, df_e, on=base_cols)

        # Load the Voronov collisional ionization rates (CI) data into a pandas DataFrame
        self.ci_df = pd.read_csv(
        os.path.join(self.path, "Collisions/Voronov/clist"),
        sep=r'\s+',
        names=["dE", "P", "A", "X", "K"]
        )

# === Public Methods ===
#############################################

# === Return name of an atom with atomic number Z ===
    def get_element_name(self, Z):
        return self.element_names[Z-1]

# === Return atomic number of a given element ===
    def get_atomic_number(self, element):
        element = element.capitalize()
        try:
            return self.element_names.index(element) + 1
        except ValueError:
            raise ValueError(f"Element '{element}' not recognized. Allowed: {self.element_names}")
    
# === Build a chemical network file, to be read by KROME ===
    def build_krome_network(self, elements: dict[str, int], output_path: str) -> None:
        """
        Build a KROME chemical network file including recombination, collisional ionization,
        and photoionization for a set of elements and maximum ionization stages.
    
        Args:
            elements (dict): Dictionary mapping element symbols (e.g. 'C') to max ionization state (int).
            output_path (str): Path to write the resulting KROME network file.
        """
        network = ''
        count = 1
    
        # Header
        network += '@var: T = Tgas\n\n'
        # Recombinations
        network += '# RECOMBINATION from Badnell website\n#http://amdpp.phys.strath.ac.uk/tamoc/RR\n'
        network += '#http://amdpp.phys.strath.ac.uk/tamoc/DR\n\n@format:idx,R,R,P,rate\n'

        lines, count = self._add_dummy_recombinations(count)
        network += '\n'.join(lines) + '\n\n'

        ## Looping on the elements
        for symbol, max_ion in elements.items():
            Z = self.element_names.index(symbol) + 1  # Atomic number
            # Looping on the ionization stages
            for ion_stage in range(max_ion):
                left = symbol + '+' * (ion_stage + 1)
                right = symbol + '+' * ion_stage
                try:
                    #Z - ion_stage - 1 is the number of bound electrons before recombination
                    rr_rate = self._get_RR_rate(Z, Z - ion_stage - 1)
                    dr_rate = self._get_DR_rate(Z, Z - ion_stage - 1)
                    rate = f"{rr_rate} + {dr_rate}"
                except Exception:
                    rate = "auto"
                network += f"{count},{left},E,{right}, {rate}\n"
                count += 1
            network += '\n'
    
        # Collisional Ionization
        network += '# COLLISIONS from Voronov et al. 1997\n@format:idx,R,R,P,P,P,rate\n'
        lines, count = self._add_dummy_collisions(count)
        network += '\n'.join(lines) + '\n\n'
        for symbol, max_ion in elements.items():
            Z = self.element_names.index(symbol) + 1
            for ion_stage in range(max_ion):
                left = symbol + '+' * ion_stage
                right = symbol + '+' * (ion_stage + 1)
                try:
                    ci_rate = self._get_CI_rate(Z, Z - ion_stage)
                except Exception:
                    ci_rate = "auto"
                network += f"{count},{left},E,{right},E,E, {ci_rate}\n"
                count += 1
            network += '\n'
    
        # Photoionization, using rates implemented in KROME ("auto")
        network += '# PHOTOIONIZATION\n@photo_start\n@format:idx,R,P,P,rate\n'
        for symbol, max_ion in elements.items():
            for ion_stage in range(max_ion):
                left = symbol + '+' * ion_stage
                right = symbol + '+' * (ion_stage + 1)
                network += f"{count},{left},{right},E, auto\n"
                count += 1
        network += '@photo_stop\n'
    
        # Write to file
        with open(output_path, 'w') as f:
            f.write(network)
        print(f"KROME network successfully written to '{os.path.abspath(output_path)}'")



# === Internal Helper Methods ===
#############################################


    def _add_dummy_recombinations(self, count):
        """Add H and He recombination reactions with 0 rate for KROME."""
        lines = []
        for Z in [1, 2]:  # H and He
            element = self.element_names[Z - 1]
            for ion_stage in range(Z):
                ion = element + '+' * (ion_stage + 1)
                product = element + '+' * ion_stage
                lines.append(f"{count},{ion},E,{product}, 0d0")
                count += 1
        return lines, count

    def _add_dummy_collisions(self, count):
        """Add dummy electron collision ionization reactions for H and He with 0 rate."""
        lines = []
        for Z in [1, 2]:  # H and He
            element = self.element_names[Z - 1]
            for ion_stage in range(Z):
                reactant = element + '+' * ion_stage
                product = element + '+' * (ion_stage + 1)
                lines.append(f"{count},{reactant},E,{product},E,E, 0d0")
                count += 1
        return lines, count


# === Gives the correct raw in the data file for recombinations ===
    def _get_row_recombination(self, df, Z, N, M):
        """
        Internal helper to select the best-fit row for given Z, N, M from a DataFrame.
        """
        if N < 0 or N > Z-1:
            raise ValueError(f"Impossible number of bound electrons (before recombination) N={N} for element with Z={Z}")
    
        subset = df[(df.Z == Z) & (df.N == N)]
    
        if subset.empty:
            raise ValueError(f"No data found for Z={Z}, N={N}.")
    
        M_max = subset.M.max()
        if M > M_max:
            M = M_max
    
        return subset[subset.M == M].iloc[0]


# === Radiative Recombination ===
    def _get_RR_rate(self, Z, N, M=1):
        """
        Returns the radiative recombination rate expression for a given ion, for a KROME chemical network.
        From http://amdpp.phys.strath.ac.uk/tamoc/RR
    
        Args:
            Z (int): Atomic number of the element (e.g., 14 for silicon).
            N (int): Number of bound electrons before recombination (e.g. N=4 and Z=6 means that C++ recombines to C+)
            M (int, optional): Fit model index if multiple fits exist. Defaults to 1.
    
        Returns:
            rr_rate_str: A Fortran-style string expression for the recombination rate as a function of temperature T.
        """

        row = self._get_row_recombination(self.rr_df, Z, N, M)
    
        ## Now building the string, in fortran-format, that will appear in the Krome chemical network.
        rr_rate_str = f"{row.A:.4e}*( sqrt(T/{row.T0:.4e}) * (1+sqrt(T/{row.T0:.4e}))**("
    
        exp_minus = f"{row.B:.4f}"
        if row.C != 0.0:
            exp_minus += f"+{row.C:.4f}*exp(-{row.T2:.4e}/T)"
        
        rr_rate_str += f"1-({exp_minus}))"
        rr_rate_str += f" * (1+sqrt(T/{row.T1:.4e}))**("
        
        exp_plus = f"{row.B:.4f}"
        if row.C != 0.0:
            exp_plus += f"+{row.C:.4f}*exp(-{row.T2:.4e}/T)"
        
        rr_rate_str += f"1+({exp_plus})) )**(-1)"

         # Convert to Fortran 'd' notation
        rr_rate_str = rr_rate_str.replace("e-", "d-").replace("e+", "d+")
        return rr_rate_str


# === Dielectronic Recombination ===
    def _get_DR_rate(self, Z, N, M=1):
        """
        Returns the dielectronic recombination rate expression for a given ion, for a KROME chemical network.
        From http://amdpp.phys.strath.ac.uk/tamoc/DR
    
        Args:
            Z (int): Atomic number of the element (e.g., 14 for silicon).
            N (int): Number of bound electrons before recombination (e.g. N=4 and Z=6 means that C++ recombines to C+)
            M (int, optional): Fit model index if multiple fits exist. Defaults to 1.
    
        Returns:
            str: A Fortran-style string expression for the recombination rate as a function of temperature T.
        """

        row = self._get_row_recombination(self.dr_df, Z, N, M)
    
         # Start building the expression
        dr_rate_str = "T**(-1.5)*("
    
        # Loop over C1–C9 and E1–E9
        terms = []
        for i in range(1, 10):
            C = row[f"C{i}"]
            E = row[f"E{i}"]
            if C == 0.0:
                break
            terms.append(f"{C:.4e}*exp(-{E:.4e}/T)")
    
        dr_rate_str += " + ".join(terms) + ")"
    
        # Convert to Fortran 'd' notation
        dr_rate_str = dr_rate_str.replace("e-", "d-").replace("e+", "d+")
    
        return dr_rate_str


# === Collisional ionization ===
    def _get_CI_rate(self, Z, N):
        """
        Returns the collisional ionization rate expression for a given ion (from Voronov et al. 1997).
    
        Args:
            Z (int): Atomic number of the element (e.g., 14 for silicon).
            N (int): Number of bound electrons before the collision (e.g. N=11 means Si⁺³ getting ionized to Si⁺⁴)
    
        Returns:
            str: Fortran-style expression for the collisional ionization rate as a function of temperature T,
                 or 'auto' if KROME should use its internal rates.
        """
    
        if N <= 0 or N > Z:
            raise ValueError(f"In get_CI_rate, invalid N={N} for Z={Z} ({self.element_names[Z-1]})")
    
        if Z < 3:
            return "auto"
    
        # Compute index based on Voronov table layout
        idx = (Z - 1) * Z // 2 + Z - N - 3
    
        row = self.ci_df.iloc[idx]
    
        # Temperature in eV (KROME uses T in K, so conversion factor 1.16045e4 is applied)
        col_rate = (
            f"{row.A:.4e}*(1 + {row.P:.4f}*sqrt({row.dE:.4e}*1.16045d4/T)) / "
            f"({row.X:.4f}+{row.dE:.4e}*1.16045d4/T) * "
            f"({row.dE:.4e}*1.16045d4/T)**{row.K:.4f} * exp(-{row.dE:.4e}*1.16045d4/T)"
        )
    
        return col_rate.replace("e-", "d-").replace("e+", "d+")
    

