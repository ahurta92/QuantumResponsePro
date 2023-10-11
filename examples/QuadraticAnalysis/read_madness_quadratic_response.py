import json
import pandas as pd
import re
from pathlib import Path


class ReadDaltonQuadraticResponse:
    # input quadratic response file
    # from file read in data lines
    # save into pandas dataframe
    def __init__(self, outfile: Path):
        self.outfile = outfile

        rows = []

        with open(self.outfile, "r") as file:
            for line in file:
                pattern = "  Results from quadratic response calculation"
                match = re.search(pattern, line)
                if match:
                    line = file.readline()
                    line = file.readline()
                    for line in file:
                        if line.startswith("@ B-freq"):
                            # Create a regex pattern that captures three letters (X, Y, or Z) within parentheses,
                            # separated by semicolons, and the floating point numbers in the string
                            pattern = r"B-freq = ([\d\.]+)\s+C-freq = ([\d\.]+)\s+beta\((X|Y|Z);(X|Y|Z),(X|Y|Z)\) = (.*)"

                            # Use the search function from the re module to get the match object
                            match = re.search(pattern, line)

                            if match:
                                b_freq, c_freq, letter1, letter2, letter3, beta_value = match.groups()
                                # Check if beta_value is another beta function call
                                if "beta" in beta_value:
                                    # If it is, replace it with the original letters
                                    beta_value = f"beta({letter1};{letter2},{letter3})"
                                a_freq = -(float(b_freq) + float(c_freq))
                                # Append the data to the DataFrame
                                rows.append(pd.Series({
                                    "A-freq": a_freq,
                                    "B-freq": float(b_freq),
                                    "C-freq": float(c_freq),
                                    "A": letter1,
                                    "B": letter2,
                                    "C": letter3,
                                    "Beta Value": beta_value if "," in beta_value else float(
                                        beta_value)
                                }))
                        else:
                            break
            self.df = pd.concat(rows, axis=1).transpose()


dalton_file = Path(
    "/mnt/data/madness_data/fd_compare/dalton/hf/Be/dipole/hyper/hyp_CO-aug-cc-pZ.out")
dalton_file_2 = Path("/mnt/data/madness_data/development/dalton/hf/H2O/dipole/quad_H2O-d-aug-cc"
                     "-pVTZ.out")
dalQuad = ReadDaltonQuadraticResponse(dalton_file_2)

df_original = dalQuad.df
# df = df_original.apply(lambda row: process_beta(row, df_original), axis=1)
# remove row if beta value is a string
df = df_original[df_original["Beta Value"].apply(lambda x: isinstance(x, float))]

df.rename(columns={"A-freq": "Afreq", "B-freq": "Bfreq", "C-freq": "Cfreq"}, inplace=True)
print(df.query("Afreq == 0.0 & Bfreq==0.0"))

madness_file = Path("/mnt/data/madness_data/development/hf/H2O/beta.json")

loaded_json = json.loads(madness_file.read_text())
df = pd.DataFrame(loaded_json)
# remove beta < 1e-3
df = df[df["Beta"] > 1e-3]
df.rename(columns={"A-freq": "Afreq", "B-freq": "Bfreq", "C-freq": "Cfreq"}, inplace=True)
print(df.query("Afreq == 0.0 & Bfreq==0.0"))
