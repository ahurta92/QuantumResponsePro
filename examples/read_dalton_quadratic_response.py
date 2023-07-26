from pathlib import Path
import re

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# Write a class that reads a Dalton quadratic response output file
# and stores the data in a pandas dataframe
# Below is a sample of a data file with the store values
# As you can see each line starts with B-freq and C-freq which indicate the two input frequencies
# follow by beta(Z;Z,Z) = val  the First parameter is the output direction, and the two second
# two are the directions of the two input frequencies.  My class should input a file name
# read the file and store the data

#  Results from quadratic response calculation
# --------------------------------------------
#
# @ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Z,Z) =     31.90375429
# @ B-freq = 0.011116  C-freq = 0.000000     beta(Z;Z,Z) =     31.93779191
# @ B-freq = 0.022233  C-freq = 0.000000     beta(Z;Z,Z) =     32.04026782
# @ B-freq = 0.033349  C-freq = 0.000000     beta(Z;Z,Z) =     32.21228015
# @ B-freq = 0.044465  C-freq = 0.000000     beta(Z;Z,Z) =     32.45568909
# @ B-freq = 0.055582  C-freq = 0.000000     beta(Z;Z,Z) =     32.77316340
# @ B-freq = 0.066698  C-freq = 0.000000     beta(Z;Z,Z) =     33.16824826
# @ B-freq = 0.077814  C-freq = 0.000000     beta(Z;Z,Z) =     33.64545758
# @ B-freq = 0.088930  C-freq = 0.000000     beta(Z;Z,Z) =     34.21039515
# @ B-freq = 0.000000  C-freq = 0.011116     beta(Z;Z,Z) = beta(Z,Z,Z)
# @ B-freq = 0.011116  C-freq = 0.011116     beta(Z;Z,Z) =     32.00604856
# @ B-freq = 0.022233  C-freq = 0.011116     beta(Z;Z,Z) =     32.14329168
# @ B-freq = 0.033349  C-freq = 0.011116     beta(Z;Z,Z) =     32.35099886
# @ B-freq = 0.044465  C-freq = 0.011116     beta(Z;Z,Z) =     32.63143213
# @ B-freq = 0.055582  C-freq = 0.011116     beta(Z;Z,Z) =     32.98769495
# @ B-freq = 0.066698  C-freq = 0.011116     beta(Z;Z,Z) =     33.42381206
# @ B-freq = 0.077814  C-freq = 0.011116     beta(Z;Z,Z) =     33.94483598
# @ B-freq = 0.088930  C-freq = 0.011116     beta(Z;Z,Z) =     34.55698503
# @ B-freq = 0.000000  C-freq = 0.022233     beta(Z;Z,Z) = beta(Z,Z,Z)
# @ B-freq = 0.011116  C-freq = 0.022233     beta(Z;Z,Z) = beta(Z,Z,Z)
# @ B-freq = 0.022233  C-freq = 0.022233     beta(Z;Z,Z) =     32.31622711
# @ B-freq = 0.033349  C-freq = 0.022233     beta(Z;Z,Z) =     32.56094956
# @ B-freq = 0.044465  C-freq = 0.022233     beta(Z;Z,Z) =     32.88014983
# @ B-freq = 0.055582  C-freq = 0.022233     beta(Z;Z,Z) =     33.27740279
# @ B-freq = 0.066698  C-freq = 0.022233     beta(Z;Z,Z) =     33.75726077
# @ B-freq = 0.077814  C-freq = 0.022233     beta(Z;Z,Z) =     34.32537631
# @ B-freq = 0.088930  C-freq = 0.022233     beta(Z;Z,Z) =     34.98865982
# @ B-freq = 0.000000  C-freq = 0.033349     beta(Z;Z,Z) = beta(Z,Z,Z)
# @ B-freq = 0.011116  C-freq = 0.033349     beta(Z;Z,Z) = beta(Z,Z,Z)
# @ B-freq = 0.022233  C-freq = 0.033349     beta(Z;Z,Z) = beta(Z,Z,Z)
# @ B-freq = 0.033349  C-freq = 0.033349     beta(Z;Z,Z) =     32.84443162
# @ B-freq = 0.044465  C-freq = 0.033349     beta(Z;Z,Z) =     33.20460035 @ B-freq = 0.055582  C-freq = 0.033349     beta(Z;Z,Z) =     33.64554879
# @ B-freq = 0.066698  C-freq = 0.033349     beta(Z;Z,Z) =     34.17241678
# @ B-freq = 0.077814  C-freq = 0.033349     beta(Z;Z,Z) =     34.79153248
# @ B-freq = 0.088930  C-freq = 0.033349     beta(Z;Z,Z) =     35.51059350
# @ B-freq = 0.000000  C-freq = 0.044465     beta(Z;Z,Z) = beta(Z,Z,Z)
# @ B-freq = 0.011116  C-freq = 0.044465     beta(Z;Z,Z) = beta(Z,Z,Z)
# @ B-freq = 0.022233  C-freq = 0.044465     beta(Z;Z,Z) = beta(Z,Z,Z)
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
    "/mnt/data/madness_data/fd_compare/dalton/hf/Be/dipole/hyper/hyp_CO-aug-cc-pCVDZ.out")
dalton_file_2 = Path("/mnt/data/madness_data/post_watoc/august/dalton/hf/H2O/dipole/hyper/hyp_H2O"
                     "-aug-cc-pVDZ.out")
dalQuad = ReadDaltonQuadraticResponse(dalton_file_2)

print(dalQuad.df)


def process_beta(row, df):
    # Check if the beta value is not a float
    if isinstance(row["Beta Value"], str):
        # Find the row with reversed B-freq and C-freq
        match_row = df[(df["B-freq"] == row["C-freq"]) & (df["C-freq"] == row["B-freq"])]
        # If a matching row is found, replace the beta value
        if not match_row.empty:
            row["Beta Value"] = match_row["Beta Value"].values[0]
    return row


df_original = dalQuad.df
df = df_original.apply(lambda row: process_beta(row, df_original), axis=1)

print(df)


def plot_beta(df):
    # Unique values of Letter1
    letter1_values = df["A"].unique()

    # Loop through each unique value of Letter1
    for letter1 in letter1_values:
        # Create a new figure for the current Letter1
        fig = plt.figure(figsize=(15, 15))

        # Filter the DataFrame for the current value of Letter1
        df_letter1 = df[df["A"] == letter1]

        # Unique combinations of Letter2 and Letter3 for the current Letter1
        directions = df_letter1[["B", "C"]].drop_duplicates()
        directions = ["X", "Y", "Z"]

        # Loop through each unique combination of Letter2 and Letter3
        j = 1
        for d1 in directions:
            for d2 in directions:
                ax = fig.add_subplot(3, 3, j, projection='3d')

                # Filter the DataFrame for the current combination of Letter2 and Letter3
                df_dir = df_letter1[
                    (df_letter1["B"] == d1) & (df_letter1["C"] == d2)]

                # Plot data on the current subplot
                ax.scatter(df_dir["B-freq"], df_dir["C-freq"], df_dir["Beta Value"])
                ax.set_title(f'Direction: {letter1},{d1},{d2}')
                ax.set_xlabel("B-freq")
                ax.set_ylabel("C-freq")
                ax.set_zlabel("Beta Value")
                j = j + 1

        plt.tight_layout()
        plt.show()


min_frequency = 0.0
max_frequency = 1

divisions = 8

plot_beta(df)
freq = np.linspace(min_frequency, max_frequency, divisions + 1)
dir = ["Z", "Y", "X"]

print(freq)
row = []
for d1 in dir:
    for d2 in dir:
        for d3 in dir:
            for i, f in enumerate(freq):
                for (j, f1) in enumerate(freq):
                    if f + f1 <= max_frequency:
                        row.append(pd.Series(
                            {"A-freq": -(f + f1), "B-freq": f1, "C-freq": f, "A": d1, "B": d2,
                             "C": d3}))

freq_def = pd.concat(row, axis=1).transpose()
print(freq_def)
