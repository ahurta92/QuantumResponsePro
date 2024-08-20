import os
import pandas as pd
from pathlib import Path


# make this into a class that way I can define class members like
# the periodic table and frequency json
# PROOT = os.getcwd()
# DALROOT = os.path.join(PROOT, os.pardir)
# DALROOT += "/dalton/"
# if not os.path.exists("dalton"):
#    os.mkdir("dalton")

# pt = pd.read_csv(PROOT+'/periodic_table.csv')
# with open(PROOT+'/molecules/frequency.json') as json_file:
#    freq_json = json.loads(json_file.read())


class madnessToDalton:
    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self.DALROOT = self.base_dir.joinpath("dalton").absolute()
        # print('madnessToDalton base dir', self.base_dir)
        # print('madnessToDalton dalroot dir', self.DALROOT)
        if not os.path.exists("dalton"):
            os.mkdir("dalton")

        self.periodic_table = pd.read_csv(self.base_dir.joinpath("periodic_table.csv"))

        # with open(self.base_dir.joinpath('json_data/frequency.json')) as json_file:
        #    self.freq_json = json.loads(json_file.read())

    def madmol_to_dalmol(self, madmol_f, basis, ncharge=0):
        print(madmol_f)

        dalton_inp = []
        dalton_inp.append("BASIS")
        dalton_inp.append(basis + "")
        param_keys={}
        with open(madmol_f) as file:  # Use file to refer to the file object
            lines = file.readlines()
            atom_dict = {}
            for line in lines:
                split = line.strip().split(" ")
                while "" in split:
                    split.remove("")
                if split[0] == "units":
                    units = split[1]
                skeys = [
                    "geometry",
                    "eprec",
                    "units",
                    "end",
                    "field",
                    "no_orient",
                    "psp_calc",
                    "pure_ae",
                    "symtol",
                    "core_type",
                ]
                skey = split[0]

                if skey in skeys:
                    # and not geometry or end
                    if skey == "end":
                        break
                    if skey == "geometry":
                        continue
                    param_keys[skey] = split[1]
                

                # if not a key word then it's probably an atom
                if not skey in skeys:
                    geometry = " ".join(split[1:])
                    if skey not in atom_dict:
                        atom_dict[skey] = []
                        atom_dict[skey].append(geometry)
                    else:
                        atom_dict[skey].append(geometry)

        
        # we validate that the keys are atoms using the dataframe defined from periodic_table.csv
        # Error if we don't find the key in the table
        atom_types = len(atom_dict)
        dalton_inp.append("blah")
        dalton_inp.append("blah")


        general_line = (
            "Atomtype="
            + str(atom_types)
            + " "
            + units.capitalize()
            + " "
            + "Charge="
            + str(ncharge)
        )
        # if param['norient'] == 'true':
        if "no_orient" in param_keys and param_keys["no_orient"] == "true":
            general_line +=" Nosymmetry"
        dalton_inp.append(general_line)

        for atom in atom_dict.keys():
            # Get the charge of the atom
            charge = self.periodic_table[self.periodic_table["Symbol"] == atom][
                "NumberofProtons"
            ].reset_index()

            charge = charge["NumberofProtons"][0]
            # number of atoms of type atom
            num_atoms = len(atom_dict[atom])
            charge_atom_line = "Charge=" + str(charge) + " Atoms=" + str(num_atoms)
            dalton_inp.append(charge_atom_line)
            lh = "a"
            for s in atom_dict[atom]:
                g_line = atom + "_" + lh + " " + s
                lh = chr(ord(lh) + 1)
                dalton_inp.append(g_line)

        dalton_inp = "\n".join(dalton_inp)
        return dalton_inp

    # Problem
    #
    # Some basis sets have characters(**) which make it difficult
    # to standardize the naming of the molfile.
    #
    # Solution?
    #
    # Create a dictionary of molnames which map to standard name that
    # do not create issues
    #
    # Or just replace * with S
    # 6-31G** becomes 6-31GSS
    #
    # assumption xc is hf or valid dft functional

    def dalton_polar_input(self, madmol, xc, operator, basis):
        molname = madmol.split(".")[0]
        dalton_inp = []
        dalton_inp.append("**DALTON INPUT")
        dalton_inp.append(".RUN RESPONSE")
        dalton_inp.append("**WAVE FUNCTIONS")
        if xc == "hf":
            dalton_inp.append(".HF")
        else:
            dalton_inp.append(".DFT")
            dalton_inp.append(xc.capitalize())
        dalton_inp.append("**RESPONSE")
        dalton_inp.append("*LINEAR")
        # looks like it's the only option for a response calculation
        if operator == "dipole":
            dalton_inp.append(".DIPLEN")
            freq = self.freq_json[molname][xc][operator]
            num_freq = len(freq)
            dalton_inp.append(".FREQUENCIES")
            dalton_inp.append(str(num_freq))

            freq_s = []
            for f in freq:
                freq_s.append(str(f))
            dalton_inp.append(" ".join(freq_s))

        dalton_inp.append("**END OF DALTON INPUT")
        dalton_inp = "\n".join(dalton_inp)
        run_dir = self.base_dir + "/dalton/" + xc + "/" + molname + "/" + operator
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        madmolfile = self.base_dir + "/molecules/" + madmol + ".mol"
        mol_input = self.madmol_to_dalmol(madmolfile, basis)
        dal_run_file = run_dir + "/freq.dal"
        with open(dal_run_file, "w") as file:  # Use file to refer to the file object
            file.write(dalton_inp)
        mol_file = run_dir + "/" + molname + "-" + basis.replace("*", "S") + ".mol"
        with open(mol_file, "w") as file:  # Use file to refer to the file object
            file.write(mol_input)
        return run_dir, dal_run_file.split("/")[-1], mol_file.split("/")[-1]
