import subprocess

import numpy as np
import shutil
import os
import pandas as pd
import json

from .daltonToJson import daltonToJson
from ..madness_to_dalton import madnessToDalton
from pathlib import Path


class DaltonRunner:

    @classmethod
    def __init__(self, base_dir: Path, run_new):
        self.run = run_new
        self.base_dir = base_dir  # what is my base directory?
        self.dalton_dir = self.base_dir.joinpath('dalton')

        print('DaltonRunner base_dir: ', self.base_dir)
        print('DaltonRunner dal_dir: ', self.dalton_dir)

        if shutil.which("mpirun") is not None:
            self.use_mpi = True
            self.Np = int(os.cpu_count() / 8)
        else:
            self.use_mpi = False
            self.Np = 1

        # where ever I run I can assume that the dalton directory will be one above cwd
        if not os.path.exists("dalton"):
            os.mkdir("dalton")
        try:
            with open(self.base_dir.joinpath("json_data/frequency.json")) as json_file:
                self.freq_json = json.loads(json_file.read())
        except FileNotFoundError as f_error:
            print("No frequency.json found. Run generate_data.py to generate it.")
            print("only excited state calculations will be available")
            pass

        # with open(DALROOT + '/dalton-dipole.json') as json_file:
        #    self.dipole_json = json.loads(json_file.read())
        # with open(DALROOT + '/dalton-excited.json') as json_file:
        #    self.excited_json = json.loads(json_file.read())

    @staticmethod
    def __write_polar_input(self, madness_molecule, xc, operator, basis):
        """writes the polar input to folder"""
        # DALTON INPUT
        molecule_input = madness_molecule.split(".")[0]
        dalton_inp = []
        dalton_inp.append("**DALTON INPUT")
        dalton_inp.append(".RUN RESPONSE")
        dalton_inp.append(".DIRECT")
        if basis.split("-")[-1] == "uc":
            dalton_inp.append("*MOLBAS ")
            dalton_inp.append(".UNCONT ")
        dalton_inp.append("**WAVE FUNCTIONS")
        # HF or DFT
        if xc == "hf":
            dalton_inp.append(".HF")
        else:
            dalton_inp.append(".DFT")
            dalton_inp.append(xc.capitalize())
        # RESPONSE
        dalton_inp.append("**RESPONSE")
        # LINEAR
        dalton_inp.append("*LINEAR")
        # looks like it's the only option for a response calculation
        if operator == "dipole":
            dalton_inp.append(".DIPLEN")
            freq = self.freq_json[molecule_input][xc][operator]
            num_freq = len(freq)
            dalton_inp.append(".FREQUENCIES")
            dalton_inp.append(str(num_freq))

            freq_s = []
            for f in freq:
                freq_s.append(str(f))
            dalton_inp.append(" ".join(freq_s))

        dalton_inp.append("**END OF DALTON INPUT")
        dalton_inp = "\n".join(dalton_inp)
        # print(dalton_inp)
        # print(molecule_input)
        # print(operator)
        # print("self.dalton", self.dalton_dir)

        run_dir = self.dalton_dir.joinpath(xc).joinpath(molecule_input).joinpath(operator)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        # Here I read the madness mol file from the molecules directory
        madness_molecule_file = self.base_dir.joinpath('molecules').joinpath(
            madness_molecule + '.mol')
        mad_to_dal = madnessToDalton(self.base_dir)
        if basis.split("-")[-1] == "uc":
            mol_input = mad_to_dal.madmol_to_dalmol(madness_molecule_file,
                                                    "-".join(basis.split("-")[:-1]))
        else:
            mol_input = mad_to_dal.madmol_to_dalmol(madness_molecule_file, basis)
        dal_run_file = run_dir.joinpath('freq.dal')
        with open(dal_run_file, "w") as file:  # Use file to refer to the file object
            file.write(dalton_inp)
        dalton_molecule_file = run_dir.joinpath(
            molecule_input + '-' + basis.replace("*", "S") + ".mol")
        with open(dalton_molecule_file, "w") as file:  # Use file to refer to the file object
            file.write(mol_input)
        molecule_input = dalton_molecule_file.stem
        dalton_input = dal_run_file.stem
        return run_dir, dalton_input, molecule_input

    @staticmethod
    def __write_quadratic_input(self, madness_molecule, xc, operator, basis):
        """writes the polar input to folder"""
        # DALTON INPUT
        molecule_input = madness_molecule.split(".")[0]
        dalton_inp = []
        dalton_inp.append("**DALTON INPUT")
        dalton_inp.append(".RUN RESPONSE")
        dalton_inp.append(".DIRECT")
        if basis.split("-")[-1] == "uc":
            dalton_inp.append("*MOLBAS ")
            dalton_inp.append(".UNCONT ")
        dalton_inp.append("**WAVE FUNCTIONS")
        # HF or DFT
        if xc == "hf":
            dalton_inp.append(".HF")
        else:
            dalton_inp.append(".DFT")
            dalton_inp.append(xc.capitalize())
        # RESPONSE
        dalton_inp.append("**RESPONSE")
        # LINEAR
        dalton_inp.append("*QUADRA")
        # looks like it's the only option for a response calculation
        if operator == "dipole":
            dalton_inp.append(".DIPLEN")
            freq = self.freq_json[molecule_input][xc][operator]
            num_freq = len(freq)
            dalton_inp.append(".FREQUENCIES")
            dalton_inp.append(str(num_freq))

            freq_s = []
            for f in freq:
                freq_s.append(str(f))
            dalton_inp.append(" ".join(freq_s))

        dalton_inp.append("**END OF DALTON INPUT")
        dalton_inp = "\n".join(dalton_inp)
        print(dalton_inp)
        print(molecule_input)
        print(operator)
        print("self.dalton", self.dalton_dir)

        run_dir = self.dalton_dir.joinpath(xc).joinpath(molecule_input).joinpath(operator)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        # Here I read the madness mol file from the molecules directory
        madness_molecule_file = self.base_dir.joinpath('molecules').joinpath(
            madness_molecule + '.mol')
        mad_to_dal = madnessToDalton(self.base_dir)
        if basis.split("-")[-1] == "uc":
            mol_input = mad_to_dal.madmol_to_dalmol(madness_molecule_file,
                                                    "-".join(basis.split("-")[:-1]))
        else:
            mol_input = mad_to_dal.madmol_to_dalmol(madness_molecule_file, basis)
        dal_run_file = run_dir.joinpath('quad.dal')
        with open(dal_run_file, "w") as file:  # Use file to refer to the file object
            file.write(dalton_inp)
        dalton_molecule_file = run_dir.joinpath(
            molecule_input + '-' + basis.replace("*", "S") + ".mol")
        with open(dalton_molecule_file, "w") as file:  # Use file to refer to the file object
            file.write(mol_input)
        molecule_input = dalton_molecule_file.stem
        dalton_input = dal_run_file.stem
        return run_dir, dalton_input, molecule_input

    def __run_dalton(self, rdir, dfile, mfile):
        dalton = shutil.which('dalton')
        # Change to run directory
        os.chdir(rdir)
        # dalton [.dal] [.mol]
        if self.use_mpi:
            daltonCommand = (
                    dalton + " -N " + str(self.Np) + " " + dfile + " " + mfile
            )
            print(daltonCommand)
        else:
            daltonCommand = "dalton " + dfile + " " + mfile
            print(daltonCommand)
        process = subprocess.Popen(daltonCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        os.chdir(self.base_dir)
        print("Changed Directory to ", self.base_dir)
        return output, error

    @staticmethod
    def __write_excited_input(self, madmol, xc, basis, num_states):
        # Given a molecule, exchange correlation functional, basis an number of states
        # generates a dalton .dal file and writes it in corresponding directory
        # /dalton/[xc]/[madmol]/excited-state

        molname = madmol.split(".")[0]
        # print(molname)
        dalton_inp = []
        dalton_inp.append("**DALTON INPUT")
        dalton_inp.append(".RUN PROPERTIES")
        dalton_inp.append(".DIRECT")
        if basis.split("-")[-1] == "uc":
            dalton_inp.append("*MOLBAS ")
            dalton_inp.append(".UNCONT ")
        dalton_inp.append("**WAVE FUNCTIONS")
        if xc == "hf":
            dalton_inp.append(".HF")
        else:
            dalton_inp.append(".DFT")
            dalton_inp.append(xc.capitalize())
        dalton_inp.append("**PROPERTIES")
        dalton_inp.append(".EXCITA")
        dalton_inp.append("*EXCITA")
        dalton_inp.append(".NEXCIT")

        states = []
        for i in range(10):
            states.append(str(num_states))
        dalton_inp.append("   " + "  ".join(states))

        dalton_inp.append("**END OF DALTON INPUT")
        dalton_inp = "\n".join(dalton_inp)
        run_dir = self.dalton_dir.joinpath(xc).joinpath(molname).joinpath('excited-state')
        print(run_dir)
        mad_to_dal = madnessToDalton(self.base_dir.absolute())
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        madmolfile = self.base_dir.joinpath('molecules').joinpath(madmol + '.mol')
        if basis.split("-")[-1] == "uc":
            mol_input = mad_to_dal.madmol_to_dalmol(madmolfile, "-".join(basis.split("-")[:-1]))
        else:
            mol_input = mad_to_dal.madmol_to_dalmol(madmolfile, basis)
        dal_run_file = run_dir.joinpath('excited.dal')
        with open(dal_run_file, "w") as file:  # Use file to refer to the file object
            file.write(dalton_inp)
        mol_file = run_dir.joinpath(molname + "-" + basis.replace("*", "S") + ".mol")
        with open(mol_file, "w") as file:  # Use file to refer to the file object
            file.write(mol_input)
        molname = mol_file.name.split("/")[-1]
        dalname = dal_run_file.name.split("/")[-1]
        print('molname', molname)
        print('dalname', dalname)
        return run_dir, dalname, molname

    @staticmethod
    def __create_frequency_json(output_json, basis):

        rdata = {
            "xx": [],
            "xy": [],
            "xz": [],
            "yx": [],
            "yy": [],
            "yz": [],
            "zx": [],
            "zy": [],
            "zz": [],
        }

        calcs = {}
        for cals in output_json['simulation']['calculations']:
            calc_type = cals['calculationType']
            calcs[calc_type] = cals
            # print(calc_type)
        # print(calcs)

        r_dict = {}
        r_dict["frequencies"] = []

        e_data = calcs['energyCalculation']
        dipole_data = calcs['Dipole']
        r_data = calcs['LinearResponse']

        p_data = r_data["calculationResults"]
        f_data = r_data["calculationSetup"]["frequencies"]

        r_dict["frequencies"] = f_data
        r_dict["values"] = p_data
        r_dict["calculationTime"] = r_data["calculationTime"]

        return {basis: {"ground": e_data, "response": r_dict, "dipole": dipole_data}}

    @staticmethod
    def __create_quadratic_json(output_json, basis):
        print(output_json)

        rdata = {
            "xx": [],
            "xy": [],
            "xz": [],
            "yx": [],
            "yy": [],
            "yz": [],
            "zx": [],
            "zy": [],
            "zz": [],
        }

        calcs = {}
        for cals in output_json['simulation']['calculations']:
            calc_type = cals['calculationType']
            calcs[calc_type] = cals
            # print(calc_type)
        # print(calcs)

        r_dict = {}
        r_dict["frequencies"] = []

        e_data = calcs['energyCalculation']
        dipole_data = calcs['Dipole']

        return {basis: {"ground": e_data, "response": r_dict, "dipole": dipole_data}}

    @staticmethod
    def __create_excited_json(output_json, basis):
        # generate tables given name of output files and basis_list used to generate output files
        data = {"totalEnergy": [], "nuclearRepulsionEnergy": [], "electronEnergy": []}
        s_dict = {}

        calcs = {}
        for cals in output_json['simulation']['calculations']:
            calc_type = cals['calculationType']
            calcs[calc_type] = cals
            print(calc_type)
        # print(calcs)

        e_data = calcs['energyCalculation']
        dipole_data = calcs['Dipole']
        r_data = calcs['SingletExcitationEnergy']

        for dkeys in data.keys():
            data[dkeys].append(float(e_data["calculationResults"][dkeys]["value"]))
        # print(r_data)

        # sort the frequencies
        freq = r_data["calculationResults"]["frequencies"]
        freq = [f for l in freq for f in l]
        freq = np.array(freq)
        sort_f = np.argsort(freq)
        freq = freq[sort_f]

        # Sort by frequency
        Mode = r_data["calculationResults"]["Mode"]
        Mode = [m for mo in Mode for m in mo]
        Mode = np.array(Mode)
        Mode = Mode[sort_f]
        Sym = r_data["calculationResults"]["Sym"]
        Sym = [s for so in Sym for s in so]
        Sym = np.array(Sym)
        Sym = Sym[sort_f]

        s_dict["Sym"] = Sym.tolist()
        s_dict["Mode"] = Mode.tolist()
        s_dict["freq"] = freq.tolist()
        s_dict["calculationTime"] = r_data["calculationTime"]

        return {basis: {"ground": e_data, "response": s_dict}}

    def polar_json_exists(self, mol, xc, operator, basis):

        run_directory, dal_input, mol_input = self.__write_polar_input(self,
                                                                       mol, xc, operator, basis
                                                                       )
        outfile = "/freq_" + "-".join([mol, basis]) + ".out"
        outfile = run_directory + outfile
        try:
            with open(outfile, "r") as daltonOutput:
                dipole_j = self.get_polar_json(mol, xc, operator, basis)
                if dipole_j is None:
                    raise TypeError(
                        'polar json does not exist for ' + mol + ' ' + xc + ' ' + operator + ' ' + basis)
                return True
        except (FileNotFoundError, KeyError, IndexError, TypeError) as e:
            print(e)
            return False

    def get_polar_json(self, mol, xc, operator, basis):

        run_directory, dal_input, mol_input = self.__write_polar_input(self,
                                                                       mol, xc, operator, basis
                                                                       )
        output_stem = "freq_" + "-".join([mol, basis])
        output_file = run_directory.joinpath(output_stem + ".out")
        output_json = run_directory.joinpath(output_stem + ".json")

        d_out, d_error = None, None
        data = None
        try:
            with open(output_file, "r") as daltonOutput:
                dj = daltonToJson()
                dalton_json = json.loads(dj.convert(daltonOutput))
                with(open(output_json, "w")) as f:
                    f.write(json.dumps(dalton_json, indent=4))
                data = self.__create_frequency_json(dalton_json, basis)
        except (FileNotFoundError, IndexError) as e:
            if self.run:
                print("Trying to run ", output_stem, " in ", run_directory)
                print('dal_input:', dal_input)
                print('mol_input:', mol_input)
                print('run_directory:', run_directory)
                print('output_file:', output_file)
                print('output_json:', output_json)
                d_out, d_error = self.__run_dalton(run_directory, dal_input, mol_input)
                # print(d_out, d_error)
                print("Finished running  ", mol, " in ", run_directory)
                print("Trying to open ", output_file)
                with open(output_file, "r") as daltonOutput:
                    dj = daltonToJson()
                    dalton_json = json.loads(dj.convert(daltonOutput))
                    with(open(output_json, "w")) as f:
                        f.write(json.dumps(dalton_json, indent=4))
                    data = self.__create_frequency_json(dalton_json, basis)
                pass
            else:
                print("Did not find ", basis, " data for", mol, "and dalton is not set to run")
                pass
        except KeyError as e:
            print(d_out, d_error)
            print("KeyError: ", e)
            print(output_stem, " in ", run_directory, " did not run correctly")
            print('most likely the basis set is not available')
            pass

        return data

    def get_quad_json(self, mol, xc, operator, basis):

        run_directory, dal_input, mol_input = self.__write_quadratic_input(self,
                                                                           mol, xc, operator, basis
                                                                           )
        output_stem = "quad_" + "-".join([mol, basis])
        output_file = run_directory.joinpath(output_stem + ".out")
        output_json = run_directory.joinpath(output_stem + ".json")

        d_out, d_error = None, None
        data = None
        try:
            with open(output_file, "r") as daltonOutput:
                dj = daltonToJson()
                dalton_json = json.loads(dj.convert(daltonOutput))
                dalton_json['Quad'] = dj.readQuadResponse(output_file).to_dict()
                with(open(output_json, "w")) as f:
                    f.write(json.dumps(dalton_json, indent=4))
                data = self.__create_quadratic_json(dalton_json, basis)
                data['Quad'] = dj.readQuadResponse(output_file)
        except (FileNotFoundError, IndexError) as e:
            if self.run:
                print("Trying to run ", output_stem, " in ", run_directory)
                print('dal_input:', dal_input)
                print('mol_input:', mol_input)
                print('run_directory:', run_directory)
                print('output_file:', output_file)
                print('output_json:', output_json)
                d_out, d_error = self.__run_dalton(run_directory, dal_input, mol_input)
                # print(d_out, d_error)
                print("Finished running  ", mol, " in ", run_directory)
                print("Trying to open ", output_file)
                with open(output_file, "r") as daltonOutput:
                    dj = daltonToJson()
                    dalton_json = json.loads(dj.convert(daltonOutput))
                    with(open(output_json, "w")) as f:
                        dalton_json['Quad'] = dj.readQuadResponse(output_file)
                        f.write(json.dumps(dalton_json, indent=4))
                    data = self.__create_quadratic_json(dalton_json, basis)
                    data['Quad'] = dj.readQuadResponse(output_file)
                pass
            else:
                print("Did not find ", basis, " data for", mol, "and dalton is not set to run")
                pass
        except KeyError as e:
            print(d_out, d_error)
            print("KeyError: ", e)
            print(output_stem, " in ", run_directory, " did not run correctly")
            print('most likely the basis set is not available')
            pass

        return data

    def get_excited_json(self, mol, xc, basis, num_states):
        """
        Get the excited state json data for a molecule
        :param mol:
        :param xc:
        :param basis: str
        :param run: bool

        :param run:
        """
        run_directory, dal_input, mol_input = self.__write_excited_input(self,
                                                                         mol, xc, basis, num_states
                                                                         )

        print('dal_input:', dal_input)
        print('mol_input:', mol_input)
        print('run_directory:', run_directory)
        run_name = "excited_" + "-".join([mol, basis]) + ".out"
        json_name = "excited_" + "-".join([mol, basis]) + ".json"
        # First look for the output file and try and convert it to a json
        outfile = run_directory.joinpath(run_name)
        outJSON = run_directory.joinpath(json_name)
        print('outfile:', outfile)
        data = None
        try:
            # open the output file
            with open(outfile, "r") as daltonOutput:
                dj = daltonToJson()
                dalton_json = json.loads(dj.convert(daltonOutput))
                with(open(outJSON, "w")) as f:
                    f.write(json.dumps(dalton_json, indent=4))
            data = self.__create_excited_json(dalton_json, basis)
        except FileNotFoundError as e:
            print("did not find output file", e)
            if self.run:
                print("Try and run molecule ", mol)
                d_out, d_error = self.__run_dalton(run_directory, dal_input, mol_input)
                # print(d_out, d_error)
                print("Finished running  ", mol, " in ", run_directory)
                print("Trying to open ", outfile)
                with open(outfile, "r") as daltonOutput:
                    dj = daltonToJson()
                    dalton_json = json.loads(dj.convert(daltonOutput))
                    with(open(outJSON, "w")) as f:
                        f.write(json.dumps(dalton_json, indent=4))
                    data = self.__create_excited_json(dalton_json, basis)
                pass
            else:
                print("Not trying to run dalton for ", mol)
                pass
        return data

    def get_excited_result(self, mol, xc, basis, run):
        excited_j = self.get_excited_json(mol, xc, basis, 4)
        if excited_j is not None:

            time = excited_j[basis]["ground"]["calculationTime"]
            results = excited_j[basis]["ground"]["calculationResults"]
            gR = {}
            gR["basis"] = basis
            # results
            gkeys = ["totalEnergy", "nuclearRepulsionEnergy", "electronEnergy"]
            for g in gkeys:
                gR[g] = float(results[g]["value"])
            # timings
            tkeys = ["cpuTime", "wallTime"]

            for t in tkeys:
                gR["g" + t] = float(time[t])
            rtime = excited_j[basis]["response"]["calculationTime"]
            for t in tkeys:
                gR["r" + t] = float(rtime[t])

            # number of electrons
            skeys = ["numberOfElectrons"]
            setup = excited_j[basis]["ground"]["calculationSetup"]
            for s in skeys:
                gR[s] = setup[s]

            gSeries = pd.Series(gR)
            rresults = excited_j[basis]["response"]
            ekeys = ["Sym", "Mode", "freq"]

            rR = {}
            for e in ekeys:
                rR[e] = rresults[e]
            rDf = pd.DataFrame.from_dict(rR)

            return gSeries, rDf
        else:
            return None

    def get_frequency_result(self, mol, xc, operator, basis):
        dipole_j = self.get_polar_json(mol, xc, operator, basis)[basis]
        time = dipole_j["ground"]["calculationTime"]
        results = dipole_j["ground"]["calculationResults"]
        ground_dipole = dipole_j["dipole"]["calculationResults"]
        gR = {}
        gR["basis"] = basis
        gR["dipole"] = pd.Series(ground_dipole)
        # results
        gkeys = ["totalEnergy", "nuclearRepulsionEnergy", "electronEnergy"]
        for g in gkeys:
            gR[g] = float(results[g]["value"])
        # timings
        tkeys = ["cpuTime", "wallTime"]
        for t in tkeys:
            gR["g" + t] = float(time[t])
        rtime = dipole_j["response"]["calculationTime"]
        for t in tkeys:
            gR["r" + t] = float(rtime[t])
        # number of electrons
        skeys = ["numberOfElectrons"]
        setup = dipole_j["ground"]["calculationSetup"]
        for s in skeys:
            gR[s] = setup[s]
        # ground results
        gSeries = pd.Series(gR)

        # response results
        rresults = dipole_j["response"]
        rkeys = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]
        rdict = {}
        rdict["frequencies"] = rresults["frequencies"]

        for r in rkeys:
            rdict[r] = rresults["values"][r]
        rdf = pd.DataFrame(rdict)

        return gSeries, rdf
