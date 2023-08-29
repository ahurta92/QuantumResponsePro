import subprocess
import seaborn as sns
import numpy as np

import matplotlib.pyplot as plt
from setuptools import glob

from .madnessReader import FrequencyData
from ..dalton.daltonrunner import DaltonRunner
from ..madness_to_dalton import *


def tensor_to_numpy(j):
    array = np.empty(j["size"])
    array[:] = j["vals"]
    return np.reshape(array, tuple(j["dims"]))


class MadnessReader:
    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.molecule_dir = data_dir.joinpath('molecules')
        # TODO: Need to make this json_data/frequency.json
        with open(self.data_dir.joinpath('json_data/frequency.json')) as \
                json_file:
            self.freq_json = json.loads(json_file.read())

    def __read_protocol_polarizability_data(self, protocol_data: json):
        polar_dfs = []
        protos = []
        kprotos = []
        iters = []
        iter_p = 0
        num_iters_per_protocol = []
        for proto in protocol_data:
            protos.append(proto["proto"])
            kprotos.append(proto["k"])
            num_iters = proto["iter_data"].__len__()
            num_iters_per_protocol.append(num_iters)
            proto_array = np.ones((num_iters, 1)) * proto["proto"]
            kproto_array = np.ones((num_iters, 1)) * proto["k"]
            polar_data = np.empty((num_iters, 9))
            i = 0
            for iter in proto["property_data"]:
                alpha = tensor_to_numpy(iter["polar"]).flatten()
                polar_data[i, :] = alpha
                i += 1
                iters.append(iter_p)
                iter_p += 1
            kproto_df = pd.DataFrame(kproto_array, columns=["k"])
            proto_df = pd.DataFrame(proto_array, columns=["thresh"])
            polar_df = pd.DataFrame(
                polar_data,
                columns=["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"],
            )
            polar_df = pd.concat([kproto_df, proto_df, polar_df], axis=1)
            polar_dfs.append(polar_df)

        iters_df = pd.DataFrame(iters, columns=["iterations"])
        final_polar = pd.concat(polar_dfs, ignore_index=True)
        final_polar = pd.concat([iters_df, final_polar], axis=1)

        return final_polar

    def __read_response_protocol_data(self, protocol_data: json, num_states, num_orbitals):
        num_protocols = protocol_data.__len__()
        # print("Number of protocols ", num_protocols)
        x_keys = []
        ax_keys = []
        rx_keys = []
        d_keys = []
        ad_keys = []
        for i in range(num_states):
            x_keys.append("X" + str(i))
            ax_keys.append("abs_X" + str(i))
            rx_keys.append("rel_X" + str(i))

            d_keys.append("D" + str(i))
            ad_keys.append("abs_D" + str(i))
        axij_keys = []
        for i in range(num_states):
            for j in range(num_orbitals):
                axij_keys.append('abs_x' + str(i) + str(j))
        for i in range(num_states):
            for j in range(num_orbitals):
                axij_keys.append('abs_y' + str(i) + str(j))
        x_norm_dfs = []
        x_abs_error_dfs = []
        d_norm_dfs = []
        d_abs_error_dfs = []
        protos = []
        kprotos = []
        iters = []
        iter_p = 0
        num_iters_per_protocol = []
        for kk in range(num_protocols):
            proto = protocol_data[kk]
            # print("proto_iter_data", kk, proto["iter_data"])
            num_iters = proto["iter_data"].__len__()
            # print(" number of iterations in this stage: ", num_iters)
            protos.append(proto["proto"])
            kprotos.append(proto["k"])
            num_iters_per_protocol.append(num_iters)
            proto_array = np.ones((num_iters, 1)) * proto["proto"]
            kproto_array = np.ones((num_iters, 1)) * proto["k"]
            x_norms = np.empty((num_iters, num_states))
            x_abs_error = np.empty((num_iters, num_states))
            d_norms = np.empty((num_iters, num_states))
            d_abs_error = np.empty((num_iters, num_states))
            i = 0
            for iter in proto["iter_data"]:
                x_norms[i, :] = tensor_to_numpy(iter["x_norms"]).flatten()
                x_abs_error[i, :] = tensor_to_numpy(iter["x_abs_error"]).flatten()
                d_norms[i, :] = tensor_to_numpy(iter["rho_norms"]).flatten()
                d_abs_error[i, :] = tensor_to_numpy(iter["rho_abs_error"]).flatten()
                i += 1
                iters.append(iter_p)
                iter_p += 1
            x_norm_dfs.append(pd.DataFrame(x_norms, columns=x_keys))
            x_abs_error_dfs.append(pd.DataFrame(x_abs_error, columns=ax_keys))
            d_norm_dfs.append(pd.DataFrame(d_norms, columns=d_keys))
            d_abs_error_dfs.append(pd.DataFrame(d_abs_error, columns=ad_keys))

        for j in range(1, num_iters_per_protocol.__len__()):
            num_iters_per_protocol[j] = (num_iters_per_protocol[j] + num_iters_per_protocol[j - 1])

        x1 = pd.concat(x_norm_dfs)
        xa = pd.concat(x_abs_error_dfs)

        d1 = pd.concat(d_norm_dfs)
        da = pd.concat(d_abs_error_dfs)

        iters_df = pd.Series(iters)
        iters_df.name = "iterations"
        full = pd.concat([x1, xa, d1, da], axis=1)
        full = pd.concat([iters_df, full.reset_index(drop=True)], axis=1)
        full.index += 1

        return num_iters_per_protocol, full

    def __open_ground_json(self, mol, xc):
        moldir = self.data_dir.joinpath(xc, mol)
        path = moldir.joinpath("moldft.calc_info.json")
        with open(path) as json_file:
            response_j = json.loads(json_file.read())
        return response_j

    def __open_ground_scf_json(self, mol, xc):
        moldir = self.data_dir.joinpath(xc, mol)
        path = moldir.joinpath("moldft.scf_info.json")
        try:
            with open(path) as json_file:
                scf_j = json.loads(json_file.read())
        except FileNotFoundError as f:
            # print(f)
            scf_j = None

        return scf_j

    def get_ground_scf_data(self, mol, xc):

        j = self.__open_ground_json(mol, xc)
        scf_j = self.__open_ground_scf_json(mol, xc)
        if scf_j:
            j.update(scf_j)
        params = j["parameters"]
        scf_e_data = j["scf_e_data"]
        timing = j["wall_time"]

        return params, scf_e_data, timing, j

    def __open_excited_rbj(self, mol, xc, num_states):

        moldir = self.data_dir + "/" + xc + "/" + mol
        dfile = "excited-" + str(num_states)
        jsonf = "response_base.json"

        path = "/".join([moldir, dfile, jsonf])

        with open(path) as json_file:
            response_j = json.loads(json_file.read())

        return response_j

    def get_excited_data(self, mol, xc):
        num_states = self.freq_json[mol][xc]["excited-state"]
        rbasej = self.__open_excited_rbj(mol, xc, num_states)
        num_orbitals = rbasej["parameters"]["num_orbitals"]
        converged = rbasej["converged"]
        num_iters_per_protocol, function_data = self.__read_response_protocol_data(
            rbasej["protocol_data"], num_states, num_orbitals
        )
        omega = self.__read_protocol_excited_state_data(
            rbasej["protocol_data"], num_states, num_orbitals
        )

        params = rbasej["parameters"]
        wall_time = rbasej["wall_time"]

        return params, wall_time, converged, num_iters_per_protocol, rbasej, function_data, omega

    def __open_frequency_rbj(self, mol, xc, operator, freq):

        sfreq = "%f" % freq
        # first number before decimal
        f1 = sfreq.split(".")[0]
        # second number after decimal
        f2 = sfreq.split(".")[1]

        moldir = self.data_dir.joinpath(xc, mol)
        dfile = operator + "_" + xc + "_" + f1 + "-" + f2
        jsonf = "response_base.json"

        path = moldir.joinpath(dfile, jsonf)
        # print(path)

        with open(path) as json_file:
            response_j = json.loads(json_file.read())

        return response_j

    def __get_polar_data(self, rbase_j):
        params = rbase_j["parameters"]
        num_orbitals = params["num_orbitals"]
        num_states = params["states"]

        (
            num_iters_per_protocol,
            full_function_data,
        ) = self.__read_response_protocol_data(rbase_j["protocol_data"], num_states, num_orbitals)
        polarizability_data = self.__read_protocol_polarizability_data(rbase_j['protocol_data'])
        polarizability_data.index += 1
        return params, num_iters_per_protocol, full_function_data, polarizability_data

    def get_polar_result(self, mol, xc, operator):
        freq = self.freq_json[mol][xc][operator]
        polar_data = {}
        last_polar_freq = {}
        function_data = {}
        time_data = {}
        converged = {}
        num_iter_proto = {}
        full_params = {}
        full_response_base = {}
        for f in freq:
            try:
                response_base_json = self.__open_frequency_rbj(mol, xc, operator, f)
                full_response_base[str(f)] = response_base_json
                converged_f = response_base_json["converged"]
                # print(converged_f)
                params, num_iters_per_protocol, full_function_data, polarizability_data = self.__get_polar_data(
                    response_base_json)
                full_params[str(f)] = params
                polar_data[str(f)] = polarizability_data
                # print(polarizability_data)
                last_polar_freq[str(f)] = polarizability_data.iloc[-1, :]

                num_iter_proto[str(f)] = num_iters_per_protocol
                function_data[str(f)] = pd.DataFrame(full_function_data)
                # fdata[str(f)] = full_function_data.iloc[-1, :]
                converged[str(f)] = converged_f

                time_data[str(f)] = response_base_json["time_data"]

            except FileNotFoundError as not_found:
                # print(f, " not found:", not_found)
                pass

        return (
            full_params,
            time_data,
            pd.Series(converged, dtype=bool),
            num_iter_proto,
            full_response_base,
            function_data,
            polar_data,
            pd.DataFrame(last_polar_freq).T
        )


class MadnessResponse:
    def __init__(self, mol, xc, operator, data_dir):
        self.data_dir = data_dir
        self.ground_info = None
        self.mol = mol
        self.xc = xc
        self.operator = operator
        self.moldir = self.data_dir.joinpath(self.xc).joinpath(self.mol)

        mad_reader = MadnessReader(self.data_dir)
        (
            self.ground_params,
            self.ground_scf_data,
            self.ground_timing,
            self.ground_info,
        ) = mad_reader.get_ground_scf_data(mol, xc)
        e_name_list = ["e_coulomb", "e_kinetic", "e_local", "e_nrep", "e_tot"]
        self.ground_e = {}
        for e_name in e_name_list:
            self.ground_e[e_name] = self.ground_scf_data[e_name][-1]

        (
            self.params,
            self.time_data,
            self.converged,
            self.num_iter_proto,
            self.response_base,
            self.function_data,
            self.full_polar_data,
            self.polar_data,
        ) = mad_reader.get_polar_result(mol, xc, operator)
        self.quad_data = self.get_quad_data(self.moldir)

        self.num_states = self.params['0.0']["states"]
        self.num_orbitals = self.params['0.0']["num_orbitals"]
        self.data = {"ground": self.__get_ground_data()}
        fdata, pdata = self.__get_response_data()
        self.frequencies = list(self.response_base.keys())
        self.data["response"] = {}
        self.data["response"]["function"] = fdata
        self.data["response"]["polarizability"] = pdata

        mad_reader = MadnessReader(self.data_dir)
        self.freq_json = mad_reader.freq_json
        self.calc_dirs = {}
        freq = self.freq_json[self.mol][self.xc][self.operator]
        for i, f in enumerate(freq):
            sfreq = "%f" % f
            f1 = sfreq.split(".")[0]
            f2 = sfreq.split(".")[1]
            dfile = operator + "_" + xc + "_" + f1 + "-" + f2
            freq_i_path = self.moldir.joinpath(dfile)
            self.calc_dirs[i] = freq_i_path
        try:
            self.data["convergence"] = self.get_convergence_dict()
        except KeyError:
            print(KeyError, "Convergence data not found for frequency ", mol, " : ", i, ",", f)
            pass

    def get_quad_data(self, moldir):
        beta_json_path = moldir.joinpath('beta.json')
        loaded_json = json.loads(beta_json_path.read_text())
        # add a column for the molecule in the beginning
        beta_json = pd.DataFrame(loaded_json)
        beta_json['molecule'] = self.mol
        beta_json['basis'] = 'MRA'
        # drop the dash from A-freq and B-freq and C-freq columns names
        beta_json.columns = beta_json.columns.str.replace('-', '')

        # set the A B C columns to categorical type
        print(beta_json)
        return beta_json

    def get_response_calc_data_dict(self, thresh, data_k):
        pdx = self.polar_data.keys()[3:]
        data_r = {}

        for k, val in data_k.items():
            try:
                if k == 'r_x':
                    df = pd.DataFrame(tensor_to_numpy(val[0]), columns=pdx).pow(1 / 2.)
                    data_r[k] = df.dropna(axis=1)
                if k == 'x_relative_residuals' or k == 'density_norms' or k == 'd' or k == 'r_d' or k == "density_residuals":
                    df = pd.DataFrame(tensor_to_numpy(val[0]), columns=['x', 'y', 'z'])
                    # if k == 'r_d':
                    #    print(k, df)
                    data_r[k] = df.dropna(axis=1)
                elif k == 'alpha' or k == 'r_alpha':
                    df = pd.DataFrame(tensor_to_numpy(val[0]), columns=pdx)
                    data_r[k] = df.dropna(axis=1)
                else:
                    df = pd.DataFrame(np.sqrt(np.absolute(tensor_to_numpy(val[0]))), columns=pdx)
                    df = df[df >= thresh]
                    data_r[k] = df.dropna(axis=1)
            except (TypeError, ValueError) as te:
                # print(te,k, val)
                pass
        return data_r

    def get_convergence_dict(self):
        convergence_dict = {}
        for freq_i in self.response_base.keys():
            data_k = self.response_base[freq_i]['response_data']['data']
            convergence_dict[freq_i] = self.get_response_calc_data_dict(0.00, data_k)
        return convergence_dict

    def get_function_keys(self):
        x_keys = []
        ax_keys = []
        d_keys = []
        ad_keys = []

        for i in range(self.num_states):
            x_keys.append("X" + str(i))
            ax_keys.append("abs_X" + str(i))

            d_keys.append("D" + str(i))
            ad_keys.append("abs_D" + str(i))
        axij_keys = []

        for i in range(self.num_states):
            for j in range(self.num_orbitals):
                axij_keys.append('abs_x' + str(i) + str(j))

        for i in range(self.num_states):
            for j in range(self.num_orbitals):
                axij_keys.append('abs_y' + str(i) + str(j))

        return {"x_norms": x_keys, "x_abs_error": ax_keys, "d_norms": d_keys,
                "d_abs_error": ad_keys, "xij_abs_error": axij_keys}

    def __get_ground_precision(self):
        gprec = self.ground_info["precision"]
        ground_data = {}
        for key, val in gprec.items():
            if key == "dconv":
                ground_data["g-dconv"] = val
            elif key == "thresh":
                ground_data["g-thresh"] = val
            elif key == "k":
                ground_data["g-k"] = val
            elif key == "eprec":
                ground_data["g-eprec"] = val
        return pd.Series(ground_data)

    def __get_response_precision(self, omega):
        p_dict = self.response_base[omega]
        r_prec = {}
        for key, val in p_dict.items():
            if key == "dconv":
                r_prec["r-dconv"] = val
            elif key == "thresh":
                r_prec["r-thresh"] = val
            elif key == "k":
                r_prec["r-k"] = val
        return pd.Series(r_prec, dtype=float)

    def __get_response_data(self):
        polar_keys = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']
        ff = []
        pp = []
        g_precision = self.__get_ground_precision()

        for om, data in self.response_base.items():
            if data["converged"]:
                fdo = self.function_data[om].iloc[-1, 1:]
                pdo = self.polar_data.loc[om, polar_keys]
                fd = {'frequency': om}
                freq = pd.Series(fd)
                r_prec = self.__get_response_precision(om)
                ff.append(pd.concat([freq, g_precision, r_prec, fdo]))
                pp.append(pd.concat([freq, g_precision, r_prec, pdo]))
        fdata = pd.DataFrame(ff)
        pdata = pd.DataFrame(pp)
        return fdata, pdata

    def __get_ground_data(self):
        ground_data = self.__get_ground_precision()
        ground_data.update(self.ground_e)
        e_data = pd.Series(ground_data)
        return e_data

    def plot_residuals(self, frequency=0):
        figsize = (10, 8 * 1)
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, constrained_layout=False)
        freq_key = self.frequencies[frequency]
        rd = self.data['convergence'][freq_key]['density_residuals'] \
            .rename(columns={'x': 'dx', 'y': 'dy', 'z': 'dz'})
        rx = self.data['convergence'][freq_key]['x_relative_residuals'] \
                 .loc[:, ['x', 'y', 'z']] \
            .rename(columns={'x': 'rx', 'y': 'ry', 'z': 'rz'})

        rd.plot(logy=True, ax=ax, colormap='Accent', markersize=12, kind='line', style='.-')
        rx.plot(logy=True, ax=ax, colormap='Accent', grid=True, markersize=12, kind='line',
                style='*-')
        iters = self.num_iter_proto[self.frequencies[frequency]]

        threshold = self.response_base[freq_key]["response_data"]['thresh']
        density_target = self.response_base[freq_key]['response_data']['density_target']
        bsh_target = self.response_base[freq_key]["response_data"]['bsh_target']

        for n, pc in enumerate(iters):

            ax.axvline(x=pc - 1, ymin=0, ymax=1, c="black", linestyle="dashed", alpha=0.5)
            if n == 0:
                ax.axhline(y=threshold[n], xmin=0, xmax=iters[-1], c="black", linestyle="dashed",
                           alpha=0.5,
                           label="Threshold")
                ax.axhline(y=density_target[n], xmin=0, xmax=iters[-1], c="red", linestyle="dashed",
                           alpha=0.5,
                           label="Density Target")
                ax.axhline(y=bsh_target[n], xmin=0, xmax=iters[-1], c="blue", linestyle="dashed",
                           alpha=0.5,
                           label="BSH Target")
            else:
                ax.axhline(y=threshold[n], xmin=0, xmax=iters[-1], c="black", linestyle="dashed",
                           alpha=0.5,
                           )
                ax.axhline(y=density_target[n], xmin=0, xmax=iters[-1], c="red", linestyle="dashed",
                           alpha=0.5,
                           )
                ax.axhline(y=bsh_target[n], xmin=0, xmax=iters[-1], c="blue", linestyle="dashed",
                           alpha=0.5,
                           )

        ax.grid(which="both")
        ax.minorticks_on()
        ax.tick_params(which="both", top="on", left="on", right="on", bottom="on", )
        ax.legend(loc='center', bbox_to_anchor=(0.5, 1.15), fancybox=True, shadow=False, ncol=3, )

        return fig, ax

    def get_line_plots(self):
        x_plots = {}
        for i, dir in self.calc_dirs.items():
            x_plot_path = dir + "/plots/.dat"

        pass

    def compare_to_mra(self, database_compare, op):

        new_x = self.polar_data[op].copy()
        nC = {}
        oC = {}
        for o in op:
            nC[o] = 'self-' + o
            oC[o] = 'compare-' + o
        new_x.rename(columns=nC, inplace=True)
        # read from two different versions of the code
        try:
            mad_compare = FrequencyData(self.mol, self.xc, self.operator, database_compare)
            old_x = mad_compare.polar_df[op].copy()

        except:
            mad_compare = MadnessResponse(self.mol, self.xc, self.operator, database_compare)
            old_x = mad_compare.polar_data[op].copy()

        old_x.rename(columns=oC, inplace=True)
        p = pd.concat([old_x, new_x], axis=1).reset_index(drop=True)
        p.plot()
        return p

    def compare_mra_convergence(self, database_compare, frequency, value, ij):
        try:
            this_data = self.data['convergence'][frequency][value][ij].copy()
            mad_compare = MadnessResponse(self.mol, self.xc, self.operator, database_compare)
            compare_data = mad_compare.data['convergence'][frequency][value][ij].copy()
        except TypeError as t:
            print(t)
            return
        except KeyError as t:
            print(self.data['convergence'][frequency][value])
            print(t)
            return
        except FileNotFoundError as f:
            print(f)
            print('file not found', database_compare)
            mad_compare = FrequencyData(self.mol, self.xc, self.operator, database_compare)
            return
        nC = {}
        oC = {}
        for o in ij:
            nC[o] = 'self-' + o
            oC[o] = 'compare-' + o
        this_data.rename(columns=nC, inplace=True)
        compare_data.rename(columns=oC, inplace=True)
        # read from two different versions of the code
        p = pd.concat([this_data, compare_data], axis=1).reset_index(drop=True)

        p.plot(title='compare: ' + value + ' at frequency ' + frequency)
        return p

    def compare_to_dalton(self, database_compare, op, basis_sets):

        new_x = self.polar_data[op].copy()
        nC = {}
        for o in op:
            nC[o] = 'mra-' + o
        new_x.rename(columns=nC, inplace=True)
        # read from two different versions of the code
        d = DaltonRunner(database_compare, False)
        bD = []
        for basis in basis_sets:
            try:
                bC = {}
                for o in op:
                    bC[o] = basis + '-' + o
                ground, response = d.get_frequency_result(self.mol, self.xc, self.operator, basis)
                basis_polar_df = response[op].copy()
                basis_polar_df.rename(columns=bC, inplace=True)
                bD.append(basis_polar_df)
            except TypeError as t:
                print("did not find basis set data")
                print(t)
                continue

        basis_polar_df = pd.concat(bD, axis=1).reset_index(drop=True)
        p = pd.concat([basis_polar_df, new_x.reset_index(drop=True)], axis=1).reset_index(drop=True)
        p.plot()
        plt.show()
        return p


class ExcitedData:
    def __init__(self, mol, xc):
        self.mol = mol
        self.xc = xc
        mad_reader = MadnessReader()
        (
            self.ground_params,
            self.ground_scf_data,
            self.ground_timing,
            self.ground_info
        ) = mad_reader.get_ground_scf_data(mol, xc)
        e_name_list = ["e_coulomb", "e_kinetic", "e_local", "e_nrep", "e_tot"]
        self.ground_e = {}
        for e_name in e_name_list:
            self.ground_e[e_name] = self.ground_scf_data[e_name][-1]
        (
            self.params,
            self.wall_time,
            self.converged,
            self.num_iter_proto,
            self.response_base,
            self.function_data,
            self.omega
        ) = mad_reader.get_excited_data(mol, xc)
        self.num_states = self.params["states"]
        self.num_orbitals = self.params["num_orbitals"]


class MadRunner:
    def run_response(self, mol, xc, operator, prec):

        if operator == "dipole":
            mad_command = "mad-freq"
        elif operator == "excited-state":
            mad_command = "mad-excited"
        else:
            print("not implemented yet")
            return 1
        madnessCommand = " ".join([mad_command, mol, xc, prec])

        process = subprocess.Popen(madnessCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error

    def run_madness_ground(self, mol, xc):
        mad_command = "database-moldft"
        madnessCommand = " ".join([mad_command, mol, xc])
        process = subprocess.Popen(madnessCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
