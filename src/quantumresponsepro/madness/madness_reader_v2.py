import subprocess
from src.quantumresponsepro.madness.madnessReader import FrequencyData
import seaborn as sns

import matplotlib.pyplot as plt
from setuptools import glob

from src.quantumresponsepro.dalton.dalton import Dalton

import numpy as np

from src.quantumresponsepro.madness_to_dalton import *


def tensor_to_numpy(j):
    array = np.empty(j["size"])
    array[:] = j["vals"]
    return np.reshape(array, tuple(j["dims"]))


class MadnessReader:
    def __init__(self, data_dir):
        self.data_dir = data_dir
        with open(self.data_dir + "/molecules/frequency.json") as json_file:
            self.freq_json = json.loads(json_file.read())

    def __read_protocol_excited_state_data(self, protocol_data: json, num_states, num_orbitals):
        num_protocols = protocol_data.__len__()
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
            omega_data = np.empty((num_iters, num_states))
            i = 0
            for iter in proto["property_data"]:
                # diagonalize the polarizability
                alpha = tensor_to_numpy(iter["omega"]).flatten()
                # alpha=.5*(alpha+alpha.transpose())
                # w, v = LA.eig(alpha)
                # print("alpha : ",alpha)
                omega_data[i, :] = alpha
                i += 1
                iters.append(iter_p)
                iter_p += 1
            kproto_df = pd.DataFrame(kproto_array, columns=["k"])
            proto_df = pd.DataFrame(proto_array, columns=["thresh"])
            polar_df = pd.DataFrame(
                omega_data
            )
            polar_df = pd.concat([kproto_df, proto_df, polar_df], axis=1)
            polar_dfs.append(polar_df)
        iters_df = pd.DataFrame(iters, columns=["iterations"])
        final_polar = pd.concat(polar_dfs, ignore_index=True)
        final_polar = pd.concat([iters_df, final_polar], axis=1)

        return final_polar

    def __read_protocol_polarizability_data(self, protocol_data: json, num_states, num_orbitals):
        num_protocols = protocol_data.__len__()
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
        moldir = self.data_dir + "/" + xc + "/" + mol
        jsonf = "moldft.calc_info.json"
        path = "/".join([moldir, jsonf])
        with open(path) as json_file:
            response_j = json.loads(json_file.read())
        return response_j

    def __open_ground_scf_json(self, mol, xc):
        moldir = self.data_dir + "/" + xc + "/" + mol
        jsonf = "moldft.scf_info.json"
        path = "/".join([moldir, jsonf])
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

        moldir = self.data_dir + "/" + xc + "/" + mol
        dfile = operator + "_" + xc + "_" + f1 + "-" + f2
        jsonf = "response_base.json"

        path = "/".join([moldir, dfile, jsonf])
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
        polarizability_data = self.__read_protocol_polarizability_data(rbase_j['protocol_data'], num_states,
                                                                       num_orbitals)
        polarizability_data.index += 1
        return params, num_iters_per_protocol, full_function_data, polarizability_data

    # TODO get the ground data
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


def get_function_keys(num_states, num_orbitals):
    x_keys = []
    ax_keys = []

    d_keys = []
    ad_keys = []

    for i in range(num_states):
        x_keys.append("X" + str(i))
        ax_keys.append("abs_X" + str(i))

        d_keys.append("D" + str(i))
        ad_keys.append("abs_D" + str(i))

    xij_keys = []
    axij_keys = []

    for i in range(num_states):
        for j in range(num_orbitals):
            xij_keys.append('x' + str(i) + str(j))
            axij_keys.append('abs_x' + str(i) + str(j))

    for i in range(num_states):
        for j in range(num_orbitals):
            xij_keys.append('y' + str(i) + str(j))
            axij_keys.append('abs_y' + str(i) + str(j))

    return {"x_norms": x_keys, "x_abs_error": ax_keys, "d_norms": d_keys,
            "d_abs_error": ad_keys, "xij_norms": xij_keys, "xij_abs_error": axij_keys}


class ResponseCalc:
    def __init__(self, mol, xc, operator, data_dir):
        self.data_dir = data_dir
        self.ground_info = None
        self.mol = mol
        self.xc = xc
        self.operator = operator
        self.moldir = self.data_dir + "/" + self.xc + "/" + self.mol

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

        self.num_states = self.params['0.0']["states"]
        self.num_orbitals = self.params['0.0']["num_orbitals"]
        self.function_keys = self.get_function_keys()
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
            freq_i_path = "/".join([self.moldir, dfile])
            self.calc_dirs[i] = freq_i_path
        try:
            self.data["convergence"] = self.get_convergence_dict()
        except KeyError:
            print(KeyError)
            pass

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
                    #if k == 'r_d':
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
        rx_keys = []

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

    def __plot_norm_and_residual_freq(self, num_i, ax):
        fkeys = get_function_keys(self.num_states, self.num_orbitals)
        abs_keys = fkeys["x_abs_error"]
        chi_norms_keys = fkeys["x_norms"]
        f_key = list(self.function_data.keys())[num_i]
        dconv = self.params[f_key]["dconv"]

        chi_norms_i = self.function_data[f_key][chi_norms_keys]
        bsh_residuals_i = self.function_data[f_key][abs_keys]
        chi_norms_i.plot(ax=ax[0], logy=False, legend=True, colormap='Accent', title='Chi Norms', marker='*',
                         markersize=12, grid=True)
        bsh_residuals_i.plot(ax=ax[1], logy=True, legend=True, colormap='Accent', title='Absolute BSH Residuals',
                             marker='*', markersize=12,
                             grid=True)

        iters = self.num_iter_proto[f_key]

        for pc in iters:
            for i in range(2):
                ax[i].axvline(x=pc, ymin=0, ymax=1, c="black", linestyle="dashed")
        for i in range(2):
            if i != 0:
                ax[i].axhline(y=dconv, xmin=0, xmax=iters[-1], c="black", linestyle="dashed", )
            ax[i].grid(which="both")
            ax[i].minorticks_on()
            ax[i].tick_params(which="both", top="on", left="on", right="on", bottom="on", )

    def plot_freq_norm_and_residual(self, save, only_static=True):
        xkeys = []
        ykeys = []
        for i in range(self.num_states):
            xkeys.append("x" + str(i))
            ykeys.append("y" + str(i))
        freq = list(self.response_base.keys())
        num_ran = len(self.converged)
        num_freqs = len(freq)
        num_plot = num_ran
        if only_static:
            num_plot = 1
        figsize = (15, 8 * num_plot)
        fig, ax = plt.subplots(nrows=num_plot, ncols=2, figsize=figsize, constrained_layout=False)
        for i in range(num_plot):
            if num_plot > 1:
                axi = ax[i]
            else:
                axi = ax
            title = 'Polarizability Convergence: ' + self.mol + r'  $\omega({}/{})$'.format(i, num_freqs - 1)
            self.__plot_norm_and_residual_freq(i, axi)
            plotname = 'freq_{}'.format(i) + ".svg"
            if save:
                if not os.path.exists("convergence"):
                    os.mkdir("convergence")
                if not os.path.exists('convergence/' + self.mol):
                    os.mkdir("convergence/" + self.mol)
                if not os.path.exists('convergence/' + self.mol + '/' + self.xc):
                    os.mkdir("convergence/" + self.mol + '/' + self.xc)
                plt.savefig("convergence/" + self.mol + '/' + self.xc + '/' + plotname)
        print(self.mol + "\n converged: ", self.converged)
        return fig, ax

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
        rx.plot(logy=True, ax=ax, colormap='Accent', grid=True, markersize=12, kind='line', style='*-')
        iters = self.num_iter_proto[self.frequencies[frequency]]

        threshold = self.response_base[freq_key]["response_data"]['thresh']
        density_target = self.response_base[freq_key]['response_data']['density_target']
        bsh_target = self.response_base[freq_key]["response_data"]['bsh_target']

        for n, pc in enumerate(iters):

            ax.axvline(x=pc - 1, ymin=0, ymax=1, c="black", linestyle="dashed", alpha=0.5)
            if n == 0:
                ax.axhline(y=threshold[n], xmin=0, xmax=iters[-1], c="black", linestyle="dashed", alpha=0.5,
                           label="Threshold")
                ax.axhline(y=density_target[n], xmin=0, xmax=iters[-1], c="red", linestyle="dashed", alpha=0.5,
                           label="Density Target")
                ax.axhline(y=bsh_target[n], xmin=0, xmax=iters[-1], c="blue", linestyle="dashed", alpha=0.5,
                           label="BSH Target")
            else:
                ax.axhline(y=threshold[n], xmin=0, xmax=iters[-1], c="black", linestyle="dashed", alpha=0.5,
                           )
                ax.axhline(y=density_target[n], xmin=0, xmax=iters[-1], c="red", linestyle="dashed", alpha=0.5,
                           )
                ax.axhline(y=bsh_target[n], xmin=0, xmax=iters[-1], c="blue", linestyle="dashed", alpha=0.5,
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
            mad_compare = ResponseCalc(self.mol, self.xc, self.operator, database_compare)
            old_x = mad_compare.polar_data[op].copy()

        old_x.rename(columns=oC, inplace=True)
        p = pd.concat([old_x, new_x], axis=1).reset_index(drop=True)
        p.plot()
        return p

    def compare_mra_convergence(self, database_compare, frequency, value, ij):
        try:
            this_data = self.data['convergence'][frequency][value][ij].copy()
            mad_compare = ResponseCalc(self.mol, self.xc, self.operator, database_compare)
            compare_data = mad_compare.data['convergence'][frequency][value][ij].copy()
        except TypeError as t:
            print(t)
            return
        except KeyError as t:
            print(self.data['convergence'][frequency][value])
            print(t)
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
        d = Dalton(database_compare, False)
        bD = []
        for basis in basis_sets:
            bC = {}
            for o in op:
                bC[o] = basis + '-' + o
            ground, response = d.get_frequency_result(self.mol, self.xc, self.operator, basis)
            basis_polar_df = response[op].copy()
            basis_polar_df.rename(columns=bC, inplace=True)
            bD.append(basis_polar_df)

        basis_polar_df = pd.concat(bD, axis=1).reset_index(drop=True)
        p = pd.concat([basis_polar_df, new_x.reset_index(drop=True)], axis=1).reset_index(drop=True)
        p.plot()
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

    def compare_dalton(self, basis, base_dir):
        dalton_reader = Dalton(base_dir)
        ground_dalton, response_dalton = dalton_reader.get_excited_result(self.mol, self.xc, basis, True)

        ground_compare = pd.concat(
            [
                ground_dalton,
                pd.Series(self.ground_timing, index=["wall_time"]),
                pd.Series(self.ground_e),
            ]
        )
        omega_df = response_dalton.iloc[0: self.num_states]
        omega_df.loc[:, "mad-omega"] = self.omega
        omega_df.loc[:, "delta-omega"] = (
                omega_df.loc[:, "freq"] - omega_df.loc[:, "mad-omega"]
        )
        omega_df.loc[:, "d-residual"] = self.d_residuals.iloc[-1, :].reset_index(
            drop=True
        )
        omega_df.loc[:, "bshx-residual"] = self.bsh_residuals.iloc[
                                           -1, 0: self.num_states
                                           ].reset_index(drop=True)
        omega_df.loc[:, "bshy-residual"] = self.bsh_residuals.iloc[
                                           -1, self.num_states::
                                           ].reset_index(drop=True)

        return ground_compare, omega_df


# input response_info json and returns a dict of response paramters
# and a list of dicts of numpy arrays holding response data


# Plotting definitions


def create_polar_table(mol, xc, basis_list, xx, database_dir):
    dalton_reader = Dalton(database_dir)
    ground_dalton, response_dalton = dalton_reader.get_frequency_result(
        mol, xc, "dipole", basis_list[0]
    )
    freq = response_dalton["frequencies"]
    g_data = {}
    xx_data = []
    for i in range(len(freq)):
        xx_data.append({})
    for basis in basis_list:
        ground_dalton, response_dalton = dalton_reader.get_frequency_result(
            mol, xc, "dipole", basis
        )
        for i in range(len(freq)):
            xx_data[i][basis] = response_dalton[xx].iloc[i]
        g_data[basis] = ground_dalton["totalEnergy"]
    g_df = pd.Series(g_data)
    g_df.name = "Total HF Energy"
    names = []
    for f in freq:
        raw_f = r"{}".format(str(f))
        # names.append(r'$$\alpha_{xx}('+raw_f+r')$$')
        names.append("a(" + "{:.3f}".format(f) + ")")
    # print(xx_data)
    r_dfs = []
    for i in range(len(freq)):
        r_dfs.append(pd.Series(xx_data[i]))
        r_dfs[i].name = names[i]
    dalton_df = pd.concat([g_df] + r_dfs, axis=1)

    moldata = ResponseCalc(mol, "hf", "dipole")
    mad_data_e = {}
    mad_data_r = {}
    mad_data_e["Total HF Energy"] = moldata.ground_e["e_tot"]

    for i in range(len(names)):
        mad_data_r[names[i]] = moldata.polar_data[xx].iloc[i]

    mad_data_e = pd.Series(mad_data_e)
    mad_data_r = pd.Series(mad_data_r)

    mad_data = pd.concat([mad_data_e, mad_data_r], axis=0)
    mad_data.name = "MRA"
    mad_data.key = ["MRA"]
    data = pd.concat([dalton_df.T, mad_data.T], axis=1)
    return data.T


def create_data(mol, xc, basis_list, database_dir):
    xx = ["xx", "yy", "zz"]
    data = []
    for x in xx:
        data.append(create_polar_table(mol, xc, basis_list, x, database_dir))
    average = (data[0] + data[1] + data[2]) / 3

    diff_data = average - average.loc["MRA"]
    diff_data = diff_data.drop(index="MRA")

    polar_diff = diff_data.drop("Total HF Energy", axis=1)

    average.name = "Average Polarizability"
    energy_diff = diff_data["Total HF Energy"]

    return average, diff_data, energy_diff, polar_diff


def polar_plot(mol, xc, basis_list, ax, database_dir):
    data, diff_data, energy_diff, polar_diff = create_data(mol, xc, basis_list, database_dir)
    num_freq = len(list(polar_diff.keys()))
    sns.set_theme(style="darkgrid")
    sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})

    btitle = basis_list[0].replace('D', 'X')

    polar_diff.iloc[:, :].plot(marker="v", linestyle="solid", markersize=12, linewidth=4,
                               colormap='magma', ax=ax, title=btitle)

    legend = []
    for i in range(num_freq):
        legend.append(r'$\omega_{}$'.format(i))

    ax.axhline(linewidth=2, ls="--", color="k")
    ax.legend(legend)
    # ax.set_xticks(rotation=20)
    yl = r" $\Delta\alpha=[\alpha($BASIS$) -\alpha($MRA$)]$"
    ax.set_ylabel(yl)
    ax.set_xlabel(ax.get_xlabel(), rotation=45)


def create_polar_diff_subplot(mol, xc, blist, dlist, database_dir):
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(25, 9), constrained_layout=True)
    title = 'Polarizability Convergence: ' + mol
    fig.suptitle(title)

    polar_plot(mol, xc, blist, ax[0], database_dir)
    polar_plot(mol, xc, dlist, ax[1], database_dir)

    btitle = blist[0].replace('D', 'X')
    save = mol + "-" + btitle
    if not os.path.exists("acs_mols"):
        os.mkdir("acs_mols")
    save = "acs_mols/" + save + ".svg"
    fig.savefig(save)
    return fig


def get_excited_mol_series(mol, basis):
    d = ExcitedData(mol, "hf")
    g, e = d.compare_dalton(basis)

    mad_series = e["mad-omega"]
    dal_series = e["freq"]
    md = []
    dd = []
    for m in range(mad_series.size):
        md.append("mra-{v}".format(v=m))
        dd.append("{b}-{v}".format(b=basis, v=m))
    mad_series.index = md
    dal_series.index = dd
    conv_series = pd.Series(d.converged)
    conv_series.index = ["Converged"]

    molseries = pd.concat([conv_series, mad_series, dal_series])
    return molseries


def create_excited_comparison_data(basis, excluded):
    data = {}
    for g in glob.glob("molecules/*.mol"):
        m = g.split("/")
        mol = m[1].split(".")[0]
        if mol not in excluded:
            data[mol] = get_excited_mol_series(mol, basis)
    excited_data = pd.DataFrame(data).round(4)
    return excited_data


def create_polar_mol_series(mol, basis):
    data = ResponseCalc(mol, "hf", "dipole")

    converged = data.converged
    freq = pd.Series(converged.keys())

    mra_keys = ["HF Energy"]
    diff_keys = [basis]
    conv_keys = []
    for f in range(freq.size):
        mra_keys.append("avg_{d}".format(d=f))
        diff_keys.append("diff_{d}".format(d=f))
        conv_keys.append("converged_{d}".format(d=f))

    xx = ["xx", "yy", "zz"]
    data = []
    for x in xx:
        data.append(create_polar_table(mol, "hf", [basis], x))
    average = (data[0] + data[1] + data[2]) / 3
    mra = average.loc["MRA"]
    basis_value = average.loc[basis]
    diff = mra - basis_value
    avg_diff = diff.mean()
    avg_diff = pd.Series(avg_diff)
    avg_diff.index = ["average diff"]
    mra.index = mra_keys
    diff.index = diff_keys
    converged.index = conv_keys
    new = pd.concat([freq, mra, pd.Series(avg_diff), converged], axis=0)

    return new


def polar_overview(basis, excluded):
    data = {}
    for g in glob.glob("molecules/*.mol"):
        m = g.split("/")
        mol = m[1].split(".")[0]
        # print(mol)
        if mol not in excluded:
            data[mol] = create_polar_mol_series(mol, basis)

    return pd.DataFrame(data)





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
