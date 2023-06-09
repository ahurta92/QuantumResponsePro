from .BasisMRADataAssembler import *
from pathlib import Path
import glob
import json


def get_polar_df(molecules, xc, op, database, basis):
    mra_data = get_mra_polar_data(molecules, xc, op, database)
    b_data = get_basis_polar_data(molecules, basis, xc, op, database)
    a_data = pd.concat([mra_data, b_data])
    a_data = a_data.reset_index(drop=True)
    return a_data


def save_compressed_polar(df, database_dir, data_file):
    compress_data_dir = database_dir + '/compress_data'
    file_name = compress_data_dir + "/" + data_file
    df.to_feather(file_name)


single = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pV5Z', 'aug-cc-pV6Z']
double = ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pV5Z', 'd-aug-cc-pV6Z']
single_polarized = ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ']
double_polarized = ['d-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']

all_basis_sets = single + double + single_polarized + double_polarized


class BasisMRADataCollection:
    """
    ResponseDataBundle: A class designed to manage and organize a collection of DataFrames related to response properties for molecules. This class simplifies the process of comparing and analyzing data from various sources, such as MADNESS and Dalton quantum chemistry packages, by consolidating the response properties information into a single, easy-to-use structure.
    """

    def __init__(self, data_dir, xc='hf', op='dipole', basis_sets=all_basis_sets):
        self.data_dir = data_dir
        self.xc = xc
        self.op = op
        self.molecules = []
        self.basis_sets = basis_sets

        feather_data = data_dir.joinpath('feather_data')
        if not feather_data.is_dir():
            feather_data.mkdir()

        all_data_path = feather_data.joinpath('all_polar_data.feather')

        if all_data_path.is_file():
            self.all_polar_data = pd.read_feather(all_data_path)
            self.molecules = list(self.all_polar_data.molecule.unique())
        else:
            for g in glob.glob('*.mol', root_dir=self.data_dir.joinpath('molecules')):
                m = g.split('/')
                mol = m[-1].split('.')[0]
                self.molecules.append(mol)
            self.all_polar_data = get_polar_df(self.molecules, xc, op, self.data_dir,
                                               self.basis_sets)
            self.all_polar_data.to_feather(all_data_path)

        iso_polar_data_path = feather_data.joinpath('iso_polar_data.feather')

        if iso_polar_data_path.is_file():
            self.iso_data = pd.read_feather(iso_polar_data_path)
        else:
            self.iso_data = get_iso_data(self.all_polar_data)
            self.iso_data.reset_index().to_feather(iso_polar_data_path)
            self.iso_data[self.iso_data['gamma'].abs() < 1e-3] = 0

        iso_diff_path = feather_data.joinpath('iso_diff_data.feather')

        if iso_diff_path.is_file():
            self.iso_diff_data = pd.read_feather(iso_diff_path)
        else:
            self.iso_diff_data = create_iso_diff_df(self.iso_data)
            self.iso_diff_data.to_feather(iso_diff_path)

        ij_df_path = feather_data.joinpath('ij_diff_data.feather')
        if ij_df_path.is_file():
            self.ij_diff = pd.read_feather(ij_df_path)
            self.ij_diff.reset_index(inplace=True)
        else:
            print('Creating ij_diff')
            self.ij_diff = create_component_diff_df(self.all_polar_data)
            self.ij_diff.reset_index().to_feather(ij_df_path)

        energy_df_path = feather_data.joinpath('energy_data.feather')

        if energy_df_path.is_file():
            self.energy_df = pd.read_feather(energy_df_path)
        else:
            self.energy_df = get_energy_data(self.molecules, self.xc, self.op, basis_sets,
                                             self.data_dir)
            self.energy_df.reset_index(inplace=True)
            self.energy_df.to_feather(energy_df_path)
            pass

        alpha_eigen_path = feather_data.joinpath('alpha_eigen_data.feather')
        if alpha_eigen_path.is_file():
            self.alpha_eigen = pd.read_feather(alpha_eigen_path)
        else:
            self.alpha_eigen = get_ij_eigen(self.all_polar_data)
            self.alpha_eigen.reset_index().to_feather(alpha_eigen_path)
        eigen_diff_path = feather_data.joinpath('eigen_diff_data.feather')
        if eigen_diff_path.is_file():
            self.eigen_diff = pd.read_feather(eigen_diff_path)
        else:
            self.eigen_diff = create_component_diff_df(self.alpha_eigen)
            self.eigen_diff.reset_index().to_feather(eigen_diff_path)

        energy_diff_path = feather_data.joinpath('energy_diff.feather')
        if energy_diff_path.is_file():
            self.energy_diff = pd.read_feather(energy_diff_path)
            # self.energy_diff.reset_index(inplace=True)

        else:
            self.energy_diff = get_energy_diff_data(self.molecules, self.xc, self.op, basis_sets,
                                                    self.data_dir)
            self.energy_diff.to_feather(energy_diff_path)
            pass

        try:

            self.detailed_iso_diff = make_detailed_df(self.iso_diff_data)
            self.detailed_ij_diff = make_detailed_df(self.ij_diff)
            self.detailed_energy_diff = make_detailed_df(self.energy_diff)
            self.detailed_eigen_diff = make_detailed_df(self.eigen_diff)
            self.detailed_alpha_eigen = make_detailed_df(self.alpha_eigen)
        except Exception as e:
            print(e)
            pass

    def __report_convergence(self):
        converged = []
        not_converged = []
        not_found = []
        type_error = []
        json_error = []
        for mol in self.molecules:
            try:
                check_mol = FrequencyData(mol, self.xc, self.op, self.database_dir)
                if check_mol.converged.all() and check_mol.converged.sum() == self.num_freq:
                    converged.append(mol)
                else:
                    not_converged.append(mol)
            except FileNotFoundError as f:
                not_found.append(mol)
            except TypeError as f:
                type_error.append(mol)
            except json.decoder.JSONDecodeError as j:
                json_error.append(mol)

        num_c = len(converged)
        num_n = len(not_converged)
        num_nf = len(not_found)
        num_json_e = len(json_error)
        num_type_e = len(type_error)
        total = num_c + num_n + num_nf + num_json_e + num_type_e
        not_converged = []
        part_converged = []
        if True:
            for mol in not_converged:
                check = FrequencyData(mol, self.xc, self.op, self.database_dir)
                if check.converged.any():
                    # print(mol,'\n',check.converged)
                    part_converged.append(mol)
                else:
                    not_converged.append(mol)
        num_not_converged = len(not_converged)
        num_part_converged = len(part_converged)
        print("converged : ", num_c)
        return converged, part_converged, not_converged, not_found, type_error, json_error


# selects outliers by change in percent error from the static to the
def select_basis_outliers(data, basis, thresh):
    om = [0, 8]
    test = data.query('omega.isin(@om) & basis==@basis')

    ma = {}
    for mol in test.molecule.unique():
        dmol = test.query('molecule==@mol')
        a0 = dmol.query('omega==0').alpha.iloc[0]
        a8 = dmol.query('omega==8').alpha.iloc[0]
        ma[mol] = (a0 - a8)
    ma = pd.Series(ma)
    out_mols = ma[ma.abs() > thresh]
    out_mols = set(out_mols.index)
