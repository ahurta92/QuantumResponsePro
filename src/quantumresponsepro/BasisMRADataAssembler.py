import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .dalton.daltonrunner import DaltonRunner
from .madness.madnessReader import FrequencyData
from .madness.madness_reader_v2 import MadnessResponse

polar_keys = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']


def get_basis_e_data(mols, basis_sets, xc, op, database):
    d = DaltonRunner(database, False)

    df = pd.DataFrame()
    for basis in basis_sets:
        basis_dict = {}
        basis_energy = {}
        for mol in mols:
            try:
                ground, response = d.get_frequency_result(mol, xc, op, basis)

                # print(ground)
                basis_energy[mol] = ground['totalEnergy']
                basis_dict[mol] = basis

            except TypeError:
                print(mol, basis)
                pass

        dd = [pd.Series(basis_dict, name='basis', dtype='category'),
              pd.Series(basis_energy, name='energy', dtype='float64')]
        df = pd.concat([df, pd.concat(dd, axis=1)], axis=0)

    return df


class MRAMolDataDoesNotExist(Exception):
    def __init__(self, message="MRA Data for molecule does not exist in database"):
        self.message = message
        super().__init__(self.message)


def get_mra_energy_data(mols, xc, op, database):
    basis = 'MRA'
    mra_e_dict = {}
    basis_dict = {}
    for mol in mols:
        try:
            mad_r = MadnessResponse(mol, xc, op, database)
            mra_e_dict[mol] = mad_r.ground_e['e_tot']
            basis_dict[mol] = 'MRA'
            # if file not found then try the old database

        except (FileNotFoundError , KeyError) as f:
            try:
                mad_r = FrequencyData(mol, xc, op, database)
                mra_e_dict[mol] = mad_r.ground_e['e_tot']
                basis_dict[mol] = 'MRA'
            except FileNotFoundError as f:
                print(f)
                print('Did not find the old or new mra data for {}'.format(mol))
                # raise a custom exception indicating that the molecule was not found to be
                pass


    dd = [pd.Series(basis_dict, name='basis'),
          pd.Series(mra_e_dict, name='energy')]
    mra_e = pd.concat(dd, axis=1)

    return mra_e


def get_energy_data(mols, xc, op, basis_list, database):
    mra = get_mra_energy_data(mols, xc, op, database)
    basis = get_basis_e_data(mols, basis_list, xc, op, database)
    df = pd.concat([mra, basis], axis=0)
    df.index.name = 'molecule'
    df.reset_index(inplace=True)
    return df


def get_energy_diff_data(mols, xc, op, basis_list, database):
    df = get_energy_data(mols, xc, op, basis_list, database)
    diff = pd.DataFrame()
    for basis in df.basis.unique():
        if basis != "MRA":
            dE = -(df.query('basis=="MRA"').set_index('molecule').energy - df.query(
                'basis ==@basis').set_index('molecule').energy)
            dE = dE.rename('error')
            print(dE)
            bS = pd.Series([basis for i in range(len(dE))], name='basis', index=dE.index)
            diff = pd.concat([diff, pd.concat([bS, dE], axis=1)])
    diff.reset_index(inplace=True)
    return diff


def get_mra_polar_data(mols, xc, op, database):
    md = []
    N = 6  # round my data
    basis = 'MRA'
    for mol in mols:
        try:
            mad_r = MadnessResponse(mol, xc, op, database)
            polar_df = mad_r.polar_data[polar_keys]
            polar_df.index.name = 'frequencies'
            polar_df = polar_df.reset_index()
        except (FileNotFoundError) as f:
            mad_r = FrequencyData(mol, xc, op, database)
            polar_df = mad_r.polar_df[polar_keys]
            polar_df.index.name = 'frequencies'
            polar_df = polar_df.reset_index()
        except (KeyError) as f:
            # did not find polar data for mol
            print("Didn't find polar data for {}".format(mol), f)
            pass
        md.append(column_polar_df(polar_df, mol, basis))
    df = pd.concat(md)
    df.loc[:, 'alpha'] = df.alpha.apply(lambda x: round(x, N - int(np.floor(np.log10(abs(x))))))
    return df


def get_mra_quad_data(mols, xc, op, database):
    df = pd.DataFrame()
    N = 6  # round my data
    basis = 'MRA'
    dfi = []
    for mol in mols:
        try:
            mad_r = MadnessResponse(mol, xc, op, database)
            beta_df = mad_r.quad_data

            dfi.append(beta_df)

        except (FileNotFoundError, KeyError) as f:
            print('did not find beta.json for {}'.format(mol))
            pass
    df = pd.concat(dfi, ignore_index=False, axis=0)

    return df


def partition_molecule_list(mol_list):
    Flist = []
    row2 = []
    rest = []
    seconds = ["Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]
    for mol_i in mol_list:
        try:
            if "F" in mol_i and not any([e2 in mol_i for e2 in seconds]):
                Flist.append(mol_i)
            elif any([e2 in mol_i for e2 in seconds]):
                row2.append(mol_i)
            else:
                rest.append(mol_i)
        except TypeError as f:
            print(f)
            pass

    return rest, row2, Flist


def isotropic_polar(polar_df: pd.DataFrame):
    frequencies = polar_df.index
    daf = {}
    for freq in frequencies:
        alpha_ij = polar_df[polar_keys].loc[freq].to_numpy().reshape((3, 3))
        xx = alpha_ij[0, 0]
        yy = alpha_ij[1, 1]
        zz = alpha_ij[2, 2]
        iso_alpha = (xx + yy + zz) / 3
        print(iso_alpha)
        daf[freq] = iso_alpha
    isoA = pd.Series(daf)
    return isoA


def anisotropic_polar(polar_df):
    frequencies = polar_df.index
    daf = {}
    for freq in frequencies:
        alpha_ij = polar_df[polar_keys].loc[freq].to_numpy().reshape((3, 3))
        xx = alpha_ij[0, 0]
        yy = alpha_ij[1, 1]
        zz = alpha_ij[2, 2]
        xy = alpha_ij[0, 1]
        yz = alpha_ij[1, 2]
        zx = alpha_ij[2, 0]
        da = 1 / np.sqrt(2) * np.sqrt(
            (xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2 + 6 * (
                    xy ** 2 + yz ** 2 + zx ** 2))
        daf[freq] = da
    deltaA = pd.Series(daf)
    return deltaA


def column_polar_df(df, mol, basis):
    om = 'omega'
    alpha = 'alpha'
    j = 0
    r = []
    for name, values in df.items():
        if name != "frequencies":
            mS = pd.Series([mol for i in range(values.size)], name='molecule', dtype='category')
            bS = pd.Series([basis for i in range(values.size)], name='basis', dtype='category')
            ij = pd.Series([name for i in range(values.size)], name='ij', dtype='category')
            fs = pd.Series([i for i in range(values.size)], name=om)
            alpha_ij = values
            alpha_ij.name = alpha
            r.append(pd.concat([mS, bS, ij, fs, alpha_ij], axis=1))
            j += 1
    c_df = pd.concat(r)
    return c_df


def compare_database(mol, xc, op, database, basis_sets):
    d = DaltonRunner(database, False)
    try:
        mad_r = MadnessResponse(mol, xc, op, database)
        polar_df = mad_r.polar_data[polar_keys]
        polar_df.index.name = 'frequencies'
        polar_df = polar_df.reset_index()
    except FileNotFoundError as f:
        mad_r = FrequencyData(mol, xc, op, database)
        polar_df = mad_r.polar_df[polar_keys]
        polar_df.index.name = 'frequencies'
        polar_df = polar_df.reset_index()
    mad_polar_df = polar_df[polar_keys]
    mad_DA = anisotropic_polar(mad_polar_df)
    mad_DA.name = 'madness'
    mad_IA = isotropic_polar(mad_polar_df)
    mad_IA.name = 'madness'
    bDA = {}
    bIA = {}

    for basis in basis_sets:
        ground, response = d.get_frequency_result(mol, 'hf', 'dipole', basis)
        basis_polar_df = response[polar_keys]
        bDA[basis] = anisotropic_polar(basis_polar_df)
        bIA[basis] = isotropic_polar(basis_polar_df)
    basisIA = pd.DataFrame(bIA)
    basisDA = pd.DataFrame(bDA)
    delta_a = pd.concat([mad_DA, basisDA], axis=1)
    iso_a = pd.concat([mad_IA, basisIA], axis=1)
    return iso_a, delta_a


def get_basis_polar_data(mols, basis_sets, xc, op, database):
    print(database)
    d = DaltonRunner(database, False)
    bd = []
    for mol in mols:
        for basis in basis_sets:
            try:
                ground, response = d.get_frequency_result(mol, xc, op, basis)
                basis_polar_df = response[polar_keys]
                bd.append(column_polar_df(basis_polar_df, mol, basis))
            except TypeError:
                print(mol, basis)
                pass
    print(bd)
    return pd.concat(bd)


def get_basis_quad_data(mols, basis_sets, xc, op, database):
    d = DaltonRunner(database, False)
    bd = []
    for mol in mols:
        for basis in basis_sets:
            try:

                quad_data = d.get_quad_json(mol, xc, op, basis)
                bd.append(quad_data)
            except TypeError:
                print(mol, basis)
                pass
            except AttributeError as e:
                print(e)
                pass
    print(bd)

    b_data = pd.concat(bd)
    # here we need to process the Beta Value coluns.  First we change the name to Beta
    b_data.rename(columns={'Beta Value': 'Beta'}, inplace=True)

    # also drop the dash from A-freq and B-freq and C-freq

    # add a basis column
    return b_data


def percent_error(df):
    return df.subtract(df['MRA'].values, axis=0).div(df['MRA'].values, axis=0) * 100


def absolute_error(df):
    return df.subtract(df['MRA'].values, axis=0)


def get_frequency_compare(iso_data, mol, basis):
    mra_mol = iso_data.query('basis==@MRA & molecule == @mol')
    mra_alpha = mra_mol.alpha.reset_index(drop=True).rename('MRA')
    mra_gamma = mra_mol.gamma.reset_index(drop=True).rename('MRA')
    b_mol = iso_data.query('basis==@basis & molecule == @mol')
    o = mra_mol.omega.reset_index(drop=True)
    a1 = pd.concat([o, mra_alpha] + \
                   [iso_data.query('basis==@ba & molecule == @mol').alpha.reset_index(
                       drop=True).rename(ba) for ba in
                    basis], axis=1)
    g1 = pd.concat(
        [o, mra_gamma] + [
            iso_data.query('basis==@ba & molecule == @mol').gamma.reset_index(drop=True).rename(
                ba)
            for ba
            in basis], axis=1)
    a1 = a1.set_index('omega')
    g1 = g1.set_index('omega')
    return [[a1, g1], [percent_error(a1), absolute_error(g1)]]


def mol_iso_convergence(iso_data, mol, basis):
    width = 2 * 5
    height = 2 * 4
    g = plt.subplots(nrows=2, ncols=2, figsize=(width, height), constrained_layout=True,
                     sharey=False)
    fig = g[0]
    ax = g[1]
    data = get_frequency_compare(iso_data, mol, basis)
    title = [r'$\alpha(\omega)$', r'$\gamma(\omega)$']

    for i, axi in enumerate(ax):
        for j, axij in enumerate(axi):
            dij = data[i][j]
            axij.set_title(title[j])
            dij.plot(ax=axij, ls='dashdot')
            axij.set_xlabel(r'$\omega_i$')
    return g


def mol_component_convergence(mol, ij_diff_detailed):
    mol_data = ij_diff_detailed.query(' molecule == @mol & ij==@ij')
    g = sns.relplot(data=mol_data,
                    x=mol_data.valence,
                    y='alpha',
                    hue='ij',
                    col='augmentation',
                    style='pC',
                    kind='line',
                    markers=True,
                    facet_kws={'sharey': False},
                    dashes=True,
                    )
    for i, ax in enumerate(g.axes):
        for j, axi in enumerate(ax):
            axi.tick_params(axis='x', rotation=0)
            axi.grid(which="both")
            axi.tick_params(which="both", top="on", left="on", right="on", bottom="on", )
            axi.minorticks_on()
    g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0).set_axis_labels("Valence",
                                                                                "Percent Error").set_titles(
        "{col_name}-cc-pV(C)nZ").tight_layout(w_pad=0)
    g.fig.suptitle(mol)
    return g


def create_iso_diff_df(iso_data):
    iso_diff = []
    blist = iso_data.copy().query("basis != 'MRA'").basis.unique()

    print("create iso diff", blist)
    for mol in iso_data.molecule.unique():
        mMRA = iso_data.query('basis=="MRA" & molecule == @mol')
        a_mra = mMRA.set_index('omega').alpha
        g_mra = mMRA.set_index('omega').gamma
        for basis in blist:
            mBASIS = iso_data.query('basis==@basis & molecule == @mol')
            a_basis = mBASIS.set_index('omega').alpha
            g_basis = mBASIS.set_index('omega').gamma
            dab = (a_basis - a_mra) / a_mra * 100
            if g_mra.mean() > 1:
                gab = (g_basis - g_mra) / g_mra * 100
            else:
                gab = (g_basis - g_mra) * 100

            dB = mBASIS.copy().set_index('omega')
            dB.gamma = gab.reset_index(drop=True)
            dB.alpha = dab.reset_index(drop=True)
            iso_diff.append(dB)
    iso_diff_df = pd.concat(iso_diff)
    # Replace infinite updated data with nan
    iso_diff_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    # iso_diff_df.dropna(inplace=True)
    iso_diff_df.reset_index(inplace=True)
    return iso_diff_df


def create_component_diff_df(a_data):
    multidex = ['ij', 'omega']
    f_ij = []
    for mol in a_data.molecule.unique():
        # get the mra data
        mMRA = a_data.query('basis=="MRA" & molecule == @mol')
        a_mra = mMRA.set_index(multidex).alpha

        # get the basis data
        basis_data = a_data.query('molecule==@mol & basis != "MRA"')
        b_mol_data = basis_data.set_index(multidex).alpha
        rep_mol_mra = pd.concat([a_mra for i in range(len(basis_data.basis.unique()))])
        # set rep_mol_mra index to match b_mol_data
        rep_mol_mra.index = b_mol_data.index

        # figure out which how many omega rep_mol_mra has
        # set b_mol_data to query the same omegas that rep_mol_mra has
        b_mol_data = b_mol_data.loc[rep_mol_mra.index]

        diff_data = pd.concat([rep_mol_mra, b_mol_data], axis=1).diff(axis=1).iloc[:, 1]
        bcol = basis_data.set_index(multidex).basis
        mcol = basis_data.set_index(multidex).molecule

        ij_diff = pd.concat([mcol, bcol, diff_data], axis=1).query('basis!="MRA"')
        f_ij.append(ij_diff)

    ij_diff = pd.concat(f_ij)
    ij_diff.reset_index(inplace=True)
    return ij_diff


def get_basis_mol_eigen(df: pd.DataFrame, mol, basis):
    ij = ['xx', 'yy', 'zz']
    ij = pd.Series(ij, name='ij')
    mol_s = pd.Series([mol for i in range(3)], name='molecule')
    b_s = pd.Series([basis for i in range(3)], name='basis')

    f_data = pd.DataFrame()
    for o in df.omega.unique():
        try:
            emra, vmra = np.linalg.eigh(np.array(df.query(
                'molecule==@mol & basis==@basis& omega==@o').alpha).reshape(3, 3))
            om = pd.Series([o for i in range(3)], name='omega')
            emra = pd.Series(emra, name='alpha')
            f_data = pd.concat([f_data, pd.concat([mol_s, b_s, ij, om, emra], axis=1)])
        except ValueError as v:
            print(v, mol, basis)
    return f_data


def get_ij_eigen(df: pd.DataFrame):
    all_eigen = pd.DataFrame()
    for mol in df.molecule.unique():
        for basis in df.basis.unique():
            all_eigen = pd.concat([all_eigen, get_basis_mol_eigen(df, mol, basis)])
    return all_eigen


def make_detailed_df(data):
    mols = list(data.molecule.unique())

    mols = [mol for mol in mols if mol is not None]

    row1, row2, flist = partition_molecule_list(mols)

    zero = ['cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pCVDZ', 'cc-pCVTZ', 'cc-pCVQZ', 'cc-pV5Z',
            'cc-pV6Z']

    single = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pCVDZ', 'aug-cc-pCVTZ',
              'aug-cc-pCVQZ', 'aug-cc-pV5Z', 'aug-cc-pV6Z']

    double = ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pCVDZ',
              'd-aug-cc-pCVTZ',
              'd-aug-cc-pCVQZ', 'd-aug-cc-pV5Z', 'd-aug-cc-pV6Z']
    zero_polarized = ['cc-pCVDZ', 'cc-pCVTZ', 'cc-pCVQZ', ]
    single_polarized = ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ', ]
    double_polarized = ['d-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ', ]

    DZ = ['cc-pVDZ', 'aug-cc-pVDZ', 'd-aug-cc-pVDZ', 'aug-cc-pCVDZ', 'd-aug-cc-pCVDZ', ]
    TZ = ['cc-pVTZ', 'aug-cc-pVTZ', 'd-aug-cc-pVTZ', 'aug-cc-pCVTZ', 'd-aug-cc-pCVTZ', ]
    QZ = ['cc-pVQZ', 'aug-cc-pVQZ', 'd-aug-cc-pVQZ', 'aug-cc-pCVQZ', 'd-aug-cc-pCVQZ', ]
    FZ = ['cc-pV5Z', 'aug-cc-pV5Z', 'd-aug-cc-pV5Z']
    SZ = ['cc-pV6Z', 'aug-cc-pV6Z', 'd-aug-cc-pV6Z']

    data = data.copy()
    data["augmentation"] = ''
    data.loc[data["basis"].isin(double), "augmentation"] = 'd-aug'
    data.loc[data["basis"].isin(single), "augmentation"] = 'aug'
    data["augmentation"] = data["augmentation"].astype("category")

    data["polarization"] = 'V'
    data.loc[data["basis"].isin(zero_polarized + single_polarized + double_polarized),
    "polarization"] = 'CV'
    data["polarization"] = data["polarization"].astype("category")
    data["mol_system"] = 'First-row'
    data.loc[data["molecule"].isin(row2), "mol_system"] = 'Second-row'
    data.loc[data["molecule"].isin(flist), "mol_system"] = 'Fluorine'
    data["mol_system"] = data["mol_system"].astype("category")
    data["valence"] = 'D'
    data.loc[data["basis"].isin(DZ), "valence"] = 'D'
    data.loc[data["basis"].isin(TZ), "valence"] = 'T'
    data.loc[data["basis"].isin(QZ), "valence"] = 'Q'
    data.loc[data["basis"].isin(FZ), "valence"] = '5'
    data.loc[data["basis"].isin(SZ), "valence"] = '6'

    data["valence"] = data["valence"].astype("category")
    try:
        valence = list(data["valence"].unique())
        print(valence)
        data["valence"] = data['valence'].cat.reorder_categories(valence)
    except ValueError:
        data["valence"] = data['valence'].cat.reorder_categories(['D', 'T', 'Q', ])

    try:
        data['polarization'].cat.reorder_categories(['V', 'CV', ])
    except ValueError as e:
        print(e)

    data['Type'] = data[['augmentation', 'polarization']].apply(
        lambda x: "-cc-p".join(x) + 'nZ',
        axis=1)
    data["Type"] = data["Type"].astype("category")

    # rename type  -cc-pVnZ and -cc-pCVnZ to cc-pVnZ and cc-pCVnZ
    data['Type'] = data['Type'].cat.rename_categories(
        {'-cc-pVnZ': 'cc-pVnZ', '-cc-pCVnZ': 'cc-pCVnZ', })

    possible_types = ['cc-pVnZ', 'cc-pCVnZ', 'aug-cc-pVnZ', 'aug-cc-pCVnZ', 'd-aug-cc-pVnZ',
                      'd-aug-cc-pCVnZ']
    # drop if not in the actual types
    actual_types = list(data['Type'].unique())
    # drop from possible types if not in actual types
    possible_types = [t for t in possible_types if t in actual_types]
    # now reorder the categories to match the actual types
    data['Type'] = data['Type'].cat.reorder_categories(possible_types, )

    return data


def get_invariant_polar(azero):
    axx = azero.query('ij=="xx"').alpha.to_numpy()
    ayy = azero.query('ij=="yy"').alpha.to_numpy()
    azz = azero.query('ij=="zz"').alpha.to_numpy()
    axy = azero.query('ij=="xy"').alpha.to_numpy()
    ayz = azero.query('ij=="yz"').alpha.to_numpy()
    azx = azero.query('ij=="zx"').alpha.to_numpy()
    a0 = (axx + ayy + azz) / 3
    g0 = 1 * np.sqrt(
        (axx - ayy) ** 2 + (ayy - azz) ** 2 + (azz - axx) ** 2 + 6 * (
                axy ** 2 + ayz ** 2 + azx ** 2)) / np.sqrt(2)
    return a0, g0


def get_iso_data(a_data):
    iso = []
    for mol in a_data.molecule.unique():

        m_data = a_data.query('molecule==@mol')
        for b in a_data.basis.unique():
            mbdata = m_data.query('basis==@b')
            for fi in a_data.omega.unique():
                azero = mbdata.query('omega == @fi')
                a0, g0 = get_invariant_polar(azero)
                iso.append(pd.concat(
                    [
                        pd.Series(mol, dtype='category', name='molecule'),
                        pd.Series(b, dtype='category', name='basis'),
                        pd.Series(fi, dtype='category', name='omega'),
                        pd.Series(a0, dtype=np.float64, name='alpha'),
                        pd.Series(g0, dtype=np.float64, name='gamma'),
                    ], axis=1))
    iso_data = pd.concat(iso)
    return iso_data
