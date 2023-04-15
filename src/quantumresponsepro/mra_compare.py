from src.quantumresponsepro.madness.madnessReader import FrequencyData
from src.quantumresponsepro.madness.madness_reader_v2 import ResponseCalc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.quantumresponsepro.Dalton.dalton import Dalton
import seaborn as sns

polar_keys = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']


def get_mra_polar_data(mols, xc, op, database):
    md = []
    N = 6  # round my data
    basis = 'MRA'
    for mol in mols:
        try:
            mad_r = ResponseCalc(mol, xc, op, database)
            polar_df = mad_r.polar_data[polar_keys]
            polar_df.index.name = 'frequencies'
            polar_df = polar_df.reset_index()
        except FileNotFoundError as f:
            mad_r = FrequencyData(mol, xc, op, database)
            polar_df = mad_r.polar_df[polar_keys]
            polar_df.index.name = 'frequencies'
            polar_df = polar_df.reset_index()
        md.append(column_polar_df(polar_df, mol, basis))
    df = pd.concat(md)
    df.loc[:, 'alpha'] = df.alpha.apply(lambda x: round(x, N - int(np.floor(np.log10(abs(x))))))
    return df


def partition_molecule_list(mol_list):
    Flist = []
    row2 = []
    rest = []
    seconds = ["Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]
    for mol_i in mol_list:
        print(mol_i)
        if "F" in mol_i and not any([e2 in mol_i for e2 in seconds]):
            Flist.append(mol_i)
        elif any([e2 in mol_i for e2 in seconds]):
            row2.append(mol_i)
        else:
            rest.append(mol_i)
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
            (xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2 + 6 * (xy ** 2 + yz ** 2 + zx ** 2))
        daf[freq] = da
    deltaA = pd.Series(daf)
    return deltaA


def column_polar_df(df, mol, basis):
    om = 'omega'
    alpha = 'alpha'
    j = 0
    r = []
    for name, values in df.iteritems():
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
    d = Dalton(database, False)
    try:
        mad_r = ResponseCalc(mol, xc, op, database)
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
    basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ',
                  'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ',
                  'd-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']
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
    d = Dalton(database, False)
    bd = []
    for mol in mols:
        for basis in basis_sets:
            try:
                ground, response = d.get_frequency_result(mol, 'hf', 'dipole', basis)
                basis_polar_df = response[polar_keys]
                bd.append(column_polar_df(basis_polar_df, mol, basis))
            except TypeError:
                print(mol, basis)
                pass
    return pd.concat(bd)


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
                   [iso_data.query('basis==@ba & molecule == @mol').alpha.reset_index(drop=True).rename(ba) for ba in
                    basis], axis=1)
    g1 = pd.concat(
        [o, mra_gamma] + [iso_data.query('basis==@ba & molecule == @mol').gamma.reset_index(drop=True).rename(ba) for ba
                          in basis], axis=1)
    a1 = a1.set_index('omega')
    g1 = g1.set_index('omega')
    return [[a1, g1], [percent_error(a1), absolute_error(g1)]]


def mol_iso_convergence(iso_data, mol, basis):
    width = 2 * 5
    height = 2 * 4
    g = plt.subplots(nrows=2, ncols=2, figsize=(width, height), constrained_layout=True, sharey=False)
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
    g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0).set_axis_labels("Valence", "Percent Error").set_titles(
        "{col_name}-cc-pV(C)nZ").tight_layout(w_pad=0)
    g.fig.suptitle(mol)
    return g


def create_iso_diff_df(iso_data):
    iso_diff = []
    blist = iso_data.basis.unique()[iso_data.basis.unique() != 'MRA']
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
    return iso_diff_df


def create_component_diff_df(a_data):
    multidex = ['ij', 'omega']
    a_diff = []
    blist = a_data.basis.unique()[a_data.basis.unique() != 'MRA']
    for mol in a_data.molecule.unique():
        mMRA = a_data.query('basis=="MRA" & molecule == @mol')
        a_mra = mMRA.set_index(multidex).alpha
        for basis in blist:
            mBASIS = a_data.query('basis==@basis & molecule == @mol')
            a_basis = mBASIS.set_index(multidex).alpha
            dab = (a_basis - a_mra) / a_mra * 100
            dB = mBASIS.copy().set_index(multidex)
            dB.alpha = dab
            a_diff.append(dB)
    ij_diff = pd.concat(a_diff)
    print(ij_diff)
    return ij_diff


def make_detailed_df(data):
    mols = list(data.molecule.unique())

    mols = [mol for mol in mols if mol is not None]

    row1, row2, flist = partition_molecule_list(mols)

    single = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ', ]

    double = ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ', ]
    single_polarized = ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ', ]
    double_polarized = ['d-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ', ]

    DZ = ['aug-cc-pVDZ', 'd-aug-cc-pVDZ', 'aug-cc-pCVDZ', 'd-aug-cc-pCVDZ', ]
    TZ = ['aug-cc-pVTZ', 'd-aug-cc-pVTZ', 'aug-cc-pCVTZ', 'd-aug-cc-pCVTZ', ]
    QZ = ['aug-cc-pVQZ', 'd-aug-cc-pVQZ', 'aug-cc-pCVQZ', 'd-aug-cc-pCVQZ', ]

    data = data.copy()
    data["augmentation"] = 'aug'
    data.loc[data["basis"].isin(double), "augmentation"] = 'd-aug'
    data.loc[data["basis"].isin(single), "augmentation"] = 'aug'
    data["augmentation"] = data["augmentation"].astype("category")

    data["polarization"] = 'V'
    data.loc[data["basis"].isin(single_polarized + double_polarized), "polarization"] = 'CV'
    data["polarization"] = data["polarization"].astype("category")
    data["mol_system"] = 'First-row'
    data.loc[data["molecule"].isin(row2), "mol_system"] = 'Second-row'
    data.loc[data["molecule"].isin(flist), "mol_system"] = 'Fluorine'
    data["mol_system"] = data["mol_system"].astype("category")
    data["valence"] = 'D'
    data.loc[data["basis"].isin(DZ), "valence"] = 'D'
    data.loc[data["basis"].isin(TZ), "valence"] = 'T'
    data.loc[data["basis"].isin(QZ), "valence"] = 'Q'
    data["valence"] = data["valence"].astype("category")
    data['valence'].cat.reorder_categories(['D', 'T', 'Q'], inplace=True)
    data['polarization'].cat.reorder_categories(['V', 'CV', ], inplace=True)

    data['Type'] = data[['augmentation', 'polarization']].apply(lambda x: "-cc-p".join(x) + 'nZ', axis=1)
    data["Type"] = data["Type"].astype("category")
    data['Type'] = data['Type'].cat.reorder_categories(
        ['aug-cc-pVnZ', 'aug-cc-pCVnZ', 'd-aug-cc-pVnZ', 'd-aug-cc-pCVnZ'])

    data["mol_size"] = 1.0
    data.loc[data["mol_system"] == 'First-Row', "mol_size"] = 1.0
    data.loc[data["mol_system"] == 'Fluorine', "mol_size"] = 10.0
    data.loc[data["mol_system"] == 'Second-Row', "mol_size"] = 100.0

    return data.reset_index()


def get_invariant_polar(azero):
    axx = azero.query('ij=="xx"').alpha.to_numpy()
    ayy = azero.query('ij=="yy"').alpha.to_numpy()
    azz = azero.query('ij=="zz"').alpha.to_numpy()
    axy = azero.query('ij=="xy"').alpha.to_numpy()
    ayz = azero.query('ij=="yz"').alpha.to_numpy()
    azx = azero.query('ij=="zx"').alpha.to_numpy()
    a0 = (axx + ayy + azz) / 3
    g0 = 1 * np.sqrt(
        (axx - ayy) ** 2 + (ayy - azz) ** 2 + (azz - axx) ** 2 + 6 * (axy ** 2 + ayz ** 2 + azx ** 2)) / np.sqrt(2)
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

    return iso_data.reset_index(drop=True)
