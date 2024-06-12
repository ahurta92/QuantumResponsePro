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

class MRAComparedBasisDF(pd.DataFrame):
    def __init__(self, polar_data, index, values: list, PercentError: bool, *args, **kwargs):
        # Use the special_parameter to modify the DataFrame or perform additional initialization
        basis_data = polar_data.query('basis!="MRA"').copy()
        basis_data = basis_data.set_index(index)

        for value in values:
            basis_data[f'{value}MRA'] = polar_data.query('basis=="MRA"').set_index(index)[
                value]
            if PercentError:
                basis_data[f'{value}E'] = ((basis_data[value] - basis_data[f'{value}MRA']) / basis_data[f'{value}MRA'] * 100)
            else:
                basis_data[f'{value}E'] = (basis_data[value] - basis_data[f'{value}MRA'])
        basis_data = basis_data.reset_index()
        # create a column of percent error in alpha
        basis_data = make_detailed_df(basis_data)
        super().__init__(basis_data, *args, **kwargs)

