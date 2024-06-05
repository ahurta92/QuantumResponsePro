
class MADMolecule:

    param_keys=['eprec','field','no_orient','psp_calc','pure_ae','symtol']

    params={
        'eprec': 1e-4,
        'field':[0.0,0.0,0.0],
        'no_orient': False,
        'psp_calc':False,
        'pure_ae': True,
        'symtol': -.01
    }
    geometry=[]
    symbols=[]

    def __init__(self,geometry,symbols,params={}):
        self.geometry=geometry
        self.symbols=symbols
        self.params.update(params)

    def __init__(self,file):

        with open(file,'r') as f:
            lines=f.readlines()

        self.geometry=[]
        self.symbols=[]
        for line in lines:
            if line[0]=='#':
                continue
            data=line.split()
            self.symbols.append(data[0])
            self.geometry.append([float(x) for x in data[1:]])
