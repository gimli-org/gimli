import pandas as pd
import pygimli.meshtools as mt
import pybert as pb

### Define all variables ###
filename = 'ERT_pos_and_scheme.xlsx'
# Including elec in geometry surface in Gmsh is complicated (Point{pnt_id} In Surface{surf_id};).
# It gives errors, therefore they are put at some depth below the surface
# and then included in the volume of the dike (Point{pnt_id} In Volume{vol_id};).
elec_depth = 0.02               # elec depth [m]


### Write electrode positions and ERT scheme in BERT .dat file ###
# Load electrode positions from Excel
pos = pd.read_excel(filename, sheet_name='elec_pos')
pos['z'] = pos['z'] - elec_depth

# Write to BERT/pyGIMLi unified data format .dat file
ne = len(pos)                   # number of electrodes
scheme = pd.read_excel(filename, sheet_name='ERT_scheme')
BERT_dat = open('BERT.dat', 'w')
BERT_dat.write(f'{ne}\n' + '# x y z\n')
for i, xyz in pos.iterrows():
    BERT_dat.write('%.3f %.3f %.3f\n' % (xyz['x'], xyz['y'], xyz['z']))
BERT_dat.write(f'{scheme.shape[0]}\n')
BERT_dat.write('# a b m n\n')
for i, abmn in scheme.iterrows():
    BERT_dat.write('%d %d %d %d\n' % (abmn['A'], abmn['B'], abmn['M'], abmn['N']))
BERT_dat.write('0')
BERT_dat.close()


### Read Gmsh .msh and model ERT ####
msh_nm = 'Dike_ERT_mod'
mesh_mod = mt.readGmsh(msh_nm + ".msh", verbose=True)

# Read electrode positions and ERT scheme from .dat file
BERT = pb.importData('BERT.dat')
# Initiate pyBERT ERT object
ert = pb.ERTManager()
# Model ERT
rhomap = [[1, 10.], [2, 10.], [3, 300.]]
data_mod = ert.simulate(mesh_mod, res=rhomap, scheme=BERT)
data_mod.save('BERT_mod.dat')





