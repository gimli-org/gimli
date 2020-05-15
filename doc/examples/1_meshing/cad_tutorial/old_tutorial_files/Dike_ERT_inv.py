import pygimli.meshtools as mt
import pybert as pb

#Read Gmsh .msh
msh_nm = 'Dike_ERT_inv'
mesh_inv = mt.readGmsh(msh_nm + ".msh", verbose=True)
# mesh_inv.save(msh_nm + '.bms')        # Save mesh in BERT .bms format

# Initiate pyBERT ERT object
ert = pb.ERTManager()
# Import the modeled data
data = pb.importData('BERT_mod.dat')
data.markInvalid(data('rhoa') < 4)
# Invert modeled data. The inversion has the following options
# Note that command line BERT has a larger number of options.
# **kwargs
#  * lam : float [20]               regularization parameter
#  * zWeight : float [0.7]          relative vertical weight
#  * maxIter : int [20]             maximum iteration number
#  * robustData : bool [False]      robust data reweighting using an L1 scheme (IRLS reweighting)
#  * blockyModel : bool [False]     blocky model constraint using L1 reweighting roughness vector
#  * startModel : array-like        starting model
#  * startModelIsReference : bool [False]   startmodel is the reference model for the inversion
ert.invert(data, mesh=mesh_inv, lam=15)

# Write inversion outcome to .vtk to visualize in paraview
paradomain = ert.paraDomain
paradomain.addData('Resistivity [Ohm.m]', ert.resistivity)
paradomain.addData('log10(sensitivity)', ert.coverageDC())
paradomain.exportVTK('invResult.vtk')
