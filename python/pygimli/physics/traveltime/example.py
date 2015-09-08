from pygimli.physics import Refraction


ra = Refraction('example_topo.sgt')
print(ra)
ra.showData()
ra.showVA()
ra.makeMesh()
ra.showMesh()
ra.run(lam=300)
ra.showResult()
