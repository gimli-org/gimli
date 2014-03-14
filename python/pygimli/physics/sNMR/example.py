#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Example for MRS inversion """

from .mrs import MRS ## change to pygimli.physics.mrs

if __name__ == "__main__":

    datafile = 'example.mrsi' # MRSmatlab inversion (data+kernel) file
    mrs = MRS(datafile) # initialize and read file
    print(mrs) # displays some parameters

    mrs.run(3,uncertainty=True) # 3 layers
    thk, wc, t2 = mrs.result()
    name = datafile.rstrip('.mrsi')
    mrs.saveResult(name+'.result')
    mrs.showResultAndFit(save=name + '.pdf', show=True)