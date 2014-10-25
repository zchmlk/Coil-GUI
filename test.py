# -*- coding: utf-8 -*-
"""
Created on Sat May 17 23:47:54 2014

@author: zchmlk
"""

from coil_sim_object import CoilSim 

print "Defining coilsim"
cs = CoilSim()

print "Reading spreadsheet for coil set"
# read the "standard" (RKT) Google spreadsheet coil set 'coil_set_spreadsheet'
cs.coil_set.read_cvs()
# this Google spreadsheet modifies the inner-most radius coils
#cs.coil_set.read_cvs(filename='coil_set_spreadsheet_2')

print "Calculating psi"
cs.calc_psi()

print "Plotting"
cs.plot()
