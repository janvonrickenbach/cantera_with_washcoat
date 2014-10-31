#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as pl
import cantera as ca
import pickle as pi
import shutil as sh
import os 


gas=ca.import_phases("cantera_input.cti",["gas"])
gas = gas[0]
i =ca.Interface("cantera_input.cti","pt_surf",[gas])

inputs = {}

inputs["wc_thickness"] = 100
inputs["h"]=4000
inputs["h_temp"]=4000
inputs["porosity"]=0.35
inputs["lambda_solid"] = 10.
inputs["tortuosity"] = 8.0
inputs["atol"] = 1E-25
inputs["rtol"]   = 1E-6
inputs["d_p"] = 1E-8
inputs["nx"]=64
inputs["x_cells"]=20
inputs["L_reactor"] = 0.005
inputs["vel"] = 2.5
inputs["A_V"] = 1000.

inputs["mintemp"] = 380
inputs["maxtemp"] = 520
inputs["trate"]= 0.1

inputs["from_file"] = False
inputs["with_energy"] = False
inputs["maxiter"] = 10000
inputs["temperature"] = 300

inputs["istorf"] = 1
inputs["wc_geo_area"] = 100
inputs["t_fac"] = 1

inputs["gas"] = gas

inputs["mf"] = 3000

inputs["dt"] = (inputs["maxtemp"]-inputs["mintemp"])/inputs["trate"] *  2 
inputs["n_output"] = 10

inputs["area_to_volume"] = 100.0/ inputs["wc_thickness"]

input_file = pi.load(open("input.pic"))

for key,item in input_file.items():
   print("Changing", key, item)
   if not key in inputs:
      print("ERROR: Key not in inputs")
      exit()
   else:
      inputs[key] = item

inputs["dt"] = (inputs["maxtemp"]-inputs["mintemp"])/inputs["trate"] *  2 


inputs["wc_thickness"] = inputs["wc_thickness"]  * 1E-6
inputs["area_to_volume"] = inputs["wc_geo_area"]/ inputs["wc_thickness"]

t_fac = inputs["t_fac"]

mf = inputs["mf"]

del inputs["wc_geo_area"]
del inputs["t_fac"]
del inputs["mf"]


i.TPX =[300,1E5,[0,0,1,0]]
i.coverages =[0,0,1,0]

temperatures = [300,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,600,700,800,900,1000]

temperatures = np.hstack([temperatures, temperatures[::-1]])

gas.TPX = [300,1E5,[0.1,0.0,0.0,0.9-mf*1E-6,mf*1E-6]]
i.advance_coverages_test(100000)

inputs["maxtemp"] = -1
inputs["mintemp"] = -1
inputs["trate"] = -1

try:
   os.mkdir("result")
except:
   pass

for idx,temp in enumerate(temperatures):

   inputs["temperature"] = temp

   cov = 1
   if idx >= len(temperatures)/2:
      cov = 0

   gas.TPX = [300,1E5,[0.1,0.0,0,0.9-mf*1E-6,mf*1E-6]]

   inputs["dt"] = 0.01

   if temp < 600:
      inputs["dt"] = 10 * t_fac
   if temp < 500:
      inputs["dt"] = 50 * t_fac

   print(inputs)
   i.advance_coverages(**inputs)

   inputs["from_file"] = True

   for point in range(inputs["x_cells"]):
      sh.copy("grid_"+str(point)+"_"+str(int(inputs["dt"]))+".dat",os.path.join("result","grid_"+str(point)+"_"+str(int(temp))+"_"+str(cov))+".dat")
      try:
        sh.copy("grid_"+str(point)+"_"+str(int(inputs["dt"]))+".dat","grid_"+str(point)+"_0.dat")
      except:
         pass

del i
del gas



