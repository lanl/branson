# This is a script to write an XML geometry for our "cubanova" calculation

import sys
import numpy
from BlockIt import block
from BlockIt import branson_run_param
from BlockIt import generate_input

# generate branson run_param

run_param = branson_run_param(t_stop=1.0e-2, dt_start=1.0e-3, photons=5000, seed=14706)

# build a symmetric cubanova
scale_size = 10**16
scale_cv = 1.0e12
smear_shock = 1.1
domain_size = 4.0*10**16/scale_size
CSM_size = 2.0*10**16/scale_size
CSM_shock_size = 2.0*10**15.831/scale_size*smear_shock
cubanova_size = 2.0*10**15.829/scale_size

block_list = []
# build up materials
void = {}
void["id"] = 1
void["density"] = 1e-6
void["cv"] = 1e16
void["opacA"] = 1e-6
void["opacB"] = 0.0
void["opacC"] = 0.0
void["opacS"] = 0.0
void["initial_T_e"] = 8.6e-4
void["initial_T_r"] = 8.6e-4

unit_cube = block([domain_size]*3, [0.0]*3, void, 1)
block_list.append(unit_cube)

CSM_material= {}
CSM_material["id"] = 2
CSM_material["density"] = 1.0e-14
CSM_material["cv"] = 1.0e-3*scale_cv
CSM_material["opacA"] = 1.0e-4*CSM_material["density"]*scale_size
CSM_material["opacB"] = 0.0
CSM_material["opacC"] = 0.0
CSM_material["opacS"] = 0.0
CSM_material["initial_T_e"] = 8.6e-4
CSM_material["initial_T_r"] = 8.6e-4

CSM_dim=[CSM_size]*3

CSM = block(CSM_dim, [0.0]*3, CSM_material, 4)
block_list.append(CSM)

CSM_Shock_material= {}
CSM_Shock_material["id"] = 3
CSM_Shock_material["density"] = 1.0e-12
CSM_Shock_material["cv"] = 1.0e-3*scale_cv
CSM_Shock_material["opacA"] = 0.3*CSM_Shock_material["density"]*scale_size/smear_shock
CSM_Shock_material["opacB"] = 0.0
CSM_Shock_material["opacC"] = 0.0
CSM_Shock_material["opacS"] = 0.0
CSM_Shock_material["initial_T_e"] = 8.6e-2
CSM_Shock_material["initial_T_r"] = 8.6e-2

CSM_Shock_dim=[CSM_shock_size]*3

CSM_Shock = block(CSM_Shock_dim, [0.0]*3, CSM_Shock_material, 10)
block_list.append(CSM_Shock)


cubanova_material= {}
cubanova_material["id"] = 4
cubanova_material["density"] = 1.0e-14
cubanova_material["cv"] = 1.0e-3*scale_cv
cubanova_material["opacA"] = 0.3*cubanova_material["density"]*scale_size
cubanova_material["opacB"] = 0.0
cubanova_material["opacC"] = 0.0
cubanova_material["opacS"] = 0.0
cubanova_material["initial_T_e"] = 8.6e-4
cubanova_material["initial_T_r"] = 8.6e-4

cubanova_dim = [cubanova_size]*3

cubanova = block(cubanova_dim, [0.0]*3, cubanova_material, 8)
block_list.append(cubanova)

generate_input(block_list, run_param)

