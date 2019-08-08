# This is a base script for generating input decks for Branson using 3D blocks

import sys
import numpy

class block:
   def  __init__(self, size, center, mat, refinement=1):
        self.size = size
        self.center = center
        self.mat = mat
        self.refinement = refinement

   def get_x_planes(self):
       return [self.center[0]-self.size[0]/2.0, self.center[0]+self.size[0]/2.0]
   def get_y_planes(self):
       return [self.center[1]-self.size[1]/2.0, self.center[1]+self.size[1]/2.0]
   def get_z_planes(self):
       return [self.center[2]-self.size[2]/2.0, self.center[2]+self.size[2]/2.0]
   def point_in_block(self, point):
       x_planes = self.get_x_planes()
       y_planes = self.get_y_planes()
       z_planes = self.get_z_planes()
       if(point[0]>=x_planes[0] and point[0]<=x_planes[1] and
          point[1]>=y_planes[0] and point[1]<=y_planes[1] and
          point[2]>=z_planes[0] and point[2]<=z_planes[1]):
           return True
       else:
           return False


class branson_run_param:
    def __init__(self, t_stop, dt_start, photons, t_start=0.0,
            t_mult=1.0, dt_max=None, seed=14706, tilt=False,
            stratified_sampling=False, use_combing=True,
            output_frequency=1, write_silo=True,
            print_verbose=False, print_mesh_info=False,
            bc_right="VACUUM", bc_left="VACUUM", bc_up="VACUUM",
            bc_down="VACUUM", bc_top="VACUUM", bc_bottom="VACUUM", method="IMC"):
        self.method = method
        self.t_stop = str(t_stop)
        self.t_start = str(t_start)
        self.dt_start = str(dt_start)
        self.t_mult = str(t_mult)
        if dt_max is None:
            self.dt_max = str(dt_start)
        else:
            self.dt_max = str(dt_max)
        self.photons = str(photons)
        self.seed = str(seed)
        self.tilt = str(tilt).upper()
        self.stratified_sampling = str(stratified_sampling).upper()
        self.use_combing = str(use_combing).upper()
        self.output_frequency = str(output_frequency)
        self.write_silo = str(write_silo).upper()
        self.print_verbose = str(print_verbose).upper()
        self.print_mesh_info = str(print_mesh_info).upper()
        self.bc_right = bc_right
        self.bc_left = bc_left
        self.bc_up = bc_up
        self.bc_down = bc_down
        self.bc_top = bc_top
        self.bc_bottom = bc_bottom

def generate_input(block_list, run_param):
    # build a list of divisions
    x_division=[]
    y_division=[]
    z_division=[]
    for my_block in block_list:
        x_division+=(my_block.get_x_planes())
        y_division+=(my_block.get_y_planes())
        z_division+=(my_block.get_z_planes())
    
    # remove duplicates
    x_division=list(set(x_division))
    y_division=list(set(y_division))
    z_division=list(set(z_division))
    
    # sort list 
    x_division.sort()
    y_division.sort()
    z_division.sort()
    
    
    # set levels
    x_levels = [1]*len(x_division[:-1])
    y_levels = [1]*len(y_division[:-1])
    z_levels = [1]*len(z_division[:-1])
    for my_block in block_list:
        for i, (x_lo, x_hi) in enumerate(zip(x_division[:-1],x_division[1:])):
            if x_lo >= my_block.get_x_planes()[0] and x_hi <= my_block.get_x_planes()[1] and  my_block.refinement > x_levels[i]:
                x_levels[i] = my_block.refinement
        for i, (y_lo, y_hi) in enumerate(zip(y_division[:-1],y_division[1:])):
            if y_lo >= my_block.get_y_planes()[0] and y_hi <= my_block.get_y_planes()[1] and  my_block.refinement > y_levels[i]:
                y_levels[i] = my_block.refinement
        for i, (z_lo, z_hi) in enumerate(zip(z_division[:-1],z_division[1:])):
            if z_lo >= my_block.get_z_planes()[0] and z_hi <= my_block.get_z_planes()[1] and  my_block.refinement > z_levels[i]:
                z_levels[i] = my_block.refinement
    
    # build up regions ids
    nx = len(x_levels)
    ny = len(y_levels)
    nz = len(z_levels)
    region_ids = [[[ -1 for z_i in range(0,nz)] for y_i in range(0,ny)] for x_i in range(0,nx)]
    # loop over all blocks and paint in the materials
    for block in block_list:
        for x_i in range(0,nx):
            for y_i in range(0,ny):
                for z_i in range(0,nz):
                    # calculate the center point of the region
                    x = x_division[x_i]+0.5*(x_division[x_i+1]-x_division[x_i])
                    y = y_division[y_i]+0.5*(y_division[y_i+1]-y_division[y_i])
                    z = z_division[z_i]+0.5*(z_division[z_i+1]-z_division[z_i])
                    point = [x, y, z]
                    # does the point reside inside the block
                    if (block.point_in_block(point)):
                        region_ids[x_i][y_i][z_i] = block.mat["id"]
    
    
    # write output
    print("<prototype>")
    print("  <common>")
    print("    <method>"+run_param.method+"</method>")
    print("    <t_start>"+run_param.t_start+"</t_start>")
    print("    <t_stop>"+run_param.t_stop+"</t_stop>")
    print("    <dt_start>"+run_param.dt_start+"</dt_start>")
    print("    <t_mult>"+run_param.t_mult+"</t_mult>")
    print("    <dt_max>"+run_param.dt_max+"</dt_max>")
    print("    <photons>"+run_param.photons+"</photons>")
    print("    <seed>"+run_param.seed+"</seed>")
    print("    <tilt>"+run_param.tilt+"</tilt>")
    print("    <stratified_sampling>"+run_param.stratified_sampling+"</stratified_sampling>")
    print("    <use_combing>"+run_param.use_combing+"</use_combing>")
    print("    <output_frequency>"+run_param.output_frequency+"</output_frequency>")
    print("    <write_silo>"+run_param.write_silo+"</write_silo>")
    print("  </common>")
    print()
    
    print("  <debug_options>")
    print("    <print_verbose>"+run_param.print_verbose+"</print_verbose>")
    print("    <print_mesh_info>"+run_param.print_mesh_info+"</print_mesh_info>")
    print("  </debug_options>")
    print()
    
    # write boundary conditions                   
    print("  <boundary>")   
    print("    <bc_right>"+run_param.bc_right.strip()+"</bc_right>")   
    print("    <bc_left>"+run_param.bc_left.strip()+"</bc_left>")   
    print("    <bc_up>"+run_param.bc_up.strip()+"</bc_up>")   
    print("    <bc_down>"+run_param.bc_down.strip()+"</bc_down>")   
    print("    <bc_top>"+run_param.bc_top.strip()+"</bc_top>")   
    print("    <bc_bottom>"+run_param.bc_bottom.strip()+"</bc_bottom>")   
    print("  </boundary>")   
    print()
    
    # write spatial and region map
    print("  <spatial>")   
    for i, (x_lo, x_hi) in enumerate(zip(x_division[:-1],x_division[1:])):
        print("    <x_division>")   
        print("       <x_start>", x_lo, "</x_start>")   
        print("       <x_end>", x_hi, "</x_end>")   
        print("       <n_x_cells>", x_levels[i], "</n_x_cells>")   
        print("    </x_division>")   
        print()
    
    for i, (y_lo, y_hi) in enumerate(zip(y_division[:-1],y_division[1:])):
        print("    <y_division>")   
        print("       <y_start>", y_lo, "</y_start>")   
        print("       <y_end>", y_hi, "</y_end>")   
        print("       <n_y_cells>", y_levels[i], "</n_y_cells>")   
        print("    </y_division>")   
        print()
    
    for i, (z_lo, z_hi) in enumerate(zip(z_division[:-1],z_division[1:])):
        print("    <z_division>")   
        print("       <z_start>", z_lo, "</z_start>")   
        print("       <z_end>", z_hi, "</z_end>")   
        print("       <n_z_cells>", z_levels[i], "</n_z_cells>")   
        print("    </z_division>")   
        print()
    
    
    for x_i in range(0,nx):
        for y_i in range(0,ny):
            for z_i in range(0,nz):
                print("    <region_map>")   
                print("       <x_div_ID>", x_i, "</x_div_ID>")   
                print("       <y_div_ID>", y_i, "</y_div_ID>")   
                print("       <z_div_ID>", z_i, "</z_div_ID>")   
                print("       <region_ID>", region_ids[x_i][y_i][z_i], "</region_ID>")   
                print("    </region_map>")   
                print()
    
    print("  </spatial>")   
    
    
    print("  <regions>")   
    for block in block_list:
        mat = block.mat
        print("    <region>")   
        print("       <ID>", mat["id"], "</ID>")   
        print("       <density>", mat["density"], "</density>")   
        print("       <CV>", mat["cv"], "</CV>")   
        print("       <opacA>", mat["opacA"], "</opacA>")   
        print("       <opacB>", mat["opacB"], "</opacB>")   
        print("       <opacC>", mat["opacC"], "</opacC>")   
        print("       <opacS>", mat["opacS"], "</opacS>")   
        print("       <initial_T_e>", mat["initial_T_e"], "</initial_T_e>")   
        print("       <initial_T_r>", mat["initial_T_r"], "</initial_T_r>")   
        print("    </region>")   
        print()
    print("  </regions>")   
    print("</prototype>")
