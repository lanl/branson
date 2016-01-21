import os
import sys
import subprocess
import re
import lxml.etree as ET
import numpy as np

################################################################################
def make_input_file(base_file, buffer_size, n_photons, n_cells):
  r_tree = ET.parse(base_file)
  r_root = r_tree.getroot()
  # set photons
  r_common = r_root.find("common")
  e_photons = r_common.find("photons")
  e_photons.text = "{0}".format(n_photons)
  # set buffer size
  r_common.find("particle_message_size").text = "{0}".format(buffer_size)
  # set cells per dimension
  e_spatial = r_root.find("spatial")
  e_spatial.find("n_x_cells").text = "{0}".format(n_cells)
  e_spatial.find("n_y_cells").text = "{0}".format(n_cells)
  e_spatial.find("n_z_cells").text = "{0}".format(n_cells)
  
  new_filename = "temp_input_n_{0}_cell_{1}.xml".format(n_photons, n_cells)
  r_tree.write(new_filename, pretty_print=True )
  return new_filename
################################################################################

if (len(sys.argv) != 2):
  print("usage: {0} <basic_input_file_name>".format(sys.argv[0]))
  sys.exit();

base_filename = sys.argv[1]

path_to_exe = "/net/scratch1/along/branson/build"
exe_name = "BRANSON"

time_r = re.compile('runtime: (.*?) $')

buffer_size = 2000

c_per_dim = [ 40, 100]
np_list = [1000000, 4000000, 16000000, 64000000, 256000000]
proc_list = [1, 2, 4, 8, 16, 32]
samples = 3

for np in np_list:
  for n_cells in c_per_dim:
    for p in proc_list:
      times = []
      temp_input_file = make_input_file(base_filename, buffer_size, np, n_cells)
      for s in range(samples):
        time = 0.0
        output_file = "temp_output.txt"
        subprocess.call(["mpirun -np {0} {1}/{2} {3} >> {4}".format( \
          p, path_to_exe, exe_name, temp_input_file, output_file) ], shell=True)
        f_out = open(output_file)
        for line in f_out:
          if (time_r.search(line)):
            time = float(time_r.findall(line)[0])
        times.append(time)
      os.remove("{0}".format(temp_input_file))
      # calculate average runtime and standard deviation
      runtime = np.average(times)
      stdev_runtime = np.std(times)
      print("{0} {1} {2} {3} {4}".format(p, c_per_dim, np, times, stdev))
