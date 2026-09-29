import xml.etree.ElementTree as ET
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import time  # <-- Add this line

parser = argparse.ArgumentParser()
parser.add_argument('type', help='aneurysm type (sym/asym)')
args = parser.parse_args()
#

# aneurysm type
type_ = args.type	# CASENAME
	
#mu_N = 5 #100
#cc_N = 5 #100

# density (points in [0,1])
dens_ = 20 #100 #25 #100 

# min degradation = 0% 
mu_max = 1.; cc_max = 1.

# max degradation
'''
if 'axis' in type_ :
	mu_min = 0.4
	cc_min = 0.6

elif 'asym' in type_: 
	mu_min = 0.1
	cc_min = 0.3

complem = False #True
# complementary
if complem and 'axis' in type_:
	mu_min = 0.
	mu_max = .39
if complem and 'asym' in type_:
	mu_min = 0.
	mu_max = .09
'''

# Max degradation
mu_min = 0.	# 100%
cc_min = 0. 	#.5 #0.

# limits cases w mass fractions < 1
'''
if 'axis' in type_ :
	cc_min = 0.85

elif 'asym' in type_: 
	cc_min = 0.63
'''


#else:
#'''
#mu_min = 0.
#cc_min = 0.
	

mu_range = mu_max - mu_min
cc_range = cc_max - cc_min

mu_N = int(mu_range*dens_)
cc_N = int(cc_range*dens_)

print('mu_N =',mu_N)
print('cc_N =',cc_N)

# 1. Define the parameter space (uniform)
cc_values = np.linspace(cc_min, cc_max, cc_N) # 
mu_values = np.linspace(mu_min, mu_max, mu_N) # 

# 1. Define the parameter space using Chebyshev spacing (denser at boundaries)
# Math: x = min + (max - min) * (1 - cos(theta)) / 2
#cc_values = cc_min + (cc_max - cc_min) * (1 - np.cos(np.linspace(0, np.pi, cc_N))) / 2
#mu_values = mu_min + (mu_max - mu_min) * (1 - np.cos(np.linspace(0, np.pi, mu_N))) / 2

# Define the folder where your FEBio file lives
work_dir = os.path.dirname(os.path.abspath(__file__)) + "/"

# input file (base)
#base_file = 'TAA-1x20x20Q-e-KG-'+type_+'-MBE-base.feb'
base_file = 'TAA-'+type_+'-MBE-BASE.feb'

# Update the base_file to use this full path
#base_file = os.path.join(work_dir, input_f)
tree = ET.parse(base_file)
root = tree.getroot()

output_files = []

# Scan the <Output><logfile> section for any data tags with a 'file' attribute
logfile_node = root.find('.//Output/logfile')
if logfile_node is not None:
	for data_tag in logfile_node:
		if 'file' in data_tag.attrib:
			output_files.append(data_tag.attrib['file'])

print(f"Found the following output files in the .feb file: {output_files}")

'''
# Ensure there is a logfile node inside the <Output> block
output_node = root.find('.//Output')
if output_node is None:
	output_node = ET.SubElement(root, 'Output')
logfile_node = output_node.find('logfile')
if logfile_node is None:
	logfile_node = ET.SubElement(output_node, 'logfile')
'''
	
# Set your custom FEBio executable path
febio_executable = __import__("os").environ.get("FEBIO", "febio4")


# NEW: List to store summary data for all simulations
summary_data_EA = []
summary_data_EB = []
summary_data_NA = []
summary_data_NB = []

total_sims = len(cc_values) * len(mu_values)
current_sim = 0

succ = 0
# 2. Loop through combinations of CC_in and MU_in
for i, cc in enumerate(cc_values):
	for j, mu in enumerate(mu_values):
		
		current_sim += 1
		percent_complete = (current_sim / total_sims) * 100
		
		# Update the material parameters in the XML tree
		mat_node = root.find('.//material[@name="Material1"]')
		if mat_node is not None:
			cc_node = mat_node.find('cc_in')
			mu_node = mat_node.find('mu_in')
			if cc_node is not None: cc_node.text = str(cc)
			if mu_node is not None: mu_node.text = str(mu)
		
		
		# Save the new .feb file for this iteration WITH the required XML header
		dummy_f = f"run_cc_var_mu_var_"+type_+".feb"
		case_filename = dummy_f
		tree.write(case_filename, xml_declaration=True)  # save file
		
		
		# 3. Run FEBio using custom build
		print('--------------------------------------------------')
		print(f"[{current_sim}/{total_sims} - {percent_complete:.2f}%]; Succ = {succ}")# ... Running FEBio for CC={cc:.2f}, MU={mu:.2f}...")
		print(f"Running sim. for CC={cc:.2f}, MU={mu:.2f} ...")
		
		# Start the timer
		start_time = time.time()
		
		subprocess.run([febio_executable, "-i", work_dir + case_filename], stdout=subprocess.DEVNULL) # remove last part for full ourput
		
		# Stop the timer and calculate the duration
		end_time = time.time()
		runtime = end_time - start_time
		
		# -------------------------------------------------------------------
		# NEW: Check for NORMAL TERMINATION in the generated log file
		log_file = dummy_f.replace('.feb', '.log')
		normal_termination = False	
		
		#if os.path.exists(log_file_path):
		with open(log_file, 'r') as log_f:
			# Leer el archivo y verificar si contiene el string exacto
			if " N O R M A L   T E R M I N A T I O N" in log_f.read():
				normal_termination = True
				succ = succ + 1
		
		# Si no terminó normalmente, saltar la lectura de datos de esta iteración
		if not normal_termination:
			#print(f"--> ERROR: Simulation CC={cc:.2f}, MU={mu:.2f} did not terminate normally. Skipping data extraction.")
			print(f"--> ERROR: Simulation did not converge. Skipping ...")
			# Opcional: summary_data.append([mu, cc, np.nan, np.nan]) para mantener el registro
			print(f"Runtime: {runtime:.2f} seconds")
			continue 
		else:
			print(f"--> SUCCESS: Normal termination.")# for CC={cc:.2f}, MU={mu:.2f}.")
			print(f"Runtime: {runtime:.2f} seconds")
		
		# -------------------------------------------------------------------
		
		# 4. Loop through and read the extracted output files
		Jmin = 1; Jmax = 1
		data0 = [mu, cc]
		for filename in output_files:
		
			# read only element solution
			# element: data="sx;sy;sz;J"
			if 'elem-' in filename:	#os.path.exists(full_log_path):
				#print(f"\nReading {filename}...")
		#try:
				# Load the text file: ignore '*' headers, split by ','
				data = np.loadtxt(filename, delimiter=',', comments='*')
			
				# Example: extracting J (column 5 / index 4) 
				if data.ndim > 1:	values0 = data[-1, :] 
				else:			values0 = np.array([data[:]])
				
				values = np.concatenate((data0,values0)) 
				
				# Extract data of single simulation
				if '-A' in filename:	# max contraction
					# Save summary (mu, cc, Jmax)
					summary_data_EA.append(values)
				
				if '-B' in filename:	# max expansion
					# Save summary (mu, cc, Jmax)
					summary_data_EB.append(values)
					
			# node: data="x;ux;y;uy"
			if 'node-' in filename:	
				#print(f"\nReading {filename}...")
				
				# Load the text file: ignore '*' headers, split by ','
				data = np.loadtxt(filename, delimiter=',', comments='*')
			
				if data.ndim > 1:	values0 = data[-1, :] 
				else:			values0 = np.array([data[:]])
				
				values = np.concatenate((data0,values0)) 
				
				# Extract data of single simulation
				if '-A' in filename:	
					summary_data_NA.append(values)
				
				if '-B' in filename:	
					summary_data_NB.append(values)
			
					

# NEW: Save data from all sims in another file
#head = '# mu_in, cc_in, J_min, J_max'
el_head = '# mu_in;	cc_in;	elem;	sx;	sy;	sz;	J; 	sed;'
summary_file_EA = os.path.join(work_dir, "summary_"+type_+"_elem_A.csv")
summary_file_EB = os.path.join(work_dir, "summary_"+type_+"_elem_B.csv")

nd_head = '# mu_in;	cc_in;	node;	x;	ux;	y;	uy'
summary_file_NA = os.path.join(work_dir, "summary_"+type_+"_node_A.csv")
summary_file_NB = os.path.join(work_dir, "summary_"+type_+"_node_B.csv")

np.savetxt(summary_file_EA, summary_data_EA, delimiter=',', header=el_head, comments='')
np.savetxt(summary_file_EB, summary_data_EB, delimiter=',', header=el_head, comments='')

np.savetxt(summary_file_NA, summary_data_NA, delimiter=',', header=nd_head, comments='')
np.savetxt(summary_file_NB, summary_data_NB, delimiter=',', header=nd_head, comments='')

print(f"All ELEMENT summary data saved to: {summary_file_EA} and {summary_file_EB}")
#
