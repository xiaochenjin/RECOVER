import numpy as np
import os
from shutil import copyfile

input_name = 'fkr-SnSn-direct-Nconfig-1000.txt'
#input_folder = '/home/xcjin/Research/other-small-things/generate-perturbed-cell/bulk-GeSn/64000-atom-cell/Sn-0.25/Random/analysis/recover-test-4/test1-64000-atom/fkr-files-original'
input_folder = '/home/xcjin/Research/other-small-things/generate-perturbed-cell/bulk-GeSn/64000-atom-cell/Sn-0.25/perlmutter/fkr-files'

#output_folder = 'fkr-files-original'

#if not os.path.exists(output_folder):
#	os.mkdir(output_folder)

mu = 0
sigma_list = ['1.0','1.2','1.4','1.6','1.8','2.0','2.3','2.5','2.7','3.0','3.3','3.5','3.7','4.0','4.3','4.5','4.7','5.0'] #,'5.5','6.5','7.0','7.5','8.0','8.5','9.0','9.5','10.0']

for sigma in sigma_list:
	directory_name = "mu-"+str(mu)+"-sgm-"+str(sigma)
	input_folder2 = os.path.join(input_folder,directory_name)
#	output_folder2 = os.path.join(output_folder,directory_name)

	output_name = 'fkr-sigma-{0}.txt'.format(sigma)

#	if not os.path.exists(output_folder2):
#		os.mkdir(output_folder2)
	
	copyfile(os.path.join(input_folder2,input_name),output_name)

	
