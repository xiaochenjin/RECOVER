import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.legend_handler import HandlerTuple

mu = 0
perturb_type = input("perturb type: ")
collect_eff = input("collect efficiency: ")
index = input("index: ")

if perturb_type == 'iso':
    sigma = input("sigma: ")
#   analysis = '../../analysis/recover-test-4/test1-64000-atom/efficiency-{0}'.format(collect_eff)
    analysis = '../simulated-structures/isotropic-perturbation/efficiency-{0}'.format(collect_eff)
    directory_name = "mu-"+str(mu)+"-sgm-"+str(sigma)

if perturb_type == 'aniso':
    sigma_z = input("sigma z: ")
    sigma_xy = input("sigma xy: ")
#   analysis = '../../analysis/recover-test-5-anisotropic/test1-64000-atom/efficiency-{0}'.format(collect_eff)
    analysis = '../simulated-structures/anisotropic-perturbation/efficiency-{0}'.format(collect_eff)
    directory_name = "mu-"+str(mu)+"-sgm_z-"+str(sigma_z)+"-sgm_xy-"+str(sigma_xy)

folder = os.path.join(analysis,directory_name)
file_name = 'displacement-list-index-{0}.txt'.format(index)

with open (os.path.join(folder,file_name)) as f1:
#with open ("displacement-list-mu-"+str(mu)+"-sgm-"+str(sigma)+".txt") as f1:
	r_list = np.loadtxt(f1, delimiter=' ', usecols=(0),unpack=True)

mu = float(mu)
sigma = float(sigma)

plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']

fontproperties = {'fontweight' : 'bold', 'fontsize' : 13}
legend_properties = {'weight':'bold','size':13,'style': 'italic'}

fig, axs = plt.subplots(1,1, figsize=(6,5))#,dpi=500)

count, bins, ignored = axs.hist(r_list, 100, density=True)

axs.tick_params(which='minor',direction="in")
axs.tick_params(which='major',direction="in")

#axs.xaxis.grid(False, which='minor')
axs.xaxis.set_minor_locator(AutoMinorLocator(2))
axs.yaxis.set_minor_locator(AutoMinorLocator(2))

axs.set_xlabel("Atomic displacement",fontsize='15',weight = 'bold')
axs.set_xlabel("Atomic displacement (Å)",fontsize='15',weight = 'bold')
axs.set_ylabel("Frequency",fontsize='15',weight = 'bold')
axs.set_xticklabels(axs.get_xticks(),fontproperties)
axs.set_yticklabels(axs.get_yticks(),fontproperties)
#axs.xaxis.set_tick_params(labelcolor='none') #hide xticks labels
axs.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f')) #https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels
axs.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

axs.plot(bins, 1/(sigma * np.sqrt(2 * np.pi))*np.exp( - (bins - mu)**2 / (2 * sigma**2) ),linewidth=2, color='r')
plt.xticks(np.arange(-2*sigma, 3*sigma, sigma))
plt.xlim(-3*sigma, 3*sigma)
#plt.ylim(top=0.2)
plt.yticks([]) #hiding y axis
plt.tight_layout()
plt.show()
