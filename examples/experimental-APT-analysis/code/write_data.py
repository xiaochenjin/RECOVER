#https://github.com/oscarbranson/apt-tools/Data Import.ipynb
import numpy as np
import pandas as pd
import struct
import re

def read_pos(f):
	""" Loads an APT .pos file as a pandas dataframe.

	Columns:
		x: Reconstructed x position
		y: Reconstructed y position
		z: Reconstructed z position
		Da: mass/charge ratio of ion"""
	# read in the data
	n = len(open(f,'rb').read())//4
	d = struct.unpack('>'+'f'*n,open(f,'rb').read(4*n))
					# '>' denotes 'big-endian' byte order
	# unpack data
	pos = pd.DataFrame({'x': d[0::4],
						'y': d[1::4],
						'z': d[2::4],
						'Da': d[3::4]})
	return pos

def read_rrng(f):
	rf = open(f,'r').readlines()

	patterns = re.compile(r'Ion([0-9]+)=([A-Za-z0-9]+).*|Range([0-9]+)=(\d+.\d+) +(\d+.\d+) +Vol:(\d+.\d+) +([A-Za-z:0-9 ]+) +Color:([A-Z0-9]{6})')

	ions = []
	rrngs = []
	for line in rf:
		m = patterns.search(line)
		if m:
			if m.groups()[0] is not None:
				ions.append(m.groups()[:2])
			else:
				rrngs.append(m.groups()[2:])

	ions = pd.DataFrame(ions, columns=['number','name'])
	ions.set_index('number',inplace=True)
	rrngs = pd.DataFrame(rrngs, columns=['number','lower','upper','vol','comp','colour'])
	rrngs.set_index('number',inplace=True)
	
	rrngs[['lower','upper','vol']] = rrngs[['lower','upper','vol']].astype(float)
	rrngs[['comp','colour']] = rrngs[['comp','colour']].astype(str)
	
	return ions,rrngs

def label_ions(pos,rrngs):
	pos['comp'] = ''
	pos['colour'] = '#FFFFFF'
	
	for n,r in rrngs.iterrows():
		pos.loc[(pos.Da >= r.lower) & (pos.Da <= r.upper),['comp','colour']] = [r['comp'],'#' + r['colour']]
	
	return pos

def deconvolve(lpos):
	"""Takes a composition-labelled pos file, and deconvolves
	the complex ions. Produces a dataframe of the same input format
	with the extra columns:
	   'element': element name
	   'n': stoichiometry
	For complex ions, the location of the different components is not
	altered - i.e. xyz position will be the same for several elements."""
	  
	out = []
	pattern = re.compile(r'([A-Za-z]+):([0-9]+)')

	for g,d in lpos.groupby('comp'):
		if g is not '':
			for i in range(len(g.split(' '))):
				tmp = d.copy()
				cn = pattern.search(g.split(' ')[i]).groups()
				tmp['element'] = cn[0]
				tmp['n'] = cn[1]
				out.append(tmp.copy())
	return pd.concat(out)

#https://stackoverflow.com/questions/31247198/python-pandas-write-content-of-dataframe-into-text-file
def write(lpos):
	file_name = 'all-raw-data.txt'
	np.savetxt(file_name, lpos.values,fmt='%s')

pos = read_pos('R6009_07581_M7_Sample2.POS')
#print(pos)
ions, rrngs = read_rrng('R6009_07581_M7_Sample2.RRNG')
lpos = label_ions(pos,rrngs)
write(lpos)
