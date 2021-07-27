
###########################################
### Reads GS and MOM energy and makes	###
### plots according to R_IC cut		###
###########################################

import numpy as np
from matplotlib import pyplot as plt

nptx = 35
npty = 18
nptz = 7

### Define input and output files

inputfile_gs='energy_gs.dat'
inputfile_mom='energy_mom.dat'

### Read input

data_gs = np.loadtxt(inputfile_gs)
data_mom = np.loadtxt(inputfile_mom)

### Shift energy wrt minima

min_gs = np.min(data_gs[:,1])
min_mom = np.min(data_mom[:,1])

print('min_gs = ',min_gs)
print('min_mom = ',min_mom)

energy_min = min(min_gs,min_mom)

data_gs[:,1] -= energy_min
data_mom[:,1] -= energy_min

### Convert from au to eV

au2ev = 27.2

data_gs[:,1] *= au2ev
data_mom[:,1] *= au2ev

### Plot data

plt.rcParams['axes.linewidth'] = 1.5

for i in range(npty):

	nstart = i*nptx*nptz 
	nend = nstart + nptx*nptz 

	fig1,ax1=plt.subplots()

	lw=2.

	ax1.plot(data_gs[:,0],data_gs[:,1],'o-',c='k',label='GS')
	ax1.plot(data_mom[:,0],data_mom[:,1],'o-',c='r',label='MOM')

	### set labels

	ax1.set_title('Energy #{}'.format(i),size=18)

	sz = 14

	ax1.set_xlabel(r'file index',size=sz)

	ax1.set_ylabel(r'Energy (eV)',size=sz)

	### set axis ticks

	labelsz = 14

	ax1.tick_params(axis='both',which='major',labelsize=labelsz)

	### plot legend

	ax1.legend()

	### plot lines separating different r_CN

	for j in range(nptz-1):
		xpoint = (j+1)*nptx + nstart
		ax1.axvline(x=xpoint,ls='--',lw=1)

	### set axis limits

	x_axis=[nstart,nend]
	y_axis=[0,10]

	ax1.set_xlim(x_axis)
	ax1.set_ylim(y_axis)

	### plot
#	plt.show()

	### Save plot
#	fig1.savefig('energy_{}.pdf'.format(i),format='pdf',dpi=1200)
	fig1.savefig('energy_{:02d}.ps'.format(i),format='ps',dpi=1200)

#=====================================================================
### End program

print('\n')
print('DONE!!')
print('\n')

