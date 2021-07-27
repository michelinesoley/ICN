
###################################
### Reads energy and plots it	###
###################################

import numpy as np
from matplotlib import pyplot as plt

nptx = 35
npty = 18
nptz = 7

### Define input and output files

inputfile_gs='energy_gs.dat'
inputfile_mom='energy_mom.dat'
inputfile_mom_excl='energy_mom_excl.dat'
inputfile_mom_incl='energy_mom_incl.dat'

### Read input

data_gs = np.loadtxt(inputfile_gs)
data_mom = np.loadtxt(inputfile_mom)
data_mom_excl = np.loadtxt(inputfile_mom_excl)
data_mom_incl = np.loadtxt(inputfile_mom_incl)

### Shift energy wrt minima

min_gs = np.min(data_gs[:,1])
min_mom = np.min(data_mom[:,1])

print('min_gs = ',min_gs)
print('min_mom = ',min_mom)

energy_min = min(min_gs,min_mom)

data_gs[:,1] -= energy_min
data_mom[:,1] -= energy_min
data_mom_excl[:,1] -= energy_min
data_mom_incl[:,1] -= energy_min

### Convert from au to eV

au2ev = 27.2

data_gs[:,1] *= au2ev
data_mom[:,1] *= au2ev
data_mom_excl[:,1] *= au2ev
data_mom_incl[:,1] *= au2ev

### Plot data

plt.rcParams['axes.linewidth'] = 1.5

for i in range(npty):

	nstart = i*nptx*nptz 
	nend = nstart + nptx*nptz 

	fig1,ax1=plt.subplots()

	#fig.subplots_adjust(bottom=0.12,right=1.0,top=0.92)

	lw=2.

	ax1.plot(data_gs[:,0],data_gs[:,1],'o-',c='k',label='GS')
	ax1.plot(data_mom[:,0],data_mom[:,1],'o-',c='r',label='MOM')
	ax1.plot(data_mom_excl[:,0],data_mom_excl[:,1],'o',c='b',label='MOM excl')

	### set labels

	ax1.set_title('Energy #{}'.format(i),size=18)

	sz = 14

	ax1.set_xlabel(r'file index',size=sz)

	ax1.set_ylabel(r'Energy (eV)',size=sz)

	### set axis limits

	#x_axis=[280,305]
	y_axis=[0,10]
	#ax.set_xlim(x_axis)
	ax1.set_ylim(y_axis)

	### set axis ticks

	labelsz = 14

	ax1.tick_params(axis='both',which='major',labelsize=labelsz)

	### plot legend

	ax1.legend()

	### plot lines separating different angles

	for j in range(nptz-1):
		xpoint = (j+1)*nptx*nptz + nstart
		ax1.axvline(x=xpoint,lw=2)

	### set axis limits

	x_axis=[nstart,nend]
	#y_axis=[0,3000]

	ax1.set_xlim(x_axis)
	#ax1.set_ylim(y_axis)

	### plot

#	plt.show()

	### Save plot
	fig1.savefig('energy_selection_{:02d}.ps'.format(i),format='ps',dpi=1200)

#=====================================================================
### End program

print( '\n')
print( 'DONE!!')
print( '\n')

