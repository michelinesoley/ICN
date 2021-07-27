
###########################################################################
### Reads energy from 'energy_mom.dat' and performs a selection of	###	
### energies based on threshold 					###
###########################################################################

import numpy as np

nptx = 35
npty = 18
nptz = 7

### Define input and output files

inputfile='energy_mom.dat'

outputfile_incl=inputfile.split('.')[0]+'_incl.dat'
outputfile_excl=inputfile.split('.')[0]+'_excl.dat'

### Read input and open output

data=np.loadtxt(inputfile)

out_incl=open(outputfile_incl,'w')
out_excl=open(outputfile_excl,'w')

### Define forward and backward selection criteria

def forward_selection(j,k):

	### get first point in the series of nptx scan
	for i in range(nptx):
		ptr_prev = j*nptz*nptx + k*nptx + i
		idx_prev, = np.where(ptr_prev == data[:,0])

		if(len(idx_prev)!=0):  
			i_first = i
			idx_prev = int(idx_prev) 
			data_prev  = data[idx_prev,1] 
			break

	print('First point in scan ({},{}) is {}.'.format(j+1,k+1,i_first+1))

	### if first point is greater than nfirst_thr then exclude whole scan
	nfirst_thr = 5
	if(i_first>nfirst_thr):
		for i in range(nptx):
			ptr_now = j*nptz*nptx + k*nptx + i 
			idx_now, = np.where(ptr_now == data[:,0])

			if (len(idx_now) != 0):
				idx_now = int(idx_now) 
				data_now  = data[idx_now,1] 
#				print 'Excluding point',ptr_now
				out_excl.write('{} {}\n'.format(ptr_now,data_now))
		
	### otherwise check every point
	else:
		### loop over nptx scan and select energies
		slope = 1000
		for i in range(nptx):
			ptr_now = j*nptz*nptx + k*nptx + i 
			idx_now, = np.where(ptr_now == data[:,0])

			if (len(idx_now) != 0):
				idx_now = int(idx_now) 
				data_now  = data[idx_now,1] 

				diff = data_now - data_prev

				if(diff>energy_thr): # if energy rises more than thr exclude point
#					print 'Excluding point',ptr_now
					out_excl.write('{} {}\n'.format(ptr_now,data_now))
				else:
#					diff /= (ptr_now - ptr_prev)
					diff = np.abs(diff)

#					slope_thr=3.*slope
					slope_thr = max(energy_thr,nslope_thr*slope)

					if(diff>slope_thr): # if energy changes more than slope_thr exclude point
#						print 'Excluding point',ptr_now
						out_excl.write('{} {}\n'.format(ptr_now,data_now))
					else:
#						print 'Including point',ptr_now
						out_incl.write('{} {}\n'.format(ptr_now,data_now))

						data_prev = data_now
						ptr_prev = ptr_now

						'''
						if(i == i_first):
							slope = 1000
						else:
							slope = diff
						'''
	return

def backward_selection(j,k):

	### get last point in the series of nptx scan
	for i in range(nptx):
		ptr_prev = j*nptz*nptx + k*nptx + nptx - 1 - i
		idx_prev, = np.where(ptr_prev == data[:,0])

		if(len(idx_prev)!=0):  
			i_last = nptx - 1 - i
			idx_prev = int(idx_prev) 
			data_prev  = data[idx_prev,1] 
			break

	print('Last point in scan ({},{}) is {}.'.format(j+1,k+1,i_last+1))

	### if last point is less than nlast_thr then exclude whole scan
	nlast_thr = nptx - 20
	if(i_last<nlast_thr):
		for i in range(nptx):
			ptr_now = j*nptz*nptx + k*nptx + i 
			idx_now, = np.where(ptr_now == data[:,0])

			if (len(idx_now) != 0):
				idx_now = int(idx_now) 
				data_now  = data[idx_now,1] 
#				print 'Excluding point',ptr_now
				out_excl.write('{} {}\n'.format(ptr_now,data_now))
		
	### otherwise check every point
	else:
		### loop over nptx scan and select energies
		slope = 1000
		for i in range(nptx):
			ptr_now = j*nptz*nptx + k*nptx + nptx - 1 - i 
			idx_now, = np.where(ptr_now == data[:,0])

			if (len(idx_now) != 0):
				idx_now = int(idx_now) 
				data_now  = data[idx_now,1] 

				diff = data_now - data_prev
	
				if(-diff>energy_thr): # if energy decrease more than thr exclude point
#					print 'Excluding point',ptr_now
					out_excl.write('{} {}\n'.format(ptr_now,data_now))
				else:
#					diff /= (ptr_now - ptr_prev)
					diff = np.abs(diff)

#					slope_thr=3.*slope
					slope_thr = max(energy_thr,nslope_thr*slope)

					if(diff>slope_thr): # if energy changes more than slope_thr exclude point
#						print 'Excluding point',ptr_now
						out_excl.write('{} {}\n'.format(ptr_now,data_now))
					else:
#						print 'Including point',ptr_now
						out_incl.write('{} {}\n'.format(ptr_now,data_now))

						data_prev = data_now
						ptr_prev = ptr_now

						'''
						if(nptx - 1 - i == i_last):
							slope = 1000
						else:
							slope = diff
						'''
	return

### Select energies

au2ev = 27.2

energy_thr = 0.3
energy_thr /= au2ev ### Convert from eV to au
nslope_thr=3

for j in range(npty): 
	for k in range(nptz): 

			forward_selection(j,k)


