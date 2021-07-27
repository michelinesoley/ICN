#!/usr/bin/python

###########################################################
### Reads QCHEM outputfile and generate response matrix	###
### excluding some points selected by energy criteria	###
###########################################################

import sys,re,os
import numpy as np

excl_list_file = 'exclusion_list.dat'
excl_list = np.loadtxt(excl_list_file,unpack=True,usecols=[0])

outputfile = open('response_mat.dat','w')

nlines=33  #number of transition considered

#====================================================================================================
### Generate array with position (in au and rad)
# x: represents the variable R_jacobi (in Angstrom) 
# y: represents the variable theta_jacobi (in degree)
# z: represents the variable r_CN (in Angstrom)

ang2au = 1./0.529177 
deg2rad = np.pi/180. 

nptx = 35
npty = 18 
nptz = 7 
vol=nptx*npty*nptz

print('')
print('nptx =',nptx)
print('npty =',npty)
print('nptz =',nptz)

x_min = 2.5
x_max = 6.
dx    = (x_max-x_min)/nptx 

print('')
print('x_min =',x_min)
print('x_max =',x_max)
print('dx =',dx)

y_min = 0.0
y_max = 180.
dy    = (y_max-y_min)/npty

print('')
print('y_min =',y_min)
print('y_max =',y_max)
print('dy =',dy)

z_min = 1.0
z_max = 1.3
dz    = (z_max-z_min)/nptz

print('')
print('z_min =',z_min)
print('z_max =',z_max)
print('dz =',dz)

x=np.zeros(vol)
y=np.zeros(vol)
z=np.zeros(vol)

ncount=0
for j in range(npty):
	temp_y = y_min + j*dy
	for k in range(nptz):
		temp_z = z_min + k*dz
		for i in range(nptx):
			temp_x = x_min + i*dx

			x[ncount]=temp_x
			y[ncount]=temp_y
			z[ncount]=temp_z

			ncount+=1

x=x*ang2au
y=y*deg2rad
z=z*ang2au

#====================================================================================================
### Extract frequencies and intensities from input file

directory = '../../'
sufix='.out'

for filename in sorted(os.listdir(directory)):
	if filename.endswith(sufix):

		print('Reading ...',filename)

		inputfile = open(directory+filename, 'r')

		### check completion of run
		last_line = inputfile.readlines()[-5]
		completion = ('Thank you very much for using Q-Chem.  Have a nice day.' in last_line)

		if(not completion):
			print('Excluding filename {} due to failed completion'.format(filename))
			continue

		iflag = 0
		istate=[]
		energy=[]
		intens=[]

		inputfile.seek(0)
		for line in inputfile:

			temp=line.split()

			if(len(temp)==0):
				iflag = 0

			if (iflag==1):
				if re.search('Strength', line):
					intens.append(float(temp[-1]))
			if re.search('Excited state', line):
				temp[2]=temp[2].replace(':','')	
				istate.append(int(temp[2]))
				energy.append(float(temp[-1]))
				iflag = 1

		inputfile.close()

### Get index of file

		temp=list(filename)
		index=temp[6]+temp[7]+temp[8]+temp[9]
		index=int(index)

### Check exclusion list

		if(index in excl_list):
			print('Excluding filename {} due to exclusion list'.format(filename))
			continue
			
### Get position for this index 

		R     = x[index]
		theta = y[index]
		rCN   = z[index]

### Save data in file in the form theta,R,rCN,energy,intensity

		for i in range(nlines):
			outputfile.write('{} {} {} {} {} \n'.format(theta,R,rCN,energy[i],intens[i])) 


#====================================================================================================
### Exit

print('\n')
print('------------')
print('\nDONE!!!\n')
print('------------')
print('\n')

