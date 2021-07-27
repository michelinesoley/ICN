
#########################################################################
# Reads energy from Qchem calculation and separate contribution 	#
# from first ground-state calculation from second MOM calculation	#	 
#########################################################################

### Define input and output files

inputfile='energy.dat'

outputfile_gs=inputfile.split('.')[0]+'_gs.dat'
outputfile_mom=inputfile.split('.')[0]+'_mom.dat'

### Read input and open output

lines=open(inputfile,'r').readlines()

out_gs=open(outputfile_gs,'w')
out_mom=open(outputfile_mom,'w')

### Separate energy of GS from MOM 

num_old=1000

for i in range(len(lines)):

	line = lines[i].split('.out')[0]
	num = line.split('_')[1]

	if(num==num_old):
		out_mom.write('{} {}'.format(num,lines[i].split('=')[1]))
	else:
		out_gs.write('{} {}'.format(num,lines[i].split('=')[1]))

	num_old = num
