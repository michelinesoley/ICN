#!/bin/bash
#SBATCH -J run_script_1
#SBATCH --partition=pi_esi
#SBATCH --time 120:00:00
#SBATCH --ntasks=12 --nodes=1
#SBATCH --mem-per-cpu=5000
#SBATCH -o run_script_1.err

### Load MPI module
#module load MPI/OpenMPI/2.1.1-intel15

### Load QCHEM module
module load Q-Chem/5.2-openmpi

export QCSCRATCH=/gpfs/loomis/scratch60/fas/batista/pev4/qchem-scratch/
 
### Run program

echo 'Starting program on'
date
echo ''

qchem -nt $SLURM_NTASKS qchem_0000.in qchem_0000.out
qchem -nt $SLURM_NTASKS qchem_0001.in qchem_0001.out
qchem -nt $SLURM_NTASKS qchem_0002.in qchem_0002.out
qchem -nt $SLURM_NTASKS qchem_0003.in qchem_0003.out
qchem -nt $SLURM_NTASKS qchem_0004.in qchem_0004.out
qchem -nt $SLURM_NTASKS qchem_0005.in qchem_0005.out
qchem -nt $SLURM_NTASKS qchem_0006.in qchem_0006.out
qchem -nt $SLURM_NTASKS qchem_0007.in qchem_0007.out
qchem -nt $SLURM_NTASKS qchem_0008.in qchem_0008.out
qchem -nt $SLURM_NTASKS qchem_0009.in qchem_0009.out
qchem -nt $SLURM_NTASKS qchem_0010.in qchem_0010.out
qchem -nt $SLURM_NTASKS qchem_0011.in qchem_0011.out
qchem -nt $SLURM_NTASKS qchem_0012.in qchem_0012.out
qchem -nt $SLURM_NTASKS qchem_0013.in qchem_0013.out
qchem -nt $SLURM_NTASKS qchem_0014.in qchem_0014.out
qchem -nt $SLURM_NTASKS qchem_0015.in qchem_0015.out
qchem -nt $SLURM_NTASKS qchem_0016.in qchem_0016.out
qchem -nt $SLURM_NTASKS qchem_0017.in qchem_0017.out
qchem -nt $SLURM_NTASKS qchem_0018.in qchem_0018.out
qchem -nt $SLURM_NTASKS qchem_0019.in qchem_0019.out
qchem -nt $SLURM_NTASKS qchem_0020.in qchem_0020.out
qchem -nt $SLURM_NTASKS qchem_0021.in qchem_0021.out
qchem -nt $SLURM_NTASKS qchem_0022.in qchem_0022.out
qchem -nt $SLURM_NTASKS qchem_0023.in qchem_0023.out
qchem -nt $SLURM_NTASKS qchem_0024.in qchem_0024.out
qchem -nt $SLURM_NTASKS qchem_0025.in qchem_0025.out
qchem -nt $SLURM_NTASKS qchem_0026.in qchem_0026.out
qchem -nt $SLURM_NTASKS qchem_0027.in qchem_0027.out
qchem -nt $SLURM_NTASKS qchem_0028.in qchem_0028.out
qchem -nt $SLURM_NTASKS qchem_0029.in qchem_0029.out
qchem -nt $SLURM_NTASKS qchem_0030.in qchem_0030.out
qchem -nt $SLURM_NTASKS qchem_0031.in qchem_0031.out
qchem -nt $SLURM_NTASKS qchem_0032.in qchem_0032.out
qchem -nt $SLURM_NTASKS qchem_0033.in qchem_0033.out
qchem -nt $SLURM_NTASKS qchem_0034.in qchem_0034.out
qchem -nt $SLURM_NTASKS qchem_0035.in qchem_0035.out
qchem -nt $SLURM_NTASKS qchem_0036.in qchem_0036.out
qchem -nt $SLURM_NTASKS qchem_0037.in qchem_0037.out
qchem -nt $SLURM_NTASKS qchem_0038.in qchem_0038.out
qchem -nt $SLURM_NTASKS qchem_0039.in qchem_0039.out
qchem -nt $SLURM_NTASKS qchem_0040.in qchem_0040.out
qchem -nt $SLURM_NTASKS qchem_0041.in qchem_0041.out
qchem -nt $SLURM_NTASKS qchem_0042.in qchem_0042.out
qchem -nt $SLURM_NTASKS qchem_0043.in qchem_0043.out
qchem -nt $SLURM_NTASKS qchem_0044.in qchem_0044.out
qchem -nt $SLURM_NTASKS qchem_0045.in qchem_0045.out
qchem -nt $SLURM_NTASKS qchem_0046.in qchem_0046.out
qchem -nt $SLURM_NTASKS qchem_0047.in qchem_0047.out
qchem -nt $SLURM_NTASKS qchem_0048.in qchem_0048.out
qchem -nt $SLURM_NTASKS qchem_0049.in qchem_0049.out
qchem -nt $SLURM_NTASKS qchem_0050.in qchem_0050.out
qchem -nt $SLURM_NTASKS qchem_0051.in qchem_0051.out
qchem -nt $SLURM_NTASKS qchem_0052.in qchem_0052.out
qchem -nt $SLURM_NTASKS qchem_0053.in qchem_0053.out
qchem -nt $SLURM_NTASKS qchem_0054.in qchem_0054.out
qchem -nt $SLURM_NTASKS qchem_0055.in qchem_0055.out
qchem -nt $SLURM_NTASKS qchem_0056.in qchem_0056.out
qchem -nt $SLURM_NTASKS qchem_0057.in qchem_0057.out
qchem -nt $SLURM_NTASKS qchem_0058.in qchem_0058.out
qchem -nt $SLURM_NTASKS qchem_0059.in qchem_0059.out
qchem -nt $SLURM_NTASKS qchem_0060.in qchem_0060.out
qchem -nt $SLURM_NTASKS qchem_0061.in qchem_0061.out
qchem -nt $SLURM_NTASKS qchem_0062.in qchem_0062.out
qchem -nt $SLURM_NTASKS qchem_0063.in qchem_0063.out
qchem -nt $SLURM_NTASKS qchem_0064.in qchem_0064.out
qchem -nt $SLURM_NTASKS qchem_0065.in qchem_0065.out
qchem -nt $SLURM_NTASKS qchem_0066.in qchem_0066.out
qchem -nt $SLURM_NTASKS qchem_0067.in qchem_0067.out
qchem -nt $SLURM_NTASKS qchem_0068.in qchem_0068.out
qchem -nt $SLURM_NTASKS qchem_0069.in qchem_0069.out
qchem -nt $SLURM_NTASKS qchem_0070.in qchem_0070.out
qchem -nt $SLURM_NTASKS qchem_0071.in qchem_0071.out
qchem -nt $SLURM_NTASKS qchem_0072.in qchem_0072.out
qchem -nt $SLURM_NTASKS qchem_0073.in qchem_0073.out
qchem -nt $SLURM_NTASKS qchem_0074.in qchem_0074.out
qchem -nt $SLURM_NTASKS qchem_0075.in qchem_0075.out
qchem -nt $SLURM_NTASKS qchem_0076.in qchem_0076.out
qchem -nt $SLURM_NTASKS qchem_0077.in qchem_0077.out
qchem -nt $SLURM_NTASKS qchem_0078.in qchem_0078.out
qchem -nt $SLURM_NTASKS qchem_0079.in qchem_0079.out
qchem -nt $SLURM_NTASKS qchem_0080.in qchem_0080.out
qchem -nt $SLURM_NTASKS qchem_0081.in qchem_0081.out
qchem -nt $SLURM_NTASKS qchem_0082.in qchem_0082.out
qchem -nt $SLURM_NTASKS qchem_0083.in qchem_0083.out
qchem -nt $SLURM_NTASKS qchem_0084.in qchem_0084.out
qchem -nt $SLURM_NTASKS qchem_0085.in qchem_0085.out
qchem -nt $SLURM_NTASKS qchem_0086.in qchem_0086.out
qchem -nt $SLURM_NTASKS qchem_0087.in qchem_0087.out
qchem -nt $SLURM_NTASKS qchem_0088.in qchem_0088.out
qchem -nt $SLURM_NTASKS qchem_0089.in qchem_0089.out
qchem -nt $SLURM_NTASKS qchem_0090.in qchem_0090.out
qchem -nt $SLURM_NTASKS qchem_0091.in qchem_0091.out
qchem -nt $SLURM_NTASKS qchem_0092.in qchem_0092.out
qchem -nt $SLURM_NTASKS qchem_0093.in qchem_0093.out
qchem -nt $SLURM_NTASKS qchem_0094.in qchem_0094.out
qchem -nt $SLURM_NTASKS qchem_0095.in qchem_0095.out
qchem -nt $SLURM_NTASKS qchem_0096.in qchem_0096.out
qchem -nt $SLURM_NTASKS qchem_0097.in qchem_0097.out
qchem -nt $SLURM_NTASKS qchem_0098.in qchem_0098.out
qchem -nt $SLURM_NTASKS qchem_0099.in qchem_0099.out
qchem -nt $SLURM_NTASKS qchem_0100.in qchem_0100.out
qchem -nt $SLURM_NTASKS qchem_0101.in qchem_0101.out
qchem -nt $SLURM_NTASKS qchem_0102.in qchem_0102.out
qchem -nt $SLURM_NTASKS qchem_0103.in qchem_0103.out
qchem -nt $SLURM_NTASKS qchem_0104.in qchem_0104.out
qchem -nt $SLURM_NTASKS qchem_0105.in qchem_0105.out
qchem -nt $SLURM_NTASKS qchem_0106.in qchem_0106.out
qchem -nt $SLURM_NTASKS qchem_0107.in qchem_0107.out
qchem -nt $SLURM_NTASKS qchem_0108.in qchem_0108.out
qchem -nt $SLURM_NTASKS qchem_0109.in qchem_0109.out
qchem -nt $SLURM_NTASKS qchem_0110.in qchem_0110.out
qchem -nt $SLURM_NTASKS qchem_0111.in qchem_0111.out
qchem -nt $SLURM_NTASKS qchem_0112.in qchem_0112.out
qchem -nt $SLURM_NTASKS qchem_0113.in qchem_0113.out
qchem -nt $SLURM_NTASKS qchem_0114.in qchem_0114.out
qchem -nt $SLURM_NTASKS qchem_0115.in qchem_0115.out
qchem -nt $SLURM_NTASKS qchem_0116.in qchem_0116.out
qchem -nt $SLURM_NTASKS qchem_0117.in qchem_0117.out
qchem -nt $SLURM_NTASKS qchem_0118.in qchem_0118.out
qchem -nt $SLURM_NTASKS qchem_0119.in qchem_0119.out
qchem -nt $SLURM_NTASKS qchem_0120.in qchem_0120.out
qchem -nt $SLURM_NTASKS qchem_0121.in qchem_0121.out
qchem -nt $SLURM_NTASKS qchem_0122.in qchem_0122.out
qchem -nt $SLURM_NTASKS qchem_0123.in qchem_0123.out
qchem -nt $SLURM_NTASKS qchem_0124.in qchem_0124.out
qchem -nt $SLURM_NTASKS qchem_0125.in qchem_0125.out
qchem -nt $SLURM_NTASKS qchem_0126.in qchem_0126.out
qchem -nt $SLURM_NTASKS qchem_0127.in qchem_0127.out
qchem -nt $SLURM_NTASKS qchem_0128.in qchem_0128.out
qchem -nt $SLURM_NTASKS qchem_0129.in qchem_0129.out
qchem -nt $SLURM_NTASKS qchem_0130.in qchem_0130.out
qchem -nt $SLURM_NTASKS qchem_0131.in qchem_0131.out
qchem -nt $SLURM_NTASKS qchem_0132.in qchem_0132.out
qchem -nt $SLURM_NTASKS qchem_0133.in qchem_0133.out
qchem -nt $SLURM_NTASKS qchem_0134.in qchem_0134.out
qchem -nt $SLURM_NTASKS qchem_0135.in qchem_0135.out
qchem -nt $SLURM_NTASKS qchem_0136.in qchem_0136.out
qchem -nt $SLURM_NTASKS qchem_0137.in qchem_0137.out
qchem -nt $SLURM_NTASKS qchem_0138.in qchem_0138.out
qchem -nt $SLURM_NTASKS qchem_0139.in qchem_0139.out
qchem -nt $SLURM_NTASKS qchem_0140.in qchem_0140.out
qchem -nt $SLURM_NTASKS qchem_0141.in qchem_0141.out
qchem -nt $SLURM_NTASKS qchem_0142.in qchem_0142.out
qchem -nt $SLURM_NTASKS qchem_0143.in qchem_0143.out
qchem -nt $SLURM_NTASKS qchem_0144.in qchem_0144.out
qchem -nt $SLURM_NTASKS qchem_0145.in qchem_0145.out
qchem -nt $SLURM_NTASKS qchem_0146.in qchem_0146.out
qchem -nt $SLURM_NTASKS qchem_0147.in qchem_0147.out
qchem -nt $SLURM_NTASKS qchem_0148.in qchem_0148.out
qchem -nt $SLURM_NTASKS qchem_0149.in qchem_0149.out
qchem -nt $SLURM_NTASKS qchem_0150.in qchem_0150.out
qchem -nt $SLURM_NTASKS qchem_0151.in qchem_0151.out
qchem -nt $SLURM_NTASKS qchem_0152.in qchem_0152.out
qchem -nt $SLURM_NTASKS qchem_0153.in qchem_0153.out
qchem -nt $SLURM_NTASKS qchem_0154.in qchem_0154.out
qchem -nt $SLURM_NTASKS qchem_0155.in qchem_0155.out
qchem -nt $SLURM_NTASKS qchem_0156.in qchem_0156.out
qchem -nt $SLURM_NTASKS qchem_0157.in qchem_0157.out
qchem -nt $SLURM_NTASKS qchem_0158.in qchem_0158.out
qchem -nt $SLURM_NTASKS qchem_0159.in qchem_0159.out
qchem -nt $SLURM_NTASKS qchem_0160.in qchem_0160.out
qchem -nt $SLURM_NTASKS qchem_0161.in qchem_0161.out
qchem -nt $SLURM_NTASKS qchem_0162.in qchem_0162.out
qchem -nt $SLURM_NTASKS qchem_0163.in qchem_0163.out
qchem -nt $SLURM_NTASKS qchem_0164.in qchem_0164.out
qchem -nt $SLURM_NTASKS qchem_0165.in qchem_0165.out
qchem -nt $SLURM_NTASKS qchem_0166.in qchem_0166.out
qchem -nt $SLURM_NTASKS qchem_0167.in qchem_0167.out
qchem -nt $SLURM_NTASKS qchem_0168.in qchem_0168.out
qchem -nt $SLURM_NTASKS qchem_0169.in qchem_0169.out
qchem -nt $SLURM_NTASKS qchem_0170.in qchem_0170.out
qchem -nt $SLURM_NTASKS qchem_0171.in qchem_0171.out
qchem -nt $SLURM_NTASKS qchem_0172.in qchem_0172.out
qchem -nt $SLURM_NTASKS qchem_0173.in qchem_0173.out
qchem -nt $SLURM_NTASKS qchem_0174.in qchem_0174.out
qchem -nt $SLURM_NTASKS qchem_0175.in qchem_0175.out
qchem -nt $SLURM_NTASKS qchem_0176.in qchem_0176.out
qchem -nt $SLURM_NTASKS qchem_0177.in qchem_0177.out
qchem -nt $SLURM_NTASKS qchem_0178.in qchem_0178.out
qchem -nt $SLURM_NTASKS qchem_0179.in qchem_0179.out
qchem -nt $SLURM_NTASKS qchem_0180.in qchem_0180.out
qchem -nt $SLURM_NTASKS qchem_0181.in qchem_0181.out
qchem -nt $SLURM_NTASKS qchem_0182.in qchem_0182.out
qchem -nt $SLURM_NTASKS qchem_0183.in qchem_0183.out
qchem -nt $SLURM_NTASKS qchem_0184.in qchem_0184.out
qchem -nt $SLURM_NTASKS qchem_0185.in qchem_0185.out
qchem -nt $SLURM_NTASKS qchem_0186.in qchem_0186.out
qchem -nt $SLURM_NTASKS qchem_0187.in qchem_0187.out
qchem -nt $SLURM_NTASKS qchem_0188.in qchem_0188.out
qchem -nt $SLURM_NTASKS qchem_0189.in qchem_0189.out
qchem -nt $SLURM_NTASKS qchem_0190.in qchem_0190.out
qchem -nt $SLURM_NTASKS qchem_0191.in qchem_0191.out
qchem -nt $SLURM_NTASKS qchem_0192.in qchem_0192.out
qchem -nt $SLURM_NTASKS qchem_0193.in qchem_0193.out
qchem -nt $SLURM_NTASKS qchem_0194.in qchem_0194.out
qchem -nt $SLURM_NTASKS qchem_0195.in qchem_0195.out
qchem -nt $SLURM_NTASKS qchem_0196.in qchem_0196.out
qchem -nt $SLURM_NTASKS qchem_0197.in qchem_0197.out
qchem -nt $SLURM_NTASKS qchem_0198.in qchem_0198.out
qchem -nt $SLURM_NTASKS qchem_0199.in qchem_0199.out
qchem -nt $SLURM_NTASKS qchem_0200.in qchem_0200.out
qchem -nt $SLURM_NTASKS qchem_0201.in qchem_0201.out
qchem -nt $SLURM_NTASKS qchem_0202.in qchem_0202.out
qchem -nt $SLURM_NTASKS qchem_0203.in qchem_0203.out
qchem -nt $SLURM_NTASKS qchem_0204.in qchem_0204.out
qchem -nt $SLURM_NTASKS qchem_0205.in qchem_0205.out
qchem -nt $SLURM_NTASKS qchem_0206.in qchem_0206.out
qchem -nt $SLURM_NTASKS qchem_0207.in qchem_0207.out
qchem -nt $SLURM_NTASKS qchem_0208.in qchem_0208.out
qchem -nt $SLURM_NTASKS qchem_0209.in qchem_0209.out
qchem -nt $SLURM_NTASKS qchem_0210.in qchem_0210.out
qchem -nt $SLURM_NTASKS qchem_0211.in qchem_0211.out
qchem -nt $SLURM_NTASKS qchem_0212.in qchem_0212.out
qchem -nt $SLURM_NTASKS qchem_0213.in qchem_0213.out
qchem -nt $SLURM_NTASKS qchem_0214.in qchem_0214.out
qchem -nt $SLURM_NTASKS qchem_0215.in qchem_0215.out
qchem -nt $SLURM_NTASKS qchem_0216.in qchem_0216.out
qchem -nt $SLURM_NTASKS qchem_0217.in qchem_0217.out
qchem -nt $SLURM_NTASKS qchem_0218.in qchem_0218.out
qchem -nt $SLURM_NTASKS qchem_0219.in qchem_0219.out
qchem -nt $SLURM_NTASKS qchem_0220.in qchem_0220.out
qchem -nt $SLURM_NTASKS qchem_0221.in qchem_0221.out
qchem -nt $SLURM_NTASKS qchem_0222.in qchem_0222.out
qchem -nt $SLURM_NTASKS qchem_0223.in qchem_0223.out
qchem -nt $SLURM_NTASKS qchem_0224.in qchem_0224.out
qchem -nt $SLURM_NTASKS qchem_0225.in qchem_0225.out
qchem -nt $SLURM_NTASKS qchem_0226.in qchem_0226.out
qchem -nt $SLURM_NTASKS qchem_0227.in qchem_0227.out
qchem -nt $SLURM_NTASKS qchem_0228.in qchem_0228.out
qchem -nt $SLURM_NTASKS qchem_0229.in qchem_0229.out
qchem -nt $SLURM_NTASKS qchem_0230.in qchem_0230.out
qchem -nt $SLURM_NTASKS qchem_0231.in qchem_0231.out
qchem -nt $SLURM_NTASKS qchem_0232.in qchem_0232.out
qchem -nt $SLURM_NTASKS qchem_0233.in qchem_0233.out
qchem -nt $SLURM_NTASKS qchem_0234.in qchem_0234.out
qchem -nt $SLURM_NTASKS qchem_0235.in qchem_0235.out
qchem -nt $SLURM_NTASKS qchem_0236.in qchem_0236.out
qchem -nt $SLURM_NTASKS qchem_0237.in qchem_0237.out
qchem -nt $SLURM_NTASKS qchem_0238.in qchem_0238.out
qchem -nt $SLURM_NTASKS qchem_0239.in qchem_0239.out
qchem -nt $SLURM_NTASKS qchem_0240.in qchem_0240.out
qchem -nt $SLURM_NTASKS qchem_0241.in qchem_0241.out
qchem -nt $SLURM_NTASKS qchem_0242.in qchem_0242.out
qchem -nt $SLURM_NTASKS qchem_0243.in qchem_0243.out
qchem -nt $SLURM_NTASKS qchem_0244.in qchem_0244.out
qchem -nt $SLURM_NTASKS qchem_0245.in qchem_0245.out
qchem -nt $SLURM_NTASKS qchem_0246.in qchem_0246.out
qchem -nt $SLURM_NTASKS qchem_0247.in qchem_0247.out
qchem -nt $SLURM_NTASKS qchem_0248.in qchem_0248.out
qchem -nt $SLURM_NTASKS qchem_0249.in qchem_0249.out
qchem -nt $SLURM_NTASKS qchem_0250.in qchem_0250.out
qchem -nt $SLURM_NTASKS qchem_0251.in qchem_0251.out
qchem -nt $SLURM_NTASKS qchem_0252.in qchem_0252.out
qchem -nt $SLURM_NTASKS qchem_0253.in qchem_0253.out
qchem -nt $SLURM_NTASKS qchem_0254.in qchem_0254.out
qchem -nt $SLURM_NTASKS qchem_0255.in qchem_0255.out
qchem -nt $SLURM_NTASKS qchem_0256.in qchem_0256.out
qchem -nt $SLURM_NTASKS qchem_0257.in qchem_0257.out
qchem -nt $SLURM_NTASKS qchem_0258.in qchem_0258.out
qchem -nt $SLURM_NTASKS qchem_0259.in qchem_0259.out
qchem -nt $SLURM_NTASKS qchem_0260.in qchem_0260.out
qchem -nt $SLURM_NTASKS qchem_0261.in qchem_0261.out
qchem -nt $SLURM_NTASKS qchem_0262.in qchem_0262.out
qchem -nt $SLURM_NTASKS qchem_0263.in qchem_0263.out
qchem -nt $SLURM_NTASKS qchem_0264.in qchem_0264.out
qchem -nt $SLURM_NTASKS qchem_0265.in qchem_0265.out
qchem -nt $SLURM_NTASKS qchem_0266.in qchem_0266.out
qchem -nt $SLURM_NTASKS qchem_0267.in qchem_0267.out
qchem -nt $SLURM_NTASKS qchem_0268.in qchem_0268.out
qchem -nt $SLURM_NTASKS qchem_0269.in qchem_0269.out
qchem -nt $SLURM_NTASKS qchem_0270.in qchem_0270.out
qchem -nt $SLURM_NTASKS qchem_0271.in qchem_0271.out
qchem -nt $SLURM_NTASKS qchem_0272.in qchem_0272.out
qchem -nt $SLURM_NTASKS qchem_0273.in qchem_0273.out
qchem -nt $SLURM_NTASKS qchem_0274.in qchem_0274.out
qchem -nt $SLURM_NTASKS qchem_0275.in qchem_0275.out
qchem -nt $SLURM_NTASKS qchem_0276.in qchem_0276.out
qchem -nt $SLURM_NTASKS qchem_0277.in qchem_0277.out
qchem -nt $SLURM_NTASKS qchem_0278.in qchem_0278.out
qchem -nt $SLURM_NTASKS qchem_0279.in qchem_0279.out
qchem -nt $SLURM_NTASKS qchem_0280.in qchem_0280.out
qchem -nt $SLURM_NTASKS qchem_0281.in qchem_0281.out
qchem -nt $SLURM_NTASKS qchem_0282.in qchem_0282.out
qchem -nt $SLURM_NTASKS qchem_0283.in qchem_0283.out
qchem -nt $SLURM_NTASKS qchem_0284.in qchem_0284.out
qchem -nt $SLURM_NTASKS qchem_0285.in qchem_0285.out
qchem -nt $SLURM_NTASKS qchem_0286.in qchem_0286.out
qchem -nt $SLURM_NTASKS qchem_0287.in qchem_0287.out
qchem -nt $SLURM_NTASKS qchem_0288.in qchem_0288.out
qchem -nt $SLURM_NTASKS qchem_0289.in qchem_0289.out
qchem -nt $SLURM_NTASKS qchem_0290.in qchem_0290.out
qchem -nt $SLURM_NTASKS qchem_0291.in qchem_0291.out
qchem -nt $SLURM_NTASKS qchem_0292.in qchem_0292.out
qchem -nt $SLURM_NTASKS qchem_0293.in qchem_0293.out
qchem -nt $SLURM_NTASKS qchem_0294.in qchem_0294.out
qchem -nt $SLURM_NTASKS qchem_0295.in qchem_0295.out
qchem -nt $SLURM_NTASKS qchem_0296.in qchem_0296.out
qchem -nt $SLURM_NTASKS qchem_0297.in qchem_0297.out
qchem -nt $SLURM_NTASKS qchem_0298.in qchem_0298.out
qchem -nt $SLURM_NTASKS qchem_0299.in qchem_0299.out
qchem -nt $SLURM_NTASKS qchem_0300.in qchem_0300.out
qchem -nt $SLURM_NTASKS qchem_0301.in qchem_0301.out
qchem -nt $SLURM_NTASKS qchem_0302.in qchem_0302.out
qchem -nt $SLURM_NTASKS qchem_0303.in qchem_0303.out
qchem -nt $SLURM_NTASKS qchem_0304.in qchem_0304.out
qchem -nt $SLURM_NTASKS qchem_0305.in qchem_0305.out
qchem -nt $SLURM_NTASKS qchem_0306.in qchem_0306.out
qchem -nt $SLURM_NTASKS qchem_0307.in qchem_0307.out
qchem -nt $SLURM_NTASKS qchem_0308.in qchem_0308.out
qchem -nt $SLURM_NTASKS qchem_0309.in qchem_0309.out
qchem -nt $SLURM_NTASKS qchem_0310.in qchem_0310.out
qchem -nt $SLURM_NTASKS qchem_0311.in qchem_0311.out
qchem -nt $SLURM_NTASKS qchem_0312.in qchem_0312.out
qchem -nt $SLURM_NTASKS qchem_0313.in qchem_0313.out
qchem -nt $SLURM_NTASKS qchem_0314.in qchem_0314.out
qchem -nt $SLURM_NTASKS qchem_0315.in qchem_0315.out
qchem -nt $SLURM_NTASKS qchem_0316.in qchem_0316.out
qchem -nt $SLURM_NTASKS qchem_0317.in qchem_0317.out
qchem -nt $SLURM_NTASKS qchem_0318.in qchem_0318.out
qchem -nt $SLURM_NTASKS qchem_0319.in qchem_0319.out
qchem -nt $SLURM_NTASKS qchem_0320.in qchem_0320.out
qchem -nt $SLURM_NTASKS qchem_0321.in qchem_0321.out
qchem -nt $SLURM_NTASKS qchem_0322.in qchem_0322.out
qchem -nt $SLURM_NTASKS qchem_0323.in qchem_0323.out
qchem -nt $SLURM_NTASKS qchem_0324.in qchem_0324.out
qchem -nt $SLURM_NTASKS qchem_0325.in qchem_0325.out
qchem -nt $SLURM_NTASKS qchem_0326.in qchem_0326.out
qchem -nt $SLURM_NTASKS qchem_0327.in qchem_0327.out
qchem -nt $SLURM_NTASKS qchem_0328.in qchem_0328.out
qchem -nt $SLURM_NTASKS qchem_0329.in qchem_0329.out
qchem -nt $SLURM_NTASKS qchem_0330.in qchem_0330.out
qchem -nt $SLURM_NTASKS qchem_0331.in qchem_0331.out
qchem -nt $SLURM_NTASKS qchem_0332.in qchem_0332.out
qchem -nt $SLURM_NTASKS qchem_0333.in qchem_0333.out
qchem -nt $SLURM_NTASKS qchem_0334.in qchem_0334.out
qchem -nt $SLURM_NTASKS qchem_0335.in qchem_0335.out
qchem -nt $SLURM_NTASKS qchem_0336.in qchem_0336.out
qchem -nt $SLURM_NTASKS qchem_0337.in qchem_0337.out
qchem -nt $SLURM_NTASKS qchem_0338.in qchem_0338.out
qchem -nt $SLURM_NTASKS qchem_0339.in qchem_0339.out
qchem -nt $SLURM_NTASKS qchem_0340.in qchem_0340.out
qchem -nt $SLURM_NTASKS qchem_0341.in qchem_0341.out
qchem -nt $SLURM_NTASKS qchem_0342.in qchem_0342.out
qchem -nt $SLURM_NTASKS qchem_0343.in qchem_0343.out
qchem -nt $SLURM_NTASKS qchem_0344.in qchem_0344.out
qchem -nt $SLURM_NTASKS qchem_0345.in qchem_0345.out
qchem -nt $SLURM_NTASKS qchem_0346.in qchem_0346.out
qchem -nt $SLURM_NTASKS qchem_0347.in qchem_0347.out
qchem -nt $SLURM_NTASKS qchem_0348.in qchem_0348.out
qchem -nt $SLURM_NTASKS qchem_0349.in qchem_0349.out
qchem -nt $SLURM_NTASKS qchem_0350.in qchem_0350.out
qchem -nt $SLURM_NTASKS qchem_0351.in qchem_0351.out
qchem -nt $SLURM_NTASKS qchem_0352.in qchem_0352.out
qchem -nt $SLURM_NTASKS qchem_0353.in qchem_0353.out
qchem -nt $SLURM_NTASKS qchem_0354.in qchem_0354.out
qchem -nt $SLURM_NTASKS qchem_0355.in qchem_0355.out
qchem -nt $SLURM_NTASKS qchem_0356.in qchem_0356.out
qchem -nt $SLURM_NTASKS qchem_0357.in qchem_0357.out
qchem -nt $SLURM_NTASKS qchem_0358.in qchem_0358.out
qchem -nt $SLURM_NTASKS qchem_0359.in qchem_0359.out
qchem -nt $SLURM_NTASKS qchem_0360.in qchem_0360.out
qchem -nt $SLURM_NTASKS qchem_0361.in qchem_0361.out
qchem -nt $SLURM_NTASKS qchem_0362.in qchem_0362.out
qchem -nt $SLURM_NTASKS qchem_0363.in qchem_0363.out
qchem -nt $SLURM_NTASKS qchem_0364.in qchem_0364.out
qchem -nt $SLURM_NTASKS qchem_0365.in qchem_0365.out
qchem -nt $SLURM_NTASKS qchem_0366.in qchem_0366.out
qchem -nt $SLURM_NTASKS qchem_0367.in qchem_0367.out
qchem -nt $SLURM_NTASKS qchem_0368.in qchem_0368.out
qchem -nt $SLURM_NTASKS qchem_0369.in qchem_0369.out
qchem -nt $SLURM_NTASKS qchem_0370.in qchem_0370.out
qchem -nt $SLURM_NTASKS qchem_0371.in qchem_0371.out
qchem -nt $SLURM_NTASKS qchem_0372.in qchem_0372.out
qchem -nt $SLURM_NTASKS qchem_0373.in qchem_0373.out
qchem -nt $SLURM_NTASKS qchem_0374.in qchem_0374.out
qchem -nt $SLURM_NTASKS qchem_0375.in qchem_0375.out
qchem -nt $SLURM_NTASKS qchem_0376.in qchem_0376.out
qchem -nt $SLURM_NTASKS qchem_0377.in qchem_0377.out
qchem -nt $SLURM_NTASKS qchem_0378.in qchem_0378.out
qchem -nt $SLURM_NTASKS qchem_0379.in qchem_0379.out
qchem -nt $SLURM_NTASKS qchem_0380.in qchem_0380.out
qchem -nt $SLURM_NTASKS qchem_0381.in qchem_0381.out
qchem -nt $SLURM_NTASKS qchem_0382.in qchem_0382.out
qchem -nt $SLURM_NTASKS qchem_0383.in qchem_0383.out
qchem -nt $SLURM_NTASKS qchem_0384.in qchem_0384.out
qchem -nt $SLURM_NTASKS qchem_0385.in qchem_0385.out
qchem -nt $SLURM_NTASKS qchem_0386.in qchem_0386.out
qchem -nt $SLURM_NTASKS qchem_0387.in qchem_0387.out
qchem -nt $SLURM_NTASKS qchem_0388.in qchem_0388.out
qchem -nt $SLURM_NTASKS qchem_0389.in qchem_0389.out
qchem -nt $SLURM_NTASKS qchem_0390.in qchem_0390.out
qchem -nt $SLURM_NTASKS qchem_0391.in qchem_0391.out
qchem -nt $SLURM_NTASKS qchem_0392.in qchem_0392.out
qchem -nt $SLURM_NTASKS qchem_0393.in qchem_0393.out
qchem -nt $SLURM_NTASKS qchem_0394.in qchem_0394.out
qchem -nt $SLURM_NTASKS qchem_0395.in qchem_0395.out
qchem -nt $SLURM_NTASKS qchem_0396.in qchem_0396.out
qchem -nt $SLURM_NTASKS qchem_0397.in qchem_0397.out
qchem -nt $SLURM_NTASKS qchem_0398.in qchem_0398.out
qchem -nt $SLURM_NTASKS qchem_0399.in qchem_0399.out
qchem -nt $SLURM_NTASKS qchem_0400.in qchem_0400.out
qchem -nt $SLURM_NTASKS qchem_0401.in qchem_0401.out
qchem -nt $SLURM_NTASKS qchem_0402.in qchem_0402.out
qchem -nt $SLURM_NTASKS qchem_0403.in qchem_0403.out
qchem -nt $SLURM_NTASKS qchem_0404.in qchem_0404.out
qchem -nt $SLURM_NTASKS qchem_0405.in qchem_0405.out
qchem -nt $SLURM_NTASKS qchem_0406.in qchem_0406.out
qchem -nt $SLURM_NTASKS qchem_0407.in qchem_0407.out
qchem -nt $SLURM_NTASKS qchem_0408.in qchem_0408.out
qchem -nt $SLURM_NTASKS qchem_0409.in qchem_0409.out
qchem -nt $SLURM_NTASKS qchem_0410.in qchem_0410.out
qchem -nt $SLURM_NTASKS qchem_0411.in qchem_0411.out
qchem -nt $SLURM_NTASKS qchem_0412.in qchem_0412.out
qchem -nt $SLURM_NTASKS qchem_0413.in qchem_0413.out
qchem -nt $SLURM_NTASKS qchem_0414.in qchem_0414.out
qchem -nt $SLURM_NTASKS qchem_0415.in qchem_0415.out
qchem -nt $SLURM_NTASKS qchem_0416.in qchem_0416.out
qchem -nt $SLURM_NTASKS qchem_0417.in qchem_0417.out
qchem -nt $SLURM_NTASKS qchem_0418.in qchem_0418.out
qchem -nt $SLURM_NTASKS qchem_0419.in qchem_0419.out
qchem -nt $SLURM_NTASKS qchem_0420.in qchem_0420.out
qchem -nt $SLURM_NTASKS qchem_0421.in qchem_0421.out
qchem -nt $SLURM_NTASKS qchem_0422.in qchem_0422.out
qchem -nt $SLURM_NTASKS qchem_0423.in qchem_0423.out
qchem -nt $SLURM_NTASKS qchem_0424.in qchem_0424.out
qchem -nt $SLURM_NTASKS qchem_0425.in qchem_0425.out
qchem -nt $SLURM_NTASKS qchem_0426.in qchem_0426.out
qchem -nt $SLURM_NTASKS qchem_0427.in qchem_0427.out
qchem -nt $SLURM_NTASKS qchem_0428.in qchem_0428.out
qchem -nt $SLURM_NTASKS qchem_0429.in qchem_0429.out
qchem -nt $SLURM_NTASKS qchem_0430.in qchem_0430.out
qchem -nt $SLURM_NTASKS qchem_0431.in qchem_0431.out
qchem -nt $SLURM_NTASKS qchem_0432.in qchem_0432.out
qchem -nt $SLURM_NTASKS qchem_0433.in qchem_0433.out
qchem -nt $SLURM_NTASKS qchem_0434.in qchem_0434.out
qchem -nt $SLURM_NTASKS qchem_0435.in qchem_0435.out
qchem -nt $SLURM_NTASKS qchem_0436.in qchem_0436.out
qchem -nt $SLURM_NTASKS qchem_0437.in qchem_0437.out
qchem -nt $SLURM_NTASKS qchem_0438.in qchem_0438.out
qchem -nt $SLURM_NTASKS qchem_0439.in qchem_0439.out
qchem -nt $SLURM_NTASKS qchem_0440.in qchem_0440.out
 
echo 'Finishing program on'
date
 