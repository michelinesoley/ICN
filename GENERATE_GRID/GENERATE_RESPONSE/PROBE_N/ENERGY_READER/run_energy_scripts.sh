### reads energy from QCHEM output
./energy_reader.sh

### separate contribution of first ground-state calculation from second MOM calculation 
python energy_cleaner.py

### plots energy cuts
python energy_plot.py

### performs selection of energies based on threshold 
python energy_selector.py

### plots energy cuts
python energy_selector_plot.py

### open plots 
open energy_*.ps 

