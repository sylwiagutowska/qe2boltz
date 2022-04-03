# qe2boltz
It is modified program qe2boltz, originally written by Georgy Samsonidze, An Li, Daehyun Wee, Bosch Research (source:  https://blog.levilentz.com/boltztrap-tutorial-for-quantum-espresso/). 

It allows to convert results obtained from Quantum Espresso (QE) to format required by Bolztrap-1.2.5 (https://www.imc.tuwien.ac.at/forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap/).

Original qe2boltz script reads data from output file and requirs verbosity='high' tag in input file.

My modification allows to read data from  prefix.xml file saved automatically by QE in output_dir instead of output file, so the verbosity='high' tag is not required. It is useful if you have already done some calculations and did not use this verbosity flag.

Usage: python qe2boltz.py pw_input_file nbnd_excluded
where:
1) pw_input_file is input file used in calculations, eg. nscf.in or dos.in
2) nbnd_excluded - the bands of index < nbnd_excluded will be omitted in calculations (I usually use 0 for safety reasons :) )
