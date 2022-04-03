#! /usr/bin/python
import xml.etree.ElementTree as ET
import numpy as np
eps3 = 1.0e-3
inf6 = 1.0e6
PRECIS=6
rydberg = 13.60569253

def search_tag(tmp,tag):
 outp=''
 for i in tmp:
  if tag in i:
   i2=i.replace(' ','').split(',')
   excluded=0
   for m in i2:
    if '!' in m:
     excluded=1
     break
    if tag in m:
     break
   if not excluded: 
    m2=m.split('=')
    outp=m2[-1].replace("'","")
 return outp

def read_input(pw_input):
 try: h=open(pw_input)
 except: raise NameError('pw input file '+pw_input+' does not exist')
 tmp=h.readlines()
 h.close()
 prefix=search_tag(tmp,'prefix')
 tmp_dir=search_tag(tmp,'outdir')
 return prefix,tmp_dir 

def read_data(prefix,tmp_dir):
  tree = ET.parse(tmp_dir+'/'+prefix+'.xml')
  root = tree.getroot()

  ENE=[]
  for i in root.findall('output/band_structure/fermi_energy'):
   ef=float(i.text.split()[0])*2
  nbnd=int(root.find('output/band_structure/nbnd').text.split()[0])
  nelec=float(root.find('output/band_structure/nelec').text.split()[0])

  print('EF='+str(ef)+', NELEC='+str(nelec))


  for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('eigenvalues'):
     ENE.append([float(m)*2 for m in  actor.text.split()])

  e,at=[],[]
  #reciprocal lattice vectors
  for i in root.findall('output/basis_set/reciprocal_lattice'):
   e.append([round(float(m),PRECIS) for m in i.find('b1').text.split()])
   e.append([round(float(m),PRECIS) for m in i.find('b2').text.split()])
   e.append([round(float(m),PRECIS) for m in i.find('b3').text.split()])
  e=np.round(np.array(e),PRECIS)


  # lattice vectors
  for i in root.findall('output/atomic_structure/cell'):
   at.append([round(float(m),PRECIS) for m in i.find('a1').text.split()])
   at.append([round(float(m),PRECIS) for m in i.find('a2').text.split()])
   at.append([round(float(m),PRECIS) for m in i.find('a3').text.split()])
  at=np.array(at)


  alat=float(root.find('output/atomic_structure').attrib["alat"])
  print('alat='+str(alat))


  # kpoints
  NONEQ,NONEQ_cryst=[],[]
  einv=np.linalg.inv(e)
  for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('k_point'):
     NONEQ.append([round(float(m),PRECIS) for m in  actor.text.split()])
     NONEQ_cryst.append(np.round(np.matmul(np.array(NONEQ[-1]),einv),PRECIS))

  print ('No of nonequivalent kpoints= '+str(len(NONEQ)))


  #symmetry operations
  SYMM_crystal=[]
  for neighbor in root.iter('rotation'):
     tmp=neighbor.text.split()
     tmp2=np.array([ [ float(m) for m in tmp[0:3]], [float(m) for m in tmp[3:6]], [float(m) for m in tmp[6:9]]])
     SYMM_crystal.append(np.transpose(tmp2))
  print ('No of symm. op.='+str(len(SYMM_crystal)))

  return alat,at,ef,nelec,len(SYMM_crystal),len(NONEQ),nbnd,SYMM_crystal,NONEQ_cryst,ENE
#alat,avec,efermi,nelec,nsym,nkpt,nbnd,rot,kpoint,energy

def main(argv = None):
    if argv is None:
        argv = sys.argv
    if len(argv) !=3:
        self = '/' + argv[0]
        self = self[len(self)-self[::-1].find('/'):]
        print("")
        print("    Converts the output of Quantum Espresso")
        print("    to the input of BoltzTraP 1.2.1. Written by Georgy Samsonidze,")
        print("    An Li, Daehyun Wee, Bosch Research (October 2011).")
        print("    Modified by Sylwia Gutowska (April 2022) [now data are read from xml file instead of output file,")
        print("    works only for PW calculations (Quantum Espresso)].")
        print("")
        print("    Usage: pw_input_file nbnd_exclude ")
        print("")
        print("  * pw_input_file = name of the pw input file, eg. dos.in, nscf.in etc.")
        print("  * nbnd_exclude = number of the lowest energy bands to exclude in the output")
        print("")
        print("    Creates files BoltzTraP.def, prefix.intrans, prefix.energy, prefix.struct.")
        print("    File names and parameters written to BoltzTraP.def and prefix.intrans")
        print("    are hard-coded into the script.")
        print("")
        return 1

 
    pw_input = argv[1]
    nbnd_exclude = int(argv[2])
    if nbnd_exclude < 0:
        print("\n    Error: invalid nbnd_exclude.\n")
        return 2
    
    prefix,tmp_dir=read_input(pw_input)
    alat,avec,efermi,nelec,nsym,nkpt,nbnd,rot,kpoint,energy=read_data(prefix,tmp_dir)

   
    fname_def = 'BoltzTraP.def'
    fname_intrans = prefix + '.intrans'
    fname_energy = prefix + '.energy'
    fname_struct = prefix + '.struct'

    deltae = 0.0005
    ecut = 0.4
    lpfac = 5
    efcut = 0.15
    tmax = 800.0
    deltat = 50.0
    ecut2 = -1.0
    dosmethod = 'TETRA'

    '''
    i = 0
    efermi_scf = (inf6 + eps3) / rydberg
    avec = []
    idxsym = []
    idxbnd = []
    spin = False
    for line in f_pw:
        if 'lattice parameter (alat)  =' in line:
            alat = float(line.split()[4])
        elif ' a(' in line:
            atext = line[23:57].split()
            avec.append([float(atext[0]) * alat, float(atext[1]) * alat, float(atext[2]) * alat])
        elif 'cryst.   s' in line:
            idxsym.append(i)
        elif 'the Fermi energy is' in line:
            efermi_scf = float(line.split()[4]) / rydberg
        elif 'highest occupied, lowest unoccupied level' in line:
            efermi_scf = (float(line.split()[6]) + float(line.split()[7])) / (2.0 * rydberg)
        elif 'number of electrons' in line:
            nelec = float(line.split()[4])
        elif 'Sym.Ops.' in line or 'Sym. Ops.' in line:
            nsym = int(line.split()[0])
        elif 'No symmetry found' in line:
	    nsym = 1
        elif 'number of k points=' in line:
            nkpt = int(line.split()[4])
        elif 'number of Kohn-Sham states=' in line:
            nbnd = int(line.split()[4])
        elif ' cryst. coord.' in line:
            idxkpt = i + 1
        elif 'band energies (ev)' in line or 'bands (ev)' in line:
            idxbnd.append(i + 2)
        elif 'SPIN' in line:
	    spin = True
        i += 1

    if abs(efermi) > inf6 / rydberg and (ftype_inp == 'pw' or ftype_inp == 'bands'):
        if abs(efermi_scf) > inf6 / rydberg:
            print("\n Error: Fermi energy not found.\n")
            return 2
        else:
            efermi = efermi_scf

    if spin:
        nelec -= nbnd_exclude
    else:
        nelec -= 2 * nbnd_exclude

    rot = []
    for ir in range(nsym):
        rot.append([])
        for i in range(3):
            rtext = f_pw[idxsym[ir] + i][19:53].split()
            rot[ir].append([int(rtext[0]), int(rtext[1]), int(rtext[2])])

    kpoint = []
    for ik in range(nkpt):
        ktext = f_pw[idxkpt + ik][20:56].split()
        kpoint.append([float(ktext[0]), float(ktext[1]), float(ktext[2])])

    if ftype_inp == 'pw':
        energy = []
        ncol = 8
        nrow = nbnd / ncol
        if nbnd % ncol != 0:
            nrow += 1
        for ik in range(nkpt):
            energy.append([])
            nelem = ncol
            for ir in range(nrow):
                etext = f_pw[idxbnd[ik] + ir].split()
                if ir == nrow - 1:
                    nelem = nbnd - ncol * (nrow - 1)
                for ie in range(nelem):
                    energy[ik].append(float(etext[ie]) / rydberg)
    elif ftype_inp == 'bands':
        energy = []
        ncol = 10
        nrow = nbnd / ncol
        if nbnd % ncol != 0:
            nrow += 1
        for ik in range(nkpt):
            energy.append([])
            nelem = ncol
            for ir in range(nrow):
                etext = f_inp[ik * (nrow + 1) + ir + 2].split()
                if ir == nrow - 1:
                    nelem = nbnd - ncol * (nrow - 1)
                for ie in range(nelem):
                    energy[ik].append(float(etext[ie]) / rydberg)
    elif ftype_inp == 'inteqp':
        nhead = 2
        ntot = len(f_inp) - nhead
        bndmin = int(f_inp[nhead].split()[1]) - 1
        bndmax = int(f_inp[nhead + ntot - 1].split()[1]) - 1
        nbnd = bndmax + 1
        energy = []
        for ik in range(nkpt):
            energy.append([])
            for ib in range(bndmin):
                energy[ik].append(0.0)
            for ib in range(bndmax - bndmin + 1):
                energy[ik].append(float(f_inp[nhead + ik + ib * nkpt].split()[6]) / rydberg)
    '''
    f_def = '5, \'' + prefix + '.intrans\',      \'old\',    \'formatted\',0\n'
    f_def += '6,\'' + prefix + '.outputtrans\',      \'unknown\',    \'formatted\',0\n'
    f_def += '20,\'' + prefix + '.struct\',         \'old\',    \'formatted\',0\n'
    f_def += '10,\'' + prefix + '.energy\',         \'old\',    \'formatted\',0\n'
    f_def += '48,\'' + prefix + '.engre\',         \'unknown\',    \'unformatted\',0\n'
    f_def += '49,\'' + prefix + '.transdos\',        \'unknown\',    \'formatted\',0\n'
    f_def += '50,\'' + prefix + '.sigxx\',        \'unknown\',    \'formatted\',0\n'
    f_def += '51,\'' + prefix + '.sigxxx\',        \'unknown\',    \'formatted\',0\n'
    f_def += '21,\'' + prefix + '.trace\',           \'unknown\',    \'formatted\',0\n'
    f_def += '22,\'' + prefix + '.condtens\',           \'unknown\',    \'formatted\',0\n'
    f_def += '24,\'' + prefix + '.halltens\',           \'unknown\',    \'formatted\',0\n'
    f_def += '30,\'' + prefix + '_BZ.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '31,\'' + prefix + '_fermi.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '32,\'' + prefix + '_sigxx.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '33,\'' + prefix + '_sigyy.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '34,\'' + prefix + '_sigzz.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '35,\'' + prefix + '_band.dat\',           \'unknown\',    \'formatted\',0\n'
    f_def += '36,\'' + prefix + '_band.gpl\',           \'unknown\',    \'formatted\',0\n'
    f_def += '37,\'' + prefix + '_deriv.dat\',           \'unknown\',    \'formatted\',0\n'
    f_def += '38,\'' + prefix + '_mass.dat\',           \'unknown\',    \'formatted\',0\n'

    f = open(fname_def, 'w')
    f.write(f_def)
    f.close()

    f_intrans = 'GENE                      # Format of DOS\n'
    f_intrans += '0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap\n'
    f_intrans += str(efermi) + ' ' + str(deltae) + ' ' + str(ecut) + ' ' + str(nelec) + '    # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons\n'
    f_intrans += 'CALC                      # CALC (calculate expansion coeff), NOCALC read from file\n'
    f_intrans += str(lpfac) + '                         # lpfac, number of latt-points per k-point\n'
    f_intrans += 'BOLTZ                     # run mode (only BOLTZ is supported)\n'
    f_intrans += str(efcut) + '                      # (efcut) energy range of chemical potential\n'
    f_intrans += str(tmax) + ' ' + str(deltat) + '                # Tmax, temperature grid\n'
    f_intrans += str(ecut2) + '                      # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)\n'
    f_intrans += dosmethod + '\n'

    f = open(fname_intrans, 'w')
    f.write(f_intrans)
    f.close()

    f_energy = prefix + '\n'
    f_energy += str(nkpt) + '\n'
    for ik in range(nkpt):
        f_energy += str(kpoint[ik][0]) + ' ' + str(kpoint[ik][1]) + ' ' + str(kpoint[ik][2]) + ' ' + str(nbnd - nbnd_exclude) + '\n'
        for ib in range(nbnd_exclude, nbnd):
            f_energy += str(energy[ik][ib]) + '\n'

    f = open(fname_energy, 'w')
    f.write(f_energy)
    f.close()

    f_struct = prefix + '\n'
    for i in range(3):
        f_struct += str(avec[i][0]) + ' ' + str(avec[i][1]) + ' ' + str(avec[i][2]) + '\n'
    f_struct += str(nsym) + '\n'
    for ir in range(nsym):
        f_struct += str(rot[ir][0][0]) + ' ' + str(rot[ir][1][0]) + ' ' + str(rot[ir][2][0]) + ' '
        f_struct += str(rot[ir][0][1]) + ' ' + str(rot[ir][1][1]) + ' ' + str(rot[ir][2][1]) + ' '
        f_struct += str(rot[ir][0][2]) + ' ' + str(rot[ir][1][2]) + ' ' + str(rot[ir][2][2]) + '\n'

    f = open(fname_struct, 'w')
    f.write(f_struct)
    f.close()

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
