      INCLUDE 'penelope.f'  ! File included to simplify compilation.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      M    M    AA   TTTTTT  EEEEEE  RRRRR   IIII    AA    L          C
C      MM  MM   A  A    TT    E       R    R   II    A  A   L          C
C      M MM M  A    A   TT    E       R    R   II   A    A  L          C
C      M    M  AAAAAA   TT    EEEE    RRRRR    II   AAAAAA  L          C
C      M    M  A    A   TT    E       R   R    II   A    A  L          C
C      M    M  A    A   TT    EEEEEE  R    R  IIII  A    A  LLLLLL     C
C                                                                      C
C                                                   (version 2005).    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2005)                                     C
C  Copyright (c) 2001-2005                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided "as is" without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  This program generates material definition files for PENELOPE, which
C  contain tables of physical properties, interaction cross sections and
C  other particle transport characteristics. These data are extracted
C  from the database, which consists of the following 465 ascii files,
C  *'pdatconf.tab': Atomic ground state configurations, ionization
C      energies and central values of the one-electron shell Compton
C      profiles for the elements, from hydrogen to uranium.
C  *'pdcompos.tab': composition data for 279 different materials.
C  *'pdeflist.tab': list of materials included in the 'pdcompos.tab'
C       file, with their identification numbers.
C  *'pdrelax.tab': data on atomic relaxation, extracted from the LLNL
C       Evaluated Atomic Data Library.
C  *92 files named 'pdeelZZ.tab' with ZZ=atomic number (01-92). These
C       files contain electron and positron elastic scattering data.
C       The same grid of energies is used for all elements.
C  *92 files named 'pdebrZZ.tab' that contain electron bremsstrahlung
C       data. These files were generated from the database of Seltzer
C       and Berger. The same grid of energies for all elements.
C  *'pdbrang.tab': parameters of the intrinsic angular distribution of
C       bremsstrahlung photons. Determined by fitting the set of bench-
C       mark partial-wave shape functions tabulated by Kissel, Quarles
C       and Pratt.
C  *92 files named 'pdgppZZ.tab' with cross sections for pair production
C       in the field of neutral atoms (sum of pair and triplet contri-
C       butions), obtained from the XCOM program of Berger and Hubbell.
C       The same energy grid for all elements.
C  *92 files named 'pdgphZZ.tab', containing total atomic photoelectric
C       cross sections and partial cross sections for inner (K and L)
C       shells, generated from the EPDL97 data library of Cullen et al.
C  *92 files named 'pdeinZZ.tab' (ZZ=01-92) that contain cross sections
C       for ionization of K shells and L subshells by electron and pos-
C       itron impact.
C
C  A material is completely characterized by its chemical composition,
C  i.e. elements present and number of atoms of each element in a
C  molecule (=stoichiometric index), mass density and mean excitation
C  energy. Alloys and mixtures are treated as compounds, with stoichio-
C  metric indices equal or proportional to the percent number of atoms
C  of each element. Information about the material is supplied by the
C  user from the keyboard, following the prompts from 'material', or
C  read from the 'pdcompos.tab' file, which contains information for 279
C  different materials. In the case of compounds, 'molecular' cross
C  sections are obtained by means of the additivity rule, i.e. as the
C  sum of the atomic cross sections.
C
C  To obtain the executable file 'material.exe', compile and link the
C  source files 'material.f' and 'penelope.f'. The database files
C  and the executable code 'material.exe' must be in the same directory.
C
C  NOTE: In the output file and in the simulation program, lengths are
C  given in cm and energies in eV. Consequently, total cross sections
C  are in cm**2, stopping powers in eV/cm, etc. However, macroscopic
C  cross sections listed in the input/output files are expressed in mass
C  thickness units, 1 mtu = 1 g/cm**2.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CALL PEMATW
      END
