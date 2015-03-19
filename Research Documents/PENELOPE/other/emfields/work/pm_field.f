      INCLUDE 'penelope.f'  ! Files included to simplify compilation.
      INCLUDE 'pengeom.f'
      INCLUDE 'penvared.f'
      INCLUDE 'timer.f'
      INCLUDE 'penfield.f'

C>>>>>>>>> EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  This source file is a modified version of program PENMAIN that
C  simulates radiation transport in matter with static electro-
C  magnetic fields. The portions of code that differ from PENMAIN are
C  marked with comment lines like the ones limiting this note.
C
C  The electromagnetic field is assumed to be constant over the entire
C  geometry. To modify the field strength, edit the present file (lines
C  1788-1814).
C
C  The present example of main program is intended only to generate a
C  file with a small number of electron or positron trajectories, e.g.
C  for plotting. For systematic simulations, the user has to edit the
C  PENMAIN.F source file and modify it appropriately.
C<<<<<<<<< EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C         PPPPP   EEEEEE  N    N  M    M    AA    IIII  N    N         C
C         P    P  E       NN   N  MM  MM   A  A    II   NN   N         C
C         P    P  E       N N  N  M MM M  A    A   II   N N  N         C
C         PPPPP   EEEE    N  N N  M    M  AAAAAA   II   N  N N         C
C         P       E       N   NN  M    M  A    A   II   N   NN         C
C         P       EEEEEE  N    N  M    M  A    A  IIII  N    N         C
C                                                                      C
C                                                   (version 2005).    C
C                                                                      C
C  This program performs Monte Carlo simulation of electron-photon     C
C  showers in material structures described with the constructive      C
C  quadric geometry package 'PENGEOM.F'.                               C
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
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  In the default mode, the program assumes that primary particles of a
C  given kind, KPARP, are emitted from a point source, with coordinates
C  (SX0,SY0,SZ0), either with fixed energy SE0 or with a specified
C  (histogram-like) energy spectrum. The initial direction of the
C  primary particles is sampled uniformly within a cone of semi-aperture
C  SALPHA and central axis in the direction (STHETA, SPHI). Thus, the
C  case SALPHA=0 defines a monodirectional source and SALPHA=180 deg
C  corresponds to an isotropic source. When SX0=SY0=0 and STHETA=0 or
C  180 deg, the source is axially symmetrical about the Z-axis.
C
C  Alternatively, the program can read the initial state variables of
C  'primary' particles from a pre-calculated phase-space file. This
C  option may be used to split the simulation of complex problems into
C  several consecutive stages.
C
C  >>>>>>>> NOTE: All energies and lengths are given in eV and cm,
C                 respectively.
C
C
C  ************  Structure of the input data file.
C
C  Each line in the input data file consists of a 6-character keyword
C  (columns 1-6) followed either by numerical data (in free format) or
C  by a character string, which start at the 8th column. Keywords are
C  explicitly used/verified by the program (which is case sensitive!).
C  Notice also that the order of the data lines is important. The
C  keyword '______' (6 blanks) indicates comment lines, these can be
C  placed anywhere in the input file. The program ignores any text
C  following the first blank after the last numerical datum, or after
C  the character string, in each line (thus, in the table given below,
C  the comments in square brackets are ignored by the program). Lines
C  with some keywords (e.g., 'SPECTR', 'IPSFN') can appear an arbitrary
C  number of times, limited only by the allocated amount of memory.
C
C  The program assigns default values to many input variables; lines
C  that declare default values may be eliminated from the input file.
C
C
C  The structure of the input file is the following,
C
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C  TITLE  Title of the job, up to 120 characters.
C
C         >>>>>>>> Source definition.
C  SKPAR  KPARP    [Primary particles: 1=electron, 2=photon, 3=positron]
C  SENERG SE0              [Initial energy (monoenergetic sources only)]
C  SPECTR Ei,Pi                 [E bin: lower-end and total probability]
C  SPOSIT SX0,SY0,SZ0                        [Coordinates of the source]
C  SDIREC STHETA,SPHI               [Beam axis direction angles, in deg]
C  SAPERT SALPHA                                 [Beam aperture, in deg]
C
C         >>>>>>>> Input phase-space file (psf).
C  IPSFN  psf_filename.ext               [Input psf name, 20 characters]
C  IPSPLI NSPLIT                                      [Splitting number]
C  EPMAX  EPMAX                 [Maximum energy of particles in the psf]
C
C         >>>>>>>> Material data and simulation parameters.
C  NMAT   NMAT                   [Number of different materials, .le.10]
C  SIMPAR M,EABS(1:3,M),C1,C2,WCC,WCR   [Sim. parameters for material M]
C  PFNAME mat_filename.ext          [Material definition file, 20 chars]
C
C         >>>>>>>> Geometry definition file.
C  GEOMFN geo_filename.ext          [Geometry definition file, 20 chars]
C  DSMAX  IBODY,DSMAX(IBODY)   [IB, maximum step length (cm) in body IB]
C
C         >>>>>>>> Interaction forcing.
C  IFORCE KB,KPAR,ICOL,FORCER,WLOW,WHIG            [Interaction forcing]
C
C         >>>>>>>> Emerging particles. Energy and angular distributions.
C  NBE    EMIN,EMAX,NBE              [E-interval and no. of energy bins]
C  NBTH   NBTH                   [No. of bins for the polar angle THETA]
C  NBPH   NBPH                 [No. of bins for the azimuthal angle PHI]
C
C         >>>>>>>> Impact detectors (up to 5 different detectors).
C  IMPDET EDIL,EDIU,NCHI,IPSF  [Energy window, no. of channels and IPSF]
C  IDPSF  pm_psf_impdet_#.dat           [Output psf file name, 20 chars]
C  IDSPC  pm_spc_impdet_#.dat      [Output spectrum file name, 20 chars]
C  IDBODY KB                       [Active body; one line for each body]
C  IDKPAR KPAR               [Kind of detected particles, one line each]
C
C         >>>>>>>> Energy-deposition detectors (up to 5).
C  ENDDET EDEL,EDEU,NCHE          [Energy window and number of channels]
C  EDSPC  pm_spc_enddet_#.dat      [Output spectrum file name, 20 chars]
C  EDBODY KB                       [Active body; one line for each body]
C
C         >>>>>>>> Dose distribution.
C  GRIDX  XL,XU                [X coordinates of the enclosure vertices]
C  GRIDY  YL,YU                [Y coordinates of the enclosure vertices]
C  GRIDZ  ZL,ZU                [Z coordinates of the enclosure vertices]
C  GRIDBN NX,NY,NZ                                     [Numbers of bins]
C
C         >>>>>>>> Job properties.
C  RESUME dumpfile_1.dat          [Resume from this dump file, 20 chars]
C  DUMPTO dumpfile_2.dat             [Generate this dump file, 20 chars]
C  DUMPP  DUMPP                                 [Dumping period, in sec]
C
C  NSIMSH DSHN                     [Desired number of simulated showers]
C  RSEED  ISEED1,ISEED2           [Seeds of the random number generator]
C  TIME   TIMEA                       [Allotted simulation time, in sec]
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C
C
C  The following listing describes the function of each of the keywords,
C  the accompanying data and their default values. For clarity, blanks
C  in keywords are indicated as '_'.
C
C  TITLE_ : Title of the job (up to 120 characters).
C             DEFAULT: none (the input file must start with this line)
C
C           The TITLE string is used to mark dump files. To prevent the
C           improper use of wrong resuming files, change the title each
C           time you modify basic parameters of your problem. The code
C           will then be able to identify the inconsistency and to print
C           an error message before stopping.
C
C  >>>>>>>> Source definition.
C
C  SKPAR_ : Kind of primary particle (1=electrons, 2=photons or
C           3=positrons).
C             DEFAULT: KPARP=1
C
C  SENERG : For a monoenergetic source, initial energy SE0 of primary
C           particles.
C             DEFAULT: SE0=1.0E6
C
C  SPECTR : For a source with continuous (stepwise constant) spectrum,
C           each 'SPECTR' line gives the lower end-point of an energy
C           bin of the source spectrum (Ei) and the associated relative
C           probability (Pi), integrated over the bin. Up to NSEM=200
C           lines, in arbitrary order. The upper end of the spectrum is
C           defined by entering a line with Ei equal to the upper energy
C           value and with Pi set to a negative value.
C             DEFAULT: none
C
C  SPOSIT : Coordinates of the (point) source.
C             DEFAULT: SX0=SY0=0, SZ0=-1.0E15
C  SDIREC : Polar and azimuthal angles of the source beam axis direc-
C           tion, in deg.
C             DEFAULTS: STHETA=0.0, SPHI=0.0
C  SAPERT : Angular aperture of the source beam, in deg.
C             DEFAULT: SALPHA=0.0
C
C           --> Notice that the default source is a pencil beam that
C           moves upwards along the Z-axis.
C
C  >>>>>>>> Input phase-space file.
C
C  The initial state variables of primary particles can be read directly
C  from a set of pre-calculated phase-space files (psf). When this
C  option is active, previous definitions about the source are ignored.
C
C  IPSFN_ : Name of an input psf (up to 20 characters).
C             DEFAULT: none
C           Up to 100 psf's may be declared. They are read sequentially.
C
C  The input psf is in ASCII format. Each line defines the initial state
C  of a particle; it contains the following quantities in free format
C  (and in the order they are listed here):
C    -- KPAR, type of particle (1, electron; 2, photon; 3, positron).
C    -- E, energy (eV).
C    -- X,Y,Z, position coordinates (cm).
C    -- U,V,W, direction cosines.
C    -- WGHT, weight.
C    -- ILB(1),ILB(2),ILB(3),ILB(4), a set of indices that provide
C           information on how the particle was generated (see the file
C           'manual.txt').
C    -- NSHI, incremental shower number (difference between the shower
C           numbers of the present particle and the one preceding it
C           in the psf).
C  Phase-space files can be generated by running PENMAIN using an impact
C  detector with the flag IPSF=1 or -1 (see below).
C
C  Because of the limited size of the psf's, the results of analogue
C  simulations tend to be 'too noisy'. This can be partially corrected
C  by splitting the particles from the psf.
C
C  IPSPLI : Splitting number. Each particle in the psf's will be split
C           into NSPLIT equivalent particles, with weights equal to
C           WGHT/NSPLIT.
C             DEFAULT: NSPLIT=1 (no splitting)
C
C  --> WARNING: Notice that there is a 'latent' uncertainty in the psf,
C  which sets a limit to the accuracy that can be attained by using
C  large splitting numbers.
C
C  EPMAX_ : Maximum energy (in eV) of particles in the psf's.
C           EPMAX is the upper limit of the energy interval covered by
C           the simulation lookup tables. To minimise interpolation
C           errors, EPMAX should not be much larger than the maximum
C           energy actually occurring during the simulation.
C
C           When the initial state variables of particles are read from
C           a psf, this parameter is required to initialise PENELOPE and
C           is critical; the code crashes if it finds a particle that
C           has energy larger than EPMAX.
C             DEFAULT: EPMAX=1.0E9 (interpolation is not optimal)
C
C  >>>>>>>> Material data and simulation parameters.
C
C  NMAT__ : Number of different materials. The original programs in the
C           distribution package allow up to 10 materials. This number
C           can be increased by changing the value of the parameter
C           MAXMAT in the original source files.
C           Materials are identified by their ordering in PENELOPE's
C           input material data file.
C             DEFAULT: NMAT=1
C
C  SIMPAR : Set of simulation parameters for material M; absorption
C           energies, EABS(1:3,M), elastic scattering parameters, C1(M)
C           and C2(M), and cutoff energy losses for inelastic collisions
C           and bremsstrahlung emission, WCC(M) and WCR(M). One line for
C           each material.
C             DEFAULTS: EABS(1,M)=EABS(3,M)=0.01*EPMAX,
C                       EABS(2,M)=0.001*EPMAX
C                       C1(M)=C2(M)=0.1, WCC=EABS(1,M), WCR=EABS(2,M)
C             EPMAX is the upper limit of the energy interval covered
C             by the simulation lookup tables.
C
C  PFNAME : Name of the PENELOPE input material data file (up to 20
C           characters).
C             DEFAULT: none
C
C  >>>>>>>> Geometry.
C
C  GEOMFN : PENGEOM geometry definition file name (a string of up to
C           20 characters).
C             DEFAULT: none.
C
C           --> The geometry definition file can be debugged/visualised
C           with the viewers GVIEW2D and GVIEW3D (operable only under
C           Windows).
C
C  DSMAX_ : Maximum step length DSMAX(IB) (in cm) of electrons and
C           positrons in body IB. This parameter is important only for
C           thin bodies; it should be given a value of the order of one
C           tenth of the body thickness or less.
C             DEFAULT: DSMAX=1.0E20 (no step length control)
C
C  >>>>>>>> Interaction forcing.
C
C  IFORCE : Activates forcing of interactions of type ICOL of particles
C           KPAR in body KB. FORCER is the forcing factor, which must
C           be larger than unity. WLOW and WHIG are the lower and upper
C           limits of the weight window where interaction forcing is
C           applied.
C             DEFAULT: no interaction forcing
C
C           If the mean free path for real interactions of type ICOL is
C           MFP, the program will simulate interactions of this type
C           (real or forced) with an effective mean free path equal to
C           MFP/FORCER.
C
C           TRICK: a negative input value of FORCER, -FN, is  assumed to
C           mean that each particle should interact, on average and
C           approximately, +FN times in a path length equal to the range
C           of that kind of particle with E=EPMAX. This is very useful
C           to generate x-ray spectra from bulk samples.
C
C  The real effect of interaction forcing on the efficiency is not easy
C  to predict. Please, do tentative runs with different FORCER values
C  and check the efficiency gain (or loss!).
C
C  >>>>>>>> Energy and angular distributions of emerging particles.
C
C  NBE___ : Limits EMIN and EMAX of the interval where energy
C           distributions of emerging particles are tallied. Number of
C           energy bins.
C             DEFAULT: EMIN=0.0, EMAX=EPMAX, NBE=100
C
C  NBTH__ : Number of bins for the polar angle THETA.
C             DEFAULT: NBTH=90
C
C           WARNING: In the output files, the terms 'transmitted' and
C           'backscattered' are used to denote particles that leave the
C           material system moving upwards (THETA>0) and downwards
C           (THETA<0), respectively. Notice that this agrees with the
C           usual meaning of these terms only when primary particles
C           impinge on the system coming from below (i.e. with THETA>0).
C
C  NBPH__ : Number of bins for the azimuthal angle PHI.
C             DEFAULT: NBPH=1 (azimuthal average).
C
C  >>>>>>>> Impact detectors.
C
C  Each impact detector consists of a set of active bodies, which must
C  have been defined as parts of the geometry. The output spectrum is
C  the energy distribution of particles that entered any of the active
C  bodies coming from a body that is not active (i.e. that is not part
C  of the detector). Notice that a detected particle can re-enter the
C  detector volume and, consequently, be 'counted' several times (except
C  when the flag IPSF is set equal to -1, see below).
C
C  Active bodies cannot be void, because the geometry routines would not
C  stop particles at its limiting surfaces. In case you need to define
C  detectors outside the material system, fill them with an arbitrary
C  material of very small density to avoid perturbing the transport
C  process.
C
C  To define each impact detector, insert the following block of lines;
C
C  IMPDET : Starts the definition of a new detector. Up to 5 different
C           detectors can be considered.
C           EDIL and EDIU are the lower and upper limits of the energy
C           window covered by the impact detector.
C           NCHI is the number of channels in the output spectrum of
C           the detector (.LE. 1000).
C           The integer flag IPSF serves to activate the creation of a
C           phase-space file (psf), which contains the state variables
C           of all particles that enter the detector. Use this option
C           with care, because psf's may grow very fast.
C           IPSF=0, the psf is not created.
C           IPSF=1, a psf is created. Particles that enter the detector
C                are transported as usually.
C           IPSF=-1, a psf is created. The simulation of a detected
C                particle is discontinued when it enters the detector.
C             DEFAULTS: None
C
C  IDPSF_ : Name of the output phase-space file (up to 20 characters).
C             DEFAULT: 'pm_psf_impdet_#.dat'
C
C  IDSPC_ : Name of the output energy spectrum file (20 characters).
C             DEFAULT: 'pm_spc_impdet_#.dat'
C
C  IDBODY : Active body of the detector. One line for each active body.
C             DEFAULT: None
C           --> Notice that a body cannot be part of more than one
C           impact detector.
C
C  IDKPAR : Kind of particle that is detected (1=electrons, 2=photons or
C           3=positrons). One line for each kind.
C             DEFAULT: All particles are detected
C
C  >>>>>>>> Energy-deposition detectors.
C
C  Each energy-deposition detector consists of a set of active bodies,
C  which must have been defined as parts of the geometry. The output
C  spectrum is the distribution of absorbed energy (per primary shower)
C  in the active bodies.
C
C           *** WARNING: The energy-deposition spectrum may be strongly
C           biased when interaction forcing is applied.
C
C  To define each energy-deposition detector insert the following block
C  of lines;
C
C  ENDDET : Starts the definition of a new energy-deposition detector.
C           Up to 5 different detectors can be considered.
C           EDEL and EDEU are the lower and upper limits of the energy
C           window covered by the detector.
C           NCHE is the number of energy channels in the output spectrum
C           (.LE. 1000).
C
C  EDSPC_ : Name of the output spectrum file (up to 20 characters).
C             DEFAULT: 'pm_spc_enddet_#.dat'
C
C  EDBODY : Active body of the detector. One line for each active body.
C             DEFAULT: None
C           --> Notice that a body cannot be part of more than one
C           energy-deposition detector.
C
C  >>>>>>>> Dose map.
C
C  The program can calculate the dose distribution inside a parallele-
C  piped (dose enclosure) whose edges are parallel to the axes of the
C  laboratory frame. The enclosure is defined by giving the coordinates
C  of its vertices. The dose is tallied using a uniform orthogonal grid
C  with NX, NY and NZ (.LE. 100) bins along the directions of the
C  coordinate axes.
C
C  GRIDX_ : X-coordinates of the vertices of the dose enclosure.
C             DEFAULT: None
C  GRIDY_ : Y-coordinates of the vertices of the dose enclosure.
C             DEFAULT: None
C  GRIDZ_ : Z-coordinates of the vertices of the dose enclosure.
C             DEFAULT: None
C  GRIDBN : Numbers of bins NX, NY, and NZ in the X, Y and Z directions,
C           respectively.
C             DEFAULTS: NX=1, NY=1, NZ=1
C
C  --> The grid defined here to calculate the dose distribution can be
C  used to tally other 3D distributions (e.g. the space distribution of
C  inner-shell ionisation events, used in electron-probe microanalysis).
C  This, however, requires to edit the present source file.
C
C  >>>>>>>> Job properties.
C
C  RESUME : The program will read the dump file named `dumpfile_1.dat'
C           (20 characters) and resume the simulation from the point
C           where it was left. Use this option very, _VERY_ carefully.
C           Make sure that the input data file is fully consistent with
C           the one used to generate the dump file.
C             DEFAULT: off
C
C  DUMPTO : Generate a dump file named 'dumpfile_2.dat' (20 characters)
C           after completing the simulation run. This allows resuming
C           the simulation later on to improve statistics.
C             DEFAULT: off
C
C           NOTE: If the file 'dumpfile_2.dat' already exists, it is
C           overwritten.
C
C  DUMPP_ : When the DUMPTO option is activated, simulation results are
C           evaluated and written in the output files every DUMPP
C           seconds. This option is useful to check the progress of long
C           simulations. It also allows running the program with a long
C           execution time and stopping it when the required statistical
C           uncertainty has been reached.
C             DEFAULT: DUMPP=1.0E15
C
C  NSIMSH : Desired number of simulated showers.
C             DEFAULT: DSHN=2.0E9
C
C  RSEED_ : Seeds of the random number generator.
C             DEFAULT: ISEED1=12345; ISEED2=54321
C
C  TIME__ : Allotted simulation time, in sec.
C             DEFAULT: TIMEA=100.0E0
C
C  The program is aborted when an incorrect input datum is found. The
C  conflicting quantity usually appears in the last line of the output
C  file 'penmain.dat'. If the trouble is with arrays having dimensions
C  smaller than required, the program indicates how the problem can be
C  solved (this usually requires editing the source file, be careful).
C
C  The clock subroutine (TIMER) may have to be adapted to your specific
C  computer-compiler configuration; standard FORTRAN 77 does not provide
C  timing tools. However, the routines in module TIMER.F do work for
C  many FORTRAN compilers.
C
C  ************  Generating the executable PENMAIN and running it.
C
C  To generate the executable binary file PENMAIN.EXE, compile and link
C  the FORTRAN 77 source files PENMAIN.F, PENELOPE.F, PENGEOM.F,
C  PENVARED.F and TIMER.F. For example, if you are using the g77
C  compiler under Windows, place these five files in the same directory,
C  open a DOS window and from that directory enter the command
C    `g77 -Wall -O PENMAIN.F -o PENMAIN.EXE'
C  (The same, with file names in lowercase, should work under Linux).
C
C  To run PENMAIN, you have to generate an input data file, let's call
C  it PENMAIN.IN, and the corresponding geometry definition and material
C  data files. Place these three files and the binary file PEMAIN.EXE in
C  the same directory and, from there, issue the command
C    `PENMAIN.EXE < PENMAIN.IN'
C
C  The calculated distributions are written in separate files, whose
C  names start with 'pm_' (for PenMain) and have the extension '.dat'.
C  These files are in a format suited for direct visualisation with
C  GNUPLOT (version 4.0).
C
C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LIT
      CHARACTER*10 LCH10
      CHARACTER*20 PFILE,OFILE,PSFI,PSFDIO,SPCDIO,SPCDEO,PFILED,PFILER
      CHARACTER*23 DATE23
      CHARACTER*120 TITLE,TITLE2,BUFFER,BUF2
      CHARACTER*200 PSFREC
C
      CHARACTER*6 KWORD,
     1  KWTITL,KWKPAR,KWSENE,KWSPEC,  KWSPOS,KWSDIR,KWSAPE,KWPSFN,
     1  KWPSPL,KWEMAX,KWNMAT,KWSIMP,  KWMATF,KWGEOM,KWSMAX,KWIFOR,
     1  KWNBE ,KWNBTH,KWNBPH,KWIDET,  KWIPSF,KWISPC,KWIBOD,KWIPAR,
     1  KWEDET,KWESPC,KWEBOD,KGRDXX,  KGRDYY,KGRDZZ,KGRDBN,KWRESU,
     1  KWDUMP,KWDMPP,KWNSIM,KWRSEE,  KWTIME,KWCOMM
      PARAMETER(
     1  KWTITL='TITLE ',KWKPAR='SKPAR ',KWSENE='SENERG',KWSPEC='SPECTR',
     1  KWSPOS='SPOSIT',KWSDIR='SDIREC',KWSAPE='SAPERT',KWPSFN='IPSFN ',
     1  KWPSPL='IPSPLI',KWEMAX='EPMAX ',KWNMAT='NMAT  ',KWSIMP='SIMPAR',
     1  KWMATF='PFNAME',KWGEOM='GEOMFN',KWSMAX='DSMAX ',KWIFOR='IFORCE',
     1  KWNBE ='NBE   ',KWNBTH='NBTH  ',KWNBPH='NBPH  ',KWIDET='IMPDET',
     1  KWIPSF='IDPSF ',KWISPC='IDSPC ',KWIBOD='IDBODY',KWIPAR='IDKPAR',
     1  KWEDET='ENDDET',KWESPC='EDSPC ',KWEBOD='EDBODY',KGRDXX='GRIDX ',
     1  KGRDYY='GRIDY ',KGRDZZ='GRIDZ ',KGRDBN='GRIDBN',KWRESU='RESUME',
     1  KWDUMP='DUMPTO',KWDMPP='DUMPP ',KWNSIM='NSIMSH',KWRSEE='RSEED ',
     1  KWTIME='TIME  ',KWCOMM='      ')
C
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI,
     1  RA2DE=180.0D0/PI, DE2RA=PI/180.0D0)
C
C  ****  Source energy spectrum.
C
      PARAMETER (NSEBM=200)
      DIMENSION ES(NSEBM),PTS(NSEBM),IAS(NSEBM),FS(NSEBM),SHIST(NSEBM)
      DATA SHIST/NSEBM*0/
C
      PARAMETER (NPSFM=100)
      DIMENSION PSFI(NPSFM)
C
C  ****  Main-PENELOPE commons.
C
      PARAMETER(MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      DIMENSION RHOI(MAXMAT)
C  ****  Geometry.
      PARAMETER (NS=2500,NB=1250,NX=150)
      COMMON/QTREE/NBODY,MATER(NB),KMOTH(NB),KDGHT(NB,NX),
     1    KSURF(NB,NX),KFLAG(NB,NX),KALIAS(NS),KSLAST
      COMMON/QKDET/KDET(NB)
      DIMENSION PARINP(20)
      DIMENSION DSMAX(NB)

C>>>>>>>>> EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ****  EM field.
      COMMON/UFIELD/EFX,EFY,EFZ,BFX,BFY,BFZ
C<<<<<<<<< EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C  ****  Interaction forcing parameters.
      DIMENSION IFORCE(NB,3),WLOW(NB,3),WHIG(NB,3)
      COMMON/CFORCE/FORCE(NB,3,8)
C
C  ************  Discrete counters.
C
      DIMENSION
     1  PRIM(3),PRIM2(3),DPRIM(3),     ! Numbers of IEXIT particles.
     1  SEC(3,3),SEC2(3,3),DSEC(3,3)   ! Generated secondary particles.
      DATA PRIM,PRIM2,SEC,SEC2/24*0.0D0/
      DIMENSION WSEC(3,3),WSEC2(3,3)
C  ----  Deposited energies in various bodies.
      DIMENSION TDEBO(NB),TDEBO2(NB),DEBO(NB)
      DATA TDEBO,TDEBO2/NB*0.0D0,NB*0.0D0/
C
C  ************  Continuous distributions.
C
C  ----  Energy distributions of emerging particles.
      PARAMETER (NBEM=100,NBE24=24*NBEM)
      DIMENSION BSE(3)
      DIMENSION PDE(3,2,NBEM),PDE2(3,2,NBEM),PDEP(3,2,NBEM),
     1          LPDE(3,2,NBEM)
      DATA PDE,PDE2,PDEP,LPDE/NBE24*0.0D0/
C
C  ----  Angular distributions of emerging particles.
      PARAMETER (NBTHM=90,NBPHM=60)
      PARAMETER (NBPH12=12*NBTHM*NBPHM)
      DIMENSION PDA(3,NBTHM,NBPHM),PDA2(3,NBTHM,NBPHM),
     1          PDAP(3,NBTHM,NBPHM),LPDA(3,NBTHM,NBPHM)
      DATA PDA,PDA2,PDAP,LPDA/NBPH12*0.0D0/
C
C  ----  Impact detectors (up to NIDM different detectors).
      PARAMETER (NIDM=5,NIDCM=1000,NIDCMT=4*NIDM*NIDCM)
      DIMENSION PSFDIO(NIDM),SPCDIO(NIDM)
      DIMENSION DIT(NIDM,NIDCM),DIT2(NIDM,NIDCM),DITP(NIDM,NIDCM),
     1  LDIT(NIDM,NIDCM),TDID(NIDM),TDID2(NIDM),DEDI(NIDM),IPSF(NIDM)
      DATA DIT,DIT2,DITP,LDIT/NIDCMT*0.0D0/
      DIMENSION EDIL(NIDM),EDIU(NIDM),BDIE(NIDM),NDICH(NIDM)
      DIMENSION KBDI(NB),KKDI(NIDM,3),RLAST(NIDM),RWRITE(NIDM)
C
C  ----  Energy-deposition detectors (up to NEDM different detectors).
      PARAMETER (NEDM=5,NEDCM=1000,NEDCMT=2*NEDM*NEDCM)
      DIMENSION SPCDEO(NEDM)
      DIMENSION DET(NEDM,NEDCM),DET2(NEDM,NEDCM),TDED(NEDM),
     1  TDED2(NEDM),DEDE(NEDM)
      DATA DET,DET2/NEDCMT*0.0D0/
      DIMENSION EDEL(NEDM),EDEU(NEDM),BDEE(NEDM),NDECH(NEDM)
      DIMENSION KBDE(NB)
C
C  ----  Dose distribution (up to NDXM bins along each coord. axis).
      PARAMETER (NDXM=100,NDYM=100,NDZM=100,NDXMT=4*NDXM*NDYM*NDZM)
      DIMENSION DOSE(NDXM,NDYM,NDZM),DOSE2(NDXM,NDYM,NDZM),
     1          DOSEP(NDXM,NDYM,NDZM),LDOSE(NDXM,NDYM,NDZM)
      DATA DOSE,DOSE2,DOSEP,LDOSE/NDXMT*0.0D0/
      DIMENSION DXL(3),DXU(3),BDOSE(3),BDOSER(3),NDB(3)
C
      EXTERNAL RAND
C
C  ****  Time counter initiation.
C
      CALL TIME0
C
C  ------------------------  Read input data file.
C
      OPEN(26,FILE='penmain.dat')
      WRITE(26,1000)
 1000 FORMAT(//3X,44('*'),/3X,'**   Program PENMAIN. ',
     1 ' Input data file.   **',/3X,44('*'))
C
      CALL PDATET(DATE23)
      WRITE(26,1001) DATE23
 1001 FORMAT(/3X,'Date and time: ',A23)
C
C  ****  Title.
C
      READ(5,'(A6,1X,A120)') KWORD,TITLE
      WRITE(26,'(/3X,A120)') TITLE
C
      ISORP=0
      ISPEC=0
      NPSF=0
      NPSN=0
      NSEB=1
C
C  ************  Description of a point source.
C
   11 CONTINUE
      READ(5,'(A6,1X,A120)') KWORD,BUFFER
      IF(KWORD.EQ.KWCOMM) GO TO 11
      IF(KWORD.EQ.KWPSFN) GO TO 21
C
      WRITE(26,1100)
 1100 FORMAT(//3X,72('-'),/3X,'>>>>>>  Source description.')
C
      IF(KWORD.EQ.KWKPAR) THEN
        READ(BUFFER,*) KPARP
   12   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 12
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        KPARP=1
      ENDIF
      IF(KPARP.LT.1.OR.KPARP.GT.3) KPARP=1
      IF(KPARP.EQ.1) WRITE(26,1110)
 1110 FORMAT(3X,'Primary particles: electrons')
      IF(KPARP.EQ.2) WRITE(26,1111)
 1111 FORMAT(3X,'Primary particles: photons')
      IF(KPARP.EQ.3) WRITE(26,1112)
 1112 FORMAT(3X,'Primary particles: positrons')
C
C  ****  Monoenergetic source.
C
      IF(KWORD.EQ.KWSENE) THEN
        READ(BUFFER,*) E0
        WRITE(26,1120) E0
 1120   FORMAT(3X,'Initial energy = ',1P,E13.6,' eV')
   13   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 13
        IF(KWORD.EQ.KWPSFN) GO TO 21
C
C  ****  Continuous energy spectrum.
C
      ELSE IF(KWORD.EQ.KWSPEC) THEN
        ISPEC=1
        NSEB=0
   14   CONTINUE
        NSEB=NSEB+1
        READ(BUFFER,*) ES(NSEB),PTS(NSEB)
        PTS(NSEB)=MAX(PTS(NSEB),0.0D0)
   15   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 15
        IF(KWORD.EQ.KWSPEC) GO TO 14
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        E0=1.0D6
        WRITE(26,1120) E0
      ENDIF
      IF(ISPEC.EQ.1) THEN
        IF(NSEB.GT.1) THEN
          CALL SORT2(ES,PTS,NSEB)
          WRITE(26,1121)
 1121     FORMAT(/3X,'Spectrum:',7X,'I',4X,'E_low(eV)',4x,'E_high(eV)',
     1      5X,'P_sum(E)',/16X,45('-'))
          DO I=1,NSEB-1
            WRITE(26,'(16X,I4,1P,3E14.6)') I,ES(I),ES(I+1),PTS(I)
          ENDDO
          WRITE(26,*) '  '
          E0=ES(NSEB)
          NSEB=NSEB-1
          CALL IRND0(PTS,FS,IAS,NSEB)
        ELSE
          WRITE(26,*) 'The source energy spectrum is not defined.'
          STOP 'The source energy spectrum is not defined.'
        ENDIF
      ENDIF
      IF(E0.LT.50.0D0) THEN
        WRITE(26,*) 'The initial energy E0 is too small.'
        STOP 'The initial energy E0 is too small.'
      ENDIF
      EPMAX=E0
C  ----  Positrons eventually give annihilation gamma-rays. The maximum
C        energy of annihilation photons is .lt. 1.21*(E0+me*c**2).
      IF(KPARP.EQ.3) EPMAX=1.21D0*(E0+5.12D5)
C
C  ****  Position of the point source.
C
      IF(KWORD.EQ.KWSPOS) THEN
        READ(BUFFER,*) SX0,SY0,SZ0
   16   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 16
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        SX0=0.0D0
        SY0=0.0D0
        SZ0=-1.0D15
      ENDIF
      WRITE(26,1132) SX0,SY0,SZ0
 1132 FORMAT(3X,'Coordinates of centre:     SX0 = ',1P,E13.6,
     1  ' cm',/30X,'SY0 = ',E13.6,' cm',/30X,'SZ0 = ',E13.6,' cm')
C
C  ****  Angular distribution of primary particles.
C
      IF(KWORD.EQ.KWSDIR) THEN
        READ(BUFFER,*) STHETA,SPHI
   17   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 17
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        STHETA=0.0D0
        SPHI=0.0D0
      ENDIF
      WRITE(26,1133) STHETA,SPHI
 1133 FORMAT(3X,'Beam direction angles:   THETA = ',1P,E13.6,' deg',
     1  /30X,'PHI = ',E13.6,' deg')
C
      IF(KWORD.EQ.KWSAPE) THEN
        READ(BUFFER,*) SALPHA
   18   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 18
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        SALPHA=0.0D0
      ENDIF
      WRITE(26,1134) SALPHA
 1134 FORMAT(3X,'Beam aperture:',11X,'ALPHA = ',1P,E13.6,' deg')
      CALL GCONE0(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA)
C
C  ************  Particle state variables read from a phase-space file.
C
   21 CONTINUE
      IF(KWORD.EQ.KWPSFN) THEN
        NPSF=NPSF+1
        IF(NPSF.EQ.1) THEN
          WRITE(26,1200)
 1200     FORMAT(//3X,72('-'),/3X,'>>>>>>  Input phase-space files.'/)
        ENDIF
        IF(NPSF.GT.NPSFM) THEN
          WRITE(26,'(/3X,''Too many phase-space files.'')')
          STOP 'Too many phase-space files'
        ENDIF
        READ(BUFFER,'(A20)') PSFI(NPSF)
        WRITE(26,1201) NPSF,PSFI(NPSF)
 1201   FORMAT(3X,'Phase-space file #',I4,': ',A20)
        ISORP=1
   22   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 22
        IF(KWORD.EQ.KWPSFN) GO TO 21
        IF(KWORD.EQ.KWPSPL) THEN
          READ(BUFFER,*) NSPLIT
          IF(NSPLIT.GT.1.AND.NSPLIT.LT.50) THEN
            WRITE(26,1210) NSPLIT
 1210       FORMAT(/3X,'Particle splitting number = ',I3)
          ELSE
            NSPLIT=1
          ENDIF
        ENDIF
   23   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 23
        IF(KWORD.EQ.KWEMAX) THEN
          READ(BUFFER,*) EPMAX
          WRITE(26,1220) EPMAX
 1220     FORMAT(3X,'Maximum particle energy = ',1P,E13.6,' eV')
   24     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 24
        ELSE
          EPMAX=1.0D9
          WRITE(26,1220) EPMAX
          WRITE(26,'(3X,''WARNING: You should have specified the '',
     1      ''maximum energy EPMAX.'')')
        ENDIF
      ENDIF
C
      IF(KWORD.EQ.KWEMAX) THEN
        READ(BUFFER,*) EPMAX
        WRITE(26,1220) EPMAX
   25   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 25
      ENDIF
C
C  ************  Material data and simulation parameters.
C
      WRITE(26,1300)
 1300 FORMAT(//3X,72('-'),/
     1  3X,'>>>>>>  Material data and simulation parameters.')
C
      IF(KWORD.EQ.KWNMAT) THEN
        READ(BUFFER,*) NMAT
   31   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 31
      ELSE
        WRITE(26,*) 'You have to specify the number of materials.'
        STOP 'You have to specify the number of materials.'
      ENDIF
      WRITE(26,1310) NMAT
 1310 FORMAT(3X,'Number of different materials = ',I2)
      IF(NMAT.LT.1.OR.NMAT.GT.MAXMAT) THEN
        WRITE(26,*) 'Wrong number of materials.'
        STOP 'Wrong number of materials.'
      ENDIF
C
C  ****  Simulation parameters.
C
      DO M=1,NMAT
        EABS(1,M)=0.010D0*EPMAX
        EABS(2,M)=0.001D0*EPMAX
        EABS(3,M)=0.010D0*EPMAX
        C1(M)=0.10D0
        C2(M)=0.10D0
        WCC(M)=EABS(1,M)
        WCR(M)=EABS(2,M)
      ENDDO
      DO IB=1,NB
        DSMAX(IB)=1.0D20
      ENDDO
      WGHT0=1.0D0   ! Primary particle weight.
      INTFOR=0      ! No interaction forcing as default.
C
      IF(KWORD.EQ.KWSIMP) THEN
        READ(BUFFER,*) M
        IF(M.LT.1.OR.M.GT.NMAT) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect material number.'
          STOP 'Incorrect material number.'
        ENDIF
        READ(BUFFER,*) M,EABS(1,M),EABS(2,M),EABS(3,M),C1(M),C2(M),
     1    WCC(M),WCR(M)
   32   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 32
        IF(KWORD.EQ.KWSIMP) THEN
          READ(BUFFER,*) M
          IF(M.LT.1.OR.M.GT.NMAT) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect material number.'
            STOP 'Incorrect material number.'
          ENDIF
          READ(BUFFER,*) M,EABS(1,M),EABS(2,M),EABS(3,M),C1(M),C2(M),
     1      WCC(M),WCR(M)
          GO TO 32
        ENDIF
      ENDIF
C
      DO M=1,NMAT
        IF(M.EQ.1) LIT='st'
        IF(M.EQ.2) LIT='nd'
        IF(M.EQ.3) LIT='rd'
        IF(M.GT.3) LIT='th'
        WRITE(26,1320) M,LIT
 1320   FORMAT(/3X,'**** ',I2,A2,' material')
        IF(EABS(1,M).LT.1.0D2) EABS(1,M)=1.0D2
        IF(EABS(2,M).LT.1.0D2) EABS(2,M)=1.0D2
        IF(EABS(3,M).LT.1.0D2) EABS(3,M)=1.0D2
        WRITE(26,1321) EABS(1,M)
 1321 FORMAT(3X,'Electron absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1322) EABS(2,M)
 1322 FORMAT(3X,'  Photon absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1323) EABS(3,M)
 1323 FORMAT(3X,'Positron absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1324) C1(M),C2(M),WCC(M),WCR(M)
 1324 FORMAT(3X,'Electron-positron simulation parameters:',
     1 /4X,'C1 =',1P,E13.6,',      C2 =',E13.6,/3X,'Wcc =',E13.6,
     1  ' eV,  Wcr =',E13.6,' eV')
      ENDDO
C
C  ****  Initialisation of PENELOPE.
C
      IF(KWORD.EQ.KWMATF) THEN
        READ(BUFFER,'(A20)') PFILE
        WRITE(26,1330) PFILE
 1330   FORMAT(/3X,'PENELOPE''s material definition file: ',A20)
   33   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 33
      ELSE
        WRITE(26,*) 'You have to specify a material file.'
        STOP 'You have to specify a material file.'
      ENDIF
C
      OPEN(15,FILE=PFILE,IOSTAT=KODE)
      IF(KODE.NE.0) THEN
        WRITE(26,'(''File '',A20,'' could not be opened.'')') PFILE
        STOP
      ENDIF
      OPEN(16,FILE='pm_material.dat')
        INFO=3
        CALL PEINIT(EPMAX,NMAT,15,16,INFO)
      CLOSE(UNIT=15)
      CLOSE(UNIT=16)
C  ----  Inverse densities are used to score the local dose.
      DO M=1,NMAT
        RHOI(M)=1.0D0/RHO(M)
      ENDDO
C
C  ************  Geometry definition.
C
C  Define here the geometry parameters that are to be altered, if any.
C     PARINP(1)=
C     PARINP(2)=  ...
      NPINP=0
C
      IF(KWORD.EQ.KWGEOM) THEN
        READ(BUFFER,'(A20)') PFILE
        WRITE(26,1340) PFILE
 1340   FORMAT(/3X,'PENGEOM''s geometry definition file: ',A20)
        OPEN(15,FILE=PFILE,IOSTAT=KODE)
        IF(KODE.NE.0) THEN
          WRITE(26,'(''File '',A20,'' could not be opened.'')') PFILE
          STOP
        ENDIF
        OPEN(16,FILE='pm_geometry.rep')
        CALL GEOMIN(PARINP,NPINP,NMATG,NBODY,15,16)
        CLOSE(UNIT=15)
        CLOSE(UNIT=16)
        IF(NMATG.LT.1) THEN
          WRITE(26,*) 'NMATG must be greater than 0.'
          STOP 'NMATG must be greater than 0.'
        ENDIF
C
        IF(NBODY.GT.NB) THEN
          WRITE(26,'(/6X,''Too many bodies.'')')
          STOP 'Too many bodies.'
        ENDIF
C
        IF(NMATG.GT.NMAT) THEN
          WRITE(26,'(/6X,''Too many different materials.'')')
          STOP 'Too many different materials.'
        ENDIF
C
   34   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 34
      ELSE
        WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
        WRITE(26,*) 'You have to specify a geometry file.'
        STOP 'You have to specify a geometry file.'
      ENDIF
C
C  ****  Maximum step lengths of electrons and positrons.
C
      IF(KWORD.EQ.KWSMAX) THEN
        READ(BUFFER,*) IB
        IF(IB.LT.1.OR.IB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect body number.'
          STOP 'Incorrect body number.'
        ENDIF
        READ(BUFFER,*) IB,DSMAX(IB)
        IF(DSMAX(IB).LT.1.0D-7) DSMAX(IB)=1.0D20
   35   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 35
        IF(KWORD.EQ.KWSMAX) THEN
          READ(BUFFER,*) IB
          IF(IB.LT.1.OR.IB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body number.'
            STOP 'Incorrect body number.'
          ENDIF
          READ(BUFFER,*) IB,DSMAX(IB)
          IF(DSMAX(IB).LT.1.0D-7) DSMAX(IB)=1.0D20
          GO TO 35
        ENDIF
      ENDIF
C
      WRITE(26,1350)
 1350 FORMAT(//3X,72('-'),/3X,'>>>>>>  Maximum allowed step lengths of',
     1  ' electrons and positrons.')
      DO IB=1,NBODY
        WRITE(26,1351) IB,DSMAX(IB)
 1351   FORMAT(3X,'* Body =',I4,',   DSMAX = ',1P,E13.6,' cm')
      ENDDO
C
C  ************  Variance reduction (only interaction forcing).
C
      DO KB=1,NB
        DO ICOL=1,8
          DO KPAR=1,3
            FORCE(KB,KPAR,ICOL)=1.0D0
          ENDDO
        ENDDO
        DO KPAR=1,3
          IFORCE(KB,KPAR)=0
          WLOW(KB,KPAR)=0.0D0
          WHIG(KB,KPAR)=1.0D6
        ENDDO
      ENDDO
      IFORON=0
C
      IF(KWORD.EQ.KWIFOR) THEN
        WRITE(26,1400)
 1400   FORMAT(//3X,72('-'),/
     1    3X,'>>>>>>  Interaction forcing: FORCE(IBODY,KPAR,ICOL)')
   41   CONTINUE
        READ(BUFFER,*) KB,KPAR,ICOL,FORCER,WWLOW,WWHIG
C  ****  Negative FORCER values are re-interpreted, as described in the
C        heading comments above.
        IF(FORCER.LT.-1.0D-6) THEN
          MM=MATER(KB)
          EVENTS=MAX(ABS(FORCER),1.0D0)
          PLT=PRANGE(E0,KPAR,MM)
          RMFP=PHMFP(E0,KPAR,MM,ICOL)
          FORCER=EVENTS*RMFP/PLT
        ENDIF
        IF(WWLOW.LT.1.0D-6) WWLOW=1.0D-6
        IF(WWHIG.GT.1.0D6) WWHIG=1.0D6
        IF(KB.LT.1.OR.KB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Inconsistent KB value.'
          STOP 'Inconsistent KB value.'
        ENDIF
        IF(KPAR.LT.1.OR.KPAR.GT.3) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of KPAR.'
          STOP 'Incorrect value of KPAR.'
        ENDIF
        IF(ICOL.LT.1.OR.ICOL.GT.8) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of ICOL.'
          STOP 'Incorrect value of ICOL.'
        ENDIF
        WLOW(KB,KPAR)=MAX(WLOW(KB,KPAR),WWLOW)
        WHIG(KB,KPAR)=MIN(WHIG(KB,KPAR),WWHIG)
        IF(WLOW(KB,KPAR).GT.WHIG(KB,KPAR)) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect weight window limits.'
          STOP 'Incorrect weight window limits.'
        ENDIF
        IF(FORCER.LT.1.0D0) STOP 'FORCER must be greater than unity.'
        IFORCE(KB,KPAR)=1
        IFORON=1
        FORCE(KB,KPAR,ICOL)=FORCER
        WRITE(26,1410) KB,KPAR,ICOL,FORCER,WLOW(KB,KPAR),WHIG(KB,KPAR)
 1410   FORMAT(3X,'FORCE(',I4,',',I1,',',I1,') =',1P,E13.6,
     1    ',  weight window = (',E9.2,',',E9.2,')')
   42   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 42
        IF(KWORD.EQ.KWIFOR) GO TO 41
      ENDIF
C
C  ************  Energy and angular distributions of emerging
C                particles.
C
      WRITE(26,1500)
 1500 FORMAT(//3X,72('-'),/
     1  3X,'>>>>>>  Energy and angular distributions of emerging',
     1  ' particles.')
C
      IF(KWORD.EQ.KWNBE) THEN
        READ(BUFFER,*) EMIN,EMAX,NBE
   51   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 51
      ELSE
        EMIN=0.0D0
        EMAX=EPMAX
        NBE=NBEM
      ENDIF
      EMIN=MAX(EMIN,0.0D0)
      WRITE(26,1510) NBE,EMIN,EMAX
 1510 FORMAT(3X,'E:       NBE = ',I3,
     1  ',  EMIN =',1P,E13.6,' eV,  EMAX =',E13.6,' eV')
      EMAX=MAX(EMAX,EPMAX)
      IF(EMAX.LT.EPMAX) THEN
        WRITE(26,*) '   WARNING: EMAX is less than EPMAX.'
      ENDIF
      IF(NBE.LT.1) THEN
        WRITE(26,*) 'Wrong number of energy bins.'
        STOP 'Wrong number of energy bins.'
      ELSE IF (NBE.GT.NBEM) THEN
        WRITE(26,*) 'NBE is too large.'
        WRITE(26,*) 'Set the parameter NBEM equal to ',NBE
        STOP 'NBE is too large.'
      ENDIF
      IF(EMIN.GT.EMAX) THEN
        WRITE(26,*) 'Energy interval is too narrow.'
        STOP 'Energy interval is too narrow.'
      ENDIF
C
      IF(KWORD.EQ.KWNBTH) THEN
        READ(BUFFER,*) NBTH
   52   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 52
      ELSE
        NBTH=NBTHM
      ENDIF
      WRITE(26,1520) NBTH
 1520 FORMAT(3X,'Theta:  NBTH = ',I3)
      IF(NBTH.LT.1) THEN
        WRITE(26,*) 'Wrong number of THETA bins.'
        STOP 'Wrong number of THETA bins.'
      ELSE IF (NBTH.GT.NBTHM) THEN
        WRITE(26,*) 'NBTH is too large.'
        WRITE(26,*) 'Set the parameter NBTHM equal to ',NBTH
        STOP 'NBTH is too large.'
      ENDIF
C
      IF(KWORD.EQ.KWNBPH) THEN
        READ(BUFFER,*) NBPH
   53   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 53
      ELSE
        NBPH=NBPHM
      ENDIF
      WRITE(26,1530) NBPH
 1530 FORMAT(3X,'Phi:    NBPH = ',I3)
      IF(NBTH.LT.1) THEN
        WRITE(26,*) 'Wrong number of PHI bins.'
        STOP 'Wrong number of PHI bins.'
      ELSE IF (NBPH.GT.NBPHM) THEN
        WRITE(26,*) 'NBPH is too large.'
        WRITE(26,*) 'Set the parameter NBPHM equal to ',NBPH
        STOP 'NBPH is too large.'
      ENDIF
C
C  ****  Bin sizes.
C
C  ----  The factor 1.0000001 serves to ensure that the upper limit of
C  the tallied interval is within the last channel (otherwise, the array
C  dimensions could be exceeded).
      BSE(1)=1.0000001D0*(EMAX-EMIN)/DBLE(NBE)
      BSE(2)=1.0000001D0*(EMAX-EMIN)/DBLE(NBE)
      BSE(3)=1.0000001D0*(EMAX-EMIN)/DBLE(NBE)
      BSTH=1.0000001D0*180.0D0/DBLE(NBTH)
      BSPH=1.0000001D0*360.0D0/DBLE(NBPH)
C
C  ************  Impact detectors.
C
      DO KB=1,NBODY+1
        KBDI(KB)=0
      ENDDO
      DO KD=1,NIDM
        KKDI(KD,1)=0
        KKDI(KD,2)=0
        KKDI(KD,3)=0
        TDID(KD)=0.0D0
        TDID2(KD)=0.0D0
      ENDDO
C
      NDIDEF=0
   61 CONTINUE
      IF(KWORD.EQ.KWIDET) THEN
        NDIDEF=NDIDEF+1
        IF(NDIDEF.GT.NIDM) THEN
          WRITE(26,'(3X,''NDIDEF = '',I4)') NDIDEF
          WRITE(26,*) 'Too many detectors.'
          STOP 'Too many detectors.'
        ENDIF
        WRITE(26,1600) NDIDEF
 1600   FORMAT(//3X,72('-'),/
     1    3X,'>>>>>>  Impact detector #', I2)
        READ(BUFFER,*) EDIL(NDIDEF),EDIU(NDIDEF),NDICH(NDIDEF),
     1    IPSF(NDIDEF)
        IF(EDIL(NDIDEF).LT.0.0D0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'EDIL must be positive.'
          STOP 'EDIL must be positive.'
        ENDIF
        IF(EDIU(NDIDEF).LT.EDIL(NDIDEF)) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect energy limits.'
          STOP 'Incorrect energy limits.'
        ENDIF
        IF(NDICH(NDIDEF).LT.10.OR.NDICH(NDIDEF).GT.NIDCM) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect number of channels.'
          STOP 'Incorrect number of channels.'
        ENDIF
        WRITE(26,1610) EDIL(NDIDEF),EDIU(NDIDEF),NDICH(NDIDEF)
 1610   FORMAT(3X,'Energy window = (',1P,E12.5,',',E12.5,') eV, no.',
     1    ' of channels = ',I4)
        BDIE(NDIDEF)=1.0000001D0*(EDIU(NDIDEF)-EDIL(NDIDEF))
     1           /DBLE(NDICH(NDIDEF))
C
   62   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 62
        IF(IPSF(NDIDEF).EQ.1.OR.IPSF(NDIDEF).EQ.-1) THEN
          IF(KWORD.EQ.KWIPSF) THEN
            READ(BUFFER,'(A20)') PSFDIO(NDIDEF)
            WRITE(26,1611) PSFDIO(NDIDEF),IPSF(NDIDEF)
 1611       FORMAT(3X,'Output phase-space file: ',A20,' (IPSF=',I2,')')
   63       CONTINUE
            READ(5,'(A6,1X,A120)') KWORD,BUFFER
            IF(KWORD.EQ.KWCOMM) GO TO 63
          ELSE
            WRITE(BUF2,'(I5)') NDIDEF
            PSFDIO(NDIDEF)='pm_psf_impdet_'//BUF2(5:5)//'.dat'
            WRITE(26,1611) PSFDIO(NDIDEF),IPSF(NDIDEF)
          ENDIF
        ENDIF
C
        IF(KWORD.EQ.KWISPC) THEN
          READ(BUFFER,'(A20)') SPCDIO(NDIDEF)
          WRITE(26,1612) SPCDIO(NDIDEF)
 1612     FORMAT(3X,'Output energy spectrum: ',A20)
   64     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 64
        ELSE
          WRITE(BUF2,'(I5)') NDIDEF
          SPCDIO(NDIDEF)='pm_spc_impdet_'//BUF2(5:5)//'.dat'
          WRITE(26,1612) SPCDIO(NDIDEF)
        ENDIF
C
   65   CONTINUE
        IF(KWORD.EQ.KWIBOD) THEN
          READ(BUFFER,*) KB
          IF(KB.LT.0.OR.KB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body label.'
            STOP 'Incorrect body label.'
          ENDIF
          IF(KBDI(KB).NE.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A body cannot be part of two detectors.'
            STOP 'A body cannot be part of two detectors.'
          ENDIF
          WRITE(26,1620) KB
 1620     FORMAT(3X,'Active body = ',I4)
          KBDI(KB)=NDIDEF
          KDET(KB)=NDIDEF
   66     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 66
          IF(KWORD.EQ.KWIBOD) GO TO 65
          IF(KWORD.EQ.KWIDET) THEN
            ITST=MAX(KKDI(NDIDEF,1),KKDI(NDIDEF,2),KKDI(NDIDEF,3))
            IF(ITST.EQ.0) THEN
              KKDI(NDIDEF,1)=1
              KKDI(NDIDEF,2)=1
              KKDI(NDIDEF,3)=1
              WRITE(26,1630)
 1630         FORMAT(3X,'Detected particles = electrons, photons and ',
     1          'positrons')
            ENDIF
            GO TO 61
          ENDIF
        ENDIF
C
   67   CONTINUE
        IF(KWORD.EQ.KWIPAR) THEN
          READ(BUFFER,*) KPARD
          IF(KPARD.EQ.1) THEN
            KKDI(NDIDEF,1)=1
            WRITE(26,1631)
 1631       FORMAT(3X,'Detected particles = electrons')
          ELSE IF(KPARD.EQ.2) THEN
            KKDI(NDIDEF,2)=1
            WRITE(26,1632)
 1632       FORMAT(3X,'Detected particles = photons')
          ELSE IF(KPARD.EQ.3) THEN
            KKDI(NDIDEF,3)=1
            WRITE(26,1633)
 1633       FORMAT(3X,'Detected particles = positrons')
          ENDIF
   68     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 68
          IF(KWORD.EQ.KWIPAR) GO TO 67
          IF(KWORD.EQ.KWIBOD) GO TO 65
          IF(KWORD.EQ.KWIDET) THEN
            ITST=MAX(KKDI(NDIDEF,1),KKDI(NDIDEF,2),KKDI(NDIDEF,3))
            IF(ITST.EQ.0) THEN
              KKDI(NDIDEF,1)=1
              KKDI(NDIDEF,2)=1
              KKDI(NDIDEF,3)=1
              WRITE(26,1630)
            ENDIF
            GO TO 61
          ENDIF
        ENDIF
      ENDIF
C
      IF(NDIDEF.GT.0) THEN
        ITST=MAX(KKDI(NDIDEF,1),KKDI(NDIDEF,2),KKDI(NDIDEF,3))
        IF(ITST.EQ.0) THEN
          KKDI(NDIDEF,1)=1
          KKDI(NDIDEF,2)=1
          KKDI(NDIDEF,3)=1
          WRITE(26,1630)
        ENDIF
      ENDIF
C
C  ************  Energy-deposition detectors.
C
      DO KB=1,NBODY+1
        KBDE(KB)=0
      ENDDO
      DO KD=1,NEDM
        TDED(KD)=0.0D0
        TDED2(KD)=0.0D0
      ENDDO
C
      NDEDEF=0
   43 CONTINUE
      IF(KWORD.EQ.KWEDET) THEN
        NDEDEF=NDEDEF+1
        IF(NDEDEF.GT.NEDM) THEN
          WRITE(26,'(3X,''NDEDEF = '',I4)') NDEDEF
          WRITE(26,*) 'Too many energy-deposition detectors.'
          STOP 'Too many energy-deposition detectors.'
        ENDIF
        WRITE(26,1650) NDEDEF
 1650   FORMAT(//3X,72('-'),/
     1    3X,'>>>>>>  Energy-deposition detector #', I2)
        IF(IFORON.NE.0) THEN
          WRITE(26,'(3X,''#  WARNING: May be strongly biased when '',
     1      ''interaction forcing is used!'')')
        ENDIF
        READ(BUFFER,*) EDEL(NDEDEF),EDEU(NDEDEF),NDECH(NDEDEF)
        IF(EDEL(NDEDEF).LT.0.0D0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'EDEL must be positive.'
          STOP 'EDEL must be positive.'
        ENDIF
        IF(EDEU(NDEDEF).LT.EDEL(NDEDEF)) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect energy limits.'
          STOP 'Incorrect energy limits.'
        ENDIF
        IF(NDECH(NDEDEF).LT.10.OR.NDECH(NDEDEF).GT.NEDCM) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect number of channels.'
          STOP 'Incorrect number of channels.'
        ENDIF
        WRITE(26,1610) EDEL(NDEDEF),EDEU(NDEDEF),NDECH(NDEDEF)
        BDEE(NDEDEF)=1.0000001D0*(EDEU(NDEDEF)-EDIL(NDEDEF))
     1           /DBLE(NDECH(NDEDEF))
C
   44   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 44
        IF(KWORD.EQ.KWESPC) THEN
          READ(BUFFER,'(A20)') SPCDEO(NDEDEF)
          WRITE(26,1651) SPCDEO(NDEDEF)
 1651     FORMAT(3X,'Output spectrum: ',A20)
   45     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 45
        ELSE
          WRITE(BUF2,'(I5)') NDEDEF
          SPCDEO(NDEDEF)='pm_spc_enddet_'//BUF2(5:5)//'.dat'
          WRITE(26,1651) SPCDEO(NDEDEF)
        ENDIF
C
   46   CONTINUE
        IF(KWORD.EQ.KWEBOD) THEN
          READ(BUFFER,*) KB
          IF(KB.LT.0.OR.KB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body label.'
            STOP 'Incorrect body label.'
          ENDIF
          IF(KBDE(KB).NE.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A body cannot be part of two detectors.'
            STOP 'A body cannot be part of two detectors.'
          ENDIF
          WRITE(26,1652) KB
 1652     FORMAT(3X,'Active body = ',I4)
          KBDE(KB)=NDEDEF
        ENDIF
   47   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 47
        IF(KWORD.EQ.KWEDET) GO TO 43
        IF(KWORD.EQ.KWEBOD) GO TO 46
      ENDIF
C
C  ************  Dose distribution.
C
      IDOSE=0
      IF(KWORD.EQ.KGRDXX) THEN
        IDOSE=1
        WRITE(26,1700)
 1700   FORMAT(//3X,72('-'),/3X,'>>>>>>  3D Dose distribution.')
        READ(BUFFER,*) DXL(1),DXU(1)
        IF(DXL(1).GT.DXU(1)) THEN
          SAVE=DXL(1)
          DXL(1)=DXU(1)
          DXU(1)=SAVE
        ENDIF
        IF(DXU(1).LT.DXL(1)+1.0D-6) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'XU must be greater than XL+1.0E-6.'
          STOP 'XU must be greater than XL+1.0E-6.'
        ENDIF
        WRITE(26,1710) DXL(1),DXU(1)
 1710   FORMAT(3X,'Dose-map enclosure:  XL = ',1P,E13.6,' cm,  XU = ',
     1    E13.6,' cm')
   71   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 71
        IF(KWORD.EQ.KGRDYY) THEN
          READ(BUFFER,*) DXL(2),DXU(2)
          IF(DXL(2).GT.DXU(2)) THEN
            SAVE=DXL(2)
            DXL(2)=DXU(2)
            DXU(2)=SAVE
          ENDIF
          IF(DXU(2).LT.DXL(2)+1.0D-6) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'YU must be greater than YL+1.0E-6.'
            STOP 'YU must be greater than YL+1.0E-6.'
          ENDIF
          WRITE(26,1711) DXL(2),DXU(2)
 1711     FORMAT(24X,'YL = ',1P,E13.6,' cm,  YU = ',E13.6,' cm')
        ELSE
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Unrecognized keyword.'
          STOP 'Unrecognized keyword.'
        ENDIF
   72   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 72
        IF(KWORD.EQ.KGRDZZ) THEN
          READ(BUFFER,*) DXL(3),DXU(3)
          IF(DXL(3).GT.DXU(3)) THEN
            SAVE=DXL(3)
            DXL(3)=DXU(3)
            DXU(3)=SAVE
          ENDIF
          IF(DXU(3).LT.DXL(3)+1.0D-6) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'ZU must be greater than ZL+1.0E-6.'
            STOP 'ZU must be greater than ZL+1.0E-6.'
          ENDIF
          WRITE(26,1712) DXL(3),DXU(3)
 1712     FORMAT(24X,'ZL = ',1P,E13.6,' cm,  ZU = ',E13.6,' cm')
        ELSE
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Unrecognized keyword.'
          STOP 'Unrecognized keyword.'
        ENDIF
   73   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 73
        IF(KWORD.EQ.KGRDBN) THEN
          READ(BUFFER,*) NDB(1),NDB(2),NDB(3)
          IF(NDB(1).LT.0.OR.NDB(1).GT.NDXM) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,'(''NDB(1) must be .GT.0. and .LT.'',I4)') NDXM
            WRITE(26,*) 'Increase the value of the parameter NDXM.'
            STOP 'NDB(1) must be .GT.0. and .LE.NDXM'
          ENDIF
          IF(NDB(2).LT.0.OR.NDB(2).GT.NDYM) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,'(''NDB(2) must be .GT.0. and .LT.'',I4)') NDYM
            WRITE(26,*) 'Increase the value of the parameter NDYM.'
            STOP 'NDB(2) must be .GT.0. and .LE.NDYM'
          ENDIF
          IF(NDB(3).LT.0.OR.NDB(3).GT.NDZM) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,'(''NDB(3) must be .GT.0. and .LT.'',I4)') NDZM
            WRITE(26,*) 'Increase the value of the parameter NDZM.'
            STOP 'NDB(3) must be .GT.0. and .LE.NDZM'
          ENDIF
          WRITE(26,1713) NDB(1),NDB(2),NDB(3)
 1713     FORMAT(3X,'Numbers of bins:     NBX =',I4,', NBY =',I4,
     1      ', NBZ =',I4)
        ELSE
          NDB(1)=1
          NDB(2)=1
          NDB(3)=1
          WRITE(26,1713) NDB(1),NDB(2),NDB(3)
        ENDIF
        DO I=1,3
          BDOSE(I)=1.0000001D0*(DXU(I)-DXL(I))/DBLE(NDB(I))
          BDOSER(I)=1.0D0/BDOSE(I)
        ENDDO
   74   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 74
      ENDIF
C
C  ************  Job characteristics.
C
      WRITE(26,1800)
 1800 FORMAT(//3X,72('-'),/
     1  3X,'>>>>>>  Job characteristics.')
C
      IRESUM=0
      IF(KWORD.EQ.KWRESU) THEN
        READ(BUFFER,'(A20)') PFILER
        WRITE(26,1810) PFILER
 1810   FORMAT(3X,'Resume simulation from previous dump file: ',A20)
        IRESUM=1
   81   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 81
      ENDIF
C
      IDUMP=0
      DUMPP=1.0D15
      IF(KWORD.EQ.KWDUMP) THEN
        READ(BUFFER,'(A20)') PFILED
        WRITE(26,1820) PFILED
 1820   FORMAT(3X,'Write final counter values on the dump file: ',A20)
        IDUMP=1
   82   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 82
      ENDIF
C
      IF(KWORD.EQ.KWDMPP) THEN
        READ(BUFFER,*) DUMPP
        IF(IDUMP.EQ.1) THEN
          IF(DUMPP.LT.15.0D0) DUMPP=15.0D0
          IF(DUMPP.GT.86400.0D0) DUMPP=86400.0D0
          WRITE(26,1830) DUMPP
 1830     FORMAT(3X,'Dumping period: DUMPP =',1P,E13.6)
        ENDIF
   83   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 83
      ENDIF
C
      IF(KWORD.EQ.KWNSIM) THEN
        READ(BUFFER,*) DSHN
        IF(DSHN.LT.1.0D0) DSHN=2.0D9
   84   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 84
      ELSE
        DSHN=2.0D9
      ENDIF
      WRITE(26,1840) DSHN
 1840 FORMAT(/3X,'Number of showers to be simulated =',1P,E13.6)
C
      IF(KWORD.EQ.KWRSEE) THEN
        READ(BUFFER,*) ISEED1,ISEED2
        WRITE(26,1850) ISEED1,ISEED2
 1850   FORMAT(3X,'Random number generator seeds = ',I10,', ',I10)
   85   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 85
      ELSE
        ISEED1=12345
        ISEED2=54321
      ENDIF
C
      IF(KWORD.EQ.KWTIME) THEN
        READ(BUFFER,*) TIMEA
      ELSE
        TIMEA=100.0D0
      ENDIF
      IF(TIMEA.LT.1.0D0) TIMEA=100.0D0
      WRITE(26,1860) TIMEA
 1860 FORMAT(3X,'Computation time available = ',1P,E13.6,' sec')
C
      CALL TIMER(TSECIN)
      CPUT0=CPUTIM()
      TSECA=TIMEA+TSECIN
      TSECAD=TSECIN
      WRITE(26,1870)
 1870 FORMAT(//3X,72('-'))
C
C  ************  If 'RESUME' is active, read previously generated
C                counters...
C
      SHNA=0.0D0
      CPUTA=0.0D0
      IRETRN=0
C
      DO ID=1,NDIDEF
        RLAST(ID)=0.0D0
        RWRITE(ID)=0.0D0
      ENDDO
C
      IF(IRESUM.EQ.1) THEN
        OPEN(9,FILE=PFILER)
        READ (9,*,ERR=91,END=91) SHNA,CPUTA
        READ (9,'(A120)') TITLE2
        IF(TITLE2.NE.TITLE) THEN
          WRITE(26,*)
     1      'The dump file is corrupted (the TITLE does not match).'
          STOP 'The dump file is corrupted (the TITLE does not match).'
        ENDIF
        READ (9,*) ISEED1,ISEED2
        READ (9,*) NPSN,RLREAD
        READ (9,*) (SHIST(I),I=1,NSEB)
        READ (9,*) (PRIM(I),I=1,3),(PRIM2(I),I=1,3)
        READ (9,*) ((SEC(K,I),I=1,3),K=1,3),((SEC2(K,I),I=1,3),K=1,3)
        READ (9,*) (TDEBO(I),I=1,NBODY), (TDEBO2(I),I=1,NBODY)
        READ (9,*) (((PDE(I,J,K),K=1,NBE),J=1,2),I=1,3),
     1             (((PDE2(I,J,K),K=1,NBE),J=1,2),I=1,3)
        READ (9,*) (((PDA(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3),
     1             (((PDA2(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3)
        READ (9,*) (TDID(I),I=1,NIDM), (TDID2(I),I=1,NIDM)
        READ (9,*) (TDED(I),I=1,NEDM), (TDED2(I),I=1,NEDM)
        IF(NDIDEF.GT.0) THEN
          READ (9,*) (RLAST(ID),ID=1,NDIDEF)
          READ (9,*) (RWRITE(ID),ID=1,NDIDEF)
          DO ID=1,NDIDEF
            READ (9,*) (DIT(ID,J),J=1,NDICH(ID))
            READ (9,*) (DIT2(ID,J),J=1,NDICH(ID))
          ENDDO
        ENDIF
        IF(NDEDEF.GT.0) THEN
          DO ID=1,NDEDEF
            READ (9,*) (DET(ID,J),J=1,NDECH(ID))
            READ (9,*) (DET2(ID,J),J=1,NDECH(ID))
          ENDDO
        ENDIF
        IF(IDOSE.NE.0) THEN
          READ (9,*)
     1      (((DOSE(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1)),
     1      (((DOSE2(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        ENDIF
        CLOSE(9)
        WRITE(26,1880) PFILER
 1880   FORMAT(3X,'Simulation has been resumed from dump file: ',A20)
        GO TO 92
   91   CONTINUE
        WRITE(26,1890)
 1890   FORMAT(3X,'WARNING: Could not resume from dump file...')
        IRESUM=0
      ENDIF
   92 CONTINUE
C
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
          IF(IPSF(ID).NE.0) THEN
            IPSFU=20+ID
            OPEN(IPSFU,FILE=PSFDIO(ID),IOSTAT=KODE)
            IF(KODE.NE.0) THEN
              WRITE(26,'(''File '',A20,'' could not be opened.'')')
     1          PSFDIO(ID)
              STOP 'File could not be opened.'
            ENDIF
            RWR=0.0D0
            IF(RWRITE(ID).GT.0) THEN
   93         CONTINUE
              CALL RDPSF(IPSFU,PSFREC,KODE)
              IF(KODE.NE.0) THEN
                GO TO 94
              ELSE
                RWR=RWR+1.0D0
                IF(RWR.LT.RWRITE(ID)-0.5D0) GO TO 93
                GO TO 94
              ENDIF
            ENDIF
   94       CONTINUE
            IF(RWR.LT.0.5D0) THEN
              WRITE(IPSFU,1901) ID
 1901         FORMAT(1X,'#  Results from PENMAIN. Phase-space fi',
     1          'le of detector no.',I3,/1X,'#')
              WRITE(IPSFU,1902)
 1902         FORMAT(1X,'#/KPAR',2X,'E',12X,'X',12X,'Y',12X,'Z',12X,
     1          'U',12X,'V',12X,'W',11X,'WGHT',5X,'ILB(1:4)',7X,'NSHI',
     1          /1X,'#',125('-'))
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
      IPSFI=20
      IF(ISORP.EQ.1) THEN
        IF(IRESUM.EQ.1) THEN
          IF(NPSN.GT.NPSF) THEN
            WRITE(6,*) '   **** The simulation was already completed.'
            WRITE(26,*) '   **** The simulation was already completed.'
            SHN=SHNA
            GO TO 203
          ENDIF
          IF(NPSN.GT.1) THEN
            DO I=1,NPSN-1
              WRITE(6,1903) PSFI(I)
              WRITE(26,1903) PSFI(I)
 1903         FORMAT(/3X,'+ The phase-space file ',A20,
     1          ' was read in previous runs.')
            ENDDO
          ENDIF
   95     CONTINUE
          OPEN(IPSFI,FILE=PSFI(NPSN))
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
 1904     FORMAT(/3X,'+ The phase-space file ',A20,' is opened.')
          IF(RLREAD.GT.0.5D0) THEN
            RI=0.0D0
   96       CONTINUE
            RI=RI+1.0D0
            CALL RDPSF(IPSFI,PSFREC,KODE)
            IF(KODE.NE.0) THEN
              NPSN=NPSN+1
              IF(NPSN.GT.NPSF) THEN
                WRITE(6,*)
     1            '   **** The simulation was already completed.'
                WRITE(26,*)
     1            '   **** The simulation was already completed.'
                SHN=SHNA
                GO TO 203
              ELSE
                WRITE(6,1905) PSFI(NPSN)
                WRITE(26,1905) PSFI(NPSN)
 1905           FORMAT(/3X,'+ The phase-space file ',A20,
     1            ' was completed in the last run.')
                CLOSE(IPSFI)
                RLREAD=0.0D0
                GO TO 95
              ENDIF
            ENDIF
            IF(RI.LT.RLREAD-0.5D0) GO TO 96
          ELSE
            RLREAD=0.0D0
          ENDIF
        ELSE
          NPSN=1
          RLREAD=0.0D0
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
          OPEN(IPSFI,FILE=PSFI(NPSN),IOSTAT=KODE)
        ENDIF
      ENDIF
C
C  ************  Initialise constants.
C
      SHN=SHNA          ! Shower counter, including the dump file.
      DSHN=DSHN+SHNA
      N=MOD(SHN,2.0D9)+0.5D0
      IF(SHN.GE.DSHN) GO TO 105

C>>>>>>>>> EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C  ************  Initializing the electromagnetic field routines.
C
      OPEN(29,FILE='track')
      EFX=-15000.0D0
      EFY=-15000.0D0
      EFZ=0.0D0
      WRITE(26,'(/3X,''E-field (V/cm)  ='',1P,3E15.7)') EFX,EFY,EFZ
      WRITE(29,'('' # E-field (V/cm)  ='',1P,3E15.7)') EFX,EFY,EFZ
C      BFX=1.5D4
      BFX=0D0
      BFY=0.0D0
      BFZ=0.0D0
      WRITE(26,'(/3X,''M-field (gauss) ='',1P,3E15.7)') BFX,BFY,BFZ
      WRITE(29,'('' # M-field (gauss) ='',1P,3E15.7)') BFX,BFY,BFZ
C
      IF(KPARP.EQ.1) THEN
        WRITE(29,*) '# Electron'
      ELSE IF(KPARP.EQ.3) THEN
        WRITE(29,*) '# Positron'
      ELSE
        STOP 'Please, give me a charged particle....'
      ENDIF
      ULDV=1.0D-2
      ULDE=1.0D-2
      ULEM=1.0D-2
C<<<<<<<<< EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  ------------------------  Shower simulation starts here.
C
  101 CONTINUE
C
C  ************  Primary particle counters.
C
      DO I=1,3
        DPRIM(I)=0.0D0
        DO K=1,3
          DSEC(K,I)=0.0D0
        ENDDO
      ENDDO
      DO KB=1,NBODY
        DEBO(KB)=0.0D0  ! Energies deposited in the various bodies KB.
      ENDDO
      IEXIT=0
C
      IF(NDIDEF.GT.0) THEN
        DO KD=1,NDIDEF
          DEDI(KD)=0.0D0  ! Energy collected by each detector.
        ENDDO
      ENDIF
C
      CALL CLEANS          ! Cleans the secondary stack.
C
C  **********  Set the initial state of the primary particle.
C
      IF(ISORP.EQ.0) THEN
C  ****  Point source.
        SHN=SHN+1.0D0
        N=N+1
        IF(N.GT.2000000000) N=N-2000000000
        KPAR=KPARP
        WGHT=WGHT0
C  ----  Initial position ...
        X=SX0
        Y=SY0
        Z=SZ0
C  ----  Initial direction ...
        CALL GCONE(U,V,W)
C  ----  initial energy ...
        IF(ISPEC.EQ.0) THEN
          E=E0    ! Monoenergetic source.
          SHIST(1)=SHIST(1)+1.0D0
        ELSE      ! Continuous spectrum. E sampled by Walker's method.
          RN=RAND(4.0D0)*NSEB+1
          K=INT(RN)
          RNF=RN-K
          IF(RNF.GT.FS(K)) THEN
            KE=IAS(K)
          ELSE
            KE=K
          ENDIF
          E=ES(KE)+RAND(5.0D0)*(ES(KE+1)-ES(KE))
          SHIST(KE)=SHIST(KE)+1.0D0
        ENDIF
        ILB(1)=1  ! Identifies primary particles.
        ILB(2)=0
        ILB(3)=0
        ILB(4)=0
        ILB(5)=0
      ELSE
C  ****  Phase-space file.
  201   CONTINUE
        CALL RDPSF(IPSFI,PSFREC,KODE)
        IF(KODE.EQ.0) THEN
          READ(PSFREC,*) KPAR,E,X,Y,Z,U,V,W,WGHT,
     1      ILB(1),ILB(2),ILB(3),ILB(4),NSHI
        ELSE
          CLOSE(IPSFI)
          NPSN=NPSN+1
          RLREAD=0.0D0
          IRETRN=0
          IF(NPSN.GT.NPSF) GO TO 203
          OPEN(IPSFI,FILE=PSFI(NPSN))
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
          GO TO 201
        ENDIF
        RLREAD=RLREAD+1.0D0
        SHN=SHN+DBLE(NSHI)
        N=N+NSHI
        IF(N.GT.2000000000) N=N-2000000000
        ILB(5)=0
      ENDIF
      IF(E.LT.EABS(KPAR,MAT)) GO TO 105
C
C  ****  Check if the trajectory intersects the material system.
C
      CALL LOCATE
C
      IF(MAT.EQ.0) THEN
        IBODYL=IBODY
        CALL STEP(1.0D30,DSEF,NCROSS)
        IF(MAT.EQ.0) THEN  ! The particle does not enter the system.
          IF(W.GT.0) THEN
            IEXIT=1        ! Labels emerging 'transmitted' particles.
          ELSE
            IEXIT=2        ! Labels emerging 'backscattered' particles.
          ENDIF
          GO TO 104        ! Exit.
        ENDIF

C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Impact detectors.
        IDET=KBDI(IBODY)
        IF(IDET.NE.0) THEN
          IF(KBDI(IBODYL).NE.IDET.AND.KKDI(IDET,KPAR).EQ.1) THEN
C
            IF(IPSF(IDET).NE.0) THEN
              NSHJ=SHN-RLAST(IDET)
              CALL N2CH10(NSHJ,LCH10,NDIG)
              WRITE(20+IDET,'(I2,1P,8E13.5,I3,I2,I2,I9,1X,A)')
     1          KPAR,E,X,Y,Z,U,V,W,WGHT,ILB(1),ILB(2),ILB(3),ILB(4),
     2          LCH10(1:NDIG)
              RWRITE(IDET)=RWRITE(IDET)+1.0D0
              RLAST(IDET)=SHN
            ENDIF
C
            DEDI(IDET)=DEDI(IDET)+E*WGHT
            IE=1.0D0+(E-EDIL(IDET))/BDIE(IDET)
            IF(IE.GT.0.AND.IE.LE.NDICH(IDET)) THEN
              IF(N.NE.LDIT(IDET,IE)) THEN
                DIT(IDET,IE)=DIT(IDET,IE)+DITP(IDET,IE)
                DIT2(IDET,IE)=DIT2(IDET,IE)+DITP(IDET,IE)**2
                DITP(IDET,IE)=WGHT
                LDIT(IDET,IE)=N
              ELSE
                DITP(IDET,IE)=DITP(IDET,IE)+WGHT
              ENDIF
            ENDIF
            IF(IPSF(IDET).EQ.-1) THEN
              DEBO(IBODY)=DEBO(IBODY)+E*WGHT
              IEXIT=3
              GO TO 104
            ENDIF
          ENDIF
        ENDIF
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      ELSE
        IF(E.LT.EABS(KPAR,MAT)) THEN  ! The energy is too low.
          DEBO(IBODY)=DEBO(IBODY)+E*WGHT
          IEXIT=3
          GO TO 104
        ENDIF
      ENDIF
C  ---------------------------------------------------------------------
C  ------------------------  Track simulation begins here.
C

C-----ParticleSplitting--ParticleSplitting--ParticleSplitting-----------
C  ****  Particle splitting (only when read from the phase-space file).
      IF(ISORP.EQ.1) THEN
        IF(NSPLIT.GT.1) THEN
          CALL VSPLIT(NSPLIT)  ! Particle splitting.
C  ----  Energy is locally deposited in the material.
          DEP=(NSPLIT-1)*E*WGHT
          DEBO(IBODY)=DEBO(IBODY)+DEP
          IF(IDOSE.NE.0) THEN  ! Particle inside the dose enclosure.
            IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
     1         (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
     1         (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
              I1=1.0D0+(X-DXL(1))*BDOSER(1)
              I2=1.0D0+(Y-DXL(2))*BDOSER(2)
              I3=1.0D0+(Z-DXL(3))*BDOSER(3)
              IF(N.NE.LDOSE(I1,I2,I3)) THEN
                DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
                DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
                DOSEP(I1,I2,I3)=DEP*RHOI(MAT)
                LDOSE(I1,I2,I3)=N
              ELSE
                DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DEP*RHOI(MAT)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C-----ParticleSplitting--ParticleSplitting--ParticleSplitting-----------

  102 CONTINUE
      CALL START           ! Starts simulation in current medium.

C>>>>>>>>> EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      WRITE(29,*) '  '
      IF(KPAR.NE.2) write(29,'(1x,1p,4e15.7,i3)') x,y,z,e,kpar
C<<<<<<<<< EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C
  103 CONTINUE
      IBODYL=IBODY
C     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,i3)')
C    1    mod(n,100),kpar,ilb(1),x,y,z,w,e,ibody
C
c     IF((IFORCE(IBODY,KPAR).EQ.1).AND.((WGHT.GE.WLOW(IBODY,KPAR)).AND.
c    1  (WGHT.LE.WHIG(IBODY,KPAR)))) THEN
c       CALL JUMPF(DSMAX(IBODY),DS)  ! Interaction forcing.
c       INTFOR=1
c     ELSE
c       CALL JUMP(DSMAX(IBODY),DS)   ! Analogue simulation.
c       INTFOR=0
c     ENDIF
c     CALL STEP(DS,DSEF,NCROSS)      ! Determines step end position.

C>>>>>>>>> EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c  the original block has been commented out.
      CALL TPEMF0(ULDV,ULDE,ULEM,TMAX)
      DSMAXL=DMIN1(TMAX,DSMAX(MAT))
      IF((IFORCE(IBODY,KPAR).EQ.1).AND.((WGHT.GT.WLOW(IBODY,KPAR)).AND.
     1  (WGHT.LT.WHIG(IBODY,KPAR)))) THEN
        CALL JUMPF(DSMAXL,DS)  ! Interaction forcing.
        INTFOR=1
      ELSE
        CALL JUMP(DSMAXL,DS)   ! Analogue simulation.
        INTFOR=0
      ENDIF
      CALL TPEMF1(DS,DSEF,NCROSS)
      IF(KPAR.NE.2) write(29,'(1x,1p,4e15.7,i3)') x,y,z,e,kpar
C<<<<<<<<< EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C  ----  Check whether the particle is outside the enclosure.
      IF(MAT.EQ.0) THEN
        IF(W.GT.0) THEN
          IEXIT=1        ! Labels emerging 'transmitted' particles.
        ELSE
          IEXIT=2        ! Labels emerging 'backscattered' particles.
        ENDIF
        GO TO 104        ! Exit.
      ENDIF

C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Impact detectors.
      IDET=KBDI(IBODY)
      IF(IDET.NE.0) THEN
        IF(KBDI(IBODYL).NE.IDET.AND.KKDI(IDET,KPAR).EQ.1) THEN
C
          IF(IPSF(IDET).NE.0) THEN
            NSHJ=SHN-RLAST(IDET)
            CALL N2CH10(NSHJ,LCH10,NDIG)
            WRITE(20+IDET,'(I2,1P,8E13.5,I3,I2,I2,I9,1X,A)')
     1        KPAR,E,X,Y,Z,U,V,W,WGHT,ILB(1),ILB(2),ILB(3),ILB(4),
     2        LCH10(1:NDIG)
            RWRITE(IDET)=RWRITE(IDET)+1.0D0
            RLAST(IDET)=SHN
          ENDIF
C
          DEDI(IDET)=DEDI(IDET)+E*WGHT
          IE=1.0D0+(E-EDIL(IDET))/BDIE(IDET)
          IF(IE.GT.0.AND.IE.LE.NDICH(IDET)) THEN
            IF(N.NE.LDIT(IDET,IE)) THEN
              DIT(IDET,IE)=DIT(IDET,IE)+DITP(IDET,IE)
              DIT2(IDET,IE)=DIT2(IDET,IE)+DITP(IDET,IE)**2
              DITP(IDET,IE)=WGHT
              LDIT(IDET,IE)=N
            ELSE
              DITP(IDET,IE)=DITP(IDET,IE)+WGHT
            ENDIF
          ENDIF
          IF(IPSF(IDET).EQ.-1) THEN
            DEBO(IBODY)=DEBO(IBODY)+E*WGHT
            IEXIT=3
            GO TO 104
          ENDIF
        ENDIF
      ENDIF
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C  ----  If the particle has crossed an interface, restart the track in
C        the new material.
      IF(NCROSS.GT.0) GO TO 102    ! The particle crossed an interface.
C  ----  Simulate next event.
      IF(INTFOR.EQ.0) THEN
        CALL KNOCK(DE,ICOL)        ! Analogue simulation.
      ELSE
        CALL KNOCKF(DE,ICOL)       ! Interaction forcing is active.
      ENDIF
      DEP=DE*WGHT
C  ----  Energy is locally deposited in the material.
      DEBO(IBODY)=DEBO(IBODY)+DEP
      IF(IDOSE.NE.0) THEN  ! The particle is inside the dose enclosure.
        IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
     1     (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
     1     (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
          I1=1.0D0+(X-DXL(1))*BDOSER(1)
          I2=1.0D0+(Y-DXL(2))*BDOSER(2)
          I3=1.0D0+(Z-DXL(3))*BDOSER(3)
          IF(N.NE.LDOSE(I1,I2,I3)) THEN
            DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
            DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
            DOSEP(I1,I2,I3)=DEP*RHOI(MAT)
            LDOSE(I1,I2,I3)=N
          ELSE
            DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DEP*RHOI(MAT)
          ENDIF
        ENDIF
      ENDIF

C-----RussianRoulette-RussianRoulette-RussianRoulette-RussianRoulette---
C  ****  Russian roulette for photons moving downstream.
cxxx  IF(KPAR.EQ.2.AND.W.LT.0.0D0) THEN
cxxx    IF(WGHT.LT.5.0D0) THEN
cxxx      PKILL=0.75D0
cxxx      CALL VKILL(PKILL)
cxxx    ENDIF
cxxx  ENDIF
C-----RussianRoulette-RussianRoulette-RussianRoulette-RussianRoulette---

C
      IF(E.LT.EABS(KPAR,MAT)) THEN  ! The particle has been absorbed.
        IEXIT=3                     ! Labels absorbed particles.
        GO TO 104                   ! Exit.
      ENDIF
C
      GO TO 103
C  ------------------------  The simulation of the track ends here.
C  ---------------------------------------------------------------------
  104 CONTINUE
C
C  ************  Increment particle counters.
C
      IF(ILB(1).EQ.1) THEN
        DPRIM(IEXIT)=DPRIM(IEXIT)+WGHT
      ELSE
        DSEC(KPAR,IEXIT)=DSEC(KPAR,IEXIT)+WGHT
      ENDIF
C
      IF(IEXIT.LT.3) THEN
C  ****  Energy distribution of emerging particles.
        K=1.0D0+(E-EMIN)/BSE(KPAR)
        IF(K.GT.0.AND.K.LE.NBE) THEN
          IF(N.NE.LPDE(KPAR,IEXIT,K)) THEN
            PDE(KPAR,IEXIT,K)=PDE(KPAR,IEXIT,K)+PDEP(KPAR,IEXIT,K)
            PDE2(KPAR,IEXIT,K)=PDE2(KPAR,IEXIT,K)+PDEP(KPAR,IEXIT,K)**2
            PDEP(KPAR,IEXIT,K)=WGHT
            LPDE(KPAR,IEXIT,K)=N
          ELSE
            PDEP(KPAR,IEXIT,K)=PDEP(KPAR,IEXIT,K)+WGHT
          ENDIF
        ENDIF
C  ****  Angular distribution of emerging particles.
        THETA=ACOS(W)
        KTH=1.0D0+THETA*RA2DE/BSTH
        IF(ABS(U).GT.1.0D-16) THEN  !  Azimuthal bin number corrected.
           PHI=ATAN2(V,U)
        ELSE IF(ABS(V).GT.1.0D-16) THEN
           PHI=ATAN2(V,U)
        ELSE
           PHI=0.0D0
        ENDIF
        IF(PHI.LT.0.0D0) PHI=TWOPI+PHI
        KPH=1.0D0+PHI*RA2DE/BSPH
        IF(N.NE.LPDA(KPAR,KTH,KPH)) THEN
          PDA(KPAR,KTH,KPH)=PDA(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)
          PDA2(KPAR,KTH,KPH)=PDA2(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)**2
          PDAP(KPAR,KTH,KPH)=WGHT
          LPDA(KPAR,KTH,KPH)=N
        ELSE
          PDAP(KPAR,KTH,KPH)=PDAP(KPAR,KTH,KPH)+WGHT
        ENDIF
      ENDIF
C
C  ************  Any secondary left?
C
  202 CONTINUE
      CALL SECPAR(LEFT)
      IF(LEFT.GT.0) THEN
C     write(6,'(/''new secondary'')')
C     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,i3)')
C    1    mod(n,100),kpar,ilb(1),x,y,z,w,e,ibody

C-----RussianRoulette-RussianRoulette-RussianRoulette-RussianRoulette---
C  ****  Russian roulette for photons moving downstream.
cxxx    IF(KPAR.EQ.2.AND.W.LT.0.0D0) THEN
cxxx      IF(WGHT.LT.5.0D0) THEN
cxxx        PKILL=0.75D0
cxxx        CALL VKILL(PKILL)
cxxx      ENDIF
cxxx      IF(E.LT.EABS(KPAR,MAT)) GO TO 202
cxxx    ENDIF
C-----RussianRoulette-RussianRoulette-RussianRoulette-RussianRoulette---

        DEBO(IBODY)=DEBO(IBODY)-E*WGHT  ! Energy is removed.
        IF(IDOSE.NE.0) THEN  ! Particle inside the dose enclosure.
          IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
     1       (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
     1       (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
            I1=1.0D0+(X-DXL(1))*BDOSER(1)
            I2=1.0D0+(Y-DXL(2))*BDOSER(2)
            I3=1.0D0+(Z-DXL(3))*BDOSER(3)
            DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)-E*WGHT*RHOI(MAT)
          ENDIF
        ENDIF
        GO TO 102
      ENDIF
C
C  ------------------------  The simulation of the shower ends here.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      DO I=1,3
        PRIM(I)=PRIM(I)+DPRIM(I)
        PRIM2(I)=PRIM2(I)+DPRIM(I)**2
        DO K=1,3
          SEC(K,I)=SEC(K,I)+DSEC(K,I)
          SEC2(K,I)=SEC2(K,I)+DSEC(K,I)**2
        ENDDO
      ENDDO
C
C  ****  Energies deposited in different bodies.
C
      IF(NDEDEF.GT.0) THEN
        DO KD=1,NDEDEF
          DEDE(KD)=0.0D0
        ENDDO
      ENDIF
      DO KB=1,NBODY
        TDEBO(KB)=TDEBO(KB)+DEBO(KB)
        TDEBO2(KB)=TDEBO2(KB)+DEBO(KB)**2
        IF(NDEDEF.GT.0) THEN
          IDET=KBDE(KB)
          IF(IDET.NE.0) THEN
            DEDE(IDET)=DEDE(IDET)+DEBO(KB)
          ENDIF
        ENDIF
      ENDDO
C
C  ****  Tallying the spectra from energy-deposition detectors.
C
      IF(NDEDEF.GT.0) THEN
        DO KD=1,NDEDEF
          TDED(KD)=TDED(KD)+DEDE(KD)
          TDED2(KD)=TDED2(KD)+DEDE(KD)**2
          IF(DEDE(KD).GT.1.0D-5) THEN
            KE=1.0D0+(DEDE(KD)-EDEL(KD))/BDEE(KD)
            IF(KE.GT.0.AND.KE.LE.NDECH(KD)) THEN
              DET(KD,KE)=DET(KD,KE)+WGHT0
              DET2(KD,KE)=DET2(KD,KE)+WGHT0**2
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C  ****  Average energies 'collected' by impact detectors.
C
      IF(NDIDEF.GT.0) THEN
        DO KD=1,NDIDEF
          TDID(KD)=TDID(KD)+DEDI(KD)
          TDID2(KD)=TDID2(KD)+DEDI(KD)**2
        ENDDO
      ENDIF
C
  105 CONTINUE
      CALL TIMER(TSEC)
C
C  ----  End the simulation after the allotted time or after completing
C        DSHN showers.
C
      IRETRN=0
      IF(TSEC.LT.TSECA.AND.SHN.LT.DSHN) THEN
C  ****  Write partial results after each dumping period.
        IF(TSEC-TSECAD.GT.DUMPP) THEN
          IF(IDUMP.EQ.1) THEN
            TSECAD=TSEC
            IRETRN=1
            GO TO 203
          ENDIF
        ENDIF
        GO TO 101
      ENDIF
  203 CONTINUE
      TSIM=CPUTA+CPUTIM()-CPUT0
C
C  ************  Transfer partial counters to global counters.
C
      DO KPAR=1,3
      DO IEXIT=1,2
      DO IE=1,NBE
        PDE(KPAR,IEXIT,IE)=PDE(KPAR,IEXIT,IE)+PDEP(KPAR,IEXIT,IE)
        PDE2(KPAR,IEXIT,IE)=PDE2(KPAR,IEXIT,IE)+PDEP(KPAR,IEXIT,IE)**2
        PDEP(KPAR,IEXIT,IE)=0.0D0
        LPDE(KPAR,IEXIT,IE)=0
      ENDDO
      ENDDO
      ENDDO
C
      DO KPAR=1,3
      DO KTH=1,NBTH
      DO KPH=1,NBPH
        PDA(KPAR,KTH,KPH)=PDA(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)
        PDA2(KPAR,KTH,KPH)=PDA2(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)**2
        PDAP(KPAR,KTH,KPH)=0.0D0
        LPDA(KPAR,KTH,KPH)=0
      ENDDO
      ENDDO
      ENDDO
C
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
        DO J=1,NDICH(ID)
          DIT(ID,J)=DIT(ID,J)+DITP(ID,J)
          DIT2(ID,J)=DIT2(ID,J)+DITP(ID,J)**2
          DITP(ID,J)=0.0D0
          LDIT(ID,J)=0
        ENDDO
        ENDDO
      ENDIF
      IF(IDOSE.NE.0) THEN
        DO I3=1,NDB(3)
        DO I2=1,NDB(2)
        DO I1=1,NDB(1)
          DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
          DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
          DOSEP(I1,I2,I3)=0.0D0
          LDOSE(I1,I2,I3)=0
        ENDDO
        ENDDO
        ENDDO
      ENDIF
C
C  ************  If 'DUMPTO' is active, write counters in a dump file.
C
      IF(IDUMP.EQ.1) THEN
        OPEN(9,FILE=PFILED)
        WRITE(9,*) SHN,TSIM
        WRITE(9,'(A120)') TITLE
        WRITE(9,*) ISEED1,ISEED2
        WRITE(9,*) NPSN,RLREAD
        WRITE(9,*) (SHIST(I),I=1,NSEB)
        WRITE(9,*) (PRIM(I),I=1,3),(PRIM2(I),I=1,3)
        WRITE(9,*) ((SEC(K,I),I=1,3),K=1,3),((SEC2(K,I),I=1,3),K=1,3)
        WRITE(9,*) (TDEBO(I),I=1,NBODY), (TDEBO2(I),I=1,NBODY)
        WRITE(9,*) (((PDE(I,J,K),K=1,NBE),J=1,2),I=1,3),
     1             (((PDE2(I,J,K),K=1,NBE),J=1,2),I=1,3)
        WRITE(9,*) (((PDA(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3),
     1             (((PDA2(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3)
        WRITE(9,*) (TDID(I),I=1,NIDM), (TDID2(I),I=1,NIDM)
        WRITE(9,*) (TDED(I),I=1,NEDM), (TDED2(I),I=1,NEDM)
        IF(NDIDEF.GT.0) THEN
          WRITE(9,*) (RLAST(ID),ID=1,NDIDEF)
          WRITE(9,*) (RWRITE(ID),ID=1,NDIDEF)
          DO ID=1,NDIDEF
            WRITE(9,*) (DIT(ID,J),J=1,NDICH(ID))
            WRITE(9,*) (DIT2(ID,J),J=1,NDICH(ID))
          ENDDO
        ENDIF
        IF(NDEDEF.GT.0) THEN
          DO ID=1,NDEDEF
            WRITE(9,*) (DET(ID,J),J=1,NDECH(ID))
            WRITE(9,*) (DET2(ID,J),J=1,NDECH(ID))
          ENDDO
        ENDIF
        IF(IDOSE.NE.0) THEN
          WRITE(9,*)
     1      (((DOSE(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1)),
     1      (((DOSE2(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        ENDIF
        WRITE(9,'(/3X,''*** END ***'')')
        CLOSE(9)
      ENDIF
C
C  ------------------------  Print simulation results.
C
C     IEXIT: 1=transmitted, 2=backscattered, 3=absorbed.
C
      TOTN=SHN*WGHT0
      WRITE(26,3000)
 3000 FORMAT(/////3X,35('*')/3X,'**   Program PENMAIN. Results.   **',
     1  /3X,35('*'))
C
      WRITE(26,3001) TSIM
 3001 FORMAT(/3X,'Simulation time ......................... ',
     1  1P,E13.6,' sec')
      TAVS=TOTN/TSIM
      WRITE(26,3002) TAVS
 3002 FORMAT(3X,'Simulation speed ........................ ',
     1  1P,E13.6,' showers/sec')
      WRITE(26,3003) TOTN
 3003 FORMAT(//3X,'Simulated primary showers ............... ',
     1  1P,E13.6)
C
      WRITE(26,3004) PRIM(1)
 3004 FORMAT(/3X,'Transmitted primary particles ........... ',
     1  1P,E13.6)
      WRITE(26,3005) PRIM(2)
 3005 FORMAT(3X,'Backscattered primary particles ......... ',
     1  1P,E13.6)
      WRITE(26,3006) PRIM(3)
 3006 FORMAT(3X,'Absorbed primary particles .............. ',
     1  1P,E13.6)
C
      IF(ISORP.EQ.0) THEN
        FNT=1.0D0/TOTN
        FT=(PRIM(1)+SEC(KPARP,1))*FNT
        ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(1)-PRIM(1)**2*FNT))
        ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,1)-SEC(KPARP,1)**2*FNT))
        ERR=ERR1+ERR2
        WRITE(26,3007) FT,ERR
 3007   FORMAT(/3X,'Fractional transmission ............ ',
     1    1P,E13.6,' +-',E8.1)
        FB=(PRIM(2)+SEC(KPARP,2))*FNT
        ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(2)-PRIM(2)**2*FNT))
        ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,2)-SEC(KPARP,2)**2*FNT))
        ERR=ERR1+ERR2
        WRITE(26,3008) FB,ERR
 3008   FORMAT(3X,'Fractional backscattering .......... ',
     1    1P,E13.6,' +-',E8.1)
        FA=PRIM(3)*FNT
        ERR=3.0D0*FNT*SQRT(ABS(PRIM2(3)-PRIM(3)**2*FNT))
        WRITE(26,3009) FA,ERR
 3009   FORMAT(3X,'Fractional absorption .............. ',
     1    1P,E13.6,' +-',E8.1)
C
        DO K=1,3
          DO I=1,3
            WSEC2(K,I)=3.0D0*FNT*SQRT(ABS(SEC2(K,I)-SEC(K,I)**2*FNT))
            WSEC(K,I)=SEC(K,I)*FNT
          ENDDO
        ENDDO
        WRITE(26,3010)
     1   WSEC(1,1),WSEC(2,1),WSEC(3,1),WSEC2(1,1),WSEC2(2,1),WSEC2(3,1),
     1   WSEC(1,2),WSEC(2,2),WSEC(3,2),WSEC2(1,2),WSEC2(2,2),WSEC2(3,2),
     1   WSEC(1,3),WSEC(2,3),WSEC(3,3),WSEC2(1,3),WSEC2(2,3),WSEC2(3,3)
 3010   FORMAT(/3X,'Secondary-particle generation probabilities:',
     1    /19X,46('-'),
     1    /19X,'|  electrons   |   photons    |  positrons   |',1P,
     1    /3X,62('-')/3X,'|  transmitted  |',3(E13.6,1X,'|'),
     1    /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1    /3X,62('-')/3X,'| backscattered |',3(E13.6,1X,'|'),
     1    /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1    /3X,62('-')/3X,'|   absorbed    |',3(E13.6,1X,'|'),
     1    /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1    /3X,62('-'))
      ENDIF
C
C  ****  Average energies deposited in bodies..
C
      DF=1.0D0/TOTN
      WRITE(26,3011)
 3011 FORMAT(/3X,'Average deposited energies (bodies):')
      DO KB=1,NBODY
        IF(MATER(KB).NE.0) THEN
          QER=3.0D0*DF*SQRT(ABS(TDEBO2(KB)-TDEBO(KB)**2*DF))
          QAV=TDEBO(KB)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(26,3012) KB,QAV,QER,EFFIC
        ENDIF
      ENDDO
 3012 FORMAT(6X,'Body ',I4, ' ...... ',1P,E13.6,' +-',E8.1,' eV',4X,
     1  '(effic. =',E9.2,')')
C
C  ****  Average energies 'collected' by impact detectors.
C
      IF(NDIDEF.GT.0) THEN
        WRITE(26,3013)
 3013   FORMAT(/3X,'Average incoming energies (impact detectors):')
        DO KD=1,NDIDEF
          QER=3.0D0*DF*SQRT(ABS(TDID2(KD)-TDID(KD)**2*DF))
          QAV=TDID(KD)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(26,3014) KD,QAV,QER,EFFIC
        ENDDO
      ENDIF
 3014 FORMAT(6X,'Detector #',I2,' ... ',1P,E13.6,' +-',E8.1,' eV',4X,
     1  '(effic. =',E9.2,')')
C
C  ****  Average deposited energies in energy-deposition detectors.
C
      IF(NDEDEF.GT.0) THEN
        WRITE(26,3015)
 3015   FORMAT(/3X,'Average deposited energies (energy detectors):')
        DO KD=1,NDEDEF
          QER=3.0D0*DF*SQRT(ABS(TDED2(KD)-TDED(KD)**2*DF))
          QAV=TDED(KD)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(26,3014) KD,QAV,QER,EFFIC
        ENDDO
      ENDIF
C
C  ************  Energy spectrum of the source.
C
      IF(ISPEC.EQ.1) THEN
        OPEN(9,FILE='psource.dat')
        WRITE(9,9000)
 9000 FORMAT(
     1  1X,'#  Results from PENMAIN. ',
     1 /1X,'#  Source energy spectrum.',
     1 /1X,'#  1st column: E (eV). 2nd column: spectrum (1/eV).',
     1 /1X,'#  3rd and 4th columns: simul. pdf limits (eV/particle).',/)
        PTOT=0.0D0
        WRITE(9,'(1X,1P,4E14.6)') ES(1),PTOT,PTOT,PTOT
        DO KE=1,NSEB
          PTOT=PTOT+PTS(KE)
        ENDDO
        DO KE=1,NSEB
          YAV=SHIST(KE)*DF
          YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV)*DF))
          EINTL=ES(KE+1)-ES(KE)
          IF(EINTL.GT.1.0D-15) THEN
            FACT=1.0D0/EINTL
          ELSE
            FACT=1.0D15
          ENDIF
          WRITE(9,'(1X,1P,4E14.6)') ES(KE),PTS(KE)*FACT/PTOT,
     1      (YAV-YERR)*FACT,(YAV+YERR)*FACT
          WRITE(9,'(1X,1P,4E14.6)') ES(KE+1),PTS(KE)*FACT/PTOT,
     1      (YAV-YERR)*FACT,(YAV+YERR)*FACT
        ENDDO
        PTOT=0.0D0
        WRITE(9,'(1X,1P,4E14.6)') ES(NSEB+1),PTOT,PTOT,PTOT
        CLOSE(9)
      ENDIF
C
C  ************  Energy distributions of emerging particles.
C
C  ****  Transmitted electrons.
      OPEN(9,FILE='pm_energy_el_trans.dat')
      WRITE(9,9110)
 9110 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of transmitted electrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(1)
        YERR=3.0D0*SQRT(ABS(PDE2(1,1,K)-PDE(1,1,K)**2*DF))
        YAV=PDE(1,1,K)*DF/BSE(1)
        YERR=YERR*DF/BSE(1)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Backscattered electrons.
      OPEN(9,FILE='pm_energy_el_back.dat')
      WRITE(9,9120)
 9120 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of backscattered electrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(1)
        YERR=3.0D0*SQRT(ABS(PDE2(1,2,K)-PDE(1,2,K)**2*DF))
        YAV=PDE(1,2,K)*DF/BSE(1)
        YERR=YERR*DF/BSE(1)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Transmitted photons.
      OPEN(9,FILE='pm_energy_ph_trans.dat')
      WRITE(9,9130)
 9130 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of transmitted photons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(2)
        YERR=3.0D0*SQRT(ABS(PDE2(2,1,K)-PDE(2,1,K)**2*DF))
        YAV=PDE(2,1,K)*DF/BSE(2)
        YERR=YERR*DF/BSE(2)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Backscattered photons.
      OPEN(9,FILE='pm_energy_ph_back.dat')
      WRITE(9,9140)
 9140 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of backscattered photons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(2)
        YERR=3.0D0*SQRT(ABS(PDE2(2,2,K)-PDE(2,2,K)**2*DF))
        YAV=PDE(2,2,K)*DF/BSE(2)
        YERR=YERR*DF/BSE(2)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Transmitted positrons.
      OPEN(9,FILE='pm_energy_po_trans.dat')
      WRITE(9,9150)
 9150 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of transmitted positrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(3)
        YERR=3.0D0*SQRT(ABS(PDE2(3,1,K)-PDE(3,1,K)**2*DF))
        YAV=PDE(3,1,K)*DF/BSE(3)
        YERR=YERR*DF/BSE(3)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Backscattered positrons.
      OPEN(9,FILE='pm_energy_po_back.dat')
      WRITE(9,9160)
 9160 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of backscattered positrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(3)
        YERR=3.0D0*SQRT(ABS(PDE2(3,2,K)-PDE(3,2,K)**2*DF))
        YAV=PDE(3,2,K)*DF/BSE(3)
        YERR=YERR*DF/BSE(3)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
C  ************  Angular distributions of emerging particles.
C
      OPEN(9,FILE='pm_angle_el.dat')
      WRITE(9,9210)
 9210 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Angular distribution of emerging electrons.',
     1 /1X,'#  1st and 2nd columns: THETA and PHI (deg).',
     1 /1X,'#  3rd and 4th columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*(BSPH*DE2RA)
        DO L=1,NBPH
          YY=(L-0.5D0)*BSPH
          YERR=3.0D0*SQRT(ABS(PDA2(1,K,L)-PDA(1,K,L)**2*DF))
          YAV=PDA(1,K,L)*DF/DSANG
          YERR=YERR*DF/DSANG
          WRITE(9,'(1X,1P,6E14.6)')
     1       XX,YY,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        WRITE(9,*) '   '
      ENDDO
      CLOSE(9)
C
      OPEN(9,FILE='pm_angle_ph.dat')
      WRITE(9,9220)
 9220 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Angular distribution of emerging photons.',
     1 /1X,'#  1st and 2nd columns: THETA and PHI (deg).',
     1 /1X,'#  3rd and 4th columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*(BSPH*DE2RA)
        DO L=1,NBPH
          YY=(L-0.5D0)*BSPH
          YERR=3.0D0*SQRT(ABS(PDA2(2,K,L)-PDA(2,K,L)**2*DF))
          YAV=PDA(2,K,L)*DF/DSANG
          YERR=YERR*DF/DSANG
          WRITE(9,'(1X,1P,6E14.6)')
     1       XX,YY,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        WRITE(9,*) '   '
      ENDDO
      CLOSE(9)
C
      OPEN(9,FILE='pm_angle_po.dat')
      WRITE(9,9230)
 9230 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Angular distribution of emerging positrons.',
     1 /1X,'#  1st and 2nd columns: THETA and PHI (deg).',
     1 /1X,'#  3rd and 4th columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*(BSPH*DE2RA)
        DO L=1,NBPH
          YY=(L-0.5D0)*BSPH
          YERR=3.0D0*SQRT(ABS(PDA2(3,K,L)-PDA(3,K,L)**2*DF))
          YAV=PDA(3,K,L)*DF/DSANG
          YERR=YERR*DF/DSANG
          WRITE(9,'(1X,1P,6E14.6)')
     1       XX,YY,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        WRITE(9,*) '   '
      ENDDO
      CLOSE(9)
C
C  ************  Spectra from impact detectors.
C
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
          OPEN(9,FILE=SPCDIO(ID))
          WRITE(9,9310) ID
 9310 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from impact detector #',I3,
     1 /1X,'#  1st column: particle energy (eV).',
     1 /1X,'#  2nd column: probability density (1/(eV*particle)).',
     1 /1X,'#  3rd column: statistical uncertainty (3 sigma).',/)
          DO J=1,NDICH(ID)
            XX=EDIL(ID)+(J-0.5D0)*BDIE(ID)
            YERR=3.0D0*SQRT(ABS(DIT2(ID,J)-DIT(ID,J)**2*DF))
            YAV=DIT(ID,J)*DF/BDIE(ID)
            YERR=YERR*DF/BDIE(ID)
            WRITE(9,'(1X,1P,3E14.6)')
     1        XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
          ENDDO
          CLOSE(9)
        ENDDO
      ENDIF
C
C  ************  Spectra from energy-deposition detectors.
C
      IF(NDEDEF.GT.0) THEN
        DO ID=1,NDEDEF
          WRITE(BUF2,'(I5)') ID
          OPEN(9,FILE=SPCDEO(ID))
          WRITE(9,9320) ID
 9320 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from energy-deposition',
     1  ' detector #',I3,
     1 /1X,'#  WARNING: May be strongly biased if interaction ',
     1  'forcing is used!',
     1 /1X,'#  1st column: deposited energy (eV).',
     1 /1X,'#  2nd column: probability density (1/(eV*particle)).',
     1 /1X,'#  3rd column: statistical uncertainty (3 sigma).',/)
          DO J=1,NDECH(ID)
            XX=EDEL(ID)+(J-0.5D0)*BDEE(ID)
            YERR=3.0D0*SQRT(ABS(DET2(ID,J)-DET(ID,J)**2*DF))
            YAV=DET(ID,J)*DF/BDEE(ID)
            YERR=YERR*DF/BDEE(ID)
            WRITE(9,'(1X,1P,3E14.6)')
     1        XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
          ENDDO
          CLOSE(9)
        ENDDO
      ENDIF
C
C  ************  Dose distributions.
C
      IF(IDOSE.NE.0) THEN
C
C  ****  Depth-dose distribution.
C
        OPEN(9,FILE='pm_depth_dose.dat')
        WRITE(9,9410)
 9410   FORMAT(
     1     1X,'#  Results from PENMAIN. Depth-dose distribution.',
     1    /1X,'#  1st column: z coordinate (cm).',
     1    /1X,'#  2nd column: depth-dose (eV/(g/cm**2)).',
     1    /1X,'#  3rd column: statistical uncertainty (3 sigma).',/)
        DO I3=1,NDB(3)
          ZZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          YAV=0.0D0
          YAV2=0.0D0
          DO I1=1,NDB(1)
            DO I2=1,NDB(2)
              YAV=YAV+DOSE(I1,I2,I3)
              YAV2=YAV2+DOSE2(I1,I2,I3)
            ENDDO
          ENDDO
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF/BDOSE(3)
          YERR=YERR*DF/BDOSE(3)
          WRITE(9,'(1X,1P,3E14.6)')
     1      ZZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        CLOSE(9)
C
        VOXEL=BDOSE(1)*BDOSE(2)*BDOSE(3)
        DMAX=0.0D0
        I1M=1
        I2M=1
        I3M=1
        DO I1=1,NDB(1)
          DO I2=1,NDB(2)
            DO I3=1,NDB(3)
              YAV=DOSE(I1,I2,I3)
              YAV2=DOSE2(I1,I2,I3)
              YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
              YAV=YAV*DF/VOXEL
              YERR=YERR*DF/VOXEL
              IF(YAV.GT.DMAX) THEN
                I1M=I1
                I2M=I2
                I3M=I3
                DMAX=YAV
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
        QAV=DOSE(I1M,I2M,I3M)
        QAV2=DOSE2(I1M,I2M,I3M)
        QER=3.0D0*SQRT(ABS(QAV2-QAV**2*DF))
        QAV=QAV*DF/VOXEL
        QER=QER*DF/VOXEL
        IF(QER.GT.1.0D-10*ABS(QAV)) THEN
          EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
        ELSE
          EFFIC=0.0D0
        ENDIF
        WRITE(26,3020) QAV,QER,EFFIC
 3020   FORMAT(/6X,'Maximum dose ... ',1P,E13.6,' +-',E8.1,' eV/g',2X,
     1    '(effic. =',E9.2,')')
C
        OPEN(9,FILE='pm_3d_dose.dat')
        WRITE(9,9420)
 9420   FORMAT(1X,'#  Results from PENMAIN. 3D dose distribution.')
        WRITE(9,9421) DXL(1),DXU(1)
 9421   FORMAT(1X,'#  Dose-map enclosure:  XL = ',1P,E13.6,
     1    ' cm,  XU = ',E13.6,' cm')
        WRITE(9,9422) DXL(2),DXU(2)
 9422   FORMAT(1X,'#',23X,'YL = ',1P,E13.6,' cm,  YU = ',E13.6,' cm')
        WRITE(9,9423) DXL(3),DXU(3)
 9423   FORMAT(1X,'#',23X,'ZL = ',1P,E13.6,' cm,  ZU = ',E13.6,' cm')
        WRITE(9,9424) NDB(1),NDB(2),NDB(3)
 9424   FORMAT(1X,'#  Numbers of bins:     NBX =',I4,', NBY =',I4,
     1        ', NBZ =',I4,/1X,'#')
        WRITE(9,9425)
 9425   FORMAT(/1X,'#  columns 1 to 3: bin indices IX, IY and IZ',
     1    /1X,'#  4th column: dose (eV/g).',
     1    /1X,'#  5th column: statistical uncertainty (3 sigma).',/)
        DO I1=1,NDB(1)
          DO I2=1,NDB(2)
            DO I3=1,NDB(3)
              YAV=DOSE(I1,I2,I3)
              YAV2=DOSE2(I1,I2,I3)
              YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
              YAV=YAV*DF/VOXEL
              YERR=YERR*DF/VOXEL
              WRITE(9,9426) I1,I2,I3,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
            ENDDO
          ENDDO
        ENDDO
 9426   FORMAT(1X,3I4,1P,E13.5,E9.1)
        CLOSE(9)
C
C  ****  Dose maps on planes perpendicular to the Z-axis (used for
C        plotting).
C
        DO I3=1,NDB(3)
          ZZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          WRITE(BUF2,'(I5)') I3
          IF(I3.LT.10) THEN
            OFILE='pm_2d_dose_'//BUF2(5:5)//'.dat'
          ELSE IF(I3.LT.100) THEN
            OFILE='pm_2d_dose_'//BUF2(4:5)//'.dat'
          ELSE IF(I3.LT.1000) THEN
            OFILE='pm_2d_dose_'//BUF2(3:5)//'.dat'
          ENDIF
          OPEN(9,FILE=OFILE)
          WRITE(9,9430) I3,ZZ
 9430   FORMAT(
     1     1X,'#  Results from PENMAIN. 3D dose distribution.',
     1    /1X,'#  Plane #',I4,'  Z= ',1P,E12.5,' cm',
     1    /1X,'#  columns 1 to 3: I1, I2, I3',
     1    /1X,'#  columns 4-6: x, y, z (cm).',
     1    /1X,'#  7th column: dose (eV/g).',
     1    /1X,'#  8th column: statistical uncertainty (3 sigma).',/)
          DO I1=1,NDB(1)
            XX=DXL(1)+(I1-0.5D0)*BDOSE(1)
            DO I2=1,NDB(2)
              YY=DXL(2)+(I2-0.5D0)*BDOSE(2)
              YAV=DOSE(I1,I2,I3)
              YAV2=DOSE2(I1,I2,I3)
              YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
              YAV=YAV*DF/VOXEL
              YERR=YERR*DF/VOXEL
              WRITE(9,'(3I4,1X,1P,4E13.5,E9.1)')
     1          I1,I2,I3,XX,YY,ZZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
            ENDDO
            WRITE(9,*) '   '
          ENDDO
          CLOSE(9)
        ENDDO
      ENDIF
C
      WRITE(26,3030) ISEED1,ISEED2
 3030 FORMAT(/3X,'Last random seeds = ',I10,' , ',I10)
C
C  ****  Continue if the DUMPTO option is active and IRETRN=1.
C
      IF(IRETRN.EQ.1) THEN
        WRITE(26,'(/3X,72(''-''))')
        WRITE(6,3040) SHN
 3040   FORMAT(2X,'Number of simulated showers =',1P,E14.7)
        GO TO 101
      ENDIF
C
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
          IF(IPSF(ID).NE.0) THEN
            IPSFU=20+ID
            CLOSE(IPSFU)
          ENDIF
        ENDDO
      ENDIF
C
      WRITE(26,'(//3X,''***  END  ***'')')
      CLOSE(26)
      STOP
      END


C  *********************************************************************
C                       SUBROUTINE GCONE
C  *********************************************************************
      SUBROUTINE GCONE(UF,VF,WF)
C
C  This subroutine samples a random direction uniformly within a cone
C  with central axis in the direction (THETA,PHI) and aperture ALPHA.
C  Parameters are initialised by calling subroutine GCONE0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI)
C  ****  Parameters for sampling directions within a cone.
      COMMON/CGCONE/CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER
C
      EXTERNAL RAND
C  ****  Define a direction relative to the z-axis.
      WT=CAPER+(1.0D0-CAPER)*RAND(1.0D0)
      DF=TWOPI*RAND(2.0D0)
      SUV=SQRT(1.0D0-WT*WT)
      UT=SUV*COS(DF)
      VT=SUV*SIN(DF)
C  **** Rotate to the beam axis direction.
      UF=CPCT*UT-SPHI*VT+CPST*WT
      VF=SPCT*UT+CPHI*VT+SPST*WT
      WF=-STHE*UT+CTHE*WT
C  ****  Ensure normalisation.
      DXY=UF*UF+VF*VF
      DXYZ=DXY+WF*WF
      IF(ABS(DXYZ-1.0D0).GT.1.0D-14) THEN
        FNORM=1.0D0/SQRT(DXYZ)
        UF=FNORM*UF
        VF=FNORM*VF
        WF=FNORM*WF
       ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GCONE0
C  *********************************************************************
      SUBROUTINE GCONE0(THETA,PHI,ALPHA)
C
C  This subroutine defines the parameters for sampling random directions
C  uniformly within a cone with axis in the direction (THETA,PHI) and
C  aperture ALPHA (in rad).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Parameters for sampling directions within a cone.
      COMMON/CGCONE/CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER
C
      CPCT=COS(PHI)*COS(THETA)
      CPST=COS(PHI)*SIN(THETA)
      SPCT=SIN(PHI)*COS(THETA)
      SPST=SIN(PHI)*SIN(THETA)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      STHE=SIN(THETA)
      CTHE=COS(THETA)
      CAPER=COS(ALPHA)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE N2CH10
C  *********************************************************************
      SUBROUTINE N2CH10(N,L,NDIG)
C
C  This subroutine writes a positive integer number N in a 10-character
C  string L. The number is written at the left, followed by unused blank
C  characters. NDIG is the number of decimal digits of N.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*10 L,LT
C
      ND=MAX(N,0)
      WRITE(LT,'(I10)') ND
      DO I=1,10
        IF(LT(I:I).NE.' ') THEN
          IT=I-1
          GO TO 1
        ENDIF
      ENDDO
      IT=9
    1 CONTINUE
      L=LT(IT+1:10)
      NDIG=10-IT
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RDPSF
C  *********************************************************************
      SUBROUTINE RDPSF(IUNIT,PSFREC,KODE)
C
C  This subroutine reads the phase space file. When KODE=0, a valid
C  particle record has been read and copied into the character variable
C  PSFREC. KODE=-1 indicates that the program tried to read after the
C  end of the phase-space file. Blank lines and lines starting with the
C  pound sign (#) are considered as comments and ignored.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*200 PSFREC
C
      KODE=0
    1 CONTINUE
      READ(IUNIT,'(A200)',END=2,ERR=1) PSFREC
      READ(PSFREC,*,ERR=1,END=1) ITEST
      RETURN
    2 CONTINUE
      KODE=-1   ! End of file
      RETURN
      END
