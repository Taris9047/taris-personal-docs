      INCLUDE 'penelope.f'  ! Inserted to simplify compilation.
      INCLUDE 'penvared.f'
      INCLUDE 'timer.f'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C            PPPPP   EEEEEE  N    N   CCCC   Y    Y  L                 C
C            P    P  E       NN   N  C    C  Y    Y  L                 C
C            P    P  E       N N  N  C        Y  Y   L                 C
C            PPPPP   EEEE    N  N N  C         YY    L                 C
C            P       E       N   NN  C    C    YY    L                 C
C            P       EEEEEE  N    N   CCCC     YY    LLLLLL            C
C                                                                      C
C                                                   (version 2005).    C
C                                                                      C
C   PENELOPE's steering main program for Monte Carlo simulation of     C
C   electron-photon showers in multilayered cylindrical structures.    C
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
C  The material system consists of one or several layers of given thick-
C  nesses. Each layer contains a number of concentric homogeneous rings
C  of given compositions and radii (and thickness equal to that of the
C  layer). The layers are perpendicular to the Z-axis and the centre of
C  the rings in each layer is specified by giving its X and Y coordi-
C  nates. When all the centres are on the Z-axis, the geometrical struc-
C  ture is symmetrical about the Z-axis.
C
C  Primary particles of a given kind, KPARP, are emitted from the
C  active volume of the source, either with fixed energy SE0 or with a
C  specified (histogram-like) energy spectrum. The initial direction of
C  the primary particles is sampled uniformly inside a cone of (semi-)
C  aperture SALPHA and with central axis in the direction (STHETA,SPHI).
C  Thus, SALPHA=0 defines a monodirectional source and SALPHA=180 deg
C  corresponds to an isotropic source.
C     The program can simulate two different types of sources:
C  a) An external source with uniform activity over a cylindrical
C     volume, which is defined separately from the geometry of the
C     material system.
C  b) A set of internal sources spread over specified bodies, each one
C     with uniform activity concentration. The original position of the
C     primary particle is sampled uniformly within the active cylinder
C     or ring, which is selected randomly with probability proportional
C     to the total activity in its volume.
C
C     In the distributed form of the program, we assume that both the
C  source and the material structure are symmetrical about the Z-axis,
C  because this eliminates the dependence on the azimuthal angle PHI. It
C  is possible to consider geometries that are not axially symmetrical,
C  but then the program only delivers values averaged over PHI. To
C  obtain the dependence of the angular distributions on the azimuthal
C  angle, we need to increase the value of the parameter NBPHM (the
C  maximum number of bins for PHI, which is set equal to 1 in the
C  distributed source file) and, in the input data file, set NBPH equal
C  to NBPHM.
C
C  ************  The input data file.
C
C  Each line in the input data file consists of a 6-character keyword
C  (columns 1-6) followed either by numerical data (in free format) or
C  by a character string, which start at the 8th column. Keywords are
C  explicitly used/verified by the program (which is case sensitive!).
C  Notice also that the order of the data lines is important. The
C  keyword '______' (6 blanks) indicates comment lines, these can be
C  placed anywhere after the end of the geometry definition list. The
C  program ignores any text following the first blank after the last
C  numerical datum, or after the character string, in each line (thus,
C  in the table given below, the comments in square brackets are
C  ignored by the program). Lines with certain keywords (e.g., 'SPECTR')
C  can appear an arbitrary number of times, limited only by the allo-
C  cated amount of memory.
C
C  The program assigns default values to many input variables; lines
C  that declare default values may be eliminated from the input file.
C  Notice that the default source is a pencil beam that moves upwards
C  along the Z-axis.
C
C  The structure of the input file is the following,
C
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C  TITLE  Title of the job, up to 120 characters.
C  GSTART >>>>>>>> Beginning of the geometry definition list.
C  LAYER  ZLOW,ZHIG                               [Z_lower and Z_higher]
C  CENTRE XCEN,YCEN                              [X_centre and Y_centre]
C  CYLIND M,RIN,ROUT                     [Material, R_inner and R_outer]
C  GEND   <<<<<<<< End of the geometry definition list.
C
C         >>>>>>>> Source definition.
C  SKPAR  KPARP    [Primary particles: 1=electron, 2=photon, 3=positron]
C  SENERG SE0              [Initial energy (monoenergetic sources only)]
C  SPECTR Ei,Pi                 [E bin: lower-end and total probability]
C  SEXTND KL,KC,RELAC [Extended source in KL,KC and rel. activity conc.]
C  STHICK STHICK                                         [Source height]
C  SRADII SRIN,SROUT                      [Source inner and outer radii]
C  SPOSIT SX0,SY0,SZ0                 [Coordinates of the source centre]
C  SDIREC STHETA,SPHI               [Beam axis direction angles, in deg]
C  SAPERT SALPHA                                 [Beam aperture, in deg]
C
C         >>>>>>>> Material data and simulation parameters.
C  NMAT   NMAT                   [Number of different materials, .le.10]
C  SIMPAR M,EABS(1:3,M),C1,C2,WCC,WCR   [Sim. parameters for material M]
C  PFNAME mat_filename.ext     [Material definition file, 20 characters]
C
C         >>>>>>>> Interaction forcing.
C  IFORCE KL,KC,KPAR,ICOL,FORCER,WLOW,WHIG         [Interaction forcing]
C
C         >>>>>>>> Counter array dimensions and pdf ranges.
C  NBE    EMIN,EMAX,NBE              [Energy interval and no. of E-bins]
C  NBTH   NBTH                   [No. of bins for the polar angle THETA]
C  NBPH   NBPH                 [No. of bins for the azimuthal angle PHI]
C  NBZ    NBZ                         [No. of bins for the Z-coordinate]
C  NBR    NBR                                       [No. of radial bins]
C  NBTL   TLMIN,TLMAX,NBTL    [Track-length interval and no. of TL-bins]
C
C         >>>>>>>> Additional distributions to be tallied.
C  ABSEN  MAT           [Tally the distr. of absorbed E in material MAT]
C  DOSE2D KL,KC,NZ,NR    [Tally 2D dose and charge dists. in body KL,KC]
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
C  >>>>>>>> Geometry definition list.
C
C           This block of lines defines the geometrical structure. The
C           only allowed keywords in the geometry definition list are
C           'GSTART', 'LAYER_', 'CENTRE', 'CYLIND' and 'GEND__' (notice
C           the blanks). The definition list must begin with the
C           'GSTART' line and terminate with the 'GEND__' line. The
C           second line must be a 'LAYER_' line, which initiates the
C           definition of the layer structure. Each 'LAYER_' line is
C           followed by a 'CENTRE' line (optional) and by one or several
C           'CYLIND' lines, which define the various concentric rings in
C           the layer; a layer may be void. No blank lines are allowed
C           in the geometry definition list.
C
C           Layers must be defined in increasing order of heights, from
C           bottom to top of the structure. If the 'CENTRE' line is not
C           entered, cylinders are assumed to be centered on the Z-axis
C           (XCEN=YCEN=0.0). Cylinders have to be defined in increasing
C           radial order, from the centre to the periphery. The two
C           lengths in each 'LAYER_' and 'CYLIND' line must be entered
C           in increasing order. All numerical data are read in free
C           format.
C
C             DEFAULT: none
C           --> The geometry definition list can be debugged/visualized
C           with the code GVIEWC (operable only under Windows 9x/NT).
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
C  SEXTND : For internal extended sources, this line defines an active
C           body KL,KC (the cylinder KC in layer KL) and its relative
C           activity concentration, RELAC.
C             DEFAULT: none.
C
C           NOTE: The labels KL,KC that identify each body are defined
C           by the ordering in the input geometry list. These labels
C           are written in the output geometry report.
C
C  STHICK : For an external source, thickness (height) of the
C           active volume of the source (cylinder or ring).
C             DEFAULT: STHICK=0.0
C  SRADII : For an external source, inner and outer radii of the active
C           volume (cylinder or ring).
C             DEFAULTS: SRIN=0.0, SROUT=0.0
C  SPOSIT : For an external source, coordinates of the centre of the
C           source volume.
C             DEFAULT: SX0=SY0=0, SZ0=-1.0E15
C
C           The initial position of each primary particle is sampled
C           uniformly within the volume of the source.
C
C  SDIREC : Polar and azimuthal angles of the source beam axis direc-
C           tion, in deg.
C             DEFAULTS: STHETA=0.0, SPHI=0.0
C  SAPERT : Angular aperture of the source beam, in deg.
C             DEFAULT: SALPHA=0.0
C
C  >>>>>>>> Material data and simulation parameters.
C
C  NMAT__ : Number of different materials (up to 10 with the original
C           program dimensions). Materials are identified by their
C           ordering in PENELOPE's input material data file.
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
C             EPMAX is the maximum energy of all particles found in the
C                   simulation (depends on the source energies).
C
C  PFNAME : Name of PENELOPE's input material data file (20 characters).
C             DEFAULT: none
C
C  >>>>>>>> Counter array dimensions and tallied ranges.
C
C  NBE___ : Limits of the interval where energy distributions of
C           emerging particles are tallied and number of energy bins.
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
C  NBZ___ : Number of bins for the Z-coordinate.
C             DEFAULT: NBZ=100
C
C  NBR___ : Number of bins for the radial variable, R=SQRT(X*X+Y*Y).
C             DEFAULT: NBR=100
C
C  NBTL__ : Limits of the interval where track-length distributions of
C           primary particles are tallied. Number of track-length bins.
C             DEFAULT: TLMIN=0, TLMAX=2.5*RANGE(EPMAX,KPARP,1), NBTL=100
C
C  >>>>>>>> Additional distributions to be tallied.
C
C  ABSEN_ : Indicates a material M for which we require the code to
C           tally the distribution of absorbed energy. Up to three
C           different materials can be selected, a separate line for
C           each one.
C             DEFAULT: off
C
C  DOSE2D : The program will tally 2D, depth-radius, dose and deposit-
C           ed charge distributions in the body KL,KC (i.e. the cylinder
C           KC of layer KL). The numbers NZ and NR of Z- and R-bins have
C           to be specified by the user, they must be in the interval
C           from 1 to 100. Up to three different bodies can be selected,
C           a separate line for each body.
C             DEFAULT: off
C
C           NOTE: The labels KL,KC that identify a given cylinder are
C           defined by the ordering in the input geometry list. These
C           labels are written in the output geometry report.
C
C  >>>>>>>> Interaction forcing.
C
C  IFORCE : Activates forcing of interactions of type ICOL of particles
C           KPAR in cylinder KC of layer KL. FORCER is the forcing
C           factor, which must be larger than unity. WLOW and WHIG are
C           the lower and upper limits of the weight window where
C           interaction forcing is applied.
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
C  to predict. Please, do tentative runs with different FORCERR values
C  and check the efficiency gain (or loss!).
C
C  >>>>>>>> Job properties.
C
C  RESUME : The program will read the dump file named 'dumpfile_1.dat'
C           (20 characters) and resume the simulation from the point
C           where it was left. Use this option very, _VERY_ carefully.
C           Make sure that the input data file is fully consistent with
C           the one used to generate the dump file.
C             DEFAULT: off
C
C  DUMPTO : Generate a dump file named 'dumpfile_2.ext' (20 characters)
C           after completing the simulation run. This allows resuming
C           the simulation later on to improve statistics.
C             DEFAULT: off
C
C           NOTE: If the file 'dumpfile_2.ext' already exists, it is
C           overwritten.
C
C  DUMPP_ : When the DUMPTO option is activated, simulation results are
C           written in the output files every DUMPP seconds. This option
C           is useful to check the progress of long simulations. It also
C           allows running the program with a long execution time and
C           stopping it when the required statistical uncertainty has
C           been reached.
C             DEFAULT: DUMPP=1.0E15
C
C  NSIMSH : Desired number of simulated showers.
C             DEFAULT: DSHN=2.0E9
C
C  RSEED_ : Seeds of the random number generator.
C             DEFAULT: ISEED1=12345; ISEED2=54321
C
C  TIME__ : Allotted simulation time, in sec.
C             DEFAULT: TIMEA=100.0D0
C
C  The program is aborted when an incorrect input datum is found. The
C  conflicting quantity usually appears in the last line of the output
C  file. If the trouble is with arrays having dimensions smaller than
C  required, the program indicates how the problem can be solved (this
C  usually requires editing the source file, be careful).
C
C  The clock subroutine (TIMER) may have to be adapted to your specific
C  computer-compiler configuration; standard FORTRAN 77 does not provide
C  timing tools. However, the routines in module 'TIMER.F' do work for
C  many FORTRAN compilers.
C
C  ************  Generating the executable PENCYL and running it.
C
C  To generate the executable binary file PENCYL.EXE, compile and link
C  the FORTRAN 77 source files PENCYL.F, PENELOPE.F, PENVARED.F and
C  TIMER.F. For example, if you are using the g77 compiler under
C  Windows, place these four files in the same directory, open a DOS
C  window and from that directory enter the command
C    `g77 -Wall -O PENCYL.F -o PENCYL.EXE'
C  (The same, with file names in lowercase, should work under Linux).
C
C  To run PENCYL, you have to generate an input data file, let's call it
C  PENCYL.IN, and the corresponding material data file. Place these
C  two files and the binary file PENCYL.EXE in the same directory and,
C  from there, issue the command
C    'PENCYL.EXE < PENCYL.IN'
C
C  The calculated distributions are written in separate files, whose
C  names start with 'pc_' (for PenCyl) and have the extension '.dat'.
C  These files are in a format suited for direct visualisation with
C  GNUPLOT (version 4.0).
C
C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      CHARACTER*2 LIT
      CHARACTER*20 PFILE,OFILE,PFILED,PFILER
      CHARACTER*23 DATE23
      CHARACTER*120 TITLE,TITLE2,BUFFER
      CHARACTER*6 KWORD,
     1  KWTITL,KWKPAR,KWSENE,KWSPEC, KWSEXT,KWSHEI,KWSRAD,KWSPOS,
     1  KWSDIR,KWSAPE,KWNMAT,KWSIMP, KWPFNA,KWNBE ,KWNBTH,KWNBPH,
     1  KWNBZ ,KWNBR ,KWABSE,KWNBTL, KWDO2D,KWIFOR,KWRESU,KWDUMP,
     1  KWDMPP,KWNSIM,KWTIME,KWRSEE,KWCOMM
      PARAMETER(
     1  KWTITL='TITLE ',KWKPAR='SKPAR ',KWSENE='SENERG',KWSPEC='SPECTR',
     1  KWSEXT='SEXTND',KWSHEI='STHICK',KWSRAD='SRADII',KWSPOS='SPOSIT',
     1  KWSDIR='SDIREC',KWSAPE='SAPERT',KWNMAT='NMAT  ',KWSIMP='SIMPAR',
     1  KWPFNA='PFNAME',KWNBE ='NBE   ',KWNBTH='NBTH  ',KWNBPH='NBPH  ',
     1  KWNBZ ='NBZ   ',KWNBR ='NBR   ',KWNBTL='NBTL  ',KWABSE='ABSEN ',
     1  KWDO2D='DOSE2D',KWIFOR='IFORCE',KWRESU='RESUME',KWDUMP='DUMPTO',
     1  KWDMPP='DUMPP ',KWNSIM='NSIMSH',KWTIME='TIME  ',KWRSEE='RSEED ',
     1  KWCOMM='      ')
C
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=REV+REV)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI,
     1  RA2DE=180.0D0/PI, DE2RA=PI/180.0D0)
C
C  ************  Main-PENELOPE commons.
C
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      DIMENSION RHOI(MAXMAT)
C
C  ****  Cylindrical geometry.
C
      PARAMETER (NLAM=100,NCYM=50,NBDM=NLAM*NCYM)
      COMMON/CYLGEO/XG(NLAM),YG(NLAM),ZG(NLAM),RG(NLAM,NCYM),
     1  RG2(NLAM,NCYM),RMAX,RMAX2,IBOD(NLAM,NCYM),MATER(NLAM,NCYM),
     1  ILAY(NBDM),ICYL(NBDM),NLAY,NCYL(NLAM),NBOD
      COMMON/CYLAUX/INOUT,KLAY,KCYL
      DIMENSION DSMAX(NBDM)
C
C  ****  Source.
C
C  ****  Extended sources.
      DIMENSION RELAC(NBDM),WSB(NBDM),KLS(NBDM),KCS(NBDM)
C  ****  Source energy spectrum.
      PARAMETER (NSEBM=200)
      DIMENSION ES(NSEBM),PTS(NSEBM),IAS(NSEBM),FS(NSEBM),SHIST(NSEBM)
      DATA SHIST/NSEBM*0.0D0/
C
C  ************  Counters.
C
C  ****  Primary particle counters.
      DIMENSION
     1  PRIM(3),PRIM2(3),DPRIM(3),     ! Numbers of IEXIT particles.
     1  SEC(3,3),SEC2(3,3),DSEC(3,3),  ! Generated secondary particles.
     1  AVTL(3),AVTL2(3),     ! Total track length (in material media).
     1  AVW(2),AVW2(2),DAVW(2),        ! Final polar director cosine.
     1  AVA(2),AVA2(2),DAVA(2),        ! Final polar angle.
     1  AVE(2),AVE2(2),DAVE(2),        ! Final energy.
     1  AVI(8),AVI2(8),DAVI(8)         ! Number of interaction events.
C
      DATA PRIM,PRIM2,SEC,SEC2,AVTL,AVTL2,AVW,AVW2,AVA,AVA2,AVE,AVE2,
     1  AVI,AVI2/58*0.0D0/
C
      DIMENSION WSEC(3,3),WSEC2(3,3),WAVTL(3),WAVTL2(3),WAVW(2),
     1  WAVW2(2),WAVA(2),WAVA2(2),WAVE(2),WAVE2(2),WAVI(8),WAVI2(8)
C
C  ****  Continuous distributions (default output).
C
      PARAMETER (NBEM=100,NBTHM=90,NBPHM=1,NBZM=100,NBRM=100,NBTLM=100)
      DIMENSION BSE(3),BSZ(NLAM)
C  ---- Deposited energies in various bodies.
      DIMENSION DEBO(NBDM),TDEBO(NBDM),TDEBO2(NBDM)
      DATA TDEBO,TDEBO2/NBDM*0.0D0,NBDM*0.0D0/
C  ----  Track length distributions.
      PARAMETER (NBTL6=6*NBTLM)
      DIMENSION PDTL(3,NBTLM),PDTL2(3,NBTLM)
      DATA PDTL,PDTL2/NBTL6*0.0D0/
C  ----  Angular distributions of emerging particles.
      PARAMETER (NBPH12=12*NBTHM*NBPHM)
      DIMENSION PDA(3,NBTHM,NBPHM),PDA2(3,NBTHM,NBPHM),
     1          PDAP(3,NBTHM,NBPHM),LPDA(3,NBTHM,NBPHM)
      DATA PDA,PDA2,PDAP,LPDA/NBPH12*0.0D0/
C  ----  Energy distributions of emerging particles.
      PARAMETER (NBE12=12*NBEM)
      DIMENSION PDE(3,2,NBEM),PDE2(3,2,NBEM),PDEP(3,2,NBEM),
     1          LPDE(3,2,NBEM)
      DATA PDE,PDE2,PDEP,LPDE/NBE12*0.0D0,NBE12*0.0D0/
C  ----  Depth-distributions of deposited energy (dose) and charge.
      PARAMETER (NDZTT=4*NLAM*NBZM)
      DIMENSION DOSZ(NLAM,NBZM),DOSZ2(NLAM,NBZM),
     1          DOSZP(NLAM,NBZM),LDOSZ(NLAM,NBZM)
      DATA DOSZ,DOSZ2,DOSZP,LDOSZ/NDZTT*0.0D0/
      DIMENSION CHRZ(NLAM,NBZM),CHRZ2(NLAM,NBZM),
     1          CHRZP(NLAM,NBZM),LCHRZ(NLAM,NBZM)
      DATA CHRZ,CHRZ2,CHRZP,LCHRZ/NDZTT*0.0D0/
C
C  ****  Continuous distributions (selected by the user).
C
C  ----  Deposited energy distributions (in up to 3 selected materials).
      PARAMETER (NDET=3,NDETT=NDET*NBEM)
      DIMENSION MATAE(NDET),IAED(MAXMAT),AESUM(MAXMAT)
      DIMENSION EDEP(NDET,NBEM),EDEP2(NDET,NBEM)
      DATA EDEP,EDEP2/NDETT*0.0D0,NDETT*0.0D0/
C  ----  3D dose and deposited charge distribution (in up to 3 selected
C        bodies).
      PARAMETER (NDOT=3,NBZR4=4*NDOT*NBZM*NBRM)
      DIMENSION KBDOD(NDOT),NZDOD(NDOT),NRDOD(NDOT),IDOD(NBDM)
      DIMENSION BZDOD(NDOT),BRDOD(NDOT)
      DIMENSION DOSE(NDOT,NBZM,NBRM),DOSE2(NDOT,NBZM,NBRM),
     1          DOSEP(NDOT,NBZM,NBRM),LDOSE(NDOT,NBZM,NBRM)
      DATA DOSE,DOSE2,DOSEP,LDOSE/NBZR4*0.0D0/
C
      DIMENSION CHAR(NDOT,NBZM,NBRM),CHAR2(NDOT,NBZM,NBRM),
     1          CHARP(NDOT,NBZM,NBRM),LCHAR(NDOT,NBZM,NBRM)
      DATA CHAR,CHAR2,CHARP,LCHAR/NBZR4*0.0D0/
C  ****  Interaction forcing parameters.
      PARAMETER (NB=1250)
      DIMENSION IFORCE(NB,3),WLOW(NB,3),WHIG(NB,3)
      COMMON/CFORCE/FORCE(NB,3,8)
C
      EXTERNAL RAND
C
C  ****  Time counter initialization.
C
      CALL TIME0
C
C  ------------------------  Read input data file.
C
      OPEN(26,FILE='pencyl.dat')
      WRITE(26,1100)
 1100 FORMAT(//3X,43('*'),/3X,'**    PROGRAM PENCYL. Input data file. ',
     1  '  **',/3X,43('*'))
C
      CALL PDATET(DATE23)
      WRITE(26,1101) DATE23
 1101 FORMAT(/3X,'Date and time: ',A23)
C
C  ****  Title.
C
      READ(5,'(A6,1X,A120)') KWORD,TITLE
      WRITE(26,'(/3X,A120)') TITLE
C
C  ************  Geometry definition and initialization of tracking
C                routines.
C
      CALL GEOINC(NMATG,5,26)
      IF(NMATG.LT.1.OR.NMATG.GT.MAXMAT) THEN
        WRITE(26,*) 'Conflicting material numbers.'
        STOP 'Conflicting material numbers.'
      ENDIF
      DO KL=1,NLAY
        DO KC=1,NCYL(KL)
          KB=IBOD(KL,KC)
          DSMAX(KB)=MIN(ZG(KL+1)-ZG(KL),RG(KL,KC+1)-RG(KL,KC))/5.0D0
          DSMAX(KB)=MAX(DSMAX(KB),1.0D-8)
        ENDDO
      ENDDO
C
C  ************  Source description.
C
      WRITE(26,1200)
 1200 FORMAT(//3X,72('-'),/3X,'>>>>>>  Source description.')
C
   21 CONTINUE
      READ(5,'(A6,1X,A120)') KWORD,BUFFER
      IF(KWORD.EQ.KWCOMM) GO TO 21
C
      IF(KWORD.EQ.KWKPAR) THEN
        READ(BUFFER,*) KPARP
   22   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 22
      ELSE
        KPARP=1
      ENDIF
      IF(KPARP.LT.1.OR.KPARP.GT.3) KPARP=1
      IF(KPARP.EQ.1) WRITE(26,1210)
 1210 FORMAT(3X,'Primary particles: electrons')
      IF(KPARP.EQ.2) WRITE(26,1211)
 1211 FORMAT(3X,'Primary particles: photons')
      IF(KPARP.EQ.3) WRITE(26,1212)
 1212 FORMAT(3X,'Primary particles: positrons')
C
C  ****  Initial energy of primary particles.
C
      ISPEC=0
      IF(KWORD.EQ.KWSENE) THEN
        NSEB=1
        READ(BUFFER,*) E0
        WRITE(26,1220) E0
 1220   FORMAT(3X,'Initial energy = ',1P,E13.6,' eV')
   23   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 23
      ELSE IF(KWORD.EQ.KWSPEC) THEN
        ISPEC=1
        NSEB=0
   24   CONTINUE
        NSEB=NSEB+1
        READ(BUFFER,*) ES(NSEB),PTS(NSEB)
   25   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 25
        IF(KWORD.EQ.KWSPEC) GO TO 24
      ELSE
        E0=1.0D6
        WRITE(26,1220) E0
      ENDIF
      IF(ISPEC.EQ.1) THEN
        IF(NSEB.GT.1) THEN
          CALL SORT2(ES,PTS,NSEB)
          WRITE(26,1221)
 1221     FORMAT(/3X,'Spectrum:',7X,'I',4X,'E_low(eV)',4x,'E_high(eV)',
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
      IF(E0.LT.100.0D0) THEN
        WRITE(26,*) 'The initial energy E0 is too small.'
        STOP 'The initial energy E0 is too small.'
      ENDIF
      EPMAX=E0
C  ----  Positrons eventually give annihilation gamma-rays. The maximum
C        energy of annihilation photons is .lt. 1.21*(E0+me*c**2).
      IF(KPARP.EQ.3) EPMAX=1.21D0*(E0+5.12D5)
C
      KSOURC=1
C
C  ****  Internal sources (uniform activity within geometry bodies).
C
      NSB=0
   35 CONTINUE
      IF(KWORD.EQ.KWSEXT) THEN
        KSOURC=2
        NSB=NSB+1
        READ(BUFFER,*) KLS(NSB),KCS(NSB),RELAC(NSB)
        WRITE(26,1235) KLS(NSB),KCS(NSB),RELAC(NSB)
 1235   FORMAT(3X,'Extended source:   KL =',I4,', KC =',I4,
     1    ',  rel. activity concentration =',1P,E13.6)
   36   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWSEXT) GO TO 35
        IF(KWORD.EQ.KWCOMM) GO TO 36
      ENDIF
C
      IF(KSOURC.EQ.2) THEN
        WSUM=0.0D0
        DO I=1,NSB
          KL=KLS(I)
          KC=KCS(I)
          WSBI=PI*(RG2(KL,KC+1)-RG2(KL,KC))*(ZG(KL+1)-ZG(KL))*RELAC(I)
          WSUM=WSUM+WSBI
          WSB(I)=WSUM
        ENDDO
        DO I=1,NSB
          WSB(I)=WSB(I)/WSUM
        ENDDO
        WSB(NSB)=1.0D0
      ENDIF
C
C  ****  External (cylindrical) source.
C
      IF(KWORD.EQ.KWSHEI) THEN
        IF(KSOURC.EQ.2) THEN
          WRITE(26,*) 'An extended source has already been defined.'
          STOP 'An extended source has already been defined.'
        ENDIF
        READ(BUFFER,*) STHICK
   26   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 26
      ELSE
        STHICK=0.0D0
      ENDIF
      IF(KSOURC.EQ.1) WRITE(26,1230) STHICK
 1230 FORMAT(/3X,'Active volume:          height = ',1P,E13.6,' cm')
      IF(STHICK.LT.-1.0D-16) THEN
        WRITE(26,*) 'Negative thickness.'
        STOP 'Negative thickness.'
      ENDIF
C
      IF(KWORD.EQ.KWSRAD) THEN
        IF(KSOURC.EQ.2) THEN
          WRITE(26,*) 'An extended source has already been defined.'
          STOP 'An extended source has already been defined.'
        ENDIF
        READ(BUFFER,*) SRIN,SROUT
        KSOURC=1
   27   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 27
      ELSE
        SRIN=0.0D0
        SROUT=0.0D0
      ENDIF
      IF(KSOURC.EQ.1) WRITE(26,1231) SRIN,SROUT
 1231 FORMAT(21X,'inner radius = ',1P,E13.6,' cm',/
     1       21X,'outer radius = ',E13.6,' cm' )
      SRIN2=SRIN**2
      SROI2=SROUT**2-SRIN**2
      IF(SROI2.LT.-1.0D-35) THEN
        WRITE(26,*) 'The source radii are inconsistent.'
        STOP 'The source radii are inconsistent.'
      ENDIF
C
      IF(KWORD.EQ.KWSPOS) THEN
        IF(KSOURC.EQ.2) THEN
          WRITE(26,*) 'An extended source has already been defined.'
          STOP 'An extended source has already been defined.'
        ENDIF
        READ(BUFFER,*) SX0,SY0,SZ0
   28   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 28
      ELSE
        SX0=0.0D0
        SY0=0.0D0
        SZ0=-1.0D15
      ENDIF
      IF(KSOURC.EQ.1) WRITE(26,1232) SX0,SY0,SZ0
 1232 FORMAT(3X,'Coordinates of centre:     SX0 = ',1P,E13.6,
     1  ' cm',/30X,'SY0 = ',E13.6,' cm',/30X,'SZ0 = ',E13.6,' cm')
C
C  ****  Source beam direction and aperture.
C
      IF(KWORD.EQ.KWSDIR) THEN
        READ(BUFFER,*) STHETA,SPHI
   29   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 29
      ELSE
        STHETA=0.0D0
        SPHI=0.0D0
      ENDIF
      WRITE(26,1233) STHETA,SPHI
 1233 FORMAT(3X,'Beam direction angles:   THETA = ',1P,E13.6,' deg',/
     1  30X,'PHI = ',E13.6,' deg')
C
      IF(KWORD.EQ.KWSAPE) THEN
        READ(BUFFER,*) SALPHA
   30   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 30
      ELSE
        SALPHA=0.0D0
      ENDIF
      WRITE(26,1234) SALPHA
 1234 FORMAT(3X,'Beam aperture:',11X,'ALPHA = ',1P,E13.6,' deg')
      CALL GCONE0(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA)
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
      IF(NMAT.LT.1.OR.NMAT.GT.NMATG) THEN
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
C  ****  Initialization of PENELOPE.
C
      IF(KWORD.EQ.KWPFNA) THEN
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
      OPEN(15,FILE=PFILE)
      OPEN(16,FILE='material.dat')
        INFO=3
        CALL PEINIT(EPMAX,NMAT,15,16,INFO)
      CLOSE(UNIT=15)
      CLOSE(UNIT=16)
C  ----  Inverse densities are used to score the local dose.
      DO M=1,NMAT
        RHOI(M)=1.0D0/RHO(M)
      ENDDO
C
C  ************  Variance reduction (only interaction forcing, for the
C                time being...).
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
        WRITE(26,1600)
 1600   FORMAT(//3X,72('-'),/
     1    3X,'>>>>>>  Interaction forcing: FORCE(IBODY,KPAR,ICOL)')
   61   CONTINUE
        READ(BUFFER,*) KL,KC,KPAR,ICOL,FORCER,WWLOW,WWHIG
C  ****  Negative FORCER values are re-interpreted, as described in the
C        heading comments above.
        IF(FORCER.LT.-1.0D-6) THEN
          MM=MATER(KL,KC)
          EVENTS=MAX(ABS(FORCER),1.0D0)
          PLT=PRANGE(E0,KPAR,MM)
          RMFP=PHMFP(E0,KPAR,MM,ICOL)
          FORCER=EVENTS*RMFP/PLT
        ENDIF
        IF(WWLOW.LT.1.0D-6) WWLOW=1.0D-6
        IF(WWHIG.GT.1.0D6) WWHIG=1.0D6
        KB=IBOD(KL,KC)
        IF(KB.LT.1.OR.KB.GT.NB) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          IF(KB.LT.1) STOP 'Inconsistent KL,KC values.'
          WRITE(26,*) 'Inconsistent KL,KC values.'
          STOP 'Inconsistent KL,KC values.'
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
        FORCE(KB,KPAR,ICOL)=FORCER
        WRITE(26,1610) KB,KPAR,ICOL,FORCER,WLOW(KB,KPAR),WHIG(KB,KPAR)
 1610   FORMAT(3X,'FORCE(',I4,',',I1,',',I1,') =',1P,E13.6,
     1    ',  weight window = (',E9.2,',',E9.2,')')
        IFORON=1
   62   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 62
        IF(KWORD.EQ.KWIFOR) GO TO 61
      ENDIF
C
C  ************  Tallied distributions.
C
      WRITE(26,1400)
 1400 FORMAT(//3X,72('-'),/
     1  3X,'>>>>>>  Dimensions of counter arrays.')
C
      IF(KWORD.EQ.KWNBE) THEN
        READ(BUFFER,*) EMIN,EMAX,NBE
   41   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 41
      ELSE
        NBE=NBEM
        EMIN=0.0D0
        EMAX=EPMAX
      ENDIF
      EMIN=MAX(EMIN,0.0D0)
      EMAX=MAX(EMAX,EPMAX)
      WRITE(26,1410) NBE,EMIN,EMAX
 1410 FORMAT(3X,'E:       NBE = ',I3,
     1  ',  EMIN =',1P,E13.6,' eV,  EMAX =',E13.6,' eV')
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
   42   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 42
      ELSE
        NBTH=NBTHM
      ENDIF
      WRITE(26,1420) NBTH
 1420 FORMAT(3X,'Theta:  NBTH = ',I3)
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
   43   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 43
      ELSE
        NBPH=NBPHM
      ENDIF
      WRITE(26,1430) NBPH
 1430 FORMAT(3X,'Phi:    NBPH = ',I3)
      IF(NBTH.LT.1) THEN
        WRITE(26,*) 'Wrong number of PHI bins.'
        STOP 'Wrong number of PHI bins.'
      ELSE IF (NBPH.GT.NBPHM) THEN
        WRITE(26,*) 'NBPH is too large.'
        WRITE(26,*) 'Set the parameter NBPHM equal to ',NBPH
        STOP 'NBPH is too large.'
      ENDIF
C
      IF(KWORD.EQ.KWNBZ) THEN
        READ(BUFFER,*) NBZ
   44   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 44
      ELSE
        NBZ=NBZM
      ENDIF
      WRITE(26,1440) NBZ
 1440 FORMAT(3X,'Z:       NBZ = ',I3)
      IF(NBZ.LT.1) THEN
        WRITE(26,*) 'Wrong number of Z bins.'
        STOP 'Wrong number of Z bins.'
      ELSE IF (NBZ.GT.NBZM) THEN
        WRITE(26,*) 'NBZ is too large.'
        WRITE(26,*) 'Set the parameter NBZM equal to ',NBZ
        STOP 'NBZ is too large.'
      ENDIF
C
      IF(KWORD.EQ.KWNBR) THEN
        READ(BUFFER,*) NBR
   45   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 45
      ELSE
        NBR=NBRM
      ENDIF
      WRITE(26,1450) NBR
 1450 FORMAT(3X,'R:       NBR = ',I3)
      IF(NBR.LT.1) THEN
        WRITE(26,*) 'Wrong number of R bins.'
        STOP 'Wrong number of R bins.'
      ELSE IF (NBR.GT.NBRM) THEN
        WRITE(26,*) 'NBR is too large.'
        WRITE(26,*) 'Set the parameter NBRM equal to ',NBR
        STOP 'NBR is too large.'
      ENDIF
C
      IF(KWORD.EQ.KWNBTL) THEN
        READ(BUFFER,*) TLMIN,TLMAX,NBTL
   46   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 46
      ELSE
        NBTL=NBTLM
        TLMIN=0.0D0
        TLMAX=2.5D0*PRANGE(EPMAX,KPARP,1)
      ENDIF
      TLMIN=MAX(TLMIN,0.0D0)
      WRITE(26,1460) NBTL,TLMIN,TLMAX
 1460 FORMAT(3X,'TL:     NBTL = ',I3,
     1  ', TLMIN =',1P,E13.6,' cm, TLMAX =',E13.6,' cm')
      IF(NBTL.LT.1) THEN
        WRITE(26,*) 'Wrong number of track length bins.'
        STOP 'Wrong number of track length bins.'
      ENDIF
      IF(TLMIN.GE.TLMAX) THEN
        WRITE(26,*) 'Track length interval is too narrow.'
        STOP 'Track length interval is too narrow.'
      ENDIF
C
C  ****  Bin sizes.
C
      BSDE=EMAX/DBLE(NBE)
      BSE(1)=(EMAX-EMIN)/DBLE(NBE)
      BSE(2)=(EMAX-EMIN)/DBLE(NBE)
      BSE(3)=(EMAX-EMIN)/DBLE(NBE)
C
      BSTH=180.0D0/DBLE(NBTH)
      BSPH=360.0D0/DBLE(NBPH)
      BSTL=(TLMAX-TLMIN)/DBLE(NBTL)
C  ----  The factor 1.0000001 serves to ensure that the upper limit of
C  the tallied interval is within the last channel (otherwise, the array
C  dimensions could be exceeded).
      BSE(1)=1.0000001D0*BSE(1)
      BSE(2)=1.0000001D0*BSE(2)
      BSE(3)=1.0000001D0*BSE(3)
      BSDE=1.0000001D0*BSDE
      BSTH=1.0000001D0*BSTH
      BSPH=1.0000001D0*BSPH
      BSTL=1.0000001D0*BSTL
C
      DO KL=1,NLAY
        BSZ(KL)=1.0000001D0*(ZG(KL+1)-ZG(KL))/DBLE(NBZ)
      ENDDO
C
C  ************  Tallied distributions (selected by the user).
C
C  ****  Absorbed energy (pulse-height) distributions in selected
C        materials.
C
      NAED=0
      DO M=1,MAXMAT
        IAED(M)=0
      ENDDO
      IF(KWORD.EQ.KWABSE.OR.KWORD.EQ.KWDO2D) THEN
        WRITE(26,1500)
 1500   FORMAT(//3X,72('-'),/
     1    3X,'>>>>>>  Additional distributions to be tallied.')
      ENDIF
C
      IF(KWORD.EQ.KWABSE) THEN
   51   CONTINUE
        READ(BUFFER,*) MAT
        IF(MAT.LT.0.OR.MAT.GT.NMAT) THEN
          WRITE(26,'(/3X,A120)') BUFFER
          WRITE(26,*) 'MAT is negative or larger than NMAT.'
          STOP 'MAT is negative or larger than NMAT.'
        ENDIF
        DO I=1,NAED
          IF(MAT.EQ.MATAE(I)) THEN
            WRITE(26,'(/3X,A120)') BUFFER
            WRITE(26,*) 'The material has been selected twice.'
            STOP 'The material has been selected twice.'
          ENDIF
        ENDDO
        NAED=NAED+1
        MATAE(NAED)=MAT
        IAED(MAT)=NAED
        WRITE(26,1510) MAT
 1510   FORMAT(3X,'Deposited energy distribution in material',I3)
        IF(IFORON.NE.0) THEN
          WRITE(26,'(3X,''#  WARNING: May be strongly biased when '',
     1      ''interaction forcing is used!'')')
        ENDIF
        IF(NAED.GT.NDET) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Too many selected materials. Increase NDET.'
          STOP 'Too many selected materials.'
        ENDIF
   52   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 52
        IF(KWORD.EQ.KWABSE) GO TO 51
      ENDIF
C
C  ****  2D dose distributions in selected bodies.
C
      NDOD=0
      DO KB=1,NBDM
        IDOD(KB)=0
      ENDDO
      IF(KWORD.EQ.KWDO2D) THEN
   53   CONTINUE
        READ(BUFFER,*) KL,KC,NZ,NR
        WRITE(26,1520) KL,KC,NZ,NR
 1520   FORMAT(3X,'2D dose/charge distribs. in body  KL =',I3,
     1    ', KC =',I3,'  (NZ =',I3,', NR =',I3,')')
        IF(KL.LT.1.OR.KL.GT.NLAY) THEN
          WRITE(26,*) 'KL is not a valid layer index...'
          STOP 'KL is not a valid layer index...'
        ENDIF
        IF(KC.LT.1.OR.KC.GT.NCYL(KL)) THEN
          WRITE(26,*) 'KC is not a valid cylinder index...'
          STOP 'KC is not a valid cylinder index...'
        ENDIF
C
        KB=IBOD(KL,KC)
        IF(NDOD.GT.0) THEN
          DO I=1,NDOD
            IF(KB.EQ.KBDOD(I)) THEN
              WRITE(26,*) 'The body has been selected twice.'
              STOP 'The body has been selected twice.'
            ENDIF
          ENDDO
        ENDIF
        NDOD=NDOD+1
        IF(NDOD.GT.NDOT) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Too many selected bodies. Increase NDOT.'
          STOP 'Too many selected bodies.'
        ENDIF
        IF(NZ.GT.NBZM) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Too many Z-bins. Increase NBZM.'
          STOP 'Too many Z-bins. Increase NBZM.'
        ENDIF
        IF(NR.GT.NBRM) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Too many R-bins. Increase NBRM.'
          STOP 'Too many R-bins. Increase NBRM.'
        ENDIF
        KBDOD(NDOD)=KB
        NZDOD(NDOD)=NZ
        NRDOD(NDOD)=NR
        BZDOD(NDOD)=1.0000001D0*(ZG(KL+1)-ZG(KL))/DBLE(NZ)
        BRDOD(NDOD)=1.0000001D0*(RG(KL,KC+1)-RG(KL,KC))/DBLE(NR)
        IDOD(KB)=NDOD
   54   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 54
        IF(KWORD.EQ.KWDO2D) GO TO 53
      ENDIF
C
C  ************  Job characteristics.
C
      WRITE(26,1700)
 1700 FORMAT(//3X,72('-'),/
     1  3X,'>>>>>>  Job characteristics.')
C
      IRESUM=0
      IF(KWORD.EQ.KWRESU) THEN
        READ(BUFFER,'(A20)') PFILER
        WRITE(26,1710) PFILER
 1710   FORMAT(3X,'Resume simulation from previous dump file: ',A20)
        IRESUM=1
   71   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 71
      ENDIF
C
      IDUMP=0
      DUMPP=1.0D15
      IF(KWORD.EQ.KWDUMP) THEN
        READ(BUFFER,'(A20)') PFILED
        WRITE(26,1720) PFILED
 1720   FORMAT(3X,'Write final counter values on the dump file: ',A20)
        IDUMP=1
   72   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 72
      ENDIF
C
      IF(KWORD.EQ.KWDMPP) THEN
        READ(BUFFER,*) DUMPP
        IF(IDUMP.EQ.1) THEN
          IF(DUMPP.LT.15.0D0) DUMPP=15.0D0
          IF(DUMPP.GT.86400.0D0) DUMPP=86400.0D0
          WRITE(26,1730) DUMPP
 1730     FORMAT(3X,'Dumping period: DUMPP =',1P,E13.6)
        ENDIF
   73   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 73
      ENDIF
C
      IF(KWORD.EQ.KWNSIM) THEN
        READ(BUFFER,*) DSHN
        IF(DSHN.LT.1.0D0) DSHN=2.0D9
   74   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 74
      ELSE
        DSHN=2.0D9
      ENDIF
      WRITE(26,1740) DSHN
 1740 FORMAT(/3X,'Number of showers to be simulated =',1P,E13.6)
C
      IF(KWORD.EQ.KWRSEE) THEN
        READ(BUFFER,*) ISEED1,ISEED2
        WRITE(26,1750) ISEED1,ISEED2
 1750   FORMAT(3X,'Random number generator seeds = ',I10,', ',I10)
   75   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 75
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
      WRITE(26,1760) TIMEA
 1760 FORMAT(3X,'Computation time available = ',1P,E13.6,' sec')
C
      CALL TIMER(TSECIN)
      CPUT0=CPUTIM()
      TSECA=TIMEA+TSECIN
      TSECAD=TSECIN
      WRITE(26,1770)
 1770 FORMAT(//3X,72('-'))
C
C  ************  If 'RESUME' is active, read previously generated
C                counters...
C
      SHNA=0
      CPUTA=0.0D0
      IRETRN=0
C
      IF(IRESUM.EQ.1) THEN
        OPEN(9,FILE=PFILER)
        READ (9,*,ERR=1800,END=1800) SHNA,CPUTA
        READ (9,'(A120)') TITLE2
        IF(TITLE2.NE.TITLE) THEN
          WRITE(26,*)
     1      'The dump file is corrupted (the TITLE does not match).'
          STOP 'The dump file is corrupted (the TITLE does not match).'
        ENDIF
        READ (9,*) ISEED1,ISEED2
        READ (9,*) (SHIST(I),I=1,NSEB)
        READ (9,*) (PRIM(I),I=1,3),(PRIM2(I),I=1,3)
        READ (9,*) ((SEC(K,I),I=1,3),K=1,3),((SEC2(K,I),I=1,3),K=1,3)
        READ (9,*) (AVTL(I),I=1,3),(AVTL2(I),I=1,3)
        READ (9,*) (AVW(I),I=1,2),(AVW2(I),I=1,2)
        READ (9,*) (AVA(I),I=1,2),(AVA2(I),I=1,2)
        READ (9,*) (AVE(I),I=1,2),(AVE2(I),I=1,2)
        READ (9,*) (AVI(I),I=1,8),(AVI2(I),I=1,8)
        READ (9,*) (TDEBO(I),I=1,NBOD),(TDEBO2(I),I=1,NBOD)
        READ (9,*) ((PDTL(I,J),J=1,NBTL),I=1,3),
     1             ((PDTL2(I,J),J=1,NBTL),I=1,3)
        READ (9,*) (((PDA(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3),
     1             (((PDA2(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3)
        READ (9,*) (((PDE(I,J,K),K=1,NBE),J=1,2),I=1,3),
     1             (((PDE2(I,J,K),K=1,NBE),J=1,2),I=1,3)
        READ (9,*) ((DOSZ(I,J),J=1,NBZ),I=1,NLAY),
     1             ((DOSZ2(I,J),J=1,NBZ),I=1,NLAY)
        READ (9,*) ((CHRZ(I,J),J=1,NBZ),I=1,NLAY),
     1             ((CHRZ2(I,J),J=1,NBZ),I=1,NLAY)
C
        IF(NAED.GT.0) THEN
          READ (9,*) ((EDEP(I,J),J=1,NBE),I=1,NAED),
     1               ((EDEP2(I,J),J=1,NBE),I=1,NAED)
        ENDIF
C
        IF(NDOD.GT.0) THEN
          READ(9,*)
     1      (((DOSE(I,J,K),K=1,NRDOD(I)),J=1,NZDOD(I)),I=1,NDOD),
     1      (((DOSE2(I,J,K),K=1,NRDOD(I)),J=1,NZDOD(I)),I=1,NDOD)
          READ(9,*)
     1      (((CHAR(I,J,K),K=1,NRDOD(I)),J=1,NZDOD(I)),I=1,NDOD),
     1      (((CHAR2(I,J,K),K=1,NRDOD(I)),J=1,NZDOD(I)),I=1,NDOD)
        ENDIF
        CLOSE(9)
        GO TO 1802
 1800   CONTINUE
        WRITE(26,1801)
 1801   FORMAT(/3X,'WARNING: Could not resume from dump file...',/)
      ENDIF
 1802 CONTINUE
C
C  ************  Initialize constants.
C
      WGHT0=1.0D0   ! Primary particle weight.
      INTFOR=0      ! No interaction forcing as default.
C
      SHN=SHNA      ! Shower counter, including the dump file.
      DSHN=DSHN+SHNA
      N=MOD(SHN,2.0D9)+0.5D0
      IF(SHN.GE.DSHN) GO TO 105
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  ------------------------  Shower simulation starts here.
C
  101 CONTINUE  ! The simulation loop starts here.
C
C  **********  We set the initial state of the primary particle.
C
      SHN=SHN+1.0D0
      N=N+1
      IF(N.GT.2000000000) N=N-2000000000
      KPAR=KPARP
      WGHT=WGHT0
C  ----  Initial position ...
      IF(KSOURC.EQ.1) THEN
        Z=SZ0+(RAND(1.0D0)-0.5D0)*STHICK
        SR=SQRT(SRIN2+RAND(2.0D0)*SROI2)
        PHIR=RAND(3.0D0)*TWOPI
        X=SX0+SR*COS(PHIR)
        Y=SY0+SR*SIN(PHIR)
      ELSE
        RU=RAND(1.0D0)
        II=NSB
        DO I=1,NSB
          II=I
          IF(RU.LT.WSB(I)) GO TO 110
        ENDDO
  110   CONTINUE
        KL=KLS(II)
        KC=KCS(II)
        Z=ZG(KL)+RAND(1.0D0)*(ZG(KL+1)-ZG(KL))
        SR=SQRT(RG2(KL,KC)+RAND(2.0D0)*(RG2(KL,KC+1)-RG2(KL,KC)))
        PHIR=RAND(3.0D0)*TWOPI
        X=SX0+SR*COS(PHIR)
        Y=SY0+SR*SIN(PHIR)
      ENDIF
C  ----  Initial direction ...
      CALL GCONE(U,V,W)
C  ----  Check if the trajectory intersects the material system.
      CALL LOCATC
      IF(MAT.EQ.0) THEN
        CALL STEPC(1.0D30,DSEF,NCROSS)
        IF(MAT.EQ.0) THEN
          GO TO 105     ! The particle does not enter the system.
        ENDIF
      ENDIF
C  ----  Initial energy ...
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
      IF(E.LT.EABS(KPAR,MAT)) GO TO 105
C
C  ************  Initialization of primary particle counters.
C
      ILB(1)=1          ! Identifies primary particles.
      ILB(2)=0
      ILB(3)=0
      ILB(4)=0
      ILB(5)=0
      DO I=1,3
        DPRIM(I)=0.0D0
        DO K=1,3
          DSEC(K,I)=0.0D0
        ENDDO
      ENDDO
      DO I=1,8
        DAVI(I)=0.0D0   ! Numbers of collisions.
      ENDDO
      DO I=1,2
        DAVW(I)=0.0D0
        DAVA(I)=0.0D0
        DAVE(I)=0.0D0
      ENDDO
      DO KB=1,NBOD
        DEBO(KB)=0.0D0  ! Energies deposited in the various bodies KB.
      ENDDO
      TL=0.0D0          !  Accumulated track length.
C  ---------------------------------------------------------------------
C  ------------------------  Track simulation begins here.
C
C     EBAL=E*WGHT       ! Used to verify that energy is conserved.
C     IF(KPARP.EQ.3) EBAL=EBAL+TREV*WGHT
      CALL CLEANS       ! Cleans secondary stack.
  102 CONTINUE
      CALL START        ! Starts simulation in current medium.
C
C  ----  Free path length to the next interaction event.
C
  103 CONTINUE
      IF((IFORCE(IBODY,KPAR).EQ.1).AND.((WGHT.GE.WLOW(IBODY,KPAR)).AND.
     1  (WGHT.LE.WHIG(IBODY,KPAR)))) THEN
        CALL JUMPF(DSMAX(IBODY),DS)  ! Interaction forcing.
        INTFOR=1
      ELSE
        CALL JUMP(DSMAX(IBODY),DS)     ! Analogue simulation.
        INTFOR=0
      ENDIF
      CALL STEPC(DS,DSEF,NCROSS)  ! Determines step end position.
      TL=TL+DSEF  ! Accumulated track length (within media).
C  ----  Check whether the particle is outside the enclosure.
      IF(MAT.EQ.0) THEN
        IF(Z.GE.ZG(NLAY+1)) THEN
          IF(W.LT.0) STOP 'Transmitted with negative W?'
          IEXIT=1  ! Labels transmitted particles.
          GO TO 104
        ELSE IF(Z.LE.ZG(1)) THEN
          IF(W.GT.0) STOP 'Backscattered with positive W?'
          IEXIT=2  ! Labels backscattered particles.
          GO TO 104
        ENDIF
      ENDIF
C  ----  If the particle has crossed an interface, restart the track in
C        the new material.
      IF(NCROSS.GT.0) GO TO 102
C  ----  Simulate next event.
      IF(INTFOR.EQ.0) THEN
        CALL KNOCK(DE,ICOL)   ! Analogue simulation.
      ELSE
        CALL KNOCKF(DE,ICOL)  ! Interaction forcing is active.
      ENDIF
      IF(ILB(1).EQ.1) DAVI(ICOL)=DAVI(ICOL)+WGHT
C      WRITE(26,'(4I3,I9,3I3,1P,3E15.7,I3)')
C     1   KPAR,ILB(1),ILB(2),ILB(3),ILB(4),ILB(5),ICOL,MAT,E,Z,W,INTFOR
      DEP=DE*WGHT
C  ----  Energy is locally deposited in the material.
      DEBO(IBODY)=DEBO(IBODY)+DEP
C
C  ****  Depth distributions of dose and deposited charge.
C
      KZ=1.0D0+(Z-ZG(KLAY))/BSZ(KLAY)   ! Depth channel number.
      IF(KZ.GT.0.AND.KZ.LE.NBZ) THEN
        IF(N.NE.LDOSZ(KLAY,KZ)) THEN
          DOSZ(KLAY,KZ)=DOSZ(KLAY,KZ)+DOSZP(KLAY,KZ)
          DOSZ2(KLAY,KZ)=DOSZ2(KLAY,KZ)+DOSZP(KLAY,KZ)**2
          DOSZP(KLAY,KZ)=DEP*RHOI(MAT)
          LDOSZ(KLAY,KZ)=N
        ELSE
          DOSZP(KLAY,KZ)=DOSZP(KLAY,KZ)+DEP*RHOI(MAT)
        ENDIF
C
        IF(E.LT.EABS(KPAR,MAT).AND.KPAR.NE.2) THEN
          IF(N.NE.LCHRZ(KLAY,KZ)) THEN
            CHRZ(KLAY,KZ)=CHRZ(KLAY,KZ)+CHRZP(KLAY,KZ)
            CHRZ2(KLAY,KZ)=CHRZ2(KLAY,KZ)+CHRZP(KLAY,KZ)**2
            CHRZP(KLAY,KZ)=(KPAR-2)*WGHT
            LCHRZ(KLAY,KZ)=N
          ELSE
            CHRZP(KLAY,KZ)=CHRZP(KLAY,KZ)+(KPAR-2)*WGHT
          ENDIF
        ENDIF
      ENDIF
C
C  ****  3D Dose and deposited charge distributions.
C
      IF(IDOD(IBODY).NE.0) THEN
        JI=IDOD(IBODY)
        JZ=1.0D0+(Z-ZG(KLAY))/BZDOD(JI)   ! Depth channel number.
        JR=1.0D0+(SQRT((X-XG(KLAY))**2+(Y-YG(KLAY))**2)
     1    -RG(KLAY,KCYL))/BRDOD(JI)  ! Radial channel number.
        IF(N.NE.LDOSE(JI,JZ,JR)) THEN
          DOSE(JI,JZ,JR)=DOSE(JI,JZ,JR)+DOSEP(JI,JZ,JR)
          DOSE2(JI,JZ,JR)=DOSE2(JI,JZ,JR)+DOSEP(JI,JZ,JR)**2
          DOSEP(JI,JZ,JR)=DEP*RHOI(MAT)
          LDOSE(JI,JZ,JR)=N
        ELSE
          DOSEP(JI,JZ,JR)=DOSEP(JI,JZ,JR)+DEP*RHOI(MAT)
        ENDIF
C
        IF(E.LT.EABS(KPAR,MAT).AND.KPAR.NE.2) THEN
          IF(N.NE.LCHAR(JI,JZ,JR)) THEN
            CHAR(JI,JZ,JR)=CHAR(JI,JZ,JR)+CHARP(JI,JZ,JR)
            CHAR2(JI,JZ,JR)=CHAR2(JI,JZ,JR)+CHARP(JI,JZ,JR)**2
            CHARP(JI,JZ,JR)=(KPAR-2)*WGHT
            LCHAR(JI,JZ,JR)=N
          ELSE
            CHARP(JI,JZ,JR)=CHARP(JI,JZ,JR)+(KPAR-2)*WGHT
          ENDIF
        ENDIF
      ENDIF
C  ----  Check if the particle has been absorbed .
      IF(E.GT.EABS(KPAR,MAT)) GO TO 103
      IEXIT=3  ! Labels absorbed particles.
C  ------------------------  The simulation of the track ends here.
C  ---------------------------------------------------------------------
  104 CONTINUE
C     EBAL=EBAL-E*WGHT
C     IF(KPAR.EQ.3.AND.IEXIT.NE.3) EBAL=EBAL-TREV*WGHT
C
C  ************  Counters.
C
C  ****  Primary particles.
C
      THETA=DACOS(W)
      IF(ILB(1).EQ.1) THEN
        DPRIM(IEXIT)=DPRIM(IEXIT)+WGHT0
        AVTL(IEXIT)=AVTL(IEXIT)+TL*WGHT0
        AVTL2(IEXIT)=AVTL2(IEXIT)+(TL*WGHT0)**2
        IF(IEXIT.LT.3) THEN
          DAVW(IEXIT)=DAVW(IEXIT)+W*WGHT0
          DAVA(IEXIT)=DAVA(IEXIT)+THETA*WGHT0
          DAVE(IEXIT)=DAVE(IEXIT)+E*WGHT0
        ENDIF
C  ****  Track-length distribution.
        KT=1.0D0+(TL-TLMIN)/BSTL
        IF(KT.GT.0.AND.KT.LE.NBTL) THEN
          PDTL(IEXIT,KT)=PDTL(IEXIT,KT)+WGHT0
          PDTL2(IEXIT,KT)=PDTL2(IEXIT,KT)+WGHT0**2
        ENDIF
      ELSE
        DSEC(KPAR,IEXIT)=DSEC(KPAR,IEXIT)+WGHT
      ENDIF
C
C  ****  Other counters.
C
      IF(IEXIT.LT.3) THEN
C  ****  Angular distribution of emerging particles.
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
      ENDIF
C
C  ************  Any secondary left?
C
      CALL SECPAR(LEFT)
      IF(LEFT.GT.0) THEN
        INOUT=1
        KLAY=ILAY(IBODY)
        KCYL=ICYL(IBODY)
C  ****  Subtract E and charge from the tallied distributions to avoid
C        double-counting.
        KZ=1.0D0+(Z-ZG(KLAY))/BSZ(KLAY)   ! Depth channel number.
        IF(N.NE.LDOSZ(KLAY,KZ)) THEN
          DOSZ(KLAY,KZ)=DOSZ(KLAY,KZ)+DOSZP(KLAY,KZ)
          DOSZ2(KLAY,KZ)=DOSZ2(KLAY,KZ)+DOSZP(KLAY,KZ)**2
          DOSZP(KLAY,KZ)=-E*WGHT*RHOI(MAT)
          LDOSZ(KLAY,KZ)=N
        ELSE
          DOSZP(KLAY,KZ)=DOSZP(KLAY,KZ)-E*WGHT*RHOI(MAT)
        ENDIF
C
        IF(N.NE.LCHRZ(KLAY,KZ)) THEN
          CHRZ(KLAY,KZ)=CHRZ(KLAY,KZ)+CHRZP(KLAY,KZ)
          CHRZ2(KLAY,KZ)=CHRZ2(KLAY,KZ)+CHRZP(KLAY,KZ)**2
          CHRZP(KLAY,KZ)=(2-KPAR)*WGHT
          LCHRZ(KLAY,KZ)=N
        ELSE
          CHRZP(KLAY,KZ)=CHRZP(KLAY,KZ)+(2-KPAR)*WGHT
        ENDIF
C
        DEBO(IBODY)=DEBO(IBODY)-E*WGHT
C
        IF(IDOD(IBODY).NE.0) THEN
          JI=IDOD(IBODY)
          JZ=1.0D0+(Z-ZG(KLAY))/BZDOD(JI)   ! Depth channel number.
          JR=1.0D0+(SQRT((X-XG(KLAY))**2+(Y-YG(KLAY))**2)
     1      -RG(KLAY,KCYL))/BRDOD(JI)  ! Radial channel number.
          IF(N.NE.LDOSE(JI,JZ,JR)) THEN
            DOSE(JI,JZ,JR)=DOSE(JI,JZ,JR)+DOSEP(JI,JZ,JR)
            DOSE2(JI,JZ,JR)=DOSE2(JI,JZ,JR)+DOSEP(JI,JZ,JR)**2
            DOSEP(JI,JZ,JR)=-E*WGHT*RHOI(MAT)
            LDOSE(JI,JZ,JR)=N
          ELSE
            DOSEP(JI,JZ,JR)=DOSEP(JI,JZ,JR)-E*WGHT*RHOI(MAT)
          ENDIF
C
          IF(N.NE.LCHAR(JI,JZ,JR)) THEN
            CHAR(JI,JZ,JR)=CHAR(JI,JZ,JR)+CHARP(JI,JZ,JR)
            CHAR2(JI,JZ,JR)=CHAR2(JI,JZ,JR)+CHARP(JI,JZ,JR)**2
            CHARP(JI,JZ,JR)=(2-KPAR)*WGHT
            LCHAR(JI,JZ,JR)=N
          ELSE
            CHARP(JI,JZ,JR)=CHARP(JI,JZ,JR)+(2-KPAR)*WGHT
          ENDIF
        ENDIF
        GO TO 102
      ENDIF
C
C  ------------------------  The simulation of the shower ends here.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO I=1,3
        PRIM(I)=PRIM(I)+DPRIM(I)
        PRIM2(I)=PRIM2(I)+DPRIM(I)**2
        DO K=1,3
          SEC(K,I)=SEC(K,I)+DSEC(K,I)
          SEC2(K,I)=SEC2(K,I)+DSEC(K,I)**2
        ENDDO
      ENDDO
      DO I=1,2
        AVW(I)=AVW(I)+DAVW(I)
        AVW2(I)=AVW2(I)+DAVW(I)**2
        AVA(I)=AVA(I)+DAVA(I)
        AVA2(I)=AVA2(I)+DAVA(I)**2
        AVE(I)=AVE(I)+DAVE(I)
        AVE2(I)=AVE2(I)+DAVE(I)**2
      ENDDO
      DO I=1,8
        AVI(I)=AVI(I)+DAVI(I)
        AVI2(I)=AVI2(I)+DAVI(I)**2
      ENDDO
C
      DO KM=1,NMAT
        AESUM(KM)=0.0D0
      ENDDO
      DO KB=1,NBOD
C       EBAL=EBAL-DEBO(KB)
        TDEBO(KB)=TDEBO(KB)+DEBO(KB)
        TDEBO2(KB)=TDEBO2(KB)+DEBO(KB)**2
        MATKB=MATER(ILAY(KB),ICYL(KB))
        IF(MATKB.GT.0) THEN
          IF(IAED(MATKB).NE.0) AESUM(MATKB)=AESUM(MATKB)+DEBO(KB)
        ENDIF
      ENDDO
      DO KM=1,NMAT
        IF(IAED(KM).NE.0) THEN
          IF(AESUM(KM).GT.1.0D-5) THEN
            KE=1.0D0+AESUM(KM)/BSDE
            IF(KE.GT.0.AND.KE.LE.NBE) THEN
              EDEP(IAED(KM),KE)=EDEP(IAED(KM),KE)+WGHT0
              EDEP2(IAED(KM),KE)=EDEP2(IAED(KM),KE)+WGHT0**2
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C
C     IF(ABS(EBAL).GT.0.01D0) THEN
C       WRITE(26,'(//3X,''Energy is not conserved:'')')
C       WRITE(26,'(5X,''EBAL ='',1P,E13.6,'' eV,  N ='',I9)') EBAL,N
C       WRITE(26,'(3X,''Correct the error, or switch off the test.'')')
C       STOP 'Verify energy conservation.'
C     ENDIF
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
C  ****  Transfer contents of partial counters to global counters.
C
      DO KPAR=1,3
        DO K=1,NBTH
          DO I=1,NBPH
            PDA(KPAR,K,I)=PDA(KPAR,K,I)+PDAP(KPAR,K,I)
            PDA2(KPAR,K,I)=PDA2(KPAR,K,I)+PDAP(KPAR,K,I)**2
            PDAP(KPAR,K,I)=0.0D0
            LPDA(KPAR,K,I)=0
          ENDDO
        ENDDO
        DO K=1,NBE
          PDE(KPAR,1,K)=PDE(KPAR,1,K)+PDEP(KPAR,1,K)
          PDE2(KPAR,1,K)=PDE2(KPAR,1,K)+PDEP(KPAR,1,K)**2
          PDEP(KPAR,1,K)=0.0D0
          LPDE(KPAR,1,K)=0
          PDE(KPAR,2,K)=PDE(KPAR,2,K)+PDEP(KPAR,2,K)
          PDE2(KPAR,2,K)=PDE2(KPAR,2,K)+PDEP(KPAR,2,K)**2
          PDEP(KPAR,2,K)=0.0D0
          LPDE(KPAR,2,K)=0
        ENDDO
      ENDDO
C
      DO L=1,NLAY
        DO J=1,NBZ
          DOSZ(L,J)=DOSZ(L,J)+DOSZP(L,J)
          DOSZ2(L,J)=DOSZ2(L,J)+DOSZP(L,J)**2
          DOSZP(L,J)=0.0D0
          LDOSZ(L,J)=0
          CHRZ(L,J)=CHRZ(L,J)+CHRZP(L,J)
          CHRZ2(L,J)=CHRZ2(L,J)+CHRZP(L,J)**2
          CHRZP(L,J)=0.0D0
          LCHRZ(L,J)=0
        ENDDO
      ENDDO
C
      DO I=1,NDOT
        DO J=1,NZDOD(I)
          DO K=1,NRDOD(I)
            DOSE(I,J,K)=DOSE(I,J,K)+DOSEP(I,J,K)
            DOSE2(I,J,K)=DOSE2(I,J,K)+DOSEP(I,J,K)**2
            DOSEP(I,J,K)=0.0D0
            LDOSE(I,J,K)=0.0D0
            CHAR(I,J,K)=CHAR(I,J,K)+CHARP(I,J,K)
            CHAR2(I,J,K)=CHAR2(I,J,K)+CHARP(I,J,K)**2
            CHARP(I,J,K)=0.0D0
            LCHAR(I,J,K)=0
          ENDDO
        ENDDO
      ENDDO
C
C  ************  If 'DUMPTO' is active, write counters in a dump file.
C
      IF(IDUMP.EQ.1) THEN
        OPEN(9,FILE=PFILED)
        WRITE(9,*) SHN,TSIM
        WRITE(9,'(A120)') TITLE
        WRITE(9,*) ISEED1,ISEED2
        WRITE(9,*) (SHIST(I),I=1,NSEB)
        WRITE(9,*) (PRIM(I),I=1,3),(PRIM2(I),I=1,3)
        WRITE(9,*) ((SEC(K,I),I=1,3),K=1,3),((SEC2(K,I),I=1,3),K=1,3)
        WRITE(9,*) (AVTL(I),I=1,3),(AVTL2(I),I=1,3)
        WRITE(9,*) (AVW(I),I=1,2),(AVW2(I),I=1,2)
        WRITE(9,*) (AVA(I),I=1,2),(AVA2(I),I=1,2)
        WRITE(9,*) (AVE(I),I=1,2),(AVE2(I),I=1,2)
        WRITE(9,*) (AVI(I),I=1,8),(AVI2(I),I=1,8)
        WRITE(9,*) (TDEBO(I),I=1,NBOD),(TDEBO2(I),I=1,NBOD)
        WRITE(9,*) ((PDTL(I,J),J=1,NBTL),I=1,3),
     1             ((PDTL2(I,J),J=1,NBTL),I=1,3)
        WRITE(9,*) (((PDA(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3),
     1             (((PDA2(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3)
        WRITE(9,*) (((PDE(I,J,K),K=1,NBE),J=1,2),I=1,3),
     1             (((PDE2(I,J,K),K=1,NBE),J=1,2),I=1,3)
        WRITE(9,*) ((DOSZ(I,J),J=1,NBZ),I=1,NLAY),
     1             ((DOSZ2(I,J),J=1,NBZ),I=1,NLAY)
        WRITE(9,*) ((CHRZ(I,J),J=1,NBZ),I=1,NLAY),
     1             ((CHRZ2(I,J),J=1,NBZ),I=1,NLAY)
C
        IF(NAED.GT.0) THEN
          WRITE(9,*) ((EDEP(I,J),J=1,NBE),I=1,NAED),
     1               ((EDEP2(I,J),J=1,NBE),I=1,NAED)
        ENDIF
C
        IF(NDOD.GT.0) THEN
          WRITE(9,*)
     1      (((DOSE(I,J,K),K=1,NRDOD(I)),J=1,NZDOD(I)),I=1,NDOD),
     1      (((DOSE2(I,J,K),K=1,NRDOD(I)),J=1,NZDOD(I)),I=1,NDOD)
          WRITE(9,*)
     1      (((CHAR(I,J,K),K=1,NRDOD(I)),J=1,NZDOD(I)),I=1,NDOD),
     1      (((CHAR2(I,J,K),K=1,NRDOD(I)),J=1,NZDOD(I)),I=1,NDOD)
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
 3000 FORMAT(///3X,34('*')/3X,'**   Program PENCYL. Results.   **',
     1  /3X,34('*'))
C
      WRITE(26,3001) TSIM
 3001 FORMAT(/3X,'Simulation time ......................... ',
     1  1P,E13.6,' sec')
      TAVS=TOTN/TSIM
      WRITE(26,3002) TAVS
 3002 FORMAT(3X,'Simulation speed ........................ ',
     1  1P,E13.6,' showers/sec')
      WRITE(26,3003) TOTN
 3003 FORMAT(//3X,'Simulated primary particles ............. ',
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
      FNT=1.0D0/TOTN
      FT=(PRIM(1)+SEC(KPARP,1))*FNT
      ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(1)-PRIM(1)**2*FNT))
      ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,1)-SEC(KPARP,1)**2*FNT))
      ERR=ERR1+ERR2
      WRITE(26,3007) FT,ERR
 3007 FORMAT(/3X,'Fractional transmission ............ ',
     1  1P,E13.6,' +-',E8.1)
      FB=(PRIM(2)+SEC(KPARP,2))*FNT
      ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(2)-PRIM(2)**2*FNT))
      ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,2)-SEC(KPARP,2)**2*FNT))
      ERR=ERR1+ERR2
      WRITE(26,3008) FB,ERR
 3008 FORMAT(3X,'Fractional backscattering .......... ',
     1  1P,E13.6,' +-',E8.1)
      FA=PRIM(3)*FNT
      ERR=3.0D0*FNT*SQRT(ABS(PRIM2(3)-PRIM(3)**2*FNT))
      WRITE(26,3009) FA,ERR
 3009 FORMAT(3X,'Fractional absorption .............. ',
     1  1P,E13.6,' +-',E8.1)
C
      DO K=1,3
        DO I=1,3
          WSEC2(K,I)=3.0D0*FNT*SQRT(ABS(SEC2(K,I)-SEC(K,I)**2*FNT))
          WSEC(K,I)=SEC(K,I)*FNT
        ENDDO
      ENDDO
      WRITE(26,3010)
     1  WSEC(1,1),WSEC(2,1),WSEC(3,1),WSEC2(1,1),WSEC2(2,1),WSEC2(3,1),
     1  WSEC(1,2),WSEC(2,2),WSEC(3,2),WSEC2(1,2),WSEC2(2,2),WSEC2(3,2),
     1  WSEC(1,3),WSEC(2,3),WSEC(3,3),WSEC2(1,3),WSEC2(2,3),WSEC2(3,3)
 3010 FORMAT(/3X,'Secondary-particle generation probabilities:',
     1 /19X,46('-'),
     1 /19X,'|  electrons   |   photons    |  positrons   |',1P,
     1 /3X,62('-')/3X,'|  transmitted  |',3(E13.6,1X,'|'),
     1 /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1 /3X,62('-')/3X,'| backscattered |',3(E13.6,1X,'|'),
     1 /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1 /3X,62('-')/3X,'|   absorbed    |',3(E13.6,1X,'|'),
     1 /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1 /3X,62('-'))
C
C  ************  Average values for primary particles.
C
      DO KC=1,8
        IF(AVI(KC).GT.1.0D-16) THEN
          WAVI2(KC)=3.0D0*FNT*SQRT(ABS(AVI2(KC)-AVI(KC)**2*FNT))
          WAVI(KC)=AVI(KC)*FNT
        ELSE
          WAVI2(KC)=0.0D0
          WAVI(KC)=0.0D0
        ENDIF
      ENDDO
      DO I=1,3
        DF=1.0D0/MAX(PRIM(I),1.0D0)
        WAVTL2(I)=3.0D0*DF*SQRT(ABS(AVTL2(I)-AVTL(I)**2*DF))
        WAVTL(I)=AVTL(I)*DF
        IF(I.LT.3) THEN
          WAVE2(I)=3.0D0*DF*SQRT(ABS(AVE2(I)-AVE(I)**2*DF))
          WAVE(I)=AVE(I)*DF
          WAVW2(I)=3.0D0*DF*SQRT(ABS(AVW2(I)-AVW(I)**2*DF))
          WAVW(I)=AVW(I)*DF
          WAVA2(I)=3.0D0*DF*RA2DE*SQRT(ABS(AVA2(I)-AVA(I)**2*DF))
          WAVA(I)=AVA(I)*RA2DE*DF
        ENDIF
      ENDDO
C
      WRITE(26,3011)
 3011 FORMAT(//3X,'Mean number of events per primary track:')
      IF(KPARP.EQ.2) THEN
        WRITE(26,3013) WAVI(1),WAVI2(1)
 3013   FORMAT(6X,'Coherent (Rayleigh) ............. ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3014) WAVI(2),WAVI2(2)
 3014   FORMAT(6X,'Incoherent (Compton) ............ ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3015) WAVI(3),WAVI2(3)
 3015   FORMAT(6X,'Photoelectric absorption ........ ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3016) WAVI(4),WAVI2(4)
 3016   FORMAT(6X,'Pair production ................. ',1P,E13.6,
     1    ' +-',E8.1)
      ELSE
        WRITE(26,3017) WAVI(1),WAVI2(1)
 3017   FORMAT(6X,'Hinges (soft events) ............ ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3018) WAVI(2),WAVI2(2)
 3018   FORMAT(6X,'Hard elastic collisions ......... ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3019) WAVI(3),WAVI2(3)
 3019   FORMAT(6X,'Hard inelastic collisions ....... ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3020) WAVI(4),WAVI2(4)
 3020   FORMAT(6X,'Hard bremsstrahlung emissions ... ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3021) WAVI(5),WAVI2(5)
 3021   FORMAT(6X,'Inner-shell ionizations ......... ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3022) WAVI(7),WAVI2(7)
 3022   FORMAT(6X,'Delta interactions .............. ',1P,E13.6,
     1    ' +-',E8.1)
      ENDIF
C
      WRITE(26,3025)
 3025 FORMAT(/3X,'Average final energy:')
      WRITE(26,3026) WAVE(1),WAVE2(1)
 3026 FORMAT(6X,'Transmitted particles ........... ',1P,E13.6,
     1  ' +-',E8.1,' eV')
      WRITE(26,3027) WAVE(2),WAVE2(2)
 3027 FORMAT(6X,'Backscattered particles ......... ',1P,E13.6,
     1  ' +-',E8.1,' eV')
C
      WRITE(26,3028)
 3028 FORMAT(/3X,'Average track length:')
      WRITE(26,3029) WAVTL(1),WAVTL2(1)
 3029 FORMAT(6X,'Transmitted particles ........... ',1P,E13.6,
     1  ' +-',E8.1,' cm')
      WRITE(26,3030) WAVTL(2),WAVTL2(2)
 3030 FORMAT(6X,'Backscattered particles ......... ',1P,E13.6,
     1  ' +-',E8.1,' cm')
      WRITE(26,3031) WAVTL(3),WAVTL2(3)
 3031 FORMAT(6X,'Absorbed particles .............. ',1P,E13.6,
     1  ' +-',E8.1,' cm')
C
      WRITE(26,3032)
 3032 FORMAT(/3X,'Mean value of the polar cosine of the exit',
     1  ' direction:')
      WRITE(26,3033) WAVW(1),WAVW2(1)
 3033 FORMAT(6X,'Transmitted particles ........... ',1P,E13.6,
     1  ' +-',E8.1)
      WRITE(26,3034) WAVW(2),WAVW2(2)
 3034 FORMAT(6X,'Backscattered particles ......... ',1P,E13.6,
     1  ' +-',E8.1)
C
      WRITE(26,3035)
 3035 FORMAT(/3X,'Mean value of the polar angle of the exit dir',
     1  'ection:')
      WRITE(26,3036) WAVA(1),WAVA2(1)
 3036 FORMAT(6X,'Transmitted particles ........... ',1P,E13.6,
     1  ' +-',E8.1,' deg')
      WRITE(26,3037) WAVA(2),WAVA2(2)
 3037 FORMAT(6X,'Backscattered particles ......... ',1P,E13.6,
     1  ' +-',E8.1,' deg')
C
C  ****  Deposited energies.
C
      DF=1.0D0/TOTN
      WRITE(26,3038)
 3038 FORMAT(/3X,'Average deposited energies (per primary shower):')
      DO KB=1,NBOD
        QER=3.0D0*DF*SQRT(ABS(TDEBO2(KB)-TDEBO(KB)**2*DF))
        QAV=TDEBO(KB)*DF
        IF(QER.GT.1.0D-10*ABS(QAV)) THEN
          EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
        ELSE
          EFFIC=0.0D0
        ENDIF
        WRITE(26,3039) KB,QAV,QER,EFFIC
      ENDDO
 3039 FORMAT(6X,'Body ',I4, ' ...... ',1P,E13.6,' +-',E8.1,' eV',4X,
     1  '(effic. =',E9.2,')')
      WRITE(26,3040) ISEED1,ISEED2
 3040 FORMAT(/3X,'Last random seeds = ',I10,' , ',I10)
C
C  ------------------------  Print tallied distributions.
C
      IF(ISPEC.EQ.1) THEN
C
C  ************  Energy spectrum of the source.
C
        DF=1.0D0/TOTN
        OPEN(9,FILE='psource.dat')
        WRITE(9,9000)
 9000 FORMAT(
     1  1X,'#  Results from PENCYL. ',
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
C  ************  Depth-dose.
C
      OPEN(9,FILE='pcddose.dat')
      WRITE(9,9010)
 9010 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Depth dose distribution.',
     1 /1X,'#  1st column: Z (cm). 4th column: corresponding layer.',
     1 /1X,'#  2nd and 3rd columns: dose and STU (eV/',
     1      '((g/cm**2)*particle)).',/)
      DF=1.0D0/TOTN
      DO KL=1,NLAY
        DO K=1,NBZ
          XX=ZG(KL)+(K-0.5D0)*BSZ(KL)
          YERR=3.0D0*SQRT(ABS(DOSZ2(KL,K)-DOSZ(KL,K)**2*DF))
          YAV=DOSZ(KL,K)*DF/BSZ(KL)
          YERR=YERR*DF/BSZ(KL)
          WRITE(9,'(1X,1P,3E14.6,I5)')
     1       XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35),KL
        ENDDO
      ENDDO
      CLOSE(9)
C
C  ************  Depth distribution of deposited charge.
C
      OPEN(9,FILE='pcdchar.dat')
      WRITE(9,9020)
 9020 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Depth distribution of deposited charge.',
     1 /1X,'#  1st column: Z (cm). 4th column: corresponding layer.',
     1 /1X,'#  2nd and 3rd columns:',
     1      ' charge density and STU (e/(cm*particle)).',/)
      DF=1.0D0/TOTN
      DO KL=1,NLAY
        DO K=1,NBZ
          XX=ZG(KL)+(K-0.5D0)*BSZ(KL)
          YERR=3.0D0*SQRT(ABS(CHRZ2(KL,K)-CHRZ(KL,K)**2*DF))
          YAV=CHRZ(KL,K)*DF/BSZ(KL)
          YERR=YERR*DF/BSZ(KL)
          WRITE(9,'(1X,1P,3E14.6,I5)') XX,YAV,YERR,KL
        ENDDO
      ENDDO
      CLOSE(9)
C
C  ************  Track length distributions of primary particles.
C
C  ****  Transmitted particles.
      OPEN(9,FILE='pctltr.dat')
      WRITE(9,9030)
 9030 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Track length distribution of transmitted particles.',
     1 /1X,'#  1st column: track length (cm).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(cm*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBTL
        XX=TLMIN+(K-0.5D0)*BSTL
        YERR=3.0D0*SQRT(ABS(PDTL2(1,K)-PDTL(1,K)**2*DF))
        YAV=PDTL(1,K)*DF/BSTL
        YERR=YERR*DF/BSTL
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Backscattered particles.
      OPEN(9,FILE='pctlbk.dat')
      WRITE(9,9040)
 9040 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Track length distribution of backscattered particles.',
     1 /1X,'#  1st column: track length (cm).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(cm*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBTL
        XX=TLMIN+(K-0.5D0)*BSTL
        YERR=3.0D0*SQRT(ABS(PDTL2(2,K)-PDTL(2,K)**2*DF))
        YAV=PDTL(2,K)*DF/BSTL
        YERR=YERR*DF/BSTL
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Absorbed particles.
      OPEN(9,FILE='pctlab.dat')
      WRITE(9,9050)
 9050 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Track length distribution of absorbed particles.',
     1 /1X,'#  1st column: track length (cm).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(cm*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBTL
        XX=TLMIN+(K-0.5D0)*BSTL
        YERR=3.0D0*SQRT(ABS(PDTL2(3,K)-PDTL(3,K)**2*DF))
        YAV=PDTL(3,K)*DF/BSTL
        YERR=YERR*DF/BSTL
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
C  ************  Angular distribution of emerging particles.
C
      OPEN(9,FILE='pcanel.dat')
      WRITE(9,9060)
 9060 FORMAT(
     1  1X,'#  Results from PENCYL. ',
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
      OPEN(9,FILE='pcanga.dat')
      WRITE(9,9070)
 9070 FORMAT(
     1  1X,'#  Results from PENCYL. ',
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
      OPEN(9,FILE='pcanpo.dat')
      WRITE(9,9080)
 9080 FORMAT(
     1  1X,'#  Results from PENCYL. ',
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
C  ************  Energy distributions of emerging particles.
C
C  ****  Transmitted electrons.
      OPEN(9,FILE='pceneltr.dat')
      WRITE(9,9090)
 9090 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Energy distribution of transmitted electrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
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
      OPEN(9,FILE='pcenelbk.dat')
      WRITE(9,9100)
 9100 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Energy distribution of backscattered electrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
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
      OPEN(9,FILE='pcengatr.dat')
      WRITE(9,9110)
 9110 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Energy distribution of transmitted photons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
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
      OPEN(9,FILE='pcengabk.dat')
      WRITE(9,9120)
 9120 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Energy distribution of backscattered photons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
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
      OPEN(9,FILE='pcenpotr.dat')
      WRITE(9,9130)
 9130 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Energy distribution of transmitted positrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
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
      OPEN(9,FILE='pcenpobk.dat')
      WRITE(9,9140)
 9140 FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Energy distribution of backscattered positrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
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
C  ************  Distribution of deposited energy in selected bodies.
C
      IF(NAED.GT.0) THEN
        DO IDE=1,NAED
          WRITE(OFILE,'(A7,I1,A4)') 'pcedepm',IDE,'.dat'
          OPEN(9,FILE=OFILE)
          WRITE(9,9150) MATAE(IDE)
 9150     FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  Distribution of deposited energy in material ',I3,
     1 /1X,'#  WARNING: Strongly biased if variance reduction is used!',
     1 /1X,'#  1st column: deposited energy (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
          DF=1.0D0/TOTN
          DO J=1,NBE
            XX=(J-0.5D0)*BSDE
            YERR=3.0D0*SQRT(ABS(EDEP2(IDE,J)-EDEP(IDE,J)**2*DF))
            YAV=EDEP(IDE,J)*DF/BSDE
            YERR=YERR*DF/BSDE
            WRITE(9,'(1X,1P,3E14.6)')
     1        XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
          ENDDO
          CLOSE(9)
        ENDDO
      ENDIF
C
C  ************  3D dose and deposited charge distributions.
C
      IF(NDOD.GT.0) THEN
        DO KI=1,NDOD
          WRITE(OFILE,'(A7,I1,A4)') 'pcdosch',KI,'.dat'
          OPEN(9,FILE=OFILE)
          KL=ILAY(KBDOD(KI))
          KC=ICYL(KBDOD(KI))
          WRITE(9,9160) KL,KC
 9160     FORMAT(
     1  1X,'#  Results from PENCYL. ',
     1 /1X,'#  3D dose and charge distribs. in body KL=',I3,', KC=',I3,
     1 /1X,'#  Columns: 1=Z (cm),  2=radius (=sqrt(X**2+Y**2), in cm) ',
     1 /1X,'#           3=dose (eV/g), 4=uncertainty (eV/g).',
     1 /1X,'#           5=charge (e/cm**3), 6=uncertainty (e/cm**3).',/)
          DF=1.0D0/TOTN
          FACT=PI*BRDOD(KI)*BZDOD(KI)
          DO KR=1,NRDOD(KI)
            R=RG(KL,KC)+(KR-0.5D0)*BRDOD(KI)
            BINVOL=((2*KR-1)*BRDOD(KI)+2.0D0*RG(KL,KC))*FACT
            DO KZ=1,NZDOD(KI)
              Z=ZG(KL)+(KZ-0.5D0)*BZDOD(KI)
              DERR=3.0D0*SQRT(ABS(DOSE2(KI,KZ,KR)-DOSE(KI,KZ,KR)**2*DF))
              DAV=DOSE(KI,KZ,KR)*DF/BINVOL
              DERR=DERR*DF/BINVOL
              CERR=3.0D0*SQRT(ABS(CHAR2(KI,KZ,KR)-CHAR(KI,KZ,KR)**2*DF))
              CAV=CHAR(KI,KZ,KR)*DF/BINVOL
              CERR=CERR*DF/BINVOL
              WRITE(9,'(1X,1P,6E14.6)') Z,R,MAX(DAV,1.0D-35),
     1          MAX(DERR,1.0D-35),CAV,MAX(CERR,1.0D-35)
            ENDDO
            WRITE(9,*) '  '
          ENDDO
          CLOSE(9)
        ENDDO
      ENDIF
C
C  ****  Continue if the DUMPTO option is active and IRETRN=1.
C
      IF(IRETRN.EQ.1) THEN
        WRITE(26,'(/3X,72(''-''))')
        WRITE(6,3045) SHN
 3045   FORMAT(2X,'Number of simulated showers =',1P,E14.7)
        GO TO 101
      ENDIF
C
      WRITE(26,3050)
 3050 FORMAT(//3X,'***  END  ***')
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
C  Parameters are initialized by calling subroutine GCONE0.
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
C  ****  Ensure normalization.
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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C            CYLINDRICAL GEOMETRY PACKAGE   (version 2003)             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  The material system consists of one or several layers of given thick-
C  nesses. Each layer contains a number of concentric homogeneous rings
C  of given compositions and radii (and thickness equal to that of the
C  layer). The layers are perpendicular to the Z-axis and the centre of
C  the rings in each layer is specified by giving its X and Y coordi-
C  nates.
C
C  The geometry specification is read from the input file (UNIT=IRD).
C  The definition list consists of a number of lines with the following
C  structure.
C
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C  GSTART >>>> Beginning of the geometry definition list. Title.
C  LAYER  ZLOW,ZHIG                               [Z_lower and Z_higher]
C  CENTRE XCEN,YCEN            [X and Y coordinates of the layer centre]
C  CYLIND M,RIN,ROUT                     [Material, R_inner and R_outer]
C  GEND   <<<< End of the geometry definition list.
C
C  Each line starts with a 6-character keyword (written in columns 1-6)
C  The only allowed keywords are 'GSTART', 'LAYER_', 'CENTRE', 'CYLIND'
C  and 'GEND__' (notice the blanks). The definition list must begin with
C  the 'GSTART' line and terminate with the 'GEND__' line. The second
C  line must be a 'LAYER_' line, which initiates the definition of the
C  layer structure. Each 'LAYER_' line is followed by a 'CENTRE' line
C  (optional) and by one or several 'CYLIND' lines, which define the
C  various concentric rings in the layer; empty layers are disregarded.
C
C  Layers must be defined in increasing order of heights, from bottom to
C  top of the structure. If the 'CENTRE' line is not entered, cylinders
C  are assumed to be centered on the Z-axis (XCEN=YCEN=0.0). Cylinders
C  have to be defined in increasing radial order, from the centre to the
C  periphery. The two lengths in each 'LAYER_' and 'CYLIND' line must
C  be entered in increasing order. All numerical data are read in free
C  format.
C
C  *********************************************************************
C                       SUBROUTINE GEOINC
C  *********************************************************************
      SUBROUTINE GEOINC(NMATG,IRD,IWR)
C
C  This subroutine reads the cylindrical geometry data and initializes
C  the geometry routines for Monte Carlo simulation of particle trans-
C  port.
C
C  Input arguments:
C     IRD ....... input file unit (opened in the main program).
C     IWR ....... output file unit (opened in the main program).
C  Outut arguments:
C     NMATG ..... number of different materials in the geometry.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*6 GSTART,LAYER,CENTRE,CYLIND,GEND,KWCOMM,KWORD
      CHARACTER*120 BUFFER
      PARAMETER (GSTART='GSTART',LAYER='LAYER ',CENTRE='CENTRE',
     1   CYLIND='CYLIND',GEND='GEND  ',KWCOMM='      ')
      PARAMETER (EPS=1.0D-9,TEPS=EPS+EPS,FEPS=5.0D0*EPS)
C  ****  Cylindrical geometry.
      PARAMETER (NLAM=100,NCYM=50,NBDM=NLAM*NCYM)
      COMMON/CYLGEO/XG(NLAM),YG(NLAM),ZG(NLAM),RG(NLAM,NCYM),
     1  RG2(NLAM,NCYM),RMAX,RMAX2,IBOD(NLAM,NCYM),MATER(NLAM,NCYM),
     1  ILAY(NBDM),ICYL(NBDM),NLAY,NCYL(NLAM),NBOD
      DIMENSION MATT(NLAM)
C
      RMAX=0.0D0
      DO L=1,NLAM
        XG(L)=0.0D0
        YG(L)=0.0D0
        MATT(L)=0
      ENDDO
      NMATG=0
C
C  ****  Geometry definition.
C
      WRITE(IWR,2000)
 2000 FORMAT(//3X,72('-'),/3X,'>>>>>>  Cylindrical Geometry.',/)
C
    1 CONTINUE
      READ(IRD,'(A6,1X,A120)') KWORD,BUFFER
      IF(KWORD.EQ.KWCOMM) GO TO 1
      IF(KWORD.NE.GSTART) THEN
        WRITE(IWR,'(3X,A6,1X,A120)') KWORD,BUFFER
        WRITE(IWR,*) 'A GSTART line is expected here ...'
        STOP 'A GSTART line is expected here ...'
      ENDIF
      WRITE(IWR,2001) BUFFER
 2001 FORMAT(3X,A120,/)
C
      READ(IRD,'(A6,1X,A120)') KWORD,BUFFER
      IF(KWORD.NE.LAYER) THEN
        WRITE(IWR,'(3X,A6,1X,A120)') KWORD,BUFFER
        WRITE(IWR,*) 'A LAYER line is expected here...'
        STOP 'A LAYER line is expected here...'
      ELSE
        READ(BUFFER,*) ZLC,ZUC
      ENDIF
      KLAY=1
C
C  ****  New layer...
C
   10 CONTINUE
      IF(ZLC.GT.ZUC-FEPS) THEN
        WRITE(IWR,2100) KLAY,ZLC,ZUC
        WRITE(IWR,*) 'Layer thickness is too small.'
        STOP 'Layer thickness is too small.'
      ENDIF
 2100 FORMAT(3X,'Layer ',I2,':  ZL =',1P,E13.6,',  ZU =',E13.6)
C
      IF(KLAY.EQ.1) THEN
        WRITE(IWR,2100) KLAY,ZLC,ZUC
        ZG(1)=ZLC
        ZG(2)=ZUC
      ELSE
        IF(ZLC.LT.ZG(KLAY)-EPS) THEN
          WRITE(IWR,2100) KLAY,ZLC,ZUC
          WRITE(IWR,*) 'Overlapping layers.'
          STOP 'Overlapping layers.'
        ENDIF
        IF(ZLC.GT.ZG(KLAY)+FEPS) THEN
          IF(KLAY.GT.NLAM-1) THEN
            WRITE(IWR,2100) KLAY,ZLC,ZUC
            WRITE(IWR,*) 'Too many layers. Increase the value of NLAM.'
            STOP 'Too many layers.'
          ENDIF
          ZG(KLAY+1)=ZLC
          RG(KLAY,2)=1.0D-10
          MATER(KLAY,1)=0
          KBOD=KBOD+1
          NCYL(KLAY)=1
          WRITE(IWR,2100) KLAY,ZG(KLAY),ZLC
          KLAY=KLAY+1
        ENDIF
        WRITE(IWR,2100) KLAY,ZLC,ZUC
        ZG(KLAY+1)=ZUC
      ENDIF
C
      KCYL=0
      RG(KLAY,1)=0.0D0
C
C  ****  Structure of the layer...
C
   20 CONTINUE
      READ(IRD,'(A6,1X,A120)') KWORD,BUFFER
      IF(KWORD.EQ.GEND) THEN
        NCYL(KLAY)=KCYL
        GO TO 30
      ELSE IF(KWORD.EQ.LAYER) THEN
        NCYL(KLAY)=KCYL
        IF(MATT(KLAY).EQ.0) KLAY=KLAY-1
        READ(BUFFER,*) ZLC,ZUC
        KLAY=KLAY+1
        IF(KLAY.GT.NLAM) THEN
          WRITE(IWR,*) 'Too many layers. Increase the parameter NLAM.'
          STOP 'Too many layers.'
        ENDIF
        GO TO 10
      ELSE IF(KWORD.EQ.CENTRE) THEN
        READ(BUFFER,*) XG(KLAY),YG(KLAY)
        WRITE(IWR,2200) XG(KLAY),YG(KLAY)
 2200   FORMAT(14X,'XG =',1P,E13.6,',  YG =',E13.6)
        GO TO 20
      ELSE IF(KWORD.EQ.CYLIND) THEN
        READ(BUFFER,*) MAT,RIC,ROC
        WRITE(IWR,2300) KCYL+1,RIC,ROC,MAT
 2300     FORMAT(5X,'Cyl ',I2,':  RI =',1P,E13.6,',  RO =',E13.6,
     1      ',  MAT  =',I3)
        IF(MAT.LT.0) MAT=0
        IF(RIC.GT.ROC-FEPS) THEN
          WRITE(IWR,*) 'The ring thickness is too small.'
          STOP 'The ring thickness is too small.'
        ENDIF
        IF(KCYL.GT.NCYM-1) THEN
          WRITE(IWR,*) 'Too many rings. Increase the parameter NCYM.'
          STOP 'Too many rings.'
        ENDIF
        IF(RIC.LT.RG(KLAY,KCYL+1)-EPS) THEN
          WRITE(IWR,*) 'Overlapping rings.'
          STOP 'Overlapping rings.'
        ENDIF
        IF(RIC.GT.RG(KLAY,KCYL+1)+FEPS) THEN
          KCYL=KCYL+1
          RG(KLAY,KCYL+1)=RIC
          MATER(KLAY,KCYL)=0
        ENDIF
        IF(KCYL.GT.NCYM-1) THEN
          WRITE(IWR,*) 'Too many rings. Increase the parameter NCYM.'
          STOP 'Too many rings.'
        ENDIF
        KCYL=KCYL+1
        RG(KLAY,KCYL+1)=ROC
        RMAX=MAX(RMAX,ROC+SQRT(XG(KLAY)**2+YG(KLAY)**2))
        MATER(KLAY,KCYL)=MAT
        NMATG=MAX(NMATG,MAT)
        MATT(KLAY)=MATT(KLAY)+MAT
        GO TO 20
      ELSE
          WRITE(IWR,'(3X,A6,1X,A120)') KWORD,BUFFER
          WRITE(IWR,*) 'Unknown keyword.'
          STOP 'Unknown keyword.'
      ENDIF
C
C  ****  End of the geometry definition list.
C
   30 CONTINUE
      NLAY=KLAY
      DO KL=1,NLAY
        RC=SQRT(XG(KL)**2+YG(KL)**2)
        IF(RG(KL,NCYL(KL)+1).LT.RMAX+RC-FEPS) THEN
          IF(MATER(KL,NCYL(KL)).EQ.0) THEN
            RG(KL,NCYL(KL)+1)=RMAX+RC
          ELSE
            NCYL(KL)=NCYL(KL)+1
            RG(KL,NCYL(KL)+1)=RMAX+RC
            MATER(KL,NCYL(KL))=0
          ENDIF
        ENDIF
      ENDDO
      RMAX2=RMAX**2
C
      WRITE(IWR,2400)
 2400 FORMAT(/3X,70('-'),/3X,'| Body  KL  KC',7X,'ZL;ZU',13X,
     1  'XC;YC',11X,'RI;RO',5X,'Mat |')
      NBOD=0
      DO L=1,NLAY
        DO J=1,NCYL(L)
          RG2(L,J)=RG(L,J)**2
          NBOD=NBOD+1
          IBOD(L,J)=NBOD
          ILAY(NBOD)=L
          ICYL(NBOD)=J
          IF(J.EQ.1) THEN
            WRITE(IWR,2500) NBOD,ILAY(NBOD),ICYL(NBOD),ZG(L),XG(L),
     1        RG(L,J),ZG(L+1),YG(L),RG(L,J+1),MATER(L,J)
          ELSE
            WRITE(IWR,2501) NBOD,ILAY(NBOD),ICYL(NBOD),
     1        RG(L,J),RG(L,J+1),MATER(L,J)
          ENDIF
        ENDDO
        RG2(L,NCYL(L)+1)=RG(L,NCYL(L)+1)**2
      ENDDO
 2500 FORMAT(3X,'|',68('-'),'|',/
     1  3X,'|',I4,2X,I3,',',I3,1P,2X,'ZL=',E13.6,2X,'XC=',E13.6,
     1  1X,E13.6,5X,'|',/3X,'|',15X,'ZU=',E13.6,2X,'YC=',E13.6,1X,
     1  E13.6,1X,I3,' |')
 2501 FORMAT(3X,'|',68('-'),'|',/
     1  3X,'|',I4,2X,I3,',',I3,35X,1P,E15.6,
     1  5X,'|',/3X,'|',13X,35X,E15.6,1X,I3,' |')
C
      WRITE(IWR,2502)
 2502 FORMAT(3X,70('-'))
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE LOCATC
C  *********************************************************************
      SUBROUTINE LOCATC
C
C  This subroutine determines the layer and cylinder that contain the
C  point with coordinates (X,Y,Z).
C
C  Input values (common /TRACK/):
C     X, Y, Z ... coordinates of the particle,
C     U, V, W ... direction of movement.
C
C  Output values (common /TRACK/):
C     IBODY ..... cylindrical ring where the particle moves,
C     MAT  ...... material in IBODY. MAT=0, indicates that the particle
C                 is outside the material system.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-9,TEPS=EPS+EPS)
C  ****  Main-PENELOPE commons.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C  ****  Cylindrical geometry.
      PARAMETER (NLAM=100,NCYM=50,NBDM=NLAM*NCYM)
      COMMON/CYLGEO/XG(NLAM),YG(NLAM),ZG(NLAM),RG(NLAM,NCYM),
     1  RG2(NLAM,NCYM),RMAX,RMAX2,IBOD(NLAM,NCYM),MATER(NLAM,NCYM),
     1  ILAY(NBDM),ICYL(NBDM),NLAY,NCYL(NLAM),NBOD
      COMMON/CYLAUX/INOUT,KL,KC
C
C  ****  We first check if the particle is within the enclosure.
C
      R2=X*X+Y*Y
      IF(Z.GE.ZG(NLAY+1).OR.Z.LT.ZG(1).OR.R2.GE.RMAX2) THEN
        INOUT=0
        IBODY=0
        MAT=0
        RETURN
      ENDIF
C  ****  If the particle is inside the enclosure, we locate the layer...
      KL=1
      KL1=NLAY+1
    1 KLT=(KL+KL1)/2
      IF(Z.GE.ZG(KLT)) THEN
        KL=KLT
      ELSE
        KL1=KLT
      ENDIF
      IF(KL1-KL.GT.1) GO TO 1
C
      IF(ABS(Z-ZG(KL)).LT.EPS) THEN
        IF(W.GT.0.0D0) THEN
          Z=ZG(KL)+TEPS
        ELSE
          Z=ZG(KL)-TEPS
          KL=KL-1
          IF(KL.LT.1) THEN
            INOUT=0
            IBODY=0
            MAT=0
          ENDIF
        ENDIF
      ELSE IF(ABS(Z-ZG(KL+1)).LT.EPS) THEN
        IF(W.GT.0.0D0) THEN
          Z=ZG(KL+1)+TEPS
          KL=KL+1
          IF(KL.GT.NLAY) THEN
            INOUT=0
            IBODY=0
            MAT=0
          ENDIF
        ELSE
          Z=ZG(KL+1)-TEPS
        ENDIF
      ENDIF
C
C  ****  ... and the cylinder.
C
      R2=(X-XG(KL))**2+(Y-YG(KL))**2
      KC=1
      KC1=NCYL(KL)+1
    2 KCT=(KC+KC1)/2
      IF(R2.GE.RG2(KL,KCT)) THEN
        KC=KCT
      ELSE
        KC1=KCT
      ENDIF
      IF(KC1-KC.GT.1) GO TO 2
C
      R=SQRT(R2)
      IF(KC.GT.1.AND.ABS(R-RG(KL,KC)).LT.EPS) THEN
        TST=(X-XG(KL))*U+(Y-YG(KL))*V
        IF(TST.GT.0.0D0) THEN
          FACT=1.0D0+TEPS/RG(KL,KC)
          X=XG(KL)+(X-XG(KL))*FACT
          Y=YG(KL)+(Y-YG(KL))*FACT
        ELSE
          FACT=1.0D0-TEPS/RG(KL,KC)
          X=XG(KL)+(X-XG(KL))*FACT
          Y=YG(KL)+(Y-YG(KL))*FACT
          KC=KC-1
        ENDIF
      ELSE IF(ABS(R-RG(KL,KC+1)).LT.EPS) THEN
        TST=(X-XG(KL))*U+(Y-YG(KL))*V
        IF(TST.GT.0.0D0) THEN
          FACT=1.0D0+TEPS/RG(KL,KC+1)
          X=XG(KL)+(X-XG(KL))*FACT
          Y=YG(KL)+(Y-YG(KL))*FACT
          KC=KC+1
          IF(KC.GT.NCYL(NLAY).AND.(X*X+Y*Y).GT.RMAX2) THEN
            INOUT=0
            IBODY=0
            MAT=0
          ENDIF
        ELSE
          FACT=1.0D0-TEPS/RG(KL,KC+1)
          X=XG(KL)+(X-XG(KL))*FACT
          Y=YG(KL)+(Y-YG(KL))*FACT
        ENDIF
      ENDIF
C
      INOUT=1
      IBODY=IBOD(KL,KC)
      MAT=MATER(KL,KC)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE STEPC
C  *********************************************************************
      SUBROUTINE STEPC(DS,DSEF,NCROSS)
C
C  This subroutine handles the geometrical part of the track simulation
C  in cylindrical geometries. The particle starts from the point (X,Y,Z)
C  and travels a length DS in the direction (U,V,W) within the material
C  where it moves. When the track leaves the initial material, the par-
C  ticle is stopped just after entering the next material body (void
C  regions with MAT=0 are crossed automatically). Furthermore, when the
C  particle arrives from a void region, it is stopped just after enter-
C  ing the first material body.
C
C  Input values (common /TRACK/):
C     X, Y, Z ... coordinates of the initial point,
C     U, V, W ... direction cosines of the displacement,
C     IBODY ..... body where the initial point is located,
C     MAT ....... material in body IBODY.
C
C  Input argument:
C     DS ........ path length to travel.
C
C  Output arguments:
C     DSEF....... travelled path length before leaving the initial
C                 material or completing the jump (less than DS if the
C                 track crosses an interface),
C     NCROSS .... = 0 if the whole step is contained in the initial
C                   material,
C                 .gt.0 if the particle has crossed an interface, i.e.
C                   if it has entered a new material.
C
C  Output values (common /TRACK/):
C     X, Y, Z ... coordinates of the final position,
C     IBODY ..... body where the final point is located,
C     MAT ....... material in IBODY. The value MAT=0 indicates that the
C                 particle has escaped from the system.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-9,TEPS=EPS+EPS)
C  ****  Main-PENELOPE commons.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C  ****  Cylindrical geometry.
      PARAMETER (NLAM=100,NCYM=50,NBDM=NLAM*NCYM)
      COMMON/CYLGEO/XG(NLAM),YG(NLAM),ZG(NLAM),RG(NLAM,NCYM),
     1  RG2(NLAM,NCYM),RMAX,RMAX2,IBOD(NLAM,NCYM),MATER(NLAM,NCYM),
     1  ILAY(NBDM),ICYL(NBDM),NLAY,NCYL(NLAM),NBOD
      COMMON/CYLAUX/INOUT,KL,KC
      DIMENSION S(5),IS(5)
C
      MAT0=MAT
      DSEF=0.0D0
      DSREM=DS
      NCROSS=0
      IF(ABS(W).GT.1.0D-30) THEN
        WW=W
      ELSE
        IF(W.GT.0.0D0) THEN
          WW=W+1.0D-30
        ELSE
          WW=W-1.0D-30
        ENDIF
      ENDIF
C
C  ****  The particle starts from the outside the enclosure.
C
      IF(INOUT.EQ.0) THEN
        S(1)=(ZG(1)-Z)/WW
        IF(S(1).LT.0.0D0) S(1)=1.0D17
        IS(1)=1
C
        S(2)=(ZG(NLAY+1)-Z)/WW
        IF(S(2).LT.0.0D0) S(2)=1.0D17
        IS(2)=2
C
        R2=X*X+Y*Y
        XUYV=X*U+Y*V
        U2V2=MAX(U*U+V*V,1.0D-16)
        DISC=XUYV**2+U2V2*(RMAX2-R2)
        IF(DISC.GT.0.0D0) THEN
          DISC=SQRT(DISC)
          SP=(DISC-XUYV)/U2V2
          IF(SP.LT.0.0D0) SP=1.0D17
          SQ=(-DISC-XUYV)/U2V2
          IF(SQ.LT.0.0D0) SQ=1.0D17
          S(3)=MIN(SP,SQ)
        ELSE
          S(3)=1.0D15
        ENDIF
        IS(3)=3
C  ****  Sort surfaces according to path-length.
        DO I=1,2
          SMIN=S(I)
          IMIN=I
          DO J=I+1,3
            IF(S(J).LT.SMIN) THEN
              SMIN=S(J)
              IMIN=J
            ENDIF
          ENDDO
          IF(IMIN.NE.I) THEN
            RSAVE=S(I)
            S(I)=S(IMIN)
            S(IMIN)=RSAVE
            ISAVE=IS(I)
            IS(I)=IS(IMIN)
            IS(IMIN)=ISAVE
          ENDIF
        ENDDO
        S(4)=S(3)+100.0D0
C
        DO I=1,3
          SS=S(I)
          IF(SS.GT.1.0D16) THEN
            X=X+U*SS
            Y=Y+V*SS
            Z=Z+W*SS
            DSEF=SS
            MAT=0
            RETURN
          ELSE
            XT=X+U*SS
            YT=Y+V*SS
            ZT=Z+W*SS
            IF(IS(I).EQ.1) THEN
              ZT=ZT+TEPS
            ELSE IF(IS(I).EQ.2) THEN
              ZT=ZT-TEPS
            ELSE
              XT=XT*(1.0D0-TEPS/RMAX)
              YT=YT*(1.0D0-TEPS/RMAX)
            ENDIF
            RT2=XT*XT+YT*YT
            IF(ZT.LT.ZG(NLAY+1).AND.ZT.GE.ZG(1).AND.RT2.LT.RMAX2) THEN
              IF(SS.LT.S(I+1)-TEPS) THEN
                DSEF=SQRT((XT-X)**2+(YT-Y)**2+(ZT-Z)**2)
                X=XT
                Y=YT
                Z=ZT
                INOUT=1
                GO TO 10
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C  ****  The particle does not enter the enclosure.
        SS=1.0D17
        DSEF=SS
        X=X+U*SS
        Y=Y+V*SS
        Z=Z+W*SS
        MAT=0
        RETURN
      ELSE
         R2=(X-XG(KL))**2+(Y-YG(KL))**2
         IF(Z.GE.ZG(KL).AND.Z.LT.ZG(KL+1).AND.R2.GE.RG2(KL,KC).
     1     AND.R2.LT.RG2(KL,KC+1)) THEN
           IF((X*X+Y*Y).LT.RMAX2) THEN
             GO TO 40
           ELSE
             CALL LOCATC
             RETURN
           ENDIF
         ELSE
           CALL LOCATC
           RETURN
         ENDIF
      ENDIF
C
C  ****  At this point, the particle is inside the enclosure, but it may
C        be in a void volume.
C
C  ****  We locate the layer...
   10 CONTINUE
      KL=1
      KL1=NLAY+1
    1 KLT=(KL+KL1)/2
      IF(Z.GE.ZG(KLT)) THEN
        KL=KLT
      ELSE
        KL1=KLT
      ENDIF
      IF(KL1-KL.GT.1) GO TO 1
C
      IF(ABS(Z-ZG(KL)).LT.EPS) THEN
        IF(W.GT.0.0D0) THEN
          Z=ZG(KL)+TEPS
        ELSE
          Z=ZG(KL)-TEPS
          KL=KL-1
          IF(KL.LT.1) THEN
            X=X+U*1.0D17
            Y=Y+V*1.0D17
            Z=Z+W*1.0D17
            IF(MAT.NE.0) NCROSS=NCROSS+1
            MAT=0
            INOUT=0
            RETURN
          ENDIF
        ENDIF
      ELSE IF(ABS(Z-ZG(KL+1)).LT.EPS) THEN
        IF(W.GT.0.0D0) THEN
          Z=ZG(KL+1)+TEPS
          KL=KL+1
          IF(KL.GT.NLAY) THEN
            X=X+U*1.0D17
            Y=Y+V*1.0D17
            Z=Z+W*1.0D17
            IF(MAT.NE.0) NCROSS=NCROSS+1
            MAT=0
            INOUT=0
            RETURN
          ENDIF
        ELSE
          Z=ZG(KL+1)-TEPS
        ENDIF
      ENDIF
C  ****  ... and the cylinder.
   20 CONTINUE
      R2=(X-XG(KL))**2+(Y-YG(KL))**2
      KC=1
      KC1=NCYL(KL)+1
    2 KCT=(KC+KC1)/2
      IF(R2.GE.RG2(KL,KCT)) THEN
        KC=KCT
      ELSE
        KC1=KCT
      ENDIF
      IF(KC1-KC.GT.1) GO TO 2
C
      R=SQRT(R2)
      IF(KC.GT.1.AND.ABS(R-RG(KL,KC)).LT.EPS) THEN
        TST=(X-XG(KL))*U+(Y-YG(KL))*V
        IF(TST.GT.0.0D0) THEN
          FACT=1.0D0+TEPS/RG(KL,KC)
          X=XG(KL)+(X-XG(KL))*FACT
          Y=YG(KL)+(Y-YG(KL))*FACT
        ELSE
          FACT=1.0D0-TEPS/RG(KL,KC)
          X=XG(KL)+(X-XG(KL))*FACT
          Y=YG(KL)+(Y-YG(KL))*FACT
          KC=KC-1
        ENDIF
      ELSE IF(ABS(R-RG(KL,KC+1)).LT.EPS) THEN
        TST=(X-XG(KL))*U+(Y-YG(KL))*V
        IF(TST.GT.0.0D0) THEN
          FACT=1.0D0+TEPS/RG(KL,KC+1)
          X=XG(KL)+(X-XG(KL))*FACT
          Y=YG(KL)+(Y-YG(KL))*FACT
          KC=KC+1
          IF((X*X+Y*Y).GE.RMAX2) THEN
            X=X+U*1.0D17
            Y=Y+V*1.0D17
            Z=Z+W*1.0D17
            MAT=0
            INOUT=0
            RETURN
          ENDIF
        ELSE
          FACT=1.0D0-TEPS/RG(KL,KC+1)
          X=XG(KL)+(X-XG(KL))*FACT
          Y=YG(KL)+(Y-YG(KL))*FACT
        ENDIF
      ENDIF
C
   30 CONTINUE
      IBODY=IBOD(KL,KC)
      MATN=MATER(KL,KC)
      IF(MATN.NE.MAT) THEN
        NCROSS=NCROSS+1
        IF(MATN.EQ.0) THEN
          MAT=0
        ELSE
          MAT=MATN
          RETURN
        ENDIF
      ENDIF
C
C  ****  We determine the distances to the near surfaces.
C
   40 CONTINUE
      IF(MAT.EQ.MAT0.AND.MAT.NE.0) THEN
        R=SQRT(R2)
        TSMIN=MIN(Z-ZG(KL),ZG(KL+1)-Z,R-RG(KL,KC),RG(KL,KC+1)-R)-TEPS
        IF(DSREM.LT.TSMIN) THEN
          X=X+U*DSREM
          Y=Y+V*DSREM
          Z=Z+W*DSREM
          DSEF=DSEF+DSREM
          RETURN
        ENDIF
      ENDIF
C  ****  Lower plane.
      IT=0
      SS=(ZG(KL)-Z)/WW
      IF(SS.GE.0.0D0) THEN
        IT=IT+1
        S(IT)=SS
        IS(IT)=1
      ENDIF
C  ****  Upper plane.
      SS=(ZG(KL+1)-Z)/WW
      IF(SS.GE.0.0D0) THEN
        IT=IT+1
        S(IT)=SS
        IS(IT)=2
      ENDIF
C  ****  Outer cylinder.
      XUYV=(X-XG(KL))*U+(Y-YG(KL))*V
      U2V2=MAX(U*U+V*V,1.0D-16)
      DISC=XUYV**2+U2V2*(RG2(KL,KC+1)-R2)
      IF(DISC.GE.0.0D0) THEN
        DISC=SQRT(DISC)
        SP=(DISC-XUYV)/U2V2
        IF(SP.LT.0.0D0) SP=1.0D17
        SQ=(-DISC-XUYV)/U2V2
        IF(SQ.LT.0.0D0) SQ=1.0D17
        SS=MIN(SP,SQ)
        IF(SS.LT.1.0D16) THEN
          IT=IT+1
          S(IT)=SS
          IS(IT)=3
        ENDIF
      ENDIF
C  ****  Inner cylinder.
      IF(KC.GT.1) THEN
        DISC=XUYV**2+U2V2*(RG2(KL,KC)-R2)
        IF(DISC.GE.0.0D0) THEN
          DISC=SQRT(DISC)
          SP=(DISC-XUYV)/U2V2
          IF(SP.LT.0.0D0) SP=1.0D17
          SQ=(-DISC-XUYV)/U2V2
          IF(SQ.LT.0.0D0) SQ=1.0D17
          SS=MIN(SP,SQ)
          IF(SS.LT.1.0D16) THEN
            IT=IT+1
            S(IT)=SS
            IS(IT)=4
          ENDIF
        ENDIF
      ENDIF
C  ****  Enclosure.
      IF(KC.EQ.NCYL(KL))THEN
        XUYV=X*U+Y*V
        U2V2=MAX(U*U+V*V,1.0D-16)
        DISC=XUYV**2+U2V2*RMAX2
        IF(DISC.GE.0.0D0) THEN
          DISC=SQRT(DISC)
          SP=(DISC-XUYV)/U2V2
          IF(SP.LT.0.0D0) SP=1.0D17
          SQ=(-DISC-XUYV)/U2V2
          IF(SQ.LT.0.0D0) SQ=1.0D17
          SS=MIN(SP,SQ)
          IF(SS.LT.1.0D16) THEN
            IT=IT+1
            S(IT)=SS
            IS(IT)=5
          ENDIF
        ENDIF
      ENDIF
C  ****  Determine what interface is crossed first.
      IF(IT.GT.0) THEN
        SS=S(1)
        IMIN=1
        IF(IT.GT.1) THEN
          DO I=2,IT
            IF(S(I).LT.SS) THEN
              SS=S(I)
              IMIN=I
            ENDIF
          ENDDO
        ENDIF
      ELSE
        CALL LOCATC
        RETURN
      ENDIF
C
C  ****  Move the particle.
C
      IF(MAT.EQ.MAT0) THEN
        IF(MAT.NE.0) THEN
          IF(DSREM.LT.SS) THEN
            SS=DSREM
            DSEF=DSEF+SS
            X=X+U*SS
            Y=Y+V*SS
            Z=Z+W*SS
            RETURN
          ELSE
            DSEF=DSEF+SS
            DSREM=DSREM-SS
          ENDIF
        ELSE
          DSEF=DSEF+SS
        ENDIF
      ENDIF
C
      X=X+U*SS
      Y=Y+V*SS
      Z=Z+W*SS
C
      IF(IS(IMIN).EQ.1) THEN        ! Lower plane
        Z=ZG(KL)-TEPS
        KL=KL-1
        IF(KL.LT.1) THEN
          X=X+U*1.0D17
          Y=Y+V*1.0D17
          Z=Z+W*1.0D17
          IF(MAT.NE.0) NCROSS=NCROSS+1
          MAT=0
          INOUT=0
          RETURN
        ENDIF
        GO TO 20
      ELSE IF(IS(IMIN).EQ.2) THEN   ! Upper plane
        KL=KL+1
        Z=ZG(KL)+TEPS
        IF(KL.GT.NLAY) THEN
          X=X+U*1.0D17
          Y=Y+V*1.0D17
          Z=Z+W*1.0D17
          IF(MAT.NE.0) NCROSS=NCROSS+1
          MAT=0
          INOUT=0
          RETURN
        ENDIF
        GO TO 20
      ELSE IF(IS(IMIN).EQ.3) THEN   ! Outer cylinder
        KC=KC+1
        FACT=1.0D0+TEPS/RG(KL,KC)
        X=XG(KL)+(X-XG(KL))*FACT
        Y=YG(KL)+(Y-YG(KL))*FACT
        IF((X*X+Y*Y).GT.RMAX2) THEN
          X=X+U*1.0D17
          Y=Y+V*1.0D17
          Z=Z+W*1.0D17
          IF(MAT.NE.0) NCROSS=NCROSS+1
          MAT=0
          INOUT=0
          RETURN
        ENDIF
        IF(Z.GT.ZG(KL)+TEPS.AND.Z.LT.ZG(KL+1)-TEPS) THEN
          R2=(X-XG(KL))**2+(Y-YG(KL))**2
          GO TO 30
        ELSE
          GO TO 10
        ENDIF
      ELSE IF(IS(IMIN).EQ.4) THEN   ! Inner cylinder
        FACT=1.0D0-TEPS/RG(KL,KC)
        KC=KC-1
        X=XG(KL)+(X-XG(KL))*FACT
        Y=YG(KL)+(Y-YG(KL))*FACT
        IF(Z.GT.ZG(KL)+TEPS.AND.Z.LT.ZG(KL+1)-TEPS) THEN
          R2=(X-XG(KL))**2+(Y-YG(KL))**2
          GO TO 30
        ELSE
          GO TO 10
        ENDIF
      ELSE                          ! Enclosure.
        KC=KC+1
        X=X+U*1.0D17
        Y=Y+V*1.0D17
        Z=Z+W*1.0D17
        IF(MAT.NE.0) NCROSS=NCROSS+1
        MAT=0
        INOUT=0
        RETURN
      ENDIF
      END
