C=========================================================================
C EQCOUNT: Counts the list of species for solving the equation of state by
C          merging the default list and species present in the line list.
C
C We assume that only neutral molecules can appear in the line list.
C For atoms, all the ions present in the table of partition functions
C are added to the list. Atomic names are case sensitive, that is the first
C character must be uppercase and for 2 character names the second character
C must be lower case.
C
C Inputs:
C   ELEMEN - the names of chemical elements in the periodic table
C   SPNAME - the names of the species present in the line lists + continuous
C            absorbers
C   ION    - ionization stage (1 -neutral, 2 - first ion etc.)
C   NLINES - the length of the line list, also dimensions of arrays SPNAME,
C            ION, SPINDX
C   NLIST  - if >0 on input, indicates that the default list of species have
C            been loaded, otherwise EQLIST loads the default list to SPLIST.
C   ELESIZ - Size of all arrays related to atomic list.
C
C   Return code  0: OK
C                1: illegal species name
C               >1: SPLSIZ is too small
C
      integer function eqcount(elemen,spname,ion,nlines,nlist,
     *                         ELESIZ)
c      integer function eqcount(elemen,spname,ion,nlines,nlist,
c     *                         environment,ELESIZ)
      INCLUDE 'SIZES.EOS'

      integer nlines,nlist,ELESIZ
      character*(3) elemen(ELESIZ)
      character*2 tmp
      character*(SPCHAR) spname(nlines)
      character*(SPCHAR) tmplist(SPLSIZ),chname
      integer ion(nlines),ionmax,ionmaxx
      real a(IONSIZ)
c      character*(*) environment
      double precision b(IONSIZ)
      INCLUDE 'DEFAULT.EOS.current'
c      INCLUDE 'DEFAULT.EOS'
C
      eqcount=0
      ionmax=0
      ncount=NDEF
c      if(environment.eq.'COLD'.or.environment.eq.'cold') then
c        do ispec=1,NDEF_cold
c          tmplist(ispec)=default_cold(ispec)
c        end do
c        ionmax=2
c        ncount=NDEF_cold
c      else if(environment.eq.'WARM'.or.environment.eq.'warm') then
c        do ispec=1,NDEF_warm
c          tmplist(ispec)=default_warm(ispec)
c        end do
c        ncount=NDEF_warm
c      else if(environment.eq.'HOT'.or.environment.eq.'hot') then
c        do ispec=1,NDEF_hot
c          tmplist(ispec)=default_hot(ispec)
c        end do
c        ncount=NDEF_hot
c      else
c        do ispec=1,NDEF_cool
c          tmplist(ispec)=default_cool(ispec)
c        end do
c        ncount=NDEF_cool
c      end if
C
C Associate each species in SPNAME with an entry in SPLIST. If SPNAME
C  contains a new species not in SPLIST, then add that new species at
C  the end of SPLIST.
C
      if(nlines.gt.0) then
        do 6 ilin=1,nlines
          call mbuild(spname(ilin),ion(ilin)-1,chname)
c          write(*,*) ncount,ilin,ionmax,spname(ilin),chname
          do ispec=1,ncount
            if(tmplist(ispec).eq.chname) goto 6
          end do
c          write(*,*) ncount,ilin,chname,ionmax,spname(ilin),ion(ilin)
c       stop
C
C Look for atomic species. Negative ions (e.g. H-) are treated as molecules
C
          if((spname(ilin)(2:2).EQ.' '.OR.
     *       (spname(ilin)(3:3).EQ.' '.AND.
     *        spname(ilin)(2:2).GE.'a'.AND.
     *        spname(ilin)(2:2).LE.'z')).AND.
     *        ion(ilin).GT.0) then
            iel=0
            tmp=spname(ilin)(1:2)
            do i=1,ELESIZ
              if(tmp.eq.elemen(i)(1:2)) then
                iel=i
                goto 4
              endif
            end do
            if(iel.lt.1) then
              eqcount=1
c              return
              write(*,*) 'eqcount: Wrong species: ',spname(ilin)
              stop
            end if
   4        call XSAHA(iel,1.,1.,1.,ionmaxx,a,b,5)
            if(ionmax.gt.0) ionmaxx=ionmax
            if(ionmaxx.lt.ion(ilin)) then
              write(*,*) ilin,ion(ilin),nlines
              write(*,*) 'XSAHA has no partition function for '//chname
              stop
            endif
            tmplist(ncount+1)=elemen(iel)(1:2)
            if(ionmaxx.gt.1) then
              do i=2,ionmaxx
                ncount=ncount+1
                i1=index(tmplist(ncount),' ')
                tmplist(ncount+1)=tmplist(ncount)(1:i1-1)//'+'
              end do
            end if
            ncount=ncount+1
          else
C
C Molecules are counted here
C
            tmplist(ncount+1)=chname
            ncount=ncount+1
          end if
   6    continue
      endif
C
C All lines have been processed, add free electrons and return
C
      nlist=ncount+1
      eqcount=0
C
      return
      end

C=========================================================================
C EQLIST: Creates the list of species for solving the equation of state by
C         merging the default list and species present in the line list.
C
C We assume that only neutral molecules can appear in the line list.
C For atoms, all the ions present in the table of partition functions
C are added to the list. Atomic names are case sensitive, that is the first
C character must be uppercase and for 2 character names the second character
C must be lower case.
C
C Inputs:
C   ELEMEN - the names of chemical elements in the periodic table
C   SPNAME - the names of the species present in the line lists + continuous
C            absorbers
C   ION    - ionization stage (1 -neutral, 2 - first ion etc.)
C   NLINES - the length of the line list, also dimensions of arrays SPNAME,
C            ION, SPINDX
C   NLIST  - if >0 on input, indicates that the default list of species have
C            been loaded, otherwise EQLIST loads the default list to SPLIST.
C   SPLDIM - maximum length of the compiled lists of species SPLIST (must
C            be smaller than SPLSIZ).
C   ELESIZ - Size of all arrays related to atomic list.
C
C Outputs:
C   SPINDX - index array of size NLINES which upon return holds pointers to
C            the complete list of species SPLIST: line L is produced by
C            species SPLIST(SPINDEX(L))
C   SPLIST - upon return contains the compiled list of all species (default
C            list + species in the line list + continuous absorbers)
C   NLIST  - the size of the compiled list of species SPLIST
C
C   Return code  0: OK
C                1: illegal species name
C                2: SPLDIM is too small)
C                3: Missing ionization stage
C                4: e- is not the last item in the list
C                5: Unreasonable abundances
C
C  2006.12.27 - converted eqlist to a function for compatibility with the SME
C
C
      integer*4 function eqlist(abund,elemen,spname,ion,spindx,splist,
     &                  nlines,nlist,SPLDIM,ELESIZ)
c      integer*4 function eqlist(abund,elemen,spname,ion,spindx,splist,
c     &                  nlines,nlist,environment,SPLDIM,ELESIZ)
      INCLUDE 'SIZES.EOS'

      integer nlines,nlist,SPLDIM,ELESIZ
      character*(SPCHAR) spname(nlines),splist(SPLDIM)
      character*(3) elemen(ELESIZ)
c      character*(*) environment
      character*2 tmp
      integer ion(nlines),spindx(nlines),ionmax,ionmaxx
      dimension abund(ELESIZ)
      real a(IONSIZ)
      double precision b(IONSIZ)
C
C SPLIST should contain all the major contributors to the electron pressure,
C and all the molecules which significantly affect atomic partial pressures.
C For each call to EQSTAT, the base set of species at the beginning of SPLIST
C are supplemented by any new species that appear in SPNAME. It is common
C for some of the species in the base set (at the beginning of SPNAME) to be
C duplicated in SPNAME. This allows one to get ZETA for these species and is
C not a problem.
C
      integer splmax
      character*(SPCHAR) chname
      INCLUDE 'DEFAULT.EOS.current'
c      INCLUDE 'DEFAULT.EOS'
C
C Determine maximum allowed number of species, based on sizes of arrays
C  defined locally (using SPLSIZ) and passed by argument (using spldim).
C
      splmax=min(SPLSIZ,SPLDIM)
C
C Load base set of species (SPLIST) with default set of species (DEFAULT),
C  if passed value of NLIST is 0. Be sure to include "e-" at the end of
C  SPLIST.
C
      idef=0
      ionmax=0
      if(nlist.eq.0) then
C
C  Copy the default list and check if we have enough space first
C
        do jdef=1,NDEF
          splist(jdef)=default(jdef)
        end do
        nlist=NDEF
cC
cC  Copy the default list and check if we have enough space first
cC
c        if(environment.eq.'COLD'.or.environment.eq.'cold') then
c          do jdef=1,NDEF_cold
c            splist(jdef)=default_cold(jdef)
c          end do
c          ionmax=2
c          nlist=NDEF_cold+idef
c        else if(environment.eq.'WARM'.or.environment.eq.'warm') then
c          do jdef=1,NDEF_warm
c            splist(jdef)=default_warm(jdef)
c          end do
c          nlist=NDEF_warm+idef
c        else if(environment.eq.'HOT'.or.environment.eq.'hot') then
c          do jdef=1,NDEF_hot
c            splist(jdef)=default_hot(jdef)
c          end do
c          nlist=NDEF_hot+idef
c        else
c          do jdef=1,NDEF_cool
c            splist(jdef)=default_cool(jdef)
c          end do
c          nlist=NDEF_cool
c        end if
        idef=nlist
        if(nlist.ge.splmax) goto 900
C
C  nlines set to -1 indicates that we need to get partial pressures for all atoms
C  This mode is meant for use within VALD
C
        if(nlines.eq.-1) then
c
c Add all atoms first (the call to XSAHA is dummy,
C just to get the number of ions available in the table)
c
          do iel=1,ELESIZ
            call XSAHA(iel,1.,1.,1.,ionmaxx,a,b,5)
            if(ionmax.gt.0) ionmaxx=ionmax
            idef=idef+1
            if(idef.gt.splmax) goto 900
            splist(idef)=elemen(iel)(1:2)
            if(ionmaxx.gt.1) then
              do i=2,ionmaxx
                idef=idef+1
                if(idef.gt.splmax) goto 900
                splist(idef)=splist(idef-1)
                isp=index(splist(idef),' ')
                if(isp.le.0) then
                  write(*,*) 'eqlist: Insufficient length of splist ',
     *                       'elements to store ion',elemen(iel)(1:2),i,
     *                       idef,SPCHAR
                  eqlist=2
                  return
                endif
                splist(idef)(isp:isp)='+'
              end do
            end if
          end do
          nlist=idef
        endif
      endif
C
C Check that abundances are sensible.
C
      absum=0.0
      do ielem=1,ELESIZ
        if(abund(ielem).lt.0.0.or.abund(ielem).gt.1.0) then
          write(*,40) ielem,abund(ielem)
  40      format('eqlist: bad abundance for element',i3,':',1pe13.4)
          write(*,*) (abund(ispec),ispec=1,99)
c          stop
          eqlist=5
          return
        endif
        absum=absum+abund(ielem)
      end do
c      do ielem=1,ELESIZ
c        abund(ielem)=abund(ielem)/absum
c      end do
c      if(abs(absum-1.0).gt.1.0e-3) then
c        write(*,70) absum
c  70    format('eqlist: warning! abundances are not normalized:'
c     &           ,1pe13.5)
c      endif

C
C Associate each species in SPNAME with an entry in SPLIST. If SPNAME
C  contains a new species not in SPLIST, then add that new species at
C  the end of SPLIST.
C
      do ispec=nlist+1,splmax
        splist(ispec)='        '
      end do
      inew=nlist+1
      if(nlines.gt.0) then
        do 150 ilin=1,nlines
          call mbuild(spname(ilin),ion(ilin)-1,chname)
          do ispec=1,nlist
            if(splist(ispec).eq.chname) then
              spindx(ilin)=ispec
              goto 150
            endif
          end do
C
C Look for atomic species. Negative ions (e.g. H-) are treated as molecules
C
          if((spname(ilin)(2:2).EQ.' '.OR.
     *       (spname(ilin)(3:3).EQ.' '.AND.
     *        spname(ilin)(2:2).GE.'a'.AND.
     *        spname(ilin)(2:2).LE.'z')).AND.
     *        ion(ilin).GT.0) then
            iel=0
            tmp=spname(ilin)(1:2)
            do i=1,ELESIZ
              if(tmp.eq.elemen(i)(1:2)) iel=i
            end do
            if(iel.lt.1) then
c              write(*,*) 'eqlist: Wrong species: "'//spname(ilin)//'"'
c              stop
              eqlist=1
              return
            end if
            call XSAHA(iel,1.,1.,1.,ionmaxx,a,b,5)
            if(ionmax.gt.0) ionmaxx=ionmax
            if(ionmaxx.lt.ion(ilin)) then
              write(*,*) 'XSAHA has no partition function for '//chname
              stop
            endif
C
C  Make sure that neutral atoms are included as well as all
C  the intermediate ions
C
            do ii=0,ionmaxx-1
              if(inew.gt.splmax) goto 900
              call mbuild(spname(ilin),ii,chname)
              splist(inew)=chname
              if(ii.eq.ion(ilin)-1) spindx(ilin)=inew
              inew=inew+1
            end do
          else
c       write(*,*) 'Molecule: '//chname,inew
            if(inew.gt.splmax) goto 900
            splist(inew)=chname
            spindx(ilin)=inew
            inew=inew+1
          end if
          nlist=inew-1
 150    continue
      endif
C
C Make sure free electrons are the last species in the list.
C
      do ispec=1,nlist-1
        if(splist(ispec).eq.'e-') then
c          write(*,*) 'eqlist: "e-" may only occur at the end of the'
c     &            // ' species list (SPLIST).'
c          stop
          eqlist=4
          return
        endif
      end do
      if(splist(nlist).ne.'e-') then
        nlist=nlist+1
        if(nlist.gt.splmax) goto 900
        splist(nlist)='e-'
      endif
C
C Make sure neutral hydrogen and neutral helium are in SPLIST. These
C  species are needed for H1FRCT and HE1FRCT. Remember the locations
C  of these species in SPLIST for later use. Code is optimized for
C  the case where H and He both occur early in SPLIST list.
C
c      ih1=-1
c      do 200 ispec=1,nlist
c        if(splist(ispec).eq.'H') then
c          ih1=ispec
c          goto 210
c        endif
c 200  continue
c      write(*,*) 'eqlist: "H" must be in species list (SPLIST)'
c      stop
c 210  ihe1=-1
c      do 220 ispec=1,nlist
c        if(splist(ispec).eq.'He') then
c          ihe1=ispec
c          goto 230
c        endif
c 220  continue
c      write(*,*) 'eqlist: "He" must be in species list (SPLIST)'
c      stop
c 230  continue
C
C Sort the list
C
      call sort2(nlist,splist,nlines,spindx,elemen,ELESIZ)
c      do 250 ispec=1,nlist
c 250  write(*,*) ispec,' "',splist(ispec),'"'
c      stop
C
      eqlist=0
      return
C
C Error handlers.
C
 900  continue
c      write(*,905) spldim,splsiz
c 905  format('eqlist: species list (SPLIST) not long enough:',2i5)
c      stop
      eqlist=2
c
      return
      end

c
C=========================================================================
C EQSTAT: Determine thermodynamic quantities required for spectroscopy.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   PTOTAL [real] Total gas pressure (in dyne/cm^2), given by NTOTAL*K*T,
C     which is to be used in calculating chemical and ionization equilibrium,
C     and partial pressures.
C   PELEC [real] Electron pressure (in dyne/cm^2), given by NELEC*K*T,
C     which is to be used in calculating ionization equilibrium.
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPNAME [character*(*) array(NLINES)] Case-sensitive species name of atom
C     or molecule. The first letter of each atom name must be uppercase. The
C     second letter of each atom name, if present, must be lowercase. Each
C     atom name may optionally be followed by a multiplicity number between
C     1 and 4. If no multiplicity number is given for a particular atom, then
C     its multiplicity is assumed to be 1. All atomic and molecular species
C     in SPNAME must be neutral, with the charge state specified separately
C     in the ION input argument.
C   ION [integer array(NLINES)] Charge state for each of the atomic and
C     molecular species specified in SPNAME. ION=-1 for negative ions (e.g.
C     H minus), ION=0 for neutrals, ION=1 for singly ionized species, etc.
C   NLINES [integer] Number of valid entries in SPNAME and ION. From an
C     external perspective, each entry in SPNAME and ION will correspond to
C     a single spectral line, so some specie/charge combinations may appear
C     more than once, while others may not appear at all.
C   SPLDIM [integer] Array sizes for the arguments SPLIST and XFRACT, which
C     contain information for each species. The maximum allowed number of
C     species is SPLMAX=MIN(SPLSIZ,SPLDIM), where SPLSIZ is a parameter
C     defined in the file SIZES.SYN and used to dimension the local arrays
C     XNPF, PFUNC, and POTION. SPLMAX must be large enough to handle the
C     base set of species used when computing the molecular equilibrium and
C     also any additional species that appear only in the line list. Ideally,
C     the calling routine will <1> Include SIZES.SYN, <2> Use SPLSIZ to
C     dimension SPLIST and XFRACT, and <3> Pass SPLSIZ in place of SPLDIM.
C     However, SPLDIM is passed separately to allow for error checking in
C     the cases when this is not done (e.g. when called from IDL).
C   MODE [integer] Determines the content of the content of the the output
C                  array xfract:
C      0    - number densities/partition functions
C      1    - number densities
C      2    - partial pressures
C      3    - number density of free electrons produced by each species
C others    - the same as 0
C   10+     - the same as above but electron density is assumed to be known
C             precisely so the input value is used instead of solving for
C             Pelec
C
C Input/Output:
C   SPLIST [character*(*) array(SPLDIM)] If NLIST is nonzero upon entry,
C     then SPLIST must contain the base set of species that must be included
C     in the molecular equilibrium calculation, regardless of which species
C     are represented by lines in SPNAME. Until the code is cleaned up, the
C     species list in SPLIST must include "e-" after the NLIST element.
C     If NLIST is zero upon entry, then SPLIST is loaded with the base set
C     of species coded into EQSTAT below (in the variable DEFAULT). Again,
C     an "e-" is appended after the base set.
C     Regardless of the whether SPLIST is valid upon entry or needs to be
C     loaded with the defaults, species that are in the lines list SPNAME,
C     but are not in the base set of species will be inserted into SPLIST
C     after the "e-" entry. Currently, the extended list is not used, but
C     in the future, we may solve for the equilibrium of all species in the
C     extended SPLIST.
C   NLIST [integer] If nonzero upon entry, NLIST is the number of species
C     in the base set of species passed in SPLIST (including the mandatory
C     "e-" at the beginning of the list). If NLIST is zero upon entry, this
C     indicates that the default base set of species coded in EQSTAT should
C     be used. Upon exit, NLIST is set to the number of species in SPLIST,
C     which contains the base set plus any additional species that occur
C     in the line list.
C
C Outputs:
C   SPINDX [integer array(NLINES)] Species index assigned to each line in
C     the input line list (specified by the input arguments SPNAME and ION).
C     The species index is used to reconstruct the species name (in SPLIST)
C     or other values (e.g in XFRACT) computed for each line in the input line
C     list. For example, ZETA(SPINDX(370)) contains the zeta value for the
C     line corresponding to SPNAME(370) and ION(370).
C   XFRACT [real array(SPLDIM)] The physical meaning and units depend on the
C     value of MODE. These values are given for all atomic or molecular
C     species in the same order as in splist.
C   PFUNC  [real array(SPLDIM)] Partition functions for all species in the
C     same order as species are listed in splist.
C   POTI   [real array(SPLDIM)] ionization potential in eV for the
C     corresponding species.
C   ATWGHT [real array(SPLDIM-1)] molecular weights in AMU for the
C     corresponding species.
cC   H1FRCT [real] Number density (in cm^-3) of neutral atomic hydgrogen,
cC     used in computing damping constants (and continuous opacities?).
cC   HE1FRCT [real] Number density (in cm^-3) of neutral atomic helium,
cC     used in computing damping constants (and continuous opacities?).
C   XNe    [real scalar] number density of free electrons per cm^3 as
C     computed by the EQSTAT. For MODE>=10 XNe is simply the input Pelec/kT.
C   XNa    [real scalar] number density of all particles except for free
C     electrons per cm^3 as computed by the EQSTAT.
C   RHO    [real scalar] density in g/cm^3 as computed by the EQSTAT.
C
      subroutine eqstat(mode,temp,Pg,Pe,abund,elemen,amass,
     &                  ELESIZ,spindx,splist,xfract,pfunc,poti,atwght,
     &                  nlines,nlist,xne,xna,rho,niter)
      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'

      integer mode,ELESIZ,niter
      integer nlines,nlist
      real temp,Tk,Pg,Pe,Pgas,Pelec,xna,xne,rho
      real Pg_old,Pe_old
      character*(SPCHAR) splist(nlist)
      character*(3) elemen(ELESIZ)
      integer spindx(nlines)
      real xfract(nlist),pfunc(nlist),poti(nlist),atwght(nlist)
      real abund(ELESIZ),amass(ELESIZ)
      logical FAILED,BARKLEM

      integer Anum(4),Natm(4),maxion,Nelm,nchg,Ntotal
      real xnpf(SPLSIZ),tol,tol1,xtotal
      real potion(IONSIZ),wtmol
      double precision awt(SPLSIZ-1),fract(IONSIZ),ratiom,part,pion
      integer icharge,iter,ispec,iel,mmode

      INTEGER MAXITER
      REAL kBol
      DOUBLE PRECISION PSI,X,amu,dummy1,dummy2
      PARAMETER (kBol=1.38065E-16,amu=1.66053886D-24,MAXITER=5000)
C
C Call equation of state solver.
C
c      open(87,file='dumpb.dat',form='unformatted',status='old')
c      read(87) temp,Pgas,Pelec,abund,elemen,amass,
c     &         mmode,spindx(nlines),splist,nlines,nlist
c      close(87)
      TOL=1.E-6
      TOL1=1.E-3
      Pgas=Pg
      Pelec=Pe
      PSI=2.d0/(1.d0+SQRT(5.d0))
      do ispec=1,nlist
        xnpf(ispec)=-1.
        pfunc(ispec)=1.
      end do
      Tk=temp*kBol
      mmode=mod(mode,10)

      if(temp.gt.12000.) then
C
C Hot gas: assume no molecules and use Saha equation
C
        niter=1
        if(mode.lt.10) then
C
C Get the number of free electrons, atomic number density and
C mean molecular weight self consistently
C
          call Nelect(temp,Pgas,abund,amass,ELESIZ,
     *                xna,xne,wtmol)
          Pelec=xne*Tk
        else
C
C MODE is larger than 10. Assume the electron pressure to be given.
C Compute mean molecular weight and atom/electron number density
C
          X=0.D0
          do iel=1,ELESIZ
            X=X+abund(iel)*amass(iel)
          end do
          wtmol=X*amu
          xne=Pelec/Tk
          xna=Pgas/Tk-xne
        endif
C
C Density is simple
C
        rho=xna*wtmol
        do 2 ispec=1,nlist-1
        CALL MPARSE(elemen,splist(ispec),Nelm,Nchg,Anum,Natm,ELESIZ)
        icharge=Nchg+1
        if(Nelm.eq.1.and.Natm(1).eq.1.and.Nchg.ge.0) then
C
C Get the number of ionization stages available in XSAHA
C
          call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,5)
C
C Get the partition function for a given species
C
          call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,3)
          pfunc(ispec)=fract(icharge)
C
C Atom. Parser returns atomic number in Anum(1)
C
          if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else if(mmode.eq.2) then
C
C  MODE=2, Return partial pressures
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*kBol*temp*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else if(mmode.eq.3) then
C
C  MODE=3, Return number of free electrons produced
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*kBol*temp*xna*abund(Anum(1))*
     *                    Nchg
            poti(ispec)=potion(icharge)
          else
C
C  Any other MODE: Return number densities / partition functions
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,1)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          endif
          atwght(ispec)=amass(Anum(1))
        else
c        else if(Nchg.ge.0) then
C
C Ignore molecules
C
c          Ntotal=0
c          ratiom=0.d0
c          dummy1=1.d0
c          dummy2=1.d0
c          do iel=1,Nelm
c            Ntotal=Ntotal+Natm(iel)
c            awt(ispec)=awt(ispec)+Natm(iel)*amass(Anum(iel))
c            ratiom=ratiom+Natm(iel)*log10(amass(Anum(iel)))
c          enddo
c          CALL MOLCON(splist(ispec),temp,Ntotal,ratiom,dummy1,
c    &                dummy2,part,pion,BARKLEM)
c          poti(ispec)=pion
c          atwght(ispec)=awt(ispec)
c          pfunc(ispec)=part
c          xfract(ispec)=0.
          if(poti(ispec).lt.0.) then
            poti(ispec)=100.
            atwght(ispec)=10.
          endif
          pfunc(ispec)=1.
          xfract(ispec)=1.e-30
        endif
c        if(Temp.gt.7950.) then
c          write(*,*) ispec,temp,splist(ispec),
c     *               xfract(ispec)*pfunc(ispec),pfunc(ispec),poti(ispec)
c        endif
c        xfract(1)=7.841741E17
c        xfract(3)=6.737E11
c        pfunc(3)=1.
c        xfract(152)=2.66e14
c        pfunc(152)=125.6
c        xfract(153)=6.85d11
c        pfunc(153)=949.2
c        xfract(169)=1.67d8
c        pfunc(169)=15817.
  2     continue
C
C Electrons
C
        if(mmode.eq.1) then
          xfract(nlist)=xne
        else if(mmode.eq.2) then
          xfract(nlist)=Pelec
        else if(mmode.eq.3) then
          xfract(nlist)=1.e-30
        else
          xfract(nlist)=xne
        endif
      else
C
C Cold gas
C
        niter=0
c        write(*,*) NLINES,NLIST,temp,Pgas,Pelec,mmode
c        write(*,'(10f8.3)') log10(abund)
C
C Initioal guess for Pelec
C
        if(mode.lt.10) then
          if(temp.gt.4000.) then
            Pe_old=Pgas*0.1
          else if(temp.gt.2000.) then
            Pe_old=Pgas*0.01
          else
            Pe_old=Pgas*0.001
          endif
        else
C
C If MODE>=10 just use Pelec that is given
C
          Pe_old=Pelec
        endif
        Pg_old=Pg
c        IF(mode.ge.10) then
c          if(temp.gt.4000.) then
c            xne_old=xnatom*0.1
c          else if(temp.gt.2000.) then
c            xne_old=xnatom*0.01
c          else
c            xne_old=xnatom*0.001
c          endif
c        else
c          xne_old=xnelec
c        endif
C
C Solve the molecular/ionization equilibrium using partial pressures (GAS)
C when Pelec is not vanishingly small and log of partial pressures (lnGAS)
C otherwise.
C
  3     continue
        if(temp.lt.2000.) then
          call lnGAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
        else
          call GAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
        endif
        niter=niter+iter
C
C Check if we reached the maximum iterations
C
        Pelec=xne*Tk
        IF(niter.ge.MAXITER) THEN
          WRITE(*,*) 'T,Pg,Pgas,Pelec,Pe_in,Pe_out,NITER=',
     *                Temp,Pg,Pgas,Pe,Pe_old,Pelec,niter,FAILED
          IF(niter.gt.MAXITER*20) STOP
        END IF
C
C Check for convergence. Repeat iterations in case we are not stable yet.
C This external loop is needed because the GAS solver internally uses XSAHA
C to computes the partition functions based on the input value of Pelec.
C The effect of screening is small but it is there and thus outer loop is
C required to reach self-consistency.
C
c        IF(mode.lt.10.and.
        IF(
     *    (abs(Pgas -Pg_old)/max(1.E-20,Pgas ).gt.tol1.or.
     *     abs(Pelec-Pe_old)/max(1.E-20,Pelec).gt.tol1)) THEN
          Pe_old=Pelec
          Pg_old=Pg
          GOTO 3
        END IF
c        write(*,*) Temp,splist(169),xnpf(169),pfunc(169),poti(169)
c        if(Temp.gt.7950.) then
c          do ispec=1,nlist-1
c            write(*,*) ispec,temp,splist(ispec),xnpf(ispec),
c     *                 pfunc(ispec),poti(ispec)
c          enddo
c        endif
c      write(*,'(F10.1,13E11.4)') Temp,xnpf(1),
c     &                                xnpf(2),
c     &                                xnpf(3),
c     &                                xnpf(4),
c     &                                xnpf(5),
c     &                                xnpf(6),
c     & (Pgas-Pelec)/Tk,xna,Pelec/Tk,xne,rho
C
C Fill the return arrays.
C
        do ispec=1,nlist-1
          atwght(ispec)=awt(ispec)
        end do
C
        if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
          do  ispec=1,nlist-1
c          write(*,*) ispec,splist(ispec),xnpf(ispec),pfunc(ispec)
           xfract(ispec)=xnpf(ispec)
          end do
          xfract(nlist)=xne
        else if(mmode.eq.2) then
C
C  MODE=2, Return partial pressures
C
          do ispec=1,nlist-1
            xfract(ispec)=xnpf(ispec)*Tk
          end do
          xfract(nlist)=xne*Tk
        else if(mmode.eq.3) then
C
C  MODE=3, Return number of free electrons
C
          do ispec=1,nlist-1
            call MPARSE(elemen,splist(ispec),nelm,nchg,Anum,Natm,ELESIZ)
            xfract(ispec)=xnpf(ispec)*nchg
          end do
          xfract(nlist)=1.
        else
C
C  Any other MODE: Return number densities / partition functions
C
          do ispec=1,nlist-1
            xfract(ispec)=xnpf(ispec)/pfunc(ispec)
c            write(*,*) ispec,SPLIST(ispec),xnpf(ispec),pfunc(ispec)
          end do
          xfract(nlist)=xne
        endif
      endif
C
      return
      end


C=========================================================================
C EQSTAT_RHO: is identical to EQSTAT except that the density is used
C             instead of the pressure.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   RHO  [real] Total gas density (in g/cm^3),
C     which is to be used in calculating chemical and ionization equilibrium,
C     and partial pressures.
C   PELEC [real] Electron pressure (in dyne/cm^2), given by NELEC*K*T,
C     which is to be used in calculating ionization equilibrium.
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPNAME [character*(*) array(NLINES)] Case-sensitive species name of atom
C     or molecule. The first letter of each atom name must be uppercase. The
C     second letter of each atom name, if present, must be lowercase. Each
C     atom name may optionally be followed by a multiplicity number between
C     1 and 4. If no multiplicity number is given for a particular atom, then
C     its multiplicity is assumed to be 1. All atomic and molecular species
C     in SPNAME must be neutral, with the charge state specified separately
C     in the ION input argument.
C   ION [integer array(NLINES)] Charge state for each of the atomic and
C     molecular species specified in SPNAME. ION=-1 for negative ions (e.g.
C     H minus), ION=0 for neutrals, ION=1 for singly ionized species, etc.
C   NLINES [integer] Number of valid entries in SPNAME and ION. From an
C     external perspective, each entry in SPNAME and ION will correspond to
C     a single spectral line, so some specie/charge combinations may appear
C     more than once, while others may not appear at all.
C   SPLDIM [integer] Array sizes for the arguments SPLIST and XFRACT, which
C     contain information for each species. The maximum allowed number of
C     species is SPLMAX=MIN(SPLSIZ,SPLDIM), where SPLSIZ is a parameter
C     defined in the file SIZES.SYN and used to dimension the local arrays
C     XNPF, PFUNC, and POTION. SPLMAX must be large enough to handle the
C     base set of species used when computing the molecular equilibrium and
C     also any additional species that appear only in the line list. Ideally,
C     the calling routine will <1> Include SIZES.SYN, <2> Use SPLSIZ to
C     dimension SPLIST and XFRACT, and <3> Pass SPLSIZ in place of SPLDIM.
C     However, SPLDIM is passed separately to allow for error checking in
C     the cases when this is not done (e.g. when called from IDL).
C   MODE [integer] Determines the content of the output:
C      1    - number densities
C      2    - partition functions
C      3    - partial pressures
C      0 or others number densities/partition functions
C   10+     - the same as above but electron density is assumed to be known
C             precisely and not re-determined in the process
C
C Input/Output:
C   SPLIST [character*(*) array(SPLDIM)] If NLIST is nonzero upon entry,
C     then SPLIST must contain the base set of species that must be included
C     in the molecular equilibrium calculation, regardless of which species
C     are represented by lines in SPNAME. Until the code is cleaned up, the
C     species list in SPLIST must include "e-" after the NLIST element.
C     If NLIST is zero upon entry, then SPLIST is loaded with the base set
C     of species coded into EQSTAT below (in the variable DEFAULT). Again,
C     an "e-" is appended after the base set.
C     Regardless of the whether SPLIST is valid upon entry or needs to be
C     loaded with the defaults, species that are in the lines list SPNAME,
C     but are not in the base set of species will be inserted into SPLIST
C     after the "e-" entry. Currently, the extended list is not used, but
C     in the future, we may solve for the equilibrium of all species in the
C     extended SPLIST.
C   NLIST [integer] If nonzero upon entry, NLIST is the number of species
C     in the base set of species passed in SPLIST (including the mandatory
C     "e-" at the beginning of the list). If NLIST is zero upon entry, this
C     indicates that the default base set of species coded in EQSTAT should
C     be used. Upon exit, NLIST is set to the number of species in SPLIST,
C     which contains the base set plus any additional species that occur
C     in the line list.
C
C Outputs:
C   SPINDX [integer array(NLINES)] Species index assigned to each line in
C     the input line list (specified by the input arguments SPNAME and ION).
C     The species index is used to reconstruct the species name (in SPLIST)
C     or "zeta" value (in XFRACT) computed for each line in the input line
C     list. For example, ZETA(SPINDX(370)) contains the zeta value for the
C     line corresponding to SPNAME(370) and ION(370).
C   Pg     [real] gas (no electrons) pressure.
C   XFRACT [real array(SPLDIM)] Zeta (in cm^-3) for the atomic or molecular
C     species in the corresponding entry of SPNAME and the charge state in
C     corresponding entry of ION. Zeta is the number density divided by the
C     partition function, and is required for spectrum synthesis.
C   POTI   [real array(SPLDIM)] ionization potential in eV for the
C     corresponding species.
C   ATWGHT [real array(SPLDIM-1)] molecular weights in AMU for the
C     corresponding species.
C   H1FRCT [real] Number density (in cm^-3) of neutral atomic hydgrogen,
C     used in computing damping constants (and continuous opacities?).
C   HE1FRCT [real] Number density (in cm^-3) of neutral atomic helium,
C     used in computing damping constants (and continuous opacities?).
C   XNA, XNE [real] Number density of gas species and free electrons as
C     compute by the EOS.
C   NITER  [integer] Number of iterations needed for the EOS.
C
      subroutine eqstat_rho(mode,temp,Pg,Pe,abund,elemen,amass,
     &                  ELESIZ,spindx,splist,xfract,poti,atwght,
     &                  nlines,nlist,xne,xna,rho,niter)
      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'

      integer mode,ELESIZ,niter
      integer nlines,nlist
      real temp,Tk,Pg,Pe,Pgas,Pelec,xna,xne,rho,xntot
c      real xnatom,xnelec,xne_old,xna_old
      real Pg_old,Pe_old,rho_new
      character*(SPCHAR) splist(nlist)
      character*(3) elemen(ELESIZ)
      integer spindx(nlines)
      real xfract(nlist),poti(nlist),atwght(nlist)
      real abund(ELESIZ),amass(ELESIZ)
      logical FAILED

      integer Anum(4),Natm(4),maxion,nelm,nchg
      real xnpf(SPLSIZ),pfunc(SPLSIZ),tol,tol1,xtotal
      real potion(IONSIZ),wtmol
      double precision awt(SPLSIZ-1),fract(IONSIZ)
      integer icharge,iter,ispec,IH1,IHe1,mmode

      INTEGER MAXITER
      REAL kBol
      DOUBLE PRECISION PSI,sum,amu
      PARAMETER (kBol=1.38065E-16,MAXITER=5000,amu=1.66053886d-24)
C
C Call equation of state solver.
C
c      open(87,file='dumpb.dat',form='unformatted',status='old')
c      read(87) temp,Pgas,Pelec,abund,elemen,amass,
c     &         mmode,spindx(nlines),splist,nlines,nlist
c      close(87)
      TOL=1.E-5
      TOL1=1.E-3
      Pelec=Pe
      PSI=2.d0/(1.d0+SQRT(5.d0))
      DO ISPEC=1,NLIST
        IF(SPLIST(ISPEC).EQ.'H  ') IH1 =ISPEC
        IF(SPLIST(ISPEC).EQ.'He ') IHE1=ISPEC
        XNPF(ISPEC)=-1.
      END DO
      Tk=temp*kBol
      mmode=mod(mode,10)
C
C================================================
C Hot gas: ignore molecules and solve ionization equilibrium only
C
      if(temp.gt.14000.) then
C
C Hot gas: assume no molecules and use Saha equation
C
C
C Compute gas pressure
C Mean molecular weight:
        sum=0.d0
        do ispec=1,ELESIZ
          sum=sum+abund(ispec)*amass(ispec)
        end do
        sum=sum*amu
C
C Number of atoms/ions and gas pressure:
        xntot=rho/sum
C
C Iterate to find gas/electron pressures consistent with the given density
C
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        niter=0
        Pgas = 2.0 * xntot * tk    
 1      niter=niter+1
        if(niter .gt. 200) stop
        Pg=Pgas
        
C
C Get number density of free electrons
C     
        call Nelect(temp,Pgas,abund,amass,ELESIZ,
     *       xna,xne,wtmol)
        
        if(mode.lt.10) then
           Pelec=xne*Tk
        else
           xne=Pelec/Tk
        endif
C     
C     If the total number of particles derived from the density and the Nelect
C     are significantly discrepant recompute Pgas and iterate
C

        if(abs((xntot-xna) / xntot) .gt. TOL) then
           Pgas = Pgas + (xntot-xna)*tk
           goto 1
        endif

c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c        niter=0
c        xna=xntot*0.5
c  1     niter=niter+1
c        Pgas=xna*Tk
c        Pg=Pgas
cC
cC Get number density of free electrons
cC
c        call Nelect(temp,Pgas,abund,amass,ELESIZ,
c     *              xna,xne,wtmol)
c        if(mode.ge.10) then
c          Pelec=xne*Tk
c        else
c          xne=Pelec/Tk
c        endif
cC
cC If the total number of particles derived from the density and the Nelect
cC are significantly discrepant scale xna and iterate
cC
c        if(abs(xna+xne-xntot)/(xna+xne).gt.TOL) then
c          xna=xna*xntot/(xna+xne)
c          go to 1
c        endif
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C We found consistent values of Pgas and Pelec. Proceed with the EOS.
C
        xna=(Pgas-Pelec)/Tk

        rho=xna*wtmol
        do 2 ispec=1,nlist-1
        CALL MPARSE(elemen,splist(ispec),Nelm,Nchg,Anum,Natm,ELESIZ)
        icharge=Nchg+1
        if(Nelm.eq.1.and.Natm(1).eq.1.and.Nchg.ge.0) then
C
C Get the number of ionization stages available in XSAHA
C
          call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,5)
C
C Atom. Parser returns atomic number in Anum(1)
C
          if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else if(mmode.eq.2) then
C
C  MODE=2, Return partition functions
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,3)
            xfract(ispec)=fract(icharge)
            poti(ispec)=potion(icharge)
          else if(mmode.eq.3) then
C
C  MODE=3, Return partial pressures
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*kBol*temp*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else
C
C  Any other MODE: Return number densities / partition functions
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,1)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          endif
          atwght(ispec)=amass(Anum(1))
        else
C
C Ignore molecules
C
          poti(ispec)  =1.
          atwght(ispec)=1.
          xfract(ispec)=0.
        endif
  2     continue
C
C Electrons
C
        if(mmode.eq.1) then
          xfract(nlist)=xne
        else if(mmode.eq.2) then
          xfract(nlist)=1.
        else if(mmode.eq.3) then
          xfract(nlist)=xne*Tk
        else
          xfract(nlist)=xne
        endif
      else
C
C================================================
C Cold gas: solve molecular and ionization equilibrium
C
C
C Compute mean molecular weight
C
        sum=0.d0
        DO ispec=1,ELESIZ
          sum=sum+abund(ispec)*amass(ispec)
        END DO
        sum=sum*amu
        wtmol=sum
C
C Gas pressure as if no molecules are present
C
        Pg_old=rho/sum
        niter=0
  3     continue
c        write(*,*) NLINES,NLIST,temp,Pgas,Pelec,mmode
c        write(*,'(10f8.3)') log10(abund)
        if(temp.gt.4000.) then
          Pe_old=Pg_old*0.1
        else if(temp.gt.2000.) then
          Pe_old=Pg_old*0.01
        else
          Pe_old=Pg_old*0.001
        endif
  4     continue
        if(temp.lt.1500.) then
          call lnGAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho_new,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
        else
          call GAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho_new,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
        endif
        niter=niter+iter
        IF(niter.ge.MAXITER) THEN
          Pelec=xne*Tk
c          WRITE(*,*) 'T,Pgas,Pnew,Pelec,Pe_in,Pe_out,NITER=',
c     *                Temp,Pgas,Pg,Pe,Pe_old,Pelec,niter,FAILED
          IF(niter.gt.MAXITER*20) STOP
        END IF
C
C Adjust pressure according to the discrepancy in density 
C
        IF(abs(Pgas -Pg_old)/max(1.E-20,Pgas ).gt.tol1.or.
     *     abs(Pelec-Pe_old)/max(1.E-20,Pelec).gt.tol1) THEN
          Pe_old=Pelec
          Pg_old=Pg
          GOTO 4
        END IF
C
C The convergence for a given value of rho is achieved.
C Iterate Pg to match the density
C
        if(abs(rho-rho_new)/rho.gt.tol) then
          Pe_old=xne*Tk*rho/rho_new
          Pg_old=Pgas*rho/rho_new
          go to 3
        endif
        Pg=Pgas
        Pe=xne*Tk
c        write(*,*) 'T, P', Temp, Pg
c        do ispec=1,nlist-1
c          write(*,*) ispec,splist(ispec),xnpf(ispec)
c        enddo
c      write(*,'(F10.1,13E11.4)') Temp,xnpf(1),
c     &                                xnpf(2),
c     &                                xnpf(3),
c     &                                xnpf(4),
c     &                                xnpf(5),
c     &                                xnpf(6),
c     & (Pgas-Pelec)/Tk,xna,Pelec/Tk,xne,rho
C
C Fill return arrays.
C
        do ispec=1,nlist-1
          atwght(ispec)=awt(ispec)
        end do
C
        if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
          do ispec=1,nlist-1
c            write(*,*) ispec,splist(ispec),xnpf(ispec),pfunc(ispec)
            xfract(ispec)=xnpf(ispec)
          end do
          xfract(nlist)=xne
        else if(mmode.eq.2) then
C
C  MODE=2, Return partition functions
C
          do ispec=1,nlist-1
            xfract(ispec)=pfunc(ispec)
          end do
          xfract(nlist)=1.
        else if(mmode.eq.3) then
C
C  MODE=3, Return partial pressures
C
          do ispec=1,nlist-1
            xfract(ispec)=xnpf(ispec)*Tk
          end do
          xfract(nlist)=xne*Tk
        else
C
C  Any other MODE: Return number densities / partition functions
C
          do ispec=1,nlist-1
            xfract(ispec)=xnpf(ispec)/pfunc(ispec)
          end do
          xfract(nlist)=xne
        endif
      endif
C
      return
      end

C=========================================================================
C LLENGTH: Returns an almost unique integer for molecule "name" which
C  is assumed to include up to 4 different types of atoms.
C  For molecule A1_n1 A2_n2 A3_n3 A4_n4 Ch
C  llength = (n1 + n2 + n3 + n4)*10000 + (Z1 + Z2 + Z3 + Z4)*10 + charge
C  Charge of -1 corresponds to 9. Positive charge is limited to +8.
C
      function llength(name,elemen,ELESIZ)
C
      integer iel(4),nat(4),charge,ELESIZ
      character*(*) name
      character*3 elemen(ELESIZ)
C
      call mparse(elemen,name,nel,charge,iel,nat,ELESIZ)
      llength=0
      do 1 i=1,nel
      llength=llength+iel(i)*10+10000*nat(i)
   1  continue
      if(charge.gt.0) then
        llength=llength+charge
      else if(charge.lt.0) then
        llength=llength+9
      end if
C
      return
      end

C=========================================================================
C NELECT: Finds consistent electron number density.
C
C Inputs:
C   T [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   P [real] Total gas pressure (in dyne/cm^2), given by NTOTAL*K*T,
C     which is to be used in calculating chemical and ionization equilibrium,
C     and partial pressures.
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   AMASS [real array(ELESIZ)] atomic weights in AMU.
C Outputs:
C   XNA    [real] Atomic number density
C   XNE    [real] Electron number density
C   H1FRC  [real] Number density (in cm^-3) of neutral atomic hydgrogen,
C      used in computing damping constants.
C   HE1FRC [real] Number density (in cm^-3) of neutral atomic helium,
C      used in computing damping constants.
C   WTMOLE [real] Mean molecular weight in AMU.
C
      SUBROUTINE NELECT(T,P,ABUND,AMASS,ELESIZ,
     *                  XNA,XNE,WTMOLE)
c     *                  XNA,XNE,H1FRC,HE1FRC,WTMOLE)
C
C
C  AUTHOR: N.Piskunov
C
C  LAST UPDATE: 29 January 1993
C
      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      INTEGER ELESIZ
      REAL T,P,XNE,XNA,WTMOLE
c      REAL T,P,XNE,XNA,H1FRC,HE1FRC,WTMOLE
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      DOUBLE PRECISION kBol,amu
      PARAMETER (kBol=1.38065D-16,amu=1.66053886D-24)

      DOUBLE PRECISION FRACT(IONSIZ)
      DOUBLE PRECISION TK,XNTOT,XNENEW,X,XA,XE,ERROR
      REAL POTI(IONSIZ)
      INTEGER L,IEL,ION,MAXION
C
      TK=kBol*T
      XNTOT=P/TK
      XE=XNTOT*0.5D0
      XA=XE
      DO 4 L=1,200
        XNENEW=0.D0
        DO 2 IEL=1,ELESIZ
          X=0.D0
          XNE=XE
          XNA=XA
C
C  Get the number of known ions
C
          CALL XSAHA(IEL,T,XNE,XNA,MAXION,POTI,FRACT,5)
C
C  Get the number of electrons contributed by all ions of atom IEL
C
          CALL XSAHA(IEL,T,XNE,XNA,MAXION,POTI,FRACT,2)
c          IF(IEL.EQ.1) H1FRC =FRACT(1)
c          IF(IEL.EQ.2) HE1FRC=FRACT(1)
          DO 1 ION=1,MIN(MAXION,IEL+1)
            X=X+FRACT(ION)*(ION-1)
   1      CONTINUE
          XNENEW=XNENEW+X*XA*ABUND(IEL)
   2    CONTINUE
        XNENEW=(XNENEW+XE)*0.5D0
        ERROR=ABS((XE-XNENEW)/XNENEW)
        XE=XNENEW
        XA=XNTOT-XE
c        write(*,'('' T,XNE,XNA,ERROR='',F8.1,3E14.6)') T,XNE,XNA,ERROR
        IF(ERROR.LT.1.D-5) THEN
          X=0.D0
          DO 3 IEL=1,99
            X=X+ABUND(IEL)*AMASS(IEL)
   3      CONTINUE
          WTMOLE=X*amu
c          WTMOLE=(X-XE*5.4857990943D-4)*amu
          RETURN
        END IF
   4  CONTINUE
      WRITE(*,*) 'Can''t converge calculating electron density'
C
      STOP
      END

C=========================================================================
C SORT2: sorts two arrays in atomic element order of the first (character) array.
C Hydrogen first, Helium next etc. All atoms/ions must end up before molecules
C that contain this atoms.
C
      subroutine sort2(nlist,list1,nlines,list2,elemen,ELESIZ)
      include 'SIZES.EOS'
c
      integer nlist,nlines,ELESIZ
      character*(*) list1(nlist)
      character*3 elemen(ELESIZ)
      character*(SPCHAR) name,name1,name2
      integer list2(nlines)
c
c Go through the list (except the last item which is e-)
c
      i=0
   1  if(i.lt.nlist-2) then
c
c Set the first entry as the minimum rank in the remaining part of the list
c
        i=i+1
        imin=i
        name2=list1(imin)
        l2=llength(name2,elemen,ELESIZ)
c
c Go through other entries. Look for smaller or identical ranks.
c
        j=i
   2    if(j.lt.nlist-1) then
          j=j+1
          name1=list1(j)
          l1=llength(name1,elemen,ELESIZ)
          if(l1.lt.l2.or.(l1.eq.l2.and.name1.lt.name2)) then
c
c Found smaller rank. Store the location of the new winner.
c
            imin=j
            name2=list1(imin)
            l2=llength(name2,elemen,ELESIZ)
c            if(list1(list2(4)).eq.'e-') write(*,*) 'A',name1,name2,
c     *      imin,list1(imin),(list2(k),k=1,nlines)
          else if(name1.eq.name2) then
c
c Found more than one candidate: kill the latter and update the index vector
c
            do 3 k=j,nlist-1
            list1(k)=list1(k+1)
   3        continue
            nlist=nlist-1
            if(nlines.gt.0) then
              do 4 k=1,nlines
              if(list2(k).eq.j) list2(k)=imin
              if(list2(k).gt.j) list2(k)=list2(k)-1
   4          continue
            endif
          end if
          go to 2
        end if
c
c Put entries in the correct order and update the index vector
c
        name=list1(i)
c        if(list1(list2(4)).eq.'e-') write(*,*) 'C',name,
c     *    list1(imin),imin,list1(imin),(list2(k),k=1,nlines)
        list1(i)=list1(imin)
        list1(imin)=name
        if(nlines.gt.0) then
          do 5 k=1,nlines
          l=list2(k)
          if(l.eq.i)    list2(k)=imin
          if(l.eq.imin) list2(k)=i
   5      continue
        endif
        go to 1
      end if
c
      return
      end

C=========================================================================
C MBUILD: Build complete name from charge value and neutral species name.
C
C Inputs:
C   SPNAME [character] Name of neutral atom or molecule,
C   ICHARGE [integer] Desired charge value (-1, 0, 1 - 4) for output
C   atomic or molecular species. The charge value is interpreted as follows:
C       -1:  negative ion
C        0:  neutral species
C       +1:  singly ionized species
C       +2:  doubly ionized species, etc.
C
C     All other charge values are invalid and generate fatal errors.
C
C Outputs:
C   CHNAME [character] Complete name of species constructed from input charge
C     value and neutral species name.
C
C 96-Jun-01 Valenti  Wrote.
C 96-Dec-12 Piskunov Expanded to IONSIZ ionization stage
C
      subroutine mbuild(spname,icharge,chname)
      INCLUDE 'SIZES.EOS'

      character*(*) spname,chname
C
C Generate a fatal error if the neutral species begins with a space.
C
      if(spname(1:1).eq.' ') then
        write(*,*) 'mbuild: species name is blank'
        stop
      endif
C
C Check that requested charge value is allowed.
C
      if(icharge.lt.-1 .or. icharge.gt.IONSIZ-1) then
        write(*,200) spname,icharge
 200    format('mbuild: invalid charge value for ',a,':',i4)
        stop
      endif
C
C Initialize the output string with spaces.
C
      chname=' '
C
C Handle the simple case where a neutral charge state was requested.
C Just copy the input neutral species name up to the first space or
C   until SPCHAR characters have been copied.
C
      if(icharge.eq.0) then
        chname=spname
        return
      endif
C
C Find location of the first space, which is where the charge string will go.
C A fatal error occurs if the output requires more than SPCHAR characters.
C
      ispace=index(spname,' ')
      if(ispace.le.0.or.ispace+abs(icharge)-1.gt.len(chname)) then
        write(*,201) spname,icharge
 201    format('mbuild: no room in string "',a,'" for charge:',i4)
        stop
      end if
C
C Copy neutral species name.
C
      chname=spname
C
C Insert charge string beginning at first space.
C
      if(icharge.lt.0) then
        chname(ispace:ispace)='-'
      else if(icharge.gt.0.and.icharge.lt.IONSIZ) then
        chname(ispace:ispace+icharge-1)='++++++++++++++++++++++++++++++'
      else
        write(*,*) 'The charge is too large. Must be less than',IONSIZ,
     *             spname,icharge
        stop
      endif
C
c      write(*,*) icharge,'"',chname,'"'
      return
      end

C=========================================================================
C MPARSE: Parse molecular name. Get number and type of atomic constituents.
C
C Inputs:
C   SPNAME [character array(*)] Case-sensitive species name of molecule.
C     First letter of each atom name must be uppercase. The second letter
C     of each atom name, if present, must be lowercase. Each atom name may
C     optionally be followed by a multiplicity number between 1 and 4. If
C     no multiplicity number is given for a particular atom, then its
C     multiplicity is assumed to be 1. Finally, a non-neutral charge state
C     for the molecule may be specified with a trailing "-", "+", or "++".
C     In the absence of such a charge indicator, the molecule is assumed
C     to be neutral.
C   ELEMEN [character array(*)] Case-sensitive list of atoms participating
C     in molecule formation (periodic table).
C
C Outputs:
C   NEL [integer] Number of elements comprising molecule. Also gives the
C     maximum valid index for IEL and NAT.
C   CHARGE [integer] Charge state of the molecule (-1, 0, +1,...,+(IONSIZ-1)).
C   IEL [integer array(4)] atomic number(s) of the atomic types comprising
C     the molecule in SPNAME.
C   NAT [integer array(4)] multiplicity (up to 4) for each of the atomic
C     types in IEL.
C
      SUBROUTINE MPARSE(ELEMEN,SPNAME,NEL,CHARGE,IEL,NAT,ELESIZ)
      INCLUDE 'SIZES.EOS'
C
      INTEGER IEL(4),NAT(4),NEL,CHARGE,ELESIZ
      CHARACTER SPNAME*(SPCHAR),TMP*2
      CHARACTER*(3) ELEMEN(ELESIZ)
C
C  Set pointer I1 to beginning of first atom name.
C
c      write(*,*) LEN(ELEMEN(1))
      CHARGE=0
      I1=1
C
C  Loop through (up to four) different atoms in a molecule.
C
      DO 4 J=1,4
C
C  Set pointer I2 to the end of the next atom's name.
C
      I2=I1
      IF(ICHAR(SPNAME(I1+1:I1+1)).GE.ICHAR('a').AND.
     *   ICHAR(SPNAME(I1+1:I1+1)).LE.ICHAR('z')) I2=I1+1
C
C  Update number of atomic species in molecule.
C
      NEL=J
C
C  Find atomic the atomic number of current atom.
C
      TMP='  '
      TMP=SPNAME(I1:I2)
      DO 1 I=1,ELESIZ
      IF(TMP.EQ.ELEMEN(I)(1:2)) GO TO 2
   1  CONTINUE
C
C  Fall through to here if atom name was not in ELEMEN list.
C
c      WRITE(*,*) 'Unknown element: ',SPNAME,i1,i2,' ',SPNAME(i1:i2)
      WRITE(*,*) 'Unknown element: ',SPNAME(I1:I2),' "',SPNAME(1:I2),'"'
      STOP
C
C  Save atomic number of current atom.
C
   2  IEL(NEL)=I
C
C  Check for optional atomic multiplicity. Default is 1; maximum is 5.
C
      I1=I2+1
      NAT(NEL)=1
      IF(SPNAME(I1:I1).EQ.'1') THEN
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'2') THEN
        NAT(NEL)=2
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'3') THEN
        NAT(NEL)=3
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'4') THEN
        NAT(NEL)=4
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'5') THEN
        NAT(NEL)=5
        I1=I1+1
      END IF
C
C   Check for optional charge on molecule. Default is neutral; "-", "+",
C   "++", etc. up to IONSIZ are allowed.
C
      IF(I1.GT.SPCHAR) RETURN
      IF(SPNAME(I1:I1).EQ.' ') RETURN
      IF(SPNAME(I1:I1).EQ.'-') THEN
        CHARGE=-1
        RETURN
      ENDIF
      IF(SPNAME(I1:I1).EQ.'+') THEN
        CHARGE=1
        DO 3 IONN=1,IONSIZ-1
        IF(SPNAME(I1+IONN:I1+IONN).NE.'+') RETURN
        CHARGE=CHARGE+1
   3    CONTINUE
      END IF
C
C  Fall through if we didn't just find a charge state and return. Loop
C  back and interpret character pointed at by I1 as beginning of atom.
C
   4  CONTINUE
C
C  There were 4 different atomic types, but presumably we are done.
C
      RETURN
      END

C=========================================================================
C EQPF: Returns partition functions interpolated for given thermodynamical
C       parameters. No equilibrium solving is apllied.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   XNELEC [real] Estimated electron number density (in 1/cm^3)
C   XNATOM [real] Number density (in 1/cm^3) of all particles other than
C     electrons (i.e. atoms or molecules), used to calculate total pressure?
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPLIST [character*(*) array(NLIST)] List of species to consider in
C     solving for the molecular equilibrium, including both the base set,
C     plus any additional species that occur in the line list.
C   NLIST [integer] the number of valid species in SPLIST.
C   TOL [real] iterative solution of the population equations terminates
C     successfully when the largest fractional change in pressure for all
C     species (including electrons) drops below this "tolerance".
C Outputs:
C   XNE [real] electron number density (in 1/cm^3) determined in GAS routine.
C   AWT [real*8] atomic weights of each species
C
      SUBROUTINE EQPF(TEMP,Pgas,Pelec,ABUND,ELEMEN,AMASS,ELESIZ,
     *               SPLIST,NLIST,PFUNC)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      CHARACTER ENAME*(SPCHAR),BLANK*1
      DOUBLE PRECISION KBOL,HMASS,AMULOG
      PARAMETER (BLANK=' ',ENAME='e-',KBOL=1.38065D-16,
     *           HMASS=1.66053D-24,AMULOG=-23.779751D0)

      INTEGER ELESIZ,NLIST
      CHARACTER*(SPCHAR) SPLIST(NLIST)
      CHARACTER*(3) ELEMEN(ELESIZ)
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      CHARACTER NAMEMX*(SPCHAR),NAMET*(SPCHAR)
      INTEGER JATOM, TYPE(SPLSIZ-1),NCH(SPLSIZ-1),IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL T,TEMP,XNELEC,XNATOM,TOL,XNE,XNA,RHO,Pgas,Pelec,Pgnew,
     *  POTI(IONSIZ),PFUNC(*),POTION(SPLSIZ),XTOTAL
      DOUBLE PRECISION IT(SPLSIZ-1),KT(SPLSIZ-1)
      DOUBLE PRECISION FRACT(IONSIZ), AWT(SPLSIZ-1)

      DOUBLE PRECISION PART(SPLSIZ-1)

      DOUBLE PRECISION PE,PG,PF,PNEW,PENEW,DP,DPE,PION,PENQ,PARTN
      INTEGER NELM,NCHG,ANUM(4),NATM(4)
      INTEGER I,J,K,NP,ISPEC,IELM
c      INTEGER IPIV(ELEDIM+1),IWORK(ELEDIM+1),
c     *  INFO,REPEAT,NSP1,NELT,NQ,K,KK,IDIR,KMAX,I,J,NEQ,
c     *  IIH2,IICO,IIH2O,NGIT
      DOUBLE PRECISION RATIOM,QPRD
c      DOUBLE PRECISION RHSTOT,SCALE,FACTOR,PNOTE,PDTOT,PU,
c     *  PD,GMU,PTOT,DELP,DELPE,PQ,RCOND,myDASUM,DELMAX,PE0,
c     *  PTOTH,PHyd,PTOTC,PTOTO,WATCOR,AQUAD,BQUAD,CQUAD,DPQ,DPTOT

      LOGICAL BARKLEM

C
C Total gas and electron pressure
C
      T=TEMP
      PG=Pgas
      PE=Pelec
      XNELEC=PE/(KBOL*TEMP)
      XNATOM=PG/(KBOL*TEMP)
C
C Avoid unpleasant surprises
C
      IF(PG.GT.PE) THEN
        XNATOM=XNATOM-XNELEC
      ELSE
        XNELEC=XNATOM*0.01
      END IF
c      PG=(XNATOM+XNELEC)*KBOL*TEMP
c      PE=XNELEC*KBOL*TEMP
C
C  Calculate equilibrium constants for each species in list (except 'e-').
C
      PION=0
      JATOM=0
      NP=0

      DO 4 ISPEC=1,NLIST-1
      CALL MPARSE(ELEMEN,SPLIST(ISPEC),NELM,NCHG,ANUM,NATM,ELESIZ)
      IF(NCHG.EQ.0) NP=ISPEC
      IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.EQ.0) THEN
C
C  Neutral atom
C
        TYPE(ISPEC)=1
        JATOM=JATOM+1
        IF(JATOM.GT.ELEDIM) THEN
          write(*,'(A,2I4,A8,I5)') 'gas: too many element types,' //
     *      ' increase ELEDIM:',ELEDIM,JATOM,SPLIST(ISPEC),ISPEC
          stop
        END IF
        IATOM(JATOM)=ANUM(1)
        INDSP(JATOM)=ISPEC
        IAT(ISPEC)=JATOM
        AWT(ISPEC)=AMASS(ANUM(1))
        INDZAT(ANUM(1))=JATOM
        NTOT(ISPEC)=1
        CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        PART(ISPEC)=FRACT(1)
      ELSE IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.NE.0) THEN
C
C  Ionized atom
C
        TYPE(ISPEC)=3
        IF(NCHG.GT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PART(ISPEC)=FRACT(NCHG+1)
        ELSE IF(NCHG.LT.0) THEN
C
C Negative ions
C
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PARTN=FRACT(1)
          CALL  NEGION(ANUM(1),TEMP,PARTN,IT(ISPEC),
     *                 PART(ISPEC),POTION(ISPEC),BARKLEM)
        END IF
C
      ELSE IF(NELM.GT.1.OR.NATM(1).GT.1) THEN
C
C  Neutral or ionized molecule
C
        TYPE(ISPEC)=2
C
C  Calculate mass ratio (RATIOM) and partition function product (QPRD)
C  needed by MOLCON. See MOLCON header for decription of these quantities.
C  While we are at it, calculate the atomic weight (AWT) of the molecule
C  and the total number of atoms (NTOT) of any type in the molecule.
C
        NTOT(ISPEC)=0
        AWT(ISPEC)=0.0D0
        RATIOM=0.0D0
C
C  Fixed the partition function ratio for ionized molecules.
C  Now we start with a product of partition functions for free
C  electrons in denominator. NP 29-12-2006.
        IF(NCHG.LE.0) THEN
          QPRD=0.0D0
        ELSE
          QPRD=-NCHG*LOG10(2.0)
        ENDIF
        DO 2 IELM=1,NELM
        NTOT(ISPEC)=NTOT(ISPEC)+NATM(IELM)
        AWT(ISPEC)=AWT(ISPEC)+NATM(IELM)*AMASS(ANUM(IELM))
        RATIOM=RATIOM+NATM(IELM)*LOG10(AMASS(ANUM(IELM)))
        CALL XSAHA(ANUM(IELM),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        QPRD=QPRD+NATM(IELM)*LOG10(FRACT(1))
   2    CONTINUE
        RATIOM=RATIOM-LOG10(AWT(ISPEC))+(NTOT(ISPEC)-1)*AMULOG
C
C  Now get the molecular constants from MOLCON.
C
        CALL MOLCON(SPLIST(ISPEC),TEMP,NTOT(ISPEC),RATIOM,QPRD,
     &              KT(ISPEC),PART(ISPEC),PION,BARKLEM)
C
C  Finally, record the charge state of the molecule.
C
        IF(NCHG.GT.0.AND.BARKLEM) THEN
C
C  The first option was used with Sauval & Tatum constants.
C  JV fits to NextGen pressures needed IT(ISPEC)=1.0 for positive
C  molecular ions.
C
C
C Positively charged molecules (single charge only!)
C
          K=1
          DO IELM=2,NELM
            IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
     *          K=IELM
          ENDDO
        ELSE IF(NCHG.LT.0) THEN
C
C Negatively charged molecules (single charge only!)
C Known negatively charged molecules are:
C H2-, CH-, C2-, CN-, OH-, SiH-, HS-
C
          IF(SPLIST(ISPEC).EQ.'H2-') THEN
            PARTN=PART(INDSP(INDZAT( 1)))
            CALL NEGION( 1,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CH-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'C2-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CN-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'OH-') THEN
            PARTN=PART(INDSP(INDZAT( 8)))
            CALL NEGION( 8,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'SiH-') THEN
            PARTN=PART(INDSP(INDZAT(14)))
            CALL NEGION(14,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'HS-') THEN
            PARTN=PART(INDSP(INDZAT(16)))
            CALL NEGION(16,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE
            IT(ISPEC)=1.D0
          ENDIF
        END IF
C
C  Store ionization potential (needed e.g. for broadening calculations)
C
        IF(PION.GT.0.D0) THEN
          POTION(ISPEC)=PION
        ELSE
c
c  If ionization potential is not available use the one for TiO!
c
          POTION(ISPEC)=6.4
        ENDIF
      ELSE
C
C  Fall through to here when the molecular formula doesn't make sense.
C
        WRITE(*,*) 'Wrong formula for the species: ',splist(ISPEC)
        STOP
      END IF
C
C  Now save results of MPARSE into arrays.
C
      NEL(ISPEC)=NELM
      NCH(ISPEC)=NCHG
      DO 3 IELM=1,NELM
      ZAT(IELM,ISPEC)=ANUM(IELM)
      NAT(IELM,ISPEC)=NATM(IELM)
   3  CONTINUE
C
C  Go back for next species.
C
   4  CONTINUE
C
C Fill up the output array and set up flags
C PNOTE is the partial pressure due to everything except electrons.
C XNA is the number density of everything except electrons.
C
      DO 5 ISPEC=1,NLIST-1
      IF(PART(ISPEC).GT.0.) THEN
        PFUNC(ISPEC)=PART(ISPEC)
      ELSE
        PFUNC(ISPEC)=1.
      END IF
   5  CONTINUE
      PFUNC(NLIST)=1.0
C
      RETURN
      END



C=========================================================================
C GAS: Determines the equilibrium ionization and partial pressure for every
C      atom and molecule in the species list, assuming no other species are
C      present. Temperature, total pressure, and elemental abundances must
C      be specified, but all atomic and molecular constants are determined
C      internally.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   XNELEC [real] Estimated electron number density (in 1/cm^3)
C   XNATOM [real] Number density (in 1/cm^3) of all particles other than
C     electrons (i.e. atoms or molecules), used to calculate total pressure?
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPLIST [character*(*) array(NLIST)] List of species to consider in
C     solving for the molecular equilibrium, including both the base set,
C     plus any additional species that occur in the line list.
C   NLIST [integer] the number of valid species in SPLIST.
C   TOL [real] iterative solution of the population equations terminates
C     successfully when the largest fractional change in pressure for all
C     species (including electrons) drops below this "tolerance".
C Outputs:
C   XNE [real] electron number density (in 1/cm^3) determined in GAS routine.
C   AWT [real*8] atomic weights of each species
C
      SUBROUTINE GAS(TEMP,Pgas,Pelec,ABUND,ELEMEN,AMASS,ELESIZ,
     *               TOL,SPLIST,NLIST,XNE,XNA,RHO,Pgnew,
     *               XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
     *               FAILED)
c      SUBROUTINE GAS(TEMP,XNELEC,XNATOM,ABUND,ELEMEN,AMASS,ELESIZ,
c     *               TOL,SPLIST,NLIST,
c     *               XNE,XNA,RHO,XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
c     *               FAILED)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      CHARACTER ENAME*(SPCHAR),BLANK*1
      INTEGER MAXIT,MAXREF
      DOUBLE PRECISION KBOL,HMASS,AMULOG
      PARAMETER (BLANK=' ',ENAME='e-',KBOL=1.38065D-16,MAXIT=1000,
     *           HMASS=1.66053D-24,AMULOG=-23.779751D0,MAXREF=10)
      LOGICAL PRINT,FAILED

      INTEGER NLIST,ELESIZ
      CHARACTER*(SPCHAR) SPLIST(NLIST)
      CHARACTER*(3) ELEMEN(ELESIZ)
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      CHARACTER NAMEMX*(SPCHAR),NAMET*(SPCHAR)
      INTEGER JATOM, TYPE(SPLSIZ-1),NCH(SPLSIZ-1),IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL T,TEMP,XNELEC,XNATOM,TOL,XNE,XNA,RHO,Pgas,Pelec,Pgnew,
     *  POTI(IONSIZ),XNPF(*),PFUNC(*),POTION(*),XTOTAL
      DOUBLE PRECISION FRACT(IONSIZ),IT(SPLSIZ-1),KT(SPLSIZ-1),
     *  AWT(SPLSIZ-1)

      DOUBLE PRECISION A(ELEDIM+1,ELEDIM+1),RHS(ELEDIM+1),
     *  AA(ELEDIM+1,ELEDIM+1),
     *  B(ELEDIM+1),BB(ELEDIM+1),
     *  P(ELEDIM+1),PP(SPLSIZ-1),PP0(SPLSIZ-1),PART(SPLSIZ-1),ND

      DOUBLE PRECISION PE,PG,PF,PNEW,PENEW,DP,DPE,PION,PENQ,PARTN
      DOUBLE PRECISION RNF(ELEDIM),AL(ELEDIM+1)
      INTEGER NELM,NCHG,ANUM(4),NATM(4),IPIV(ELEDIM+1),IWORK(ELEDIM+1),
     *  INFO,REPEAT,ISPEC,NSP1,NELT,NQ,K,KK,IDIR,KMAX,I,J,NEQ,IELM,NP,
     *  IIH2,IICO,IIH2O,NGIT
      DOUBLE PRECISION RATIOM,QPRD,RHSTOT,SCALE,FACTOR,PNOTE,PDTOT,PU,
     *  PD,GMU,PTOT,DELP,DELPE,PQ,RCOND,myDASUM,DELMAX,PE0,
     *  PTOTH,PHyd,PTOTC,PTOTO,WATCOR,AQUAD,BQUAD,CQUAD,DPQ,DPTOT
c      DOUBLE PRECISION PZS,COMPZ

      DOUBLE PRECISION RSCL(ELEDIM+1),CSCL(ELEDIM+1)
      DOUBLE PRECISION FERR(1),BERR(1),WORK(5*(ELEDIM+1))
      CHARACTER*1 EQUED
      LOGICAL BARKLEM
      INTEGER JDAMAX
      EXTERNAL JDAMAX,myDASUM,myDGESVX,xDCOPY

cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real ttt(101)
c      real*8 Kttt(101)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C Initialize the Reciprocal Neutral Fraction (RNF). The RNF is used to
C adjust the initial neutral atomic partial pressures used in the linear
C solver. Originally, atomic species were assumed to be predominantly
C neutral, but at low electron pressures, this is a poor assumption for
C species with low ionization potentials.
C
      DO 1 I=1,ELEDIM
      RNF(I)=1.0D0
   1  CONTINUE
C
C Total gas and electron pressure
C
c      T=MAX(1200.,TEMP)
      T=TEMP
      PG=Pgas
      PE=Pelec
      XNELEC=PE/(KBOL*TEMP)
      XNATOM=PG/(KBOL*TEMP)
C
C Avoid unpleasant surprises
C
      IF(PG.GT.PE) THEN
        XNATOM=XNATOM-XNELEC
      ELSE
        XNELEC=XNATOM*0.01
      END IF
c      PG=(XNATOM+XNELEC)*KBOL*TEMP
c      PE=XNELEC*KBOL*TEMP
C
C  Calculate equilibrium constants for each species in list (except 'e-').
C
c      PRINT=.TRUE.
      PRINT=.FALSE.
      PION=0
      IIH2=0
      IICO=0
      IIH2O=0
      JATOM=0
      NP=0
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      open(13,file='KT_eos.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')
c      write(13) NLIST,LEN(SPLIST(1))
c      write(*,*) 'NLIST=',NLIST,splist(17)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      do 4 ISPEC=17,17
      DO 4 ISPEC=1,NLIST-1
      PP0(ISPEC)=0.D0
      CALL MPARSE(ELEMEN,SPLIST(ISPEC),NELM,NCHG,ANUM,NATM,ELESIZ)
      IF(NCHG.EQ.0) NP=ISPEC
      IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.EQ.0) THEN
C
C  Neutral atom
C
        TYPE(ISPEC)=1
        KT(ISPEC)=1.0
        IT(ISPEC)=1.0
        JATOM=JATOM+1
        IF(JATOM.GT.ELEDIM) THEN
          write(*,'(A,2I4,A8,I5)') 'gas: too many element types,' //
     *      ' increase ELEDIM:',ELEDIM,JATOM,SPLIST(ISPEC),ISPEC
          stop
        END IF
        IATOM(JATOM)=ANUM(1)
        INDSP(JATOM)=ISPEC
        IAT(ISPEC)=JATOM
        AWT(ISPEC)=AMASS(ANUM(1))
        INDZAT(ANUM(1))=JATOM
        NTOT(ISPEC)=1
        CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        PART(ISPEC)=FRACT(1)
        POTION(ISPEC)=POTI(1)
      ELSE IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.NE.0) THEN
C
C  Ionized atom
C
        TYPE(ISPEC)=3
        IF(NCHG.GT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,2)
          IT(ISPEC)=FRACT(NCHG+1)/FRACT(1)*PE**NCHG
          RNF(ANUM(1))=RNF(ANUM(1))+FRACT(NCHG+1)/FRACT(1)
c          if(ANUM(1).eq.26) write(*,*) SPLIST(ISPEC),NCHG,
c     *                      (FRACT(I),I=1,IONSIZ)
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PART(ISPEC)=FRACT(NCHG+1)
c          if(ANUM(1).eq.62) write(*,*) 'pf: ',SPLIST(ISPEC),NCHG,FRACT
          POTION(ISPEC)=POTI(NCHG+1)
          KT(ISPEC)=1.0
        ELSE IF(NCHG.LT.0) THEN
C
C Negative ions
C
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PARTN=FRACT(1)
          CALL  NEGION(ANUM(1),TEMP,PARTN,IT(ISPEC),
     *                 PART(ISPEC),POTION(ISPEC),BARKLEM)
        END IF
C
        KT(ISPEC)=1.
        AWT(ISPEC)=AMASS(ANUM(1))
        NTOT(ISPEC)=1
      ELSE IF(NELM.GT.1.OR.NATM(1).GT.1) THEN
C
C  Neutral or ionized molecule
C
        TYPE(ISPEC)=2
C
C  Calculate mass ratio (RATIOM) and partition function product (QPRD)
C  needed by MOLCON. See MOLCON header for decription of these quantities.
C  While we are at it, calculate the atomic weight (AWT) of the molecule
C  and the total number of atoms (NTOT) of any type in the molecule.
C
        NTOT(ISPEC)=0
        AWT(ISPEC)=0.0D0
        RATIOM=0.0D0
C
C  Fixed the partition function ratio for ionized molecules.
C  Now we start with a product of partition functions for free
C  electrons in denominator. NP 29-12-2006.
        IF(NCHG.LE.0) THEN
          QPRD=0.0D0
        ELSE
          QPRD=-NCHG*LOG10(2.0)
        ENDIF
        DO 2 IELM=1,NELM
        NTOT(ISPEC)=NTOT(ISPEC)+NATM(IELM)
        AWT(ISPEC)=AWT(ISPEC)+NATM(IELM)*AMASS(ANUM(IELM))
        RATIOM=RATIOM+NATM(IELM)*LOG10(AMASS(ANUM(IELM)))
        CALL XSAHA(ANUM(IELM),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        IF(SPLIST(ISPEC).EQ.'H2')  IIH2=ISPEC
        IF(SPLIST(ISPEC).EQ.'CO')  IICO=ISPEC
        IF(SPLIST(ISPEC).EQ.'H2O') IIH2O=ISPEC
        QPRD=QPRD+NATM(IELM)*LOG10(FRACT(1))
   2    CONTINUE
        RATIOM=RATIOM-LOG10(AWT(ISPEC))+(NTOT(ISPEC)-1)*AMULOG
C
C  Now get the molecular constants from MOLCON.
C
        CALL MOLCON(SPLIST(ISPEC),TEMP,NTOT(ISPEC),RATIOM,QPRD,
     &              KT(ISPEC),PART(ISPEC),PION,BARKLEM)
c       if(SPLIST(ISPEC).eq.'TiO')write(*,*) TEMP,KT(ISPEC),PART(ISPEC)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        do ittt=0,100
c          ttt(ittt+1)=20.*ittt+1000.
c          CALL MOLCON(SPLIST(ISPEC),ttt(ittt+1),NTOT(ISPEC),
c     &                RATIOM,QPRD,Kttt(ittt+1),PART(ISPEC),PION)
c        enddo
c        write(13) SPLIST(ispec),ttt,Kttt
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  Finally, record the charge state of the molecule.
C
        IT(ISPEC)=1.D0
c        write(*,*) ISPEC,SPLIST(ISPEC)
        IF(NCHG.GT.0.AND.BARKLEM) THEN
C
C  The first option was used with Sauval & Tatum constants.
C  JV fits to NextGen pressures needed IT(ISPEC)=1.0 for positive
C  molecular ions.
C
c-----------------------------------------------------------------------
c          IF(SPLIST(ISPEC).EQ.'H2+'.OR.SPLIST(ISPEC).EQ.'NO+') THEN
c            K=1
c            DO IELM=2,NELM
c              IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
c     *          K=IELM
c            ENDDO
c            IT(ISPEC)=IT(INDSP(ANUM(K))+1)
c            KT(ISPEC)=KT(ISPEC)/IT(ISPEC)
c          ENDIF
c          IT(ISPEC)=1.0
c-----------------------------------------------------------------------
C
C Positively charged molecules (single charge only!)
C
          K=1
          DO IELM=2,NELM
            IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
     *          K=IELM
          ENDDO
          IT(ISPEC)=IT(INDSP(ANUM(K))+1)
        ELSE IF(NCHG.LT.0) THEN
C
C Negatively charged molecules (single charge only!)
C Known negatively charged molecules are:
C H2-, CH-, C2-, CN-, OH-, SiH-, HS-
C
          IF(SPLIST(ISPEC).EQ.'H2-') THEN
            PARTN=PART(INDSP(INDZAT( 1)))
            CALL NEGION( 1,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CH-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'C2-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CN-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'OH-') THEN
            PARTN=PART(INDSP(INDZAT( 8)))
            CALL NEGION( 8,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'SiH-') THEN
            PARTN=PART(INDSP(INDZAT(14)))
            CALL NEGION(14,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'HS-') THEN
            PARTN=PART(INDSP(INDZAT(16)))
            CALL NEGION(16,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE
            IT(ISPEC)=1.D0
          ENDIF
          IT(ISPEC)=1.D0
        END IF
C
C  Store ionization potential (needed e.g. for broadening calculations)
C
        IF(PION.GT.0.D0) THEN
          POTION(ISPEC)=PION
        ELSE
c
c  If ionization potential is not available use the one for TiO!
c
          POTION(ISPEC)=6.4
        ENDIF
      ELSE
C
C  Fall through to here when the molecular formula doesn't make sense.
C
        WRITE(*,*) 'Wrong formula for the species: ',splist(ISPEC)
        STOP
      END IF
C
C  Now save results of MPARSE into arrays.
C
      NEL(ISPEC)=NELM
      NCH(ISPEC)=NCHG
      DO 3 IELM=1,NELM
      ZAT(IELM,ISPEC)=ANUM(IELM)
      NAT(IELM,ISPEC)=NATM(IELM)
   3  CONTINUE
C
C  Go back for next species.
C
c      write(*,'(f10.2,I4,A12,4E15.4)') T,ISPEC,SPLIST(ISPEC),
c     *     PART(ISPEC),
c     *     KT(ISPEC),IT(ISPEC),KT(ISPEC)/MAX(IT(ISPEC),1.D-30)
   4  CONTINUE
c      write(*,*) 'GAS completed',TEMP,KBOL,Pgas,Pelec,NLIST
c      stop
c      return
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      close(13)
c      stop
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NEQ=JATOM+1
C==================================
C== End of species list parsing. ==
C==================================
C
C Print diagnostic: neutral fractions.
C
c     write(*,*) 'Reciprocal Neutral Fractions'
c     do 850 i=1,JATOM/7
c       write(*,860) (jeff(iatom(j)),j=7*i-6,7*i)
c850  continue
c860  format(1p7e10.3,a)
c     if(JATOM.gt.7*(JATOM/7)) write(*,860)
c    *  (jeff(iatom(j)),j=7*(JATOM/7)+1,JATOM)
c      do 52 i=1,nlist-1
c  52  write(*,'(I4,1P2E12.4,3I3,A6,0Pf8.2,8I4)')
c     *  i,IT(i),KT(i),NCH(i),NTOT(i),NEL(i),SPLIST(i),AWT(i),
c     *  (ZAT(j,i),NAT(j,i),j=1,NEL(i))
C================================================================
C== UPDATE MAIN ARRAYS                                         ==
C================================================================
c
c Make the initial estimate of the partial pressures for neutral atoms. These
c pressures are used as input to the linear solver. When only abundances are
c considered, the largest errors occur for low ionization elements, which can
c be highly ionized at low electron pressures. Thus, we apply a correction
c to recover the neutral fraction for each atom. The neutral fraction only
c corrects for losses into ionization states included in the species list.
c When the ionization correction is included, the largest error in the inital
c guess for carbon, which has unaccounted for losses into CO. Late in the
c convergence process, nitrogen becomes the dominant source of error.
c
      DO 5 J=1,JATOM
      P(J)=PG*ABUND(IATOM(J))/RNF(IATOM(J))
      ISPEC=INDSP(J)
      PP0(ISPEC)=P(J)
   5  CONTINUE
c
c Make an initial guess at the balance between H and H2.
c Assumes pressures of species other than H, H2, He, and Ne are negligible.
c Constraints:
c   KT(IIH2)*PP(IIH2)=P(1)**2           <-- chemical equilibrium
c   P(1)+2*PP(IIH2)=ABUND(1)*(PG-PE)    <-- H particle conservation
c
      IF(IIH2.GT.0) THEN
        PHyd=0.5*(-KT(IIH2)+SQRT(KT(IIH2)**2
     &        +4.0*KT(IIH2)*(PG-PE-P(2)-P(10))))
      ELSE
        PHyd=(PG-PE)*ABUND(1)
      ENDIF
c      IF(PHyd.GT.0.) P(1)=PHyd
c
c Make an initial guess at the balance between C, O, CO, and H2O.
c Constraints:
c   KT(IICO)*PP(IICO)=P(6)*P(8)         <-- chemical equilibrium
c   KT(IIH2O)*PP(IIH2O)=P(1)**2*P(8)    <-- chemical equilibrium
c   PTOTH=P(1)+2*PP(IIH2)       <-- defines density of H nuclei
c   PTOTC=P(6)+PP(IICO)                 <-- defines density of C nuclei
c   PTOTO=P(8)+PP(IICO)+PP(IIH2O)       <-- defines density of O nuclei
c   PTOTC=PTOTH*ABUND(6)/ABUND(1)       <-- abundance constraint
c   PTOTO=PTOTH*ABUND(8)/ABUND(1)       <-- abundance constraint
c
      PTOTH=P(1)
      IF(IIH2.GT.0) PTOTH=PTOTH+2.0*P(1)**2/KT(IIH2)
      PTOTC=PTOTH*ABUND(6)/ABUND(1)
      PTOTO=PTOTH*ABUND(8)/ABUND(1)
      IF(IIH2O.GT.0) THEN
        WATCOR=1.0+P(1)**2/KT(IIH2O)
        AQUAD=1.0/WATCOR
        IF(IICO.GT.0) THEN
          BQUAD=KT(IICO)+(PTOTO-PTOTC)/WATCOR
          CQUAD=-KT(IICO)*PTOTC
c          P(6)=(-BQUAD+SQRT(BQUAD**2-4.0*AQUAD*CQUAD))/(2.0*AQUAD)
c          P(8)=(P(6)+PTOTO-PTOTC)/WATCOR
        ELSE
c          P(6)=PTOTC
c          P(8)=PTOTO
        ENDIF
      ELSE
c        P(6)=PTOTC
c        P(8)=PTOTO
      ENDIF
c      IF(P(6).LE.0.) P(6)=PTOTC
c      IF(P(8).LE.0.) P(8)=PTOTO
      PE0=PE
      NAMEMX=BLANK
      DELMAX=0.0D0
c      COMPZ=0.0D0
c      PZS=0.0D0
c      write(*,*) SPLIST(1),P(1),SPLIST(IIH2),P(IIH2),
c     *           SPLIST(IIH2+1),P(IIH2+1),
c     *           SPLIST(IIH2+2),P(IIH2+2)
c      DO 6 J=1,JATOM
c      NN=INDSP(J)
c      IF(IPR(NN).NE.2) GOTO 3
c      NNP=INDX(3,ITAB(ZAT(1,NN)),1,1,1)
c      COMPZ=COMPZ+ABUND(IATOM(J))
c      IF(PE.EQ.0.0D0) PZS= PZS + P(J)
c      IF(PE.GT.0.0D0) PZS= PZS + (1.0D0+IT(NNP)/PE)*P(J)
c   6  CONTINUE
c      do J=1,JATOM
c        write(*,*) J,P(J),ABUND(IATOM(J)),SPLIST(INDSP(J))
c      enddo
c      write(*,*) JATOM+1,PE,'e-'
c      stop
C================================================================
C== MAIN LOOP: FILL LINEARIZED COEFFICIENT MATRIX AND RHS VECTOR,
C== AND SOLVE SYSTEM FOR PARTIAL PRESSURE CORRECTIONS.         ==
C== ISOLV=1: LINEARIZE ONLY THE PARTIAL PRESSURES OF THE NEUTRAL=
C== ATOMS FOR WHICH IPR(J)=1 (MAJOR SPECIES). THE ELECTRON     ==
C== PRESSURE PE IS ASSUMED TO BE GIVEN IN THIS CASE, AND SO IS ==
C== NOT INCLUDED IN THE LINEARIZATION. THIS IS NECESSARY SINCE ==
C== MOST OF THESE ELECTRONS (AT COOL TEMPS.) ORIGINATE FROM    ==
C== ELEMENTS NOT CONSIDERED IN THE LINEARIZATION. IN ORDER TO  ==
C== OBTAIN A GOOD VALUE FOR PE IN THE FIRST PLACE, IT IS       ==
C== NECESSARY TO CALL GAS WITH ISOLV=2.                        ==
C== ISOLV=2: THIS LINEARIZES THE PARTIAL PRESSURES OF THE NEUTRAL
C== ATOMS FOR WHICH IPR(J)=1 OR 2. THIS LIST OF ELEMENTS SHOULD==
C== INCLUDE ALL THE SIGNIFICANT CONTRIBUTORS TO THE TOTAL      ==
C== PRESSURE PG, AS WELL AS THE ELECTON PRESSURE PE. ANY ELEMENT=
C== (IPR(J)=3) NOT INCLUDED IS ASSUMED TO HAVE A NEGLIGIBLE    ==
C== EFFECT ON BOTH P AND PE.                                   ==
C== IN BOTH CASES, THE PARTIAL PRESSURES OF THE NEUTRAL ATOMS  ==
C== FOR ELEMENTS NOT INCLUDED IN THE LINEARIZATION ARE         ==
C== CALCULATED DIRECTLY FROM THE NOW DETERMINED PRESSURES OF   ==
C== THE LINEARIZED ELEMENTS.                                   ==
C================================================================
      NGIT=0
      RHSTOT=1.D99
C
C Top of loop in which linearized equations are solved recursively.
C
      REPEAT=0
      KMAX=1
   7  IF(NGIT.GE.MAXIT) THEN
        WRITE(*,208)
 208    FORMAT('*** ERROR: TOO MANY ITERATIONS IN ROUTINE "GAS"')
        WRITE(*,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),RHSTOT
        write(*,*) TEMP,PG,P(1),XNATOM,XNELEC
        STOP
      END IF
      NGIT=NGIT+1
      P(NEQ)=PE

      SCALE=10.D0
      IDIR=0
   9  CALL EOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
     *  IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)

c      write(*,*) 'Pe,SCALE,B(1),Pg=',PE,SCALE,B(1),PG,NGIT

      IF(B(1).GT.1.D2) THEN
        IF(IDIR.NE.-1) THEN
          SCALE=SQRT(SCALE)
          IDIR=-1
        ENDIF
C
C Neutral atomic pressures are too high. Scale them down until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)/SCALE
        ENDDO
        GOTO 9
      ELSE IF(B(1).LT.-1.D2) THEN
        IF(IDIR.NE.1) THEN
          SCALE=SQRT(SCALE)
          IDIR=1
        ENDIF
C
C Neutral atomic pressures are too low. Scale them up until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)*SCALE
        ENDDO
        GOTO 9
      ENDIF

      CALL EOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
     *  IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
C================================================================
C== NOW SOLVE THE LINEARIZED EQUATIONS (USING ROUTINE "LINEQ") ==
C================================================================
      IF(PRINT) THEN
        WRITE(*,200) NGIT
 200    FORMAT('LOG OF COEFFICIENT MATRIX AT ITERATION #',I5//)
        KK=MIN(30,NEQ-1)
        WRITE(*,201) (SPLIST(INDSP(K)),K=1,KK-1),'e-','RHS'
 201    FORMAT(4x,31(1x,a3,2x))
        DO 21 I=1,KK-1
        DO 20 J=1,KK-1
        AL(J)=LOG10(ABS(A(J,I))+1.0D-50)
  20    CONTINUE
        AL(KK)=LOG10(ABS(A(NEQ,I))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(I))+1.0D-50)
        NAMET=SPLIST(INDSP(I))
        WRITE(*,202) NAMET,(AL(J),J=1,KK+1)
  21    CONTINUE
        DO 22 J=1,KK-1
        AL(J)=LOG10(ABS(A(J,NEQ))+1.0D-50)
  22    CONTINUE
        AL(KK)=LOG10(ABS(A(NEQ,NEQ))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(NEQ))+1.0D-50)
        NAMET='e-'
        WRITE(*,202) NAMET,(AL(J),J=1,KK+1)
 202    FORMAT(A2,31F6.1)
        WRITE(*,'(/)')
c        stop
      END IF
C
C  Save a copy of the RHS for future step refinement
C
      DO 23 I=1,NEQ
      RHS(I)=B(I)
  23  CONTINUE
      RHSTOT=myDASUM(NEQ,RHS,1)
C
C  Solve linear system for corrections
C  In order not to solve for Pelect, one should use NEQ-1 as the first
C  argument. NEQ solves the whole system including electron pressure
C
c
c  Using LAPACK routine
c
c        open(unit=4,file='dump.bin',form='UNFORMATTED')
c        write(4) NEQ
c        write(4) ((A(i,j),i=1,NEQ),j=1,NEQ)
c        write(4) (B(i),i=1,NEQ)
      CALL myDGESVX('E','N',NEQ,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
     *            RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
     *            WORK,IWORK,INFO)
c       write(4) (BB(i),i=1,NEQ)
c       stop
      CALL xDCOPY(NEQ,BB,1,B,1)
c      DO I=1,NEQ
c        B(I)=BB(I)
c      END DO
c
c  The same thing using LINEQ2 or LINEQ and BLAS 2/3
c          open(unit=4,file='dump.bin',form='UNFORMATTED')
c          write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
c          close(4)
c      CALL LINEQ(NEQ,1,A,ELEDIM+1,IPIV,B,ELEDIM+1,INFO)
      IF(INFO.NE.0) THEN
        IF(REPEAT.LT.2) THEN
          DO J=1,NEQ-1
           P(J)=P(J)*0.999D0
          END DO
          REPEAT=REPEAT+1
          GO TO 7
        ELSE IF(REPEAT.LT.4) THEN
          DO J=1,NEQ-1
           P(J)=P(J)*1.001D0
          END DO
          REPEAT=REPEAT+1
          GO TO 7
        ELSE
          WRITE(*,*) 'EOS: LINEQ failed to solved for corrections to'
          WRITE(*,*) '     the partial pressures. Matrix is degenerate'
          WRITE(*,*) '     Temp=',TEMP,', Natom=',XNATOM,', Nelec=',
     *                XNELEC
          WRITE(*,*) '     INFO=',INFO,' Iter=',NGIT,' EQUED=',EQUED
cc          open(unit=4,file='dump.bin',form='UNFORMATTED')
cc          write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
cc          close(4)
cc          write(1) 0
cc          close(1)
c          STOP
          CALL myDGESVX('E','N',NEQ-1,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,
     *                  EQUED,RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,
     *                  FERR,BERR,WORK,IWORK,INFO)
          CALL xDCOPY(NEQ,BB,1,B,1)
c          DO J=1,NEQ
c            B(J)=BB(J)
c          END DO
          PTOT=0.D0
          DO J=1,NEQ-1
            PTOT=PTOT+P(J)
          END DO
          PE=MAX(PG-PTOT,1.D-20)
        END IF
      END IF
      REPEAT=0

c
C=================================================================
C== FINALLY, UPDATE THE PARTIAL PRESSURES FOR THE MAJOR SPECIES ==
C== BY ADDING THE PRESSURE CORRECTIONS OBTAINED FOR EACH ATOM   ==
C== FROM THE LINEARIZATION PROCEDURE.                           ==
C=================================================================
      DELMAX=-1.0D0
      KMAX=1
      DO 31 K=1,NEQ
      ISPEC=INDSP(K)
C
C Compute the maximum correction in order to computer the under-relaxation factor
C
      DP=B(K)
      DELP=ABS(DP/MAX(P(K),1.D-50))
      IF(DELP.GT.DELMAX) THEN
        DELMAX=DELP
      END IF
  31  CONTINUE
C
C  Under-relaxation factor
C
      FACTOR=0.2D0/(DELMAX+0.2D0)
c      FACTOR=1.D0
C
C Apply corrections
C
      DELMAX=-1.0D0
      KMAX=1
      DO 32 K=1,JATOM
      ISPEC=INDSP(K)
C
C  Restrict the correction to avoid getting negative pressures
C
      PNEW=P(K)-B(K)*FACTOR
      IF(PNEW.LT.0.D0) PNEW=MIN(MIN(P(K),ABS(PNEW)),PG)
c      IF(PNEW.LT.0.D0) PNEW=ABS(PNEW)
      DP=PNEW-P(K)
      IF(ABS(DP).GT.1.D-15) DP=DP*MIN(1.D0,0.4D0*P(K)/ABS(DP))
      P(K)=PNEW
      DELP=ABS(DP/MAX(P(K),1.D-50))
      IF(DELP.GT.DELMAX) THEN
        NAMEMX=SPLIST(ISPEC)
        DELMAX=DELP
        KMAX=K
      END IF
  32  CONTINUE

c      PENEW=BBB(NEQ)
      PENEW=PE-B(NEQ)*FACTOR
c      write(*,*) NEQ,PE,PENEW,B(NEQ),NGIT
      IF(PENEW.LT.0.D0) PENEW=MIN(PE,ABS(PENEW))
c      IF(PENEW.LT.0.D0) PENEW=ABS(PENEW)
      DPE=PENEW-PE
      IF(ABS(DPE).GT.1.D-15) DPE=DPE*MIN(1.D0,0.4D0*PE/ABS(DPE))
      PE=PENEW
      IF(ABS(PE/PG).GE.1.0D-15) THEN
        DELPE=ABS(DPE/PE)
        IF(DELPE.GT.DELMAX) NAMEMX=ENAME
        IF(DELPE.GT.DELMAX) DELMAX=DELPE
      END IF
C================================================================
C== PRINT OUT SUMMARY LINE FOR EACH ITERATION                  ==
C================================================================
      PTOT=PE
      PQ=0.0D0
c      write(*,*) 0,'e-',PE,PTOT,PG,NGIT
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=LOG(MAX(IT(ISPEC),1.D-115))-LOG(KT(ISPEC))-
     -     LOG(MAX(PE,1.D-115))*NQ
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+LOG(MAX(P(J),1.D-115))*NAT(I,ISPEC)
        ENDDO
c        PENQ=1.0D0
c        IF(PE.GT.0.0D0.AND.NQ.NE.0) PENQ=PE**NQ
c        PP(ISPEC)=IT(ISPEC)/(KT(ISPEC)*PENQ)*PF
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG
      ENDDO
c      stop
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PE-PQ)/PG
c      write(*,*) PG,PTOT,DELMAX,DPTOT,DPQ,FACTOR
      IF(PRINT) THEN
        WRITE(*,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),
     *               PTOT/TEMP/KBOL,DPTOT,PE/TEMP/KBOL,DPQ
 203    FORMAT(I10,2X,A8,1P9E11.3)
      END IF
      IF((DPTOT.GT.TOL.OR.DPQ.GT.TOL.OR.DELMAX.GT.TOL)
     *   .AND.NGIT.LT.MAXIT) GOTO 7
C
C Bottom of the loop in which linearized equations are solved recursively.
C
C================================================================
C== CALCULATE FINAL PARTIAL PRESSURES AFTER CONVERGENCE OBTAINED=
C================================================================
      PTOT=PE
      PD=0.0D0
      PU=0.0D0
      PU=PE*0.000548597D0
      PQ=0.0D0
      DO 34 ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=1.0D0
        DO 33 I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF*P(J)**NAT(I,ISPEC)
  33    CONTINUE
        PENQ=1.0D0
        IF(PE.GT.0.0D0) PENQ=PE**NQ
        PP(ISPEC)=IT(ISPEC)/(KT(ISPEC)*PENQ)*PF
        PTOT=PTOT+PP(ISPEC)
        PD=PD+NTOT(ISPEC)*PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
        PU=PU+AWT(ISPEC)*PP(ISPEC)
  34  CONTINUE
      PP(NLIST)=PE
      PDTOT=PD+PE
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PQ-PE)/PG
      GMU=PU/PTOT
      ND=PTOT/(TEMP*KBOL)
      RHO=ND*GMU*HMASS
      XNE=PE/(TEMP*KBOL)
C================================================================
C== WRITE OUT FINAL PARTIAL PRESSURES                          ==
C================================================================
      IF(PRINT) THEN
c      IF(myDASUM(NLIST-1,PP,1)+PE.GT.PG*1.01D0) THEN
        write(*,'(''AFTER '',I3,'' iterations.   Max change of:'',G10.3,
     #      ''  in element:'',A)') NGIT,DELMAX,NAMEMX
        WRITE(*,'(''AFTER '',I3,'' ITERATIONS WITH ''/
     #            ''T='',1PE10.3,''   P='',E10.3)') NGIT,TEMP,
     #    myDASUM(NLIST-1,PP,1)+PE
        WRITE(*,'(''PDTOT='',1PE10.3,''   DPTOT='',E10.3,
     #            ''  DPQ='',E10.3,''  Nelectron='',E10.3,'' cm^3''/
     #    '' Nparticle='',1PE10.3,'' cm^3   Mean At.Wt.='',
     #    0PF7.3,''   Density='',1PE10.3,'' g/cm^3''//
     #    '' # Species   Abundance   Initial P   Final P'',
     #    ''      IT         KT         pf''/)')
     #    PDTOT,DPTOT,DPQ,XNE,ND-XNE,GMU,RHO
        NSP1=NLIST
        DO 35 ISPEC=1,NLIST-1
        IF(TYPE(ISPEC).NE.1) THEN
          WRITE(*,206) ISPEC,SPLIST(ISPEC),PP0(ISPEC),PP(ISPEC),
     #                 IT(ISPEC),KT(ISPEC),PART(ISPEC)
 206      FORMAT(I3,1X,A8,11X,1P5E11.3)
        ELSE
          J=IAT(ISPEC)
          WRITE(*,207) ISPEC,splist(ISPEC),ABUND(IATOM(J)),PP0(ISPEC),
     #                 PP(ISPEC),IT(ISPEC),KT(ISPEC),PART(ISPEC)
 207      FORMAT(I3,1X,A8,1P6E11.3)
        END IF
  35    CONTINUE
        WRITE(*,206) NSP1,ENAME,PE0,PE
        WRITE(*,*) JDAMAX(NLIST-1,PP,1),SPLIST(JDAMAX(NLIST-1,PP,1))
c        stop
      END IF
C
C Fill up the output array and set up flags
C PNOTE is the partial pressure due to everything except electrons.
C XNA is the number density of everything except electrons.
C
      PNOTE=0.D0
      DO 36 ISPEC=1,NLIST-1
      IF(PART(ISPEC).GT.0.) THEN
        IF(PP(ISPEC)/KBOL/TEMP.GE.1.D-20) THEN
c          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP*PART(ISPEC))
          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP)
        ELSE
          XNPF(ISPEC)=0.0
        END IF
        PFUNC(ISPEC)=PART(ISPEC)
      ELSE
        XNPF(ISPEC)=0.
        PFUNC(ISPEC)=1.
      END IF
      PNOTE=PNOTE+PP(ISPEC)
c      write(*,*) ISPEC,PNOTE,PP(ISPEC),SPLIST(ISPEC)
c      write(*,*) ISPEC,SPLIST(ISPEC),PFUNC(ISPEC)
  36  CONTINUE
c      write(*,*) 'e-',XNE
c      stop
      XNPF(NLIST)=XNE
      PFUNC(NLIST)=1.0
      XTOTAL=PD/(KBOL*TEMP)
      XNA=PNOTE/(KBOL*TEMP)
      Pgnew=PTOT
C
      RETURN
      END

C=========================================================================
C LOGARITHMIC version: the solution is found for the logs of ficticious
C                      partial pressures.
C GAS: Determines the equilibrium ionization and partial pressure for every
C      atom and molecule in the species list, assuming no other species are
C      present. Temperature, total pressure, and elemental abundances must
C      be specified, but all atomic and molecular constants are determined
C      internally.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   XNELEC [real] Estimated electron number density (in 1/cm^3)
C   XNATOM [real] Number density (in 1/cm^3) of all particles other than
C     electrons (i.e. atoms or molecules), used to calculate total pressure?
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPLIST [character*(*) array(NLIST)] List of species to consider in
C     solving for the molecular equilibrium, including both the base set,
C     plus any additional species that occur in the line list.
C   NLIST [integer] the number of valid species in SPLIST.
C   TOL [real] iterative solution of the population equations terminates
C     successfully when the largest fractional change in pressure for all
C     species (including electrons) drops below this "tolerance".
C Outputs:
C   XNE [real] electron number density (in 1/cm^3) determined in GAS routine.
C   AWT [real*8] atomic weights of each species
C
      SUBROUTINE lnGAS(TEMP,Pgas,Pelec,ABUND,ELEMEN,AMASS,ELESIZ,
     *                 TOL,SPLIST,NLIST,XNE,XNA,RHO,Pgnew,
     *                 XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
     *                 FAILED)
c      SUBROUTINE lnGAS(TEMP,XNELEC,XNATOM,ABUND,ELEMEN,AMASS,ELESIZ,
c     *                 TOL,SPLIST,NLIST,
c     *                 XNE,XNA,RHO,XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
c     *                 FAILED)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      CHARACTER ENAME*(SPCHAR),BLANK*1
      INTEGER MAXIT,MAXREF
      DOUBLE PRECISION KBOL,HMASS,AMULOG
      PARAMETER (BLANK=' ',ENAME='e-',KBOL=1.38065D-16,MAXIT=10000,
     *           HMASS=1.66053D-24,AMULOG=-23.779751D0,MAXREF=10)

      LOGICAL PRINT,FAILED

      INTEGER NLIST,ELESIZ
      CHARACTER*(SPCHAR) SPLIST(NLIST)
      CHARACTER*(3) ELEMEN(ELESIZ)
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      CHARACTER NAMEMX*(SPCHAR),NAMET*(SPCHAR)
      INTEGER JATOM, TYPE(SPLSIZ-1),NCH(SPLSIZ-1),IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL T,TEMP,XNELEC,XNATOM,TOL,XNE,XNA,RHO,Pgas,Pelec,Pgnew,
     *  POTI(IONSIZ),XNPF(*),PFUNC(*),POTION(*),XTOTAL
      DOUBLE PRECISION FRACT(IONSIZ),IT(SPLSIZ-1),KT(SPLSIZ-1),
     *  AWT(SPLSIZ-1)

      DOUBLE PRECISION A(ELEDIM+1,ELEDIM+1),RHS(ELEDIM+1),
     *  AA(ELEDIM+1,ELEDIM+1),
     *  B(ELEDIM+1),BB(ELEDIM+1),
     *  P(ELEDIM+1),PP(SPLSIZ-1),PP0(SPLSIZ-1),PART(SPLSIZ-1),ND

      DOUBLE PRECISION PE,PG,PF,PNEW,PENEW,DP,DPE,PION,PARTN
c      DOUBLE PRECISION AT,BT,PN,DPF(4),CRATIO,BBB(ELEDIM+1),
c     *  PENQ,DPP,DPPE
      DOUBLE PRECISION RNF(ELEDIM),AL(ELEDIM+1)
      INTEGER NELM,NCHG,ANUM(4),NATM(4),IPIV(ELEDIM+1),IWORK(ELEDIM+1),
     *  INFO,ISPEC,NSP1,NELT,NQ,K,KK,IDIR,KMAX,I,J,NEQ,IELM,NP,
     *  IIH2,IICO,IIH2O,NGIT,REPEAT
      DOUBLE PRECISION RATIOM,QPRD,RHSTOT,SCALE,FACTOR,PNOTE,PDTOT,PU,
     *  PD,GMU,PTOT,DELP,DELPE,PQ,RCOND,myDASUM,DELMAX,PE0,PTOTH,
     *  PHyd,PTOTC,PTOTO,WATCOR,AQUAD,BQUAD,CQUAD,DPQ,DPTOT,RENORM
c      DOUBLE PRECISION DUMMY,SCOLD,RHS0,RHS1,RHS2

c      DOUBLE PRECISION BOLD(ELEDIM+1),S(ELEDIM+1),GAMMA,BNORM,BOLDN
      DOUBLE PRECISION RSCL(ELEDIM+1),CSCL(ELEDIM+1)
c      DOUBLE PRECISION ROWCND,COLCND,AMX
      DOUBLE PRECISION FERR(1),BERR(1),WORK(5*(ELEDIM+1))
      CHARACTER*1 EQUED
      LOGICAL BARKLEM
      EXTERNAL myDASUM

      INTEGER NFIELDS
      PARAMETER (NFIELDS=40)
      CHARACTER*(*) FORMAT201,FORMAT202
c      CHARACTER*(*) AFIELDS
c      PARAMETER (AFIELDS=CHAR(NFIELDS/10+ICHAR('0'))//
c     *                   CHAR(MOD(NFIELDS,10)+ICHAR('0')))
c      PARAMETER (FORMAT201='(4x,'//AFIELDS//'(1X,A3,2X))')
c      PARAMETER (FORMAT202='(A2,'//AFIELDS//'F6.1)')
      PARAMETER (FORMAT201='(4x,48(1X,A3,2X))')
      PARAMETER (FORMAT202='(A2,48F6.1)')

cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real ttt(101)
c      real*8 Kttt(101)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C Initialize the Reciprocal Neutral Fraction (RNF). The RNF is used to
C adjust the initial neutral atomic partial pressures used in the linear
C solver. Originally, atomic species were assumed to be predominantly
C neutral, but at low electron pressures, this is a poor assumption for
C species with low ionization potentials.
C
      DO 1 I=1,ELEDIM
      RNF(I)=1.0D0
   1  CONTINUE
C
C Total gas and electron pressure
C
c      T=MAX(1200.,TEMP)
      T=TEMP
      PG=Pgas
      PE=Pelec
      XNELEC=PE/(KBOL*TEMP)
      XNATOM=PG/(KBOL*TEMP)
C
C Avoid unpleasant surprises
C
      if(PG.GT.PE) THEN
        XNATOM=XNATOM-XNELEC
      ELSE
        XNELEC=XNATOM*0.01
      END IF
c      PG=(XNATOM+XNELEC)*KBOL*TEMP
c      PE=XNELEC*KBOL*TEMP
C
C  Calculate equilibrium constants for each species in list (except 'e-').
C
c      PRINT=.TRUE.
      PRINT=.FALSE.
      PION=0
      IIH2=0
      IICO=0
      IIH2O=0
      JATOM=0
      NP=0
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      open(13,file='KT_eos.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')
c      write(13) NLIST,LEN(SPLIST(1))
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO 4 ISPEC=1,NLIST-1
      PP0(ISPEC)=0.D0
      CALL MPARSE(ELEMEN,SPLIST(ISPEC),NELM,NCHG,ANUM,NATM,ELESIZ)
c      write(*,*) ISPEC,'"'//SPLIST(ISPEC)//'"',NELM,NCHG,
c     *           ANUM,NATM,ELESIZ
      IF(NCHG.EQ.0) NP=ISPEC
      IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.EQ.0) THEN
C
C  Neutral atom
C
        TYPE(ISPEC)=1
        KT(ISPEC)=1.0
        IT(ISPEC)=1.0
        JATOM=JATOM+1
        IF(JATOM.GT.ELEDIM) THEN
          write(*,'(a,2i4)') 'gas: too many element types,' //
     *      ' increase ELEDIM:',ELEDIM,JATOM
          stop
        END IF
        IATOM(JATOM)=ANUM(1)
        INDSP(JATOM)=ISPEC
        IAT(ISPEC)=JATOM
        AWT(ISPEC)=AMASS(ANUM(1))
        INDZAT(ANUM(1))=JATOM
        NTOT(ISPEC)=1
        CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        PART(ISPEC)=FRACT(1)
        POTION(ISPEC)=POTI(1)
      ELSE IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.NE.0) THEN
C
C  Ionized atom
C
        TYPE(ISPEC)=3
        IF(NCHG.GT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,2)
          IT(ISPEC)=FRACT(NCHG+1)/FRACT(1)*PE**NCHG
          RNF(ANUM(1))=RNF(ANUM(1))+FRACT(NCHG+1)/FRACT(1)
c          if(ANUM(1).eq.26) write(*,*) SPLIST(ISPEC),NCHG,
c     *                      (FRACT(I),I=1,IONSIZ)
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PART(ISPEC)=FRACT(NCHG+1)
c          if(ANUM(1).eq.62) write(*,*) 'pf: ',SPLIST(ISPEC),NCHG,FRACT
          POTION(ISPEC)=POTI(NCHG+1)
          KT(ISPEC)=1.0
        ELSE IF(NCHG.LT.0) THEN
C
C Negative ions
C
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PARTN=FRACT(1)
          CALL  NEGION(ANUM(1),TEMP,PARTN,IT(ISPEC),
     *                 PART(ISPEC),POTION(ISPEC),BARKLEM)
        END IF
C
        KT(ISPEC)=1.D0
        AWT(ISPEC)=AMASS(ANUM(1))
        NTOT(ISPEC)=1
      ELSE IF(NELM.GT.1.OR.NATM(1).GT.1) THEN
C
C  Neutral or ionized molecule
C
        TYPE(ISPEC)=2
C
C  Calculate mass ratio (RATIOM) and partition function product (QPRD)
C  needed by MOLCON. See MOLCON header for decription of these quantities.
C  While we are at it, calculate the atomic weight (AWT) of the molecule
C  and the total number of atoms (NTOT) of any type in the molecule.
C
        NTOT(ISPEC)=0
        AWT(ISPEC)=0.0D0
        RATIOM=0.0D0
C
C  Fixed the partition function ratio for ionized molecules.
C  Now we start with a product of partition functions for free
C  electrons in denominator. NP 29-12-2006.
        IF(NCHG.LE.0) THEN
          QPRD=0.0D0
        ELSE
          QPRD=-NCHG*LOG10(2.0)
        ENDIF
        DO 2 IELM=1,NELM
        NTOT(ISPEC)=NTOT(ISPEC)+NATM(IELM)
        AWT(ISPEC)=AWT(ISPEC)+NATM(IELM)*AMASS(ANUM(IELM))
        RATIOM=RATIOM+NATM(IELM)*LOG10(AMASS(ANUM(IELM)))
        CALL XSAHA(ANUM(IELM),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        IF(SPLIST(ISPEC).EQ.'H2')  IIH2=ISPEC
        IF(SPLIST(ISPEC).EQ.'CO')  IICO=ISPEC
        IF(SPLIST(ISPEC).EQ.'H2O') IIH2O=ISPEC
        QPRD=QPRD+NATM(IELM)*LOG10(FRACT(1))
   2    CONTINUE
        RATIOM=RATIOM-LOG10(AWT(ISPEC))+(NTOT(ISPEC)-1)*AMULOG
C
C  Now get the molecular constants from MOLCON.
C
        CALL MOLCON(SPLIST(ISPEC),TEMP,NTOT(ISPEC),RATIOM,QPRD,
     *              KT(ISPEC),PART(ISPEC),PION,BARKLEM)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        do ittt=0,100
c          ttt(ittt+1)=20.*ittt+1000.
c          CALL MOLCON(SPLIST(ISPEC),ttt(ittt+1),NTOT(ISPEC),
c     *                RATIOM,QPRD,Kttt(ittt+1),PART(ISPEC),PION)
c        END DO
c        write(13) SPLIST(ispec),ttt,Kttt
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  Finally, record the charge state of the molecule.
C
        IT(ISPEC)=1.D0
        IF(NCHG.GT.0.AND.BARKLEM) THEN
C
C  The first option was used with Sauval & Tatum constants.
C  JV fits to NextGen pressures needed IT(ISPEC)=1.0 for positive
C  molecular ions.
C
c-----------------------------------------------------------------------
c          IF(SPLIST(ISPEC).EQ.'H2+'.OR.SPLIST(ISPEC).EQ.'NO+') THEN
c            K=1
c            DO IELM=2,NELM
c              IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
c     *          K=IELM
c            ENDDO
c            IT(ISPEC)=IT(INDSP(ANUM(K))+1)
c            KT(ISPEC)=KT(ISPEC)/IT(ISPEC)
c          ENDIF
c          IT(ISPEC)=1.0
c-----------------------------------------------------------------------
C
C Positively charged molecules (single charge only!)
C
          K=1
          DO IELM=2,NELM
            IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
     *          K=IELM
          ENDDO
          IT(ISPEC)=IT(INDSP(ANUM(K))+1)
        ELSE IF(NCHG.LT.0) THEN
C
C Negatively charged molecules (single charge only!)
C Known negatively charged molecules are:
C H2-, CH-, C2-, CN-, OH-, SiH-, HS-
C
          IF(SPLIST(ISPEC).EQ.'H2-') THEN
            PARTN=PART(INDSP(INDZAT( 1)))
            CALL NEGION( 1,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CH-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'C2-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CN-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'OH-') THEN
            PARTN=PART(INDSP(INDZAT( 8)))
            CALL NEGION( 8,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'SiH-') THEN
            PARTN=PART(INDSP(INDZAT(14)))
            CALL NEGION(14,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'HS-') THEN
            PARTN=PART(INDSP(INDZAT(16)))
            CALL NEGION(16,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE
            IT(ISPEC)=1.D0
          ENDIF
          IT(ISPEC)=1.D0
        END IF
C
C  Store ionization potential (needed e.g. for broadening calculations)
C
        IF(PION.GT.0.D0) THEN
          POTION(ISPEC)=PION
        ELSE
c
c  If ionization potential is not available use the one for TiO!
c
          POTION(ISPEC)=6.4
        ENDIF
      ELSE
C
C  Fall through to here when the molecular formula doesn't make sense.
C
        WRITE(*,*) 'Wrong formula for the species: ',splist(ISPEC)
        STOP
      END IF
C
C  Now save results of MPARSE into arrays.
C
      NEL(ISPEC)=NELM
      NCH(ISPEC)=NCHG
      DO 3 IELM=1,NELM
      ZAT(IELM,ISPEC)=ANUM(IELM)
      NAT(IELM,ISPEC)=NATM(IELM)
   3  CONTINUE
C
C  Go back for next species.
C
c      write(*,*) ISPEC,SPLIST(ISPEC),IT(ISPEC),KT(ISPEC)
c      IT(ISPEC)=MIN(MAX(1.D-250,IT(ISPEC)),1.D250)
c      KT(ISPEC)=MIN(MAX(1.D-250,KT(ISPEC)),1.D250)
c      write(*,'(f10.2,I4,A12,4E13.4)') TEMP,ISPEC,SPLIST(ISPEC),
c     *     PART(ISPEC),KT(ISPEC),IT(ISPEC)
c     *     ,KT(ISPEC)/MAX(IT(ISPEC),1.D-150)
   4  CONTINUE
c      RENORM=LOG(SQRT(myDASUM(NLIST-1,KT,1)))
c      write(*,*) RENORM
c      DO ISPEC=1,NLIST-1
c        KT(ISPEC)=LOG(KT(ISPEC))+RENORM*NTOT(ISPEC)
c      END DO
      
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      close(13)
c      stop
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      NEQ=JATOM+1
C==================================
C== End of species list parsing. ==
C==================================
C
C Print diagnostic: neutral fractions.
C
c     write(*,*) 'Reciprocal Neutral Fractions'
c     do 850 i=1,JATOM/7
c       write(*,860) (jeff(iatom(j)),j=7*i-6,7*i)
c850  continue
c860  format(1p,7e10.3,a)
c     if(JATOM.gt.7*(JATOM/7)) write(*,860)
c    *  (jeff(iatom(j)),j=7*(JATOM/7)+1,JATOM)
c      do 52 i=1,nlist-1
c  52  write(*,'(I4,1P2E12.4,3I3,A6,0Pf8.2,8I4)')
c     *  i,IT(i),KT(i),NCH(i),NTOT(i),NEL(i),SPLIST(i),AWT(i),
c     *  (ZAT(j,i),NAT(j,i),j=1,NEL(i))
C================================================================
C== UPDATE MAIN ARRAYS                                         ==
C================================================================
c
c Make the initial estimate of the partial pressures for neutral atoms. These
c pressures are used as input to the linear solver. When only abundances are
c considered, the largest errors occur for low ionization elements, which can
c be highly ionized at low electron pressures. Thus, we apply a correction
c to recover the neutral fraction for each atom. The neutral fraction only
c corrects for losses into ionization states included in the species list.
c When the ionization correction is included, the largest error in the inital
c guess for carbon, which has unaccounted for losses into CO. Late in the
c convergence process, nitrogen becomes the dominant source of error.
c
      DO 5 J=1,JATOM
      P(J)=PG*ABUND(IATOM(J))/RNF(IATOM(J))
      ISPEC=INDSP(J)
      PP0(ISPEC)=P(J)
   5  CONTINUE
c
c Make an initial guess at the balance between H and H2.
c Assumes pressures of species other than H, H2, He, and Ne are negligible.
c Constraints:
c   KT(IIH2)*PP(IIH2)=P(1)**2           <-- chemical equilibrium
c   P(1)+2*PP(IIH2)=ABUND(1)*(PG-PE)    <-- H particle conservation
c
      IF(IIH2.GT.0) THEN
        PHyd=0.5*(-KT(IIH2)+SQRT(KT(IIH2)**2
     *        +4.0*KT(IIH2)*(PG-PE-P(2)-P(10))))
      ELSE
        PHyd=(PG-PE)*ABUND(1)
      END IF
c      IF(PHyd.GT.0.0.AND.PHyd.LT.Pgas-Pelec) P(1)=PHyd
c
c Make an initial guess at the balance between C, O, CO, and H2O.
c Constraints:
c   KT(IICO)*PP(IICO)=P(6)*P(8)         <-- chemical equilibrium
c   KT(IIH2O)*PP(IIH2O)=P(1)**2*P(8)    <-- chemical equilibrium
c   PTOTH=P(1)+2*PP(IIH2)               <-- defines density of H nuclei
c   PTOTC=P(6)+PP(IICO)                 <-- defines density of C nuclei
c   PTOTO=P(8)+PP(IICO)+PP(IIH2O)       <-- defines density of O nuclei
c   PTOTC=PTOTH*ABUND(6)/ABUND(1)       <-- abundance constraint
c   PTOTO=PTOTH*ABUND(8)/ABUND(1)       <-- abundance constraint
c
      PTOTH=P(1)
      IF(IIH2.GT.0) PTOTH=PTOTH+2.0*P(1)**2/KT(IIH2)
      PTOTC=PTOTH*ABUND(6)/ABUND(1)
      PTOTO=PTOTH*ABUND(8)/ABUND(1)
      IF(IIH2O.GT.0) THEN
        WATCOR=1.0+P(1)**2/KT(IIH2O)
        AQUAD=1.0/WATCOR
        IF(IICO.GT.0) THEN
          BQUAD=KT(IICO)+(PTOTO-PTOTC)/WATCOR
          CQUAD=-KT(IICO)*PTOTC
c          P(6)=(-BQUAD+SQRT(BQUAD**2-4.0*AQUAD*CQUAD))/(2.0*AQUAD)
c          P(8)=(P(6)+PTOTO-PTOTC)/WATCOR
        ELSE
c          P(6)=PTOTC
c          P(8)=PTOTO
        END IF
      ELSE
c        P(6)=PTOTC
c        P(8)=PTOTO
      END IF
c      IF(P(6).LE.0.0.OR.P(6).GT.0.1*P(1)) P(6)=PTOTC
c      IF(P(8).LE.0.0.OR.P(8).GT.0.1*P(1)) P(8)=PTOTO
      PE0=PE
      NAMEMX=BLANK
      DELMAX=0.0D0
c      COMPZ=0.0D0
c      PZS=0.0D0
c      DO 6 J=1,JATOM
c      NN=INDSP(J)
c      IF(IPR(NN).NE.2) GOTO 3
c      NNP=INDX(3,ITAB(ZAT(1,NN)),1,1,1)
c      COMPZ=COMPZ+ABUND(IATOM(J))
c      IF(PE.EQ.0.0D0) PZS= PZS + P(J)
c      IF(PE.GT.0.0D0) PZS= PZS + (1.0D0+IT(NNP)/PE)*P(J)
c   6  CONTINUE
c      do J=1,JATOM
c        write(*,*) J,P(J),ABUND(IATOM(J)),SPLIST(INDSP(J))
c      END DO
c      write(*,*) JATOM+1,PE,'e-'
c      stop
C================================================================
C== MAIN LOOP: FILL LINEARIZED COEFFICIENT MATRIX AND RHS VECTOR,
C== AND SOLVE SYSTEM FOR PARTIAL PRESSURE CORRECTIONS.         ==
C== ISOLV=1: LINEARIZE ONLY THE PARTIAL PRESSURES OF THE NEUTRAL=
C== ATOMS FOR WHICH IPR(J)=1 (MAJOR SPECIES). THE ELECTRON     ==
C== PRESSURE PE IS ASSUMED TO BE GIVEN IN THIS CASE, AND SO IS ==
C== NOT INCLUDED IN THE LINEARIZATION. THIS IS NECESSARY SINCE ==
C== MOST OF THESE ELECTRONS (AT COOL TEMPS.) ORIGINATE FROM    ==
C== ELEMENTS NOT CONSIDERED IN THE LINEARIZATION. IN ORDER TO  ==
C== OBTAIN A GOOD VALUE FOR PE IN THE FIRST PLACE, IT IS       ==
C== NECESSARY TO CALL GAS WITH ISOLV=2.                        ==
C== ISOLV=2: THIS LINEARIZES THE PARTIAL PRESSURES OF THE NEUTRAL
C== ATOMS FOR WHICH IPR(J)=1 OR 2. THIS LIST OF ELEMENTS SHOULD==
C== INCLUDE ALL THE SIGNIFICANT CONTRIBUTORS TO THE TOTAL      ==
C== PRESSURE PG, AS WELL AS THE ELECTON PRESSURE PE. ANY ELEMENT=
C== (IPR(J)=3) NOT INCLUDED IS ASSUMED TO HAVE A NEGLIGIBLE    ==
C== EFFECT ON BOTH P AND PE.                                   ==
C== IN BOTH CASES, THE PARTIAL PRESSURES OF THE NEUTRAL ATOMS  ==
C== FOR ELEMENTS NOT INCLUDED IN THE LINEARIZATION ARE         ==
C== CALCULATED DIRECTLY FROM THE NOW DETERMINED PRESSURES OF   ==
C== THE LINEARIZED ELEMENTS.                                   ==
C================================================================
      FACTOR=1.D0
      NGIT=0
      RHSTOT=1.D99
c      goto 2222
C
C Top of loop in which linearized equations are solved recursively.
C
      KMAX=1
c      PG=PG+myDASUM(NEQ-1,P)*(RENORM-1)
      DO J=1,NEQ-1
c        P(J)=LOG(P(J))+RENORM
        P(J)=LOG(P(J))
      END DO
      PE=LOG(MAX(PE,1.D-150))
c      open(unit=4,file='dump.bin',form='UNFORMATTED')
c      write(4) NEQ
      REPEAT=0
   7  IF(NGIT.GE.MAXIT) THEN
        WRITE(*,208)
 208    FORMAT('*** ERROR: TOO MANY ITERATIONS IN ROUTINE "GAS"')
        WRITE(*,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),RHSTOT
        write(*,*) TEMP,PG,P(1),XNATOM,XNELEC
        STOP
      END IF
      NGIT=NGIT+1
      P(NEQ)=PE

c      do J=1,NEQ
c        p(J)=exp(p(j))
c      enddo
c      write(*,*) (P(J),J=1,NEQ)
c      CALL lnEOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
c     *     IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c      CALL lnEOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
c     *     IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c      do j=1,NEQ
c        SCALE=P(J)
c        P(J)=P(J)+0.1d0
c        CALL lnEOSFCN(NEQ,P,BB,A,1,PG,NCH,NLIST,
c     *    IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c        write(*,*) J,SCALE
c        write(*,'(40e10.3)')(a(i,j)-(bb(i)-b(i))/0.1d0
c     *                            ,i=1,40)
c        write(*,'(40e10.3)')(a(i,j),i=1,40)
c        write(*,'(40e10.3)')((bb(i)-b(i))/0.1d0,i=1,40)
c        write(*,'(40e10.3)')(bb(i),i=1,40)
c        P(J)=SCALE
c      enddo
c      stop

      SCALE=10.D0
      IDIR=0
c      do j=1,NEQ
c        write(*,*) J,P(J),PG
c      enddo
c      write(*,*) B(1),PG
   9  CALL lnEOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
     *              IATOM,INDSP,NAT,ZAT,NTOT,
     *              NEL,IAT,INDZAT,ABUND,KT,IT)
c      write(*,*) SCALE,B(1),PG
      IF(B(1).GT.0.001D0*PG) THEN
        IF(IDIR.NE.-1) THEN
          SCALE=SQRT(SCALE)
          IDIR=-1
        END IF
C
C Neutral atomic pressures are too high. Scale them down until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)-LOG(SCALE)
        END DO
        GOTO 9
      ELSE IF(B(1).LT.-0.001D0*PG) THEN
        IF(IDIR.NE.1) THEN
          SCALE=SQRT(SCALE)
          IDIR=1
        END IF
C
C Neutral atomic pressures are too low. Scale them up until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)+LOG(SCALE)
        END DO
        GOTO 9
      END IF

c      IF(B(1).GT.0.02D0*PG) THEN
c        IF(IDIR.NE.1) THEN
c          SCALE=SQRT(SCALE)
c          IDIR=1
c        END IF
cC
cC Neutral atomic pressures are too high. Scale them down until
cC partical conservation equation will become negative
cC
c        DO ISPEC=1,NLIST-1
c          J=0
c          DO I=1,NEL(ISPEC)
c            J=J+NAT(I,ISPEC)
c          END DO
c          write(*,*) ISPEC,SPLIST(ISPEC),J,NCH(ISPEC)
c          KT(ISPEC)=KT(ISPEC)*SCALE**J
c          IT(ISPEC)=IT(ISPEC)*SCALE**NCH(ISPEC)
c        END DO
c        GOTO 9
c      ELSE IF(B(1).LT.-0.02D0*PG) THEN
c        IF(IDIR.NE.-1) THEN
c          SCALE=SQRT(SCALE)
c          IDIR=-1
c        END IF
cC
cC Neutral atomic pressures are too low. Scale them up until
cC partical conservation equation will become negative
cC
c        DO ISPEC=1,NLIST-1
c          J=0
c          DO I=1,NEL(ISPEC)
c            J=J+NAT(I,ISPEC)
c          END DO
c          KT(ISPEC)=KT(ISPEC)/SCALE**J
c          IT(ISPEC)=IT(ISPEC)/SCALE**NCH(ISPEC)
c        END DO
c        GOTO 9
c      END IF

c      do j=1,NEQ
c        write(*,*) J,P(J),PG
c      enddo
c      write(*,*) B(1),PG
      CALL lnEOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
     *              IATOM,INDSP,NAT,ZAT,NTOT,
     *              NEL,IAT,INDZAT,ABUND,KT,IT)
c      DO I=1,NEQ-1
c      WRITE(*,FORMAT202) SPLIST(INDSP(I)),(A(I,J),J=1,NEQ-1),B(I)
c      END DO
c      stop
C
C================================================================
C== NOW SOLVE THE LINEARIZED EQUATIONS (USING ROUTINE "LINEQ") ==
C================================================================
      IF(PRINT) THEN
        WRITE(*,200) NGIT
 200    FORMAT('LOG OF COEFFICIENT MATRIX AT ITERATION #',I5/)
        KK=MIN(NFIELDS,NEQ-1)
        WRITE(*,FORMAT201) (SPLIST(INDSP(K)),K=1,KK-1),'e-','RHS'
        DO 21 I=1,KK-1
        DO 20 J=1,KK-1
        AL(J)=LOG10(ABS(A(J,I))+1.0D-50)
  20    CONTINUE
        AL(KK)=LOG10(ABS(A(NEQ,I))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(I))+1.0D-50)
        NAMET=SPLIST(INDSP(I))
        WRITE(*,FORMAT202) NAMET,(AL(J),J=1,KK+1)
  21    CONTINUE
        DO 22 J=1,KK-1
        AL(J)=LOG10(ABS(A(J,NEQ))+1.0D-50)
  22    CONTINUE
        AL(KK)=LOG10(ABS(A(NEQ,NEQ))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(NEQ))+1.0D-50)
        NAMET='e-'
        WRITE(*,FORMAT202) NAMET,(AL(J),J=1,KK+1)
        WRITE(*,'(/)')
      END IF
c      stop
C
C  Save a copy of the RHS for future step refinement
C
      DO 23 I=1,NEQ
      RHS(I)=B(I)
  23  CONTINUE
      RHSTOT=myDASUM(NEQ,RHS,1)
C
C  Solve linear system for corrections
C  In order not to solve for Pelect, one should use NEQ-1 as the first
C  argument. NEQ solves the whole system including electron pressure
C
c
c  Using LAPACK routine
c
c        open(unit=4,file='dump.bin',form='UNFORMATTED')
c        write(4) NEQ
c        write(4) ((A(i,j),i=1,NEQ),j=1,NEQ)
c        write(4) (B(i),i=1,NEQ)
c      write(4) ((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
      CALL myDGESVX('E','N',NEQ,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
     *            RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
     *            WORK,IWORK,INFO)
c      stop
      CALL xDCOPY(NEQ,BB,1,B,1)
c      DO I=1,NEQ
c        B(I)=BB(I)
c      ENDDO
c      write(4) ((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
c
c  The same thing using LINEQ2 or LINEQ and BLAS 2/3
c      CALL LINEQ(NEQ,1,A,ELEDIM+1,IPIV,B,ELEDIM+1,INFO)
      IF(INFO.NE.0) THEN
        IF(REPEAT.LT.2) THEN
          DO J=1,NEQ-1
           P(J)=P(J)-0.01D0
          END DO
          REPEAT=REPEAT+1
          GO TO 7
        ELSE IF(REPEAT.LT.4) THEN
          DO J=1,NEQ-1
           P(J)=P(J)+0.01D0
          END DO
          REPEAT=REPEAT+1
          GO TO 7
        ELSE
          WRITE(*,*) 'lnGAS: DGESVX failed to solved for corrections to'
          WRITE(*,*) '  the partial pressures. Matrix is degenerate'
          WRITE(*,*) '  Temp=',TEMP,', Natom=',XNATOM,', Nelec=',XNELEC
          IF(INFO.EQ.NEQ) THEN
            WRITE(*,*) '  Pg=',PG,', INFO=',INFO,
     *                 ', Element: e-',
     *                 ', Iter=',NGIT,' EQUED=',EQUED
          ELSE
            WRITE(*,*) '  Pg=',PG,', INFO=',INFO,
     *                 ', Element: ',SPLIST(INDSP(INFO)),
     *                 ', Iter=',NGIT,' EQUED=',EQUED
          END IF
          CALL lnEOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,IATOM,INDSP,
     *                  NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
          open(unit=4,file='dump.bin',form='UNFORMATTED')
          write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
          close(4)
          WRITE(*,*) '  Matrix and the RHS were dumped to file dump.bin'
          STOP
c          CALL myDGESVX('E','N',NEQ-1,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
c     *                RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
c     *                WORK,IWORK,INFO)
c          CALL xDCOPY(NEQ-1,BB,1,B,1)
cc          DO I=1,NEQ
cc            B(I)=BB(I)
cc          END DO
c          PTOT=0.D0
c          DO J=1,NEQ-1
c            PTOT=PTOT+exp(P(J)-B(J))
c          END DO
c          PE=MAX(PG-PTOT,1.D-20)
c          Pe=log(Pe)
        END IF
      END IF
      REPEAT=0
c      IF(INFO.NE.0) THEN
c        WRITE(*,*) 'lnEOS: LINEQ failed to solved for corrections to'
c        WRITE(*,*) '     the partial pressures. Matrix is degenerate'
c        WRITE(*,*) '     Temp=',TEMP,', Natom=',XNATOM,', Nelec=',XNELEC
c        WRITE(*,*) '     Pg=',PG,', INFO=',INFO,
c     *             ', Element: ',SPLIST(INDSP(INFO)),
c     *             ', Iter=',NGIT,' EQUED=',EQUED
cc        open(unit=4,file='dump.bin',form='UNFORMATTED')
cc        write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
cc        close(4)
cc        write(1) 0
cc        close(1)
c        IF(PRINT) THEN
cc          close(4)
c          STOP
c        END IF
cc        DO J=1,NEQ
cc          P(J)=MAX(P(J)+0.1D0,-115.d0)
cc          write(*,*) J,P(J),B(J),B(J)*FACTOR
cc        END DO
c        write(*,*) P(INFO),B(INFO),B(INFO)*FACTOR
c        P(INFO)=MAX(P(INFO)+0.1D0,-115.d0)
c        PRINT=.TRUE.
c        GO TO 9
c      END IF
c
C=================================================================
C== FINALLY, UPDATE THE PARTIAL PRESSURES FOR THE MAJOR SPECIES ==
C== BY ADDING THE PRESSURE CORRECTIONS OBTAINED FOR EACH ATOM   ==
C== FROM THE LINEARIZATION PROCEDURE.                           ==
C=================================================================
      DELMAX=-200.0D0
      KMAX=1
      DO K=1,JATOM
c        write(*,*) K,P(K),B(K)
        ISPEC=INDSP(K)
c        DP=ABS(P(K))
        DELP=ABS(B(K))
c        IF(DP.GT.1.D-10) DELP=DELP/DP
        IF(DELP.GT.DELMAX) THEN
          NAMEMX=SPLIST(ISPEC)
          DELMAX=DELP
          KMAX=K
        END IF
      END DO
c      DPE=ABS(P(NEQ))
      DELPE=ABS(B(NEQ))
c      IF(DPE.GT.1.D-10) DELPE=DELPE/DPE
      IF(DELPE.GT.DELMAX) THEN
        NAMEMX=ENAME
        DELMAX=DELPE
        KMAX=NEQ
      END IF
c      write(*,*) KMAX,EXP(P(KMAX)),EXP(B(KMAX)),P(KMAX),B(KMAX)
C
C  Under-relaxation factor
C
      FACTOR=0.2D0/(DELMAX+0.2D0)
      DO K=1,JATOM
C
C  Apply corrections
C
        DP=B(K)*FACTOR
c        DP=10.D0*DP/MAX(10.D0,ABS(DP))
        PNEW=P(K)-DP
        P(K)=MAX(PNEW,-115.D0)
      END DO
      DP=B(NEQ)*FACTOR
c      DP=10.D0*DP/MAX(10.D0,ABS(DP))
      PENEW=PE-DP
      PE=MAX(PENEW,-115.D0)
C================================================================
C== PRINT OUT SUMMARY LINE FOR EACH ITERATION                  ==
C================================================================
      PTOT=EXP(PE)
      PQ=0.0D0
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=-PE*NQ+LOG(IT(ISPEC))-LOG(KT(ISPEC))
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+P(J)*NAT(I,ISPEC)
        END DO
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG,NQ,PQ,EXP(PE)
      END DO
c      stop
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(EXP(PE)-PQ)/PG
c      write(*,*) DELMAX,DPTOT,DPQ
      IF(PRINT) THEN
        WRITE(*,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),
     *               PTOT/TEMP/KBOL,DPTOT,EXP(PE)/TEMP/KBOL,DPQ,FACTOR
 203    FORMAT(I10,2X,A8,1P,9E11.3)
      END IF
c      write(*,*) NGIT,TOL,DPTOT,DELMAX,PTOT,PG
      IF((RHSTOT.GT.TOL.OR.DPTOT.GT.TOL.OR.DELMAX.GT.TOL)
     *   .AND.NGIT.LT.MAXIT) GO TO 7
C
C Bottom of the loop in which linearized equations are solved recursively.
C
C================================================================
C== CALCULATE FINAL PARTIAL PRESSURES AFTER CONVERGENCE OBTAINED=
C================================================================
c      write(*,*) RHSTOT,DELMAX,DPTOT,DPQ,TOL
      PTOT=EXP(PE)
      PD=0.0D0
      PU=0.0D0
      PQ=0.0D0
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=-PE*NQ+LOG(IT(ISPEC))-LOG(KT(ISPEC))
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+P(J)*NAT(I,ISPEC)
        END DO
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PD=PD+NTOT(ISPEC)*PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
        PU=PU+AWT(ISPEC)*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG,NQ,PQ,EXP(PE)
      END DO
      PE=EXP(PE)
      DO J=1,JATOM
        P(J)=EXP(P(J))
      END DO
      PP(NLIST)=PE
      PDTOT=PD+PE
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PQ-PE)/PG
      GMU=PU/PTOT
      ND=PTOT/(TEMP*KBOL)
      RHO=ND*GMU*HMASS
      XNE=PE/(TEMP*KBOL)
C================================================================
C== WRITE OUT FINAL PARTIAL PRESSURES                          ==
C================================================================
      IF(PRINT) THEN
        write(*,'(''AFTER '',I3,'' iterations.   Max change of:'',G10.3,
     #      ''  in element:'',A)') NGIT,DELMAX,NAMEMX
        WRITE(*,'(''AFTER '',I3,'' ITERATIONS WITH ''/
     #            ''T='',1PE10.3,''   P='',E10.3)') NGIT,TEMP,PG
        WRITE(*,'(''PDTOT='',1PE10.3,''   DPTOT='',E10.3,
     #            ''  DPQ='',E10.3,''  Nelectron='',E10.3,'' cm^3''/
     #    '' Nparticle='',1PE10.3,'' cm^3   Mean At.Wt.='',
     #    0PF7.3,''   Density='',1PE10.3,'' g/cm^3''/
     #    '' # Species   Abundance   Initial P   Final P'',
     #    ''      IT         KT         pf''//)')
     #   PDTOT,DPTOT,DPQ,XNE,ND-XNE,GMU,RHO
        NSP1=NLIST
        DO 35 ISPEC=1,NLIST-1
        IF(TYPE(ISPEC).NE.1) THEN
          WRITE(*,206) ISPEC,SPLIST(ISPEC),PP0(ISPEC),PP(ISPEC),
     #                 IT(ISPEC),KT(ISPEC),PART(ISPEC)
 206      FORMAT(I3,1X,A8,11X,1P,5E11.3)
        ELSE
          J=IAT(ISPEC)
          WRITE(*,207) ISPEC,splist(ISPEC),ABUND(IATOM(J)),PP0(ISPEC),
     #                 PP(ISPEC),IT(ISPEC),KT(ISPEC),PART(ISPEC)
 207      FORMAT(I3,1X,A8,1P,6E11.3)
        END IF
  35    CONTINUE
        WRITE(*,206) NSP1,ENAME,PE0,EXP(PE)
      END IF
C
C Fill up the output array and set up flags
C PNOTE is the partial pressure due to everything except electrons.
C XNA is the number density of everything except electrons.
C
2222  continue
      PNOTE=0.0
      DO 36 ISPEC=1,NLIST-1
      IF(PART(ISPEC).GT.0.) THEN
        IF(PP(ISPEC)/KBOL/TEMP.GE.1.D-20) THEN
c          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP*PART(ISPEC))
          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP)
        ELSE
          XNPF(ISPEC)=0.0
        END IF
        PFUNC(ISPEC)=PART(ISPEC)
      ELSE
        XNPF(ISPEC)=0.
        PFUNC(ISPEC)=1.
      END IF
      PNOTE=PNOTE+PP(ISPEC)
c      write(*,'(I4,2E12.4,2X,A)') ISPEC,PNOTE,PP(ISPEC),SPLIST(ISPEC)
c      write(*,*) ISPEC,SPLIST(ISPEC),PFUNC(ISPEC)
  36  CONTINUE
      XNPF(NLIST)=XNE
      PFUNC(NLIST)=1.0
      XTOTAL=PD/(TEMP*KBOL)
      XNA=PNOTE/(TEMP*KBOL)
c      write(*,*) 'Pg,PD,PNOTE,PE,PNOTE+PE',Pg,PD,PTOT,PE,PNOTE+PE
      Pgnew=Ptot
C
      RETURN
      END


C=========================================================================
C MOLCON: Returns equilibrium constant and partition function for a given
C   molecule and temperature.
C
C Inputs:
C   SPNAME [character(*)] Name of molecule, chosen from SPLIST below.
C   T [real] Temperature (in K) at which EQK and PART are to be found.
C   NTOT [real] Total number of atoms in the molecule.
C   RATIOM [real] Logarithm (base 10) of mass ratio (in g^(natoms-1)):
C     ratiom = Sum{log10(Atomic Masses)} - log10(Sum{Atomic Masses})
C   QPRD [double] Logarithm of product of atomic partition functions:
C     qprd = Sum{log10(Atomic Partition Functions)}
C
C Outputs:
C   EQK [real] Equilibrium constant (in dynes/cm/cm) at temperature T,
C     calculated from dissociation energy and partition function.
C   PART [real] Partition function at temperature T, calculated from
C     expressions in the references cited below.
C
C References:
C   For diatomic molecules: Sauval & Tatum (1984, ApJS, 56, 193).
C
      SUBROUTINE MOLCON(SPNAME,T,NTOT,RATIOM,QPRD,EQK,PART,PION,
     *                  BARKLEM)
C
      INCLUDE 'SIZES.EOS'
C
      INTEGER MSPEC,NTOT
      DOUBLE PRECISION KERG,KEV
      DOUBLE PRECISION RATIOM,QPRD,PION
      PARAMETER (KERG=1.38065D-16,KEV=KERG/1.60219D-12)
      PARAMETER (CONST=25947.256)
C
      REAL T
      DOUBLE PRECISION TLIM,TH,LOGTH,EQK,PART,Qm_spln,Kp_spln
c      DOUBLE PRECISION EQK_ST
      LOGICAL BARKLEM
C
C Combine equilibrium constant coefficients into one large array.
C
      PARAMETER (MSPEC=424)
      PARAMETER (NEQCOE=7)
      DOUBLE PRECISION COEF(NEQCOE,MSPEC)
      DOUBLE PRECISION C01(NEQCOE,50),C02(NEQCOE,50),
     *                 C03(NEQCOE,50),C04(NEQCOE,50),
     *                 C05(NEQCOE,50),C06(NEQCOE,50),
     *                 C07(NEQCOE,50),C08(NEQCOE,50),
     *                 C09(NEQCOE,24)
      EQUIVALENCE (C01(1,1),COEF(1,  1)),(C02(1,1),COEF(1, 51))
      EQUIVALENCE (C03(1,1),COEF(1,101)),(C04(1,1),COEF(1,151))
      EQUIVALENCE (C05(1,1),COEF(1,201)),(C06(1,1),COEF(1,251))
      EQUIVALENCE (C07(1,1),COEF(1,301)),(C08(1,1),COEF(1,351))
      EQUIVALENCE (C09(1,1),COEF(1,401))
C
C Combine partition function coefficients into one large array.
C
      PARAMETER (NPCOEF=11)
      DOUBLE PRECISION PCOEF(NPCOEF,MSPEC)
      DOUBLE PRECISION P01(NPCOEF,50),P02(NPCOEF,50),
     *                 P03(NPCOEF,50),P04(NPCOEF,50),
     *                 P05(NPCOEF,50),P06(NPCOEF,50),
     *                 P07(NPCOEF,50),P08(NPCOEF,50),
     *                 P09(NPCOEF,24)
      EQUIVALENCE (P01(1,1),PCOEF(1,  1)),(P02(1,1),PCOEF(1, 51))
      EQUIVALENCE (P03(1,1),PCOEF(1,101)),(P04(1,1),PCOEF(1,151))
      EQUIVALENCE (P05(1,1),PCOEF(1,201)),(P06(1,1),PCOEF(1,251))
      EQUIVALENCE (P07(1,1),PCOEF(1,301)),(P08(1,1),PCOEF(1,351))
      EQUIVALENCE (P09(1,1),PCOEF(1,401))
C
      CHARACTER SPNAME*(*),SPLIST(MSPEC)*(SPCHAR)
      SAVE
C
C Molecular species list from NextGen models (Allard & Hauschildt).
C See old/eos.4.f for molecular species list from Sauval & Tatum (1984).
C
      DATA SPLIST/
     * 'H2      ','CO      ','H2O     ','OH      ','N2      ',
     * 'SiO     ','HS      ','H2S     ','NH      ','SiH     ',
     * 'CH      ','H2+     ','NO      ','MgH     ','HCl     ',
     * 'SiS     ','AlOH    ','NH2     ','AlH     ','CN      ',
     * 'CO2     ','SO      ','TiO     ','S2      ','FeH     ',
     * 'NH3     ','HCN     ','HCO     ','O2      ','CH2     ',
     * 'HF      ','H3+     ','CaH     ','Al2O    ','AlO     ',
     * 'CH3     ','SiH2    ','MgO     ','C2      ','TiO2    ',
     * 'VO2     ','NaH     ','AlCl    ','AlF     ','VO      ',
     * 'CS      ','MgOH    ','PO2     ','CaOH    ','PH2     ',
     * 'C2H     ','ScO     ','AlO2H   ','AlS     ','FeO     ',
     * 'CrO     ','CH4     ','NS      ','SO2     ','SiN     ',
     * 'OH-     ','ZrO     ','NO+     ','ZrO2    ','BO      ',
     * 'SiO2    ','HBO     ','SiC     ','YO2     ','TiS     ',
     * 'HBO2    ','C2H2    ','OCS     ','ZrO+    ','NaOH    ',
     * 'CaCl    ','AlOF    ','YO      ','NaCl    ','C2O     ',
     * 'CHP     ','HS-     ','H2-     ','TiH     ','PH3     ',
     * 'MgS     ','TiO+    ','LaO2    ','Si2     ','SiH4    ',
     * 'BH2     ','AlOCl   ','LaO     ','C2N     ','AlBO2   ',
     * 'KCl     ','SiH-    ','CaF     ','CaO2H2  ','KOH     ',
     * 'CN-     ','Al2O2   ','BaOH    ','SrOH    ','BO2     ',
     * 'SiF     ','CH-     ','C3      ','C2-     ','MgO2H2  ',
     * 'BeOH    ','HBS     ','SiC2    ','FeO2H2  ','CrO2    ',
     * 'BeH2O2  ','BH3     ','NaCN    ','BeH2    ','Si2N    ',
     * 'CaCl2   ','NaBO2   ','C3H     ','OBF     ','CS2     ',
     * 'LiOH    ','Al2     ','LiCl    ','TiOCl   ','C2H4    ',
     * 'CHCl    ','TiCl    ','AlOF2   ','KBO2    ','Si2C    ',
     * 'CHF     ','BO-     ','AlO2    ','BaO2H2  ','OTiF    ',
     * 'CS-     ','C2N2    ','SrO2H2  ','ClCN    ','AlClF   ',
     * 'KCN     ','AlCl2   ','BaCl2   ','AlF2    ','MgCl2   ',
     * 'FeO-    ','BO2H2   ','SiH3Cl  ','FeCl2   ','Si3     ',
     * 'SiH3F   ','CH3Cl   ','SrCl2   ','CaF2    ','TiF2    ',
     * 'LiBO2   ','MgClF   ','BeBO2   ','C2HCl   ','TiCl2   ',
     * 'C4      ','H3BO3   ','MgF2    ','BaClF   ','BeF2    ',
     * 'C2HF    ','BeCl2   ','TiOCl2  ','ZrCl2   ','BaF2    ',
     * 'BeC2    ','Be2O    ','SrF2    ','ZrF2    ','FeF2    ',
     * 'P4      ','SiH2F2  ','H3O+    ','C5      ','TiF3    ',
     * 'TiCl3   ','ZrCl3   ','Na2Cl2  ','Na2O2H2 ','Be3O3   ',
     * 'K2Cl2   ','K2O2H2  ','ZrCl4   ','Na2C2N2 ','ZrF4    ',
     * 'Li2O2H2 ','CrH     ','Li2     ','B2      ','F2      ',
     * 'Na2     ','Mg2     ','P2      ','Cl2     ','K2      ',
     * 'Cu2     ','As2     ','Se2     ','Sb2     ','Te2     ',
     * 'I2      ','Cs2     ','He2+    ','C2+     ','N2+     ',
     * 'O2+     ','Ne2+    ','P2+     ','S2+     ','LiH     ',
     * 'BeH     ','BH      ','PH      ','KH      ','MnH     ',
     * 'CoH     ','NiH     ','CuH     ','ZnH     ','GaH     ',
     * 'GeH     ','AsH     ','SeH     ','HBr     ','RbH     ',
     * 'SrH     ','AgH     ','CdH     ','InH     ','SnH     ',
     * 'SbH     ','TeH     ','HI      ','CsH     ','BaH     ',
     * 'YbH     ','PtH     ','AuH     ','HgH     ','TlH     ',
     * 'PbH     ','BiH     ','HeH+    ','BeH+    ','CH+     ',
     * 'NH+     ','OH+     ','HF+     ','NeH+    ','MgH+    ',
     * 'AlH+    ','SiH+    ','PH+     ','SH+     ','HCl+    ',
     * 'ZnH+    ','HBr+    ','CdH+    ','HgH+    ','CF      ',
     * 'CP      ','CCl     ','CSe     ','CBr     ','RhC     ',
     * 'IrC     ','PtC     ','CN+     ','CO+     ','BN      ',
     * 'NF      ','AlN     ','PN      ','NCl     ','TiN     ',
     * 'AsN     ','SeN     ','ZrN     ','NS+     ','LiO     ',
     * 'BeO     ','FO      ','NaO     ','PO      ','ClO     ',
     * 'KO      ','CaO     ','MnO     ','NiO     ','CuO     ',
     * 'GaO     ','GeO     ','AsO     ','SeO     ','BrO     ',
     * 'RbO     ','SrO     ','NbO     ','InO     ','SnO     ',
     * 'SbO     ','TeO     ','IO      ','BaO     ','TbO     ',
     * 'LuO     ','HfO     ','TaO     ','WO      ','PtO     ',
     * 'PbO     ','BiO     ','ThO     ','BO+     ','SiO+    ',
     * 'PO+     ','SO+     ','AsO+    ','TaO+    ','LiF     ',
     * 'BeF     ','BF      ','NaF     ','MgF     ','PF      ',
     * 'SF      ','KF      ','ScF     ','MnF     ','NiF     ',
     * 'CuF     ','ZnF     ','GaF     ','GeF     ','AsF     ',
     * 'SeF     ','BrF     ','RbF     ','SrF     ','YF      ',
     * 'AgF     ','CdF     ','InF     ','SnF     ','SbF     ',
     * 'IF      ','CsF     ','BaF     ','LaF     ','HoF     ',
     * 'YbF     ','LuF     ','HgF     ','TlF     ','PbF     ',
     * 'LiNa    ','AsP     ','SbP     ','BeS     ','BS      ',
     * 'PS      ','CaS     ','ScS     ','CrS     ','CuS     ',
     * 'GeS     ','AsS     ','SeS     ','SrS     ','YS      ',
     * 'SnS     ','TeS     ','BaS     ','LaS     ','PbS     ',
     * 'BiS     ','BeCl    ','BCl     ','MgCl    ','SiCl    ',
     * 'PCl     ','ScCl    ','MnCl    ','FeCl    ','CuCl    ',
     * 'ZnCl    ','GaCl    ','GeCl    ','AsCl    ','SeCl    ',
     * 'BrCl    ','RbCl    ','SrCl    ','YCl     ','AgCl    ',
     * 'CdCl    ','InCl    ','SnCl    ','SbCl    ','ICl     ',
     * 'CsCl    ','BaCl    ','YbCl    ','AuCl    ','HgCl    ',
     * 'TlCl    ','PbCl    ','AlSe    ','SiSe    ','GeSe    ',
     * 'KBr     ','SiTe    ','GeTe    ','KI      '/
C
C Dissociation energy (first column, in eV) and equilibrium constant
C   coefficients. See the file "atomiz.notes" for the information on the
C   origin of the dissociation energies. The polynomial fit coefficients
C   for the equilibrium constants were determined with "ng_kfit.pro" and
C   are meant to reproduce the constants used in constructing the NextGen
C   models. The NextGen equilibrium constants were fit over the temperature
C   range 1600 < T < 7730 K. The fits are likely to diverge rapidly from
C   the truth outside this temperature range.
C Equilibrium constants may be constructed from the coefficients using:
C
C     log10(Kp) = Sum{i=2,7}{COEF(i)*log10(THETA)**(i-2)} - COEF(1)*THETA
C
      DATA C01/
     *   4.4781, 12.1354, -0.7752, -0.7821,  0.1464,  0.1603, -0.0626,  H2
     *  11.0920, 13.2368, -0.8342, -0.0477, -0.2923, -0.4557,  0.6108,  CO
     *   9.6221, 24.7774, -2.3428,  1.6868, -1.2845, -2.9925,  3.6555,  H2O
     *   4.3920, 11.8016, -0.8507, -0.5193,  0.0502, -0.3409,  0.4836,  OH
     *   9.7594, 12.8868, -0.8813,  0.2639, -1.5912,  1.5866, -0.5407,  N2
     *   8.2600, 12.9252, -0.7608, -0.3541,  1.5620, -3.5952,  2.5962,  SiO
     *   3.5500, 11.4382, -0.7816, -0.4659,  0.4314, -1.2144,  0.9648,  HS
     *   7.5946, 23.8543, -0.9525, -0.8118,  0.2051, -1.0299,  1.1555,  H2S
     *   3.4700, 11.4658, -0.7258, -0.6418, -0.0442,  0.2836, -0.1618,  NH
     *   3.0600, 11.2595, -0.6962, -0.6435,  0.6663, -0.3357, -0.4151,  SiH
     *   3.4650, 11.5333, -0.5255, -0.7105,  0.2264, -0.9271,  0.9577,  CH
     *   2.6508, 15.8052, 33.7578, 34.5956, 27.3455, 16.6214,  9.9717,  H2+
     *   6.4968, 11.9347, -0.7596,  0.0953, -0.9731,  0.8265, -0.2151,  NO
     *   1.3400, 10.2911, -0.3698, -0.0655, -2.9771,  6.1325, -4.3869,  MgH
     *   4.4336, 11.9041, -0.8281, -0.6163,  0.1580, -0.5068,  0.5164,  HCl
     *   6.4200, 12.6363, -0.7355,  0.0488,  0.8442, -2.0131,  1.3603,  SiS
     *  10.1252, 25.2575, -0.6810, -0.3051, -1.5765,  2.7536, -1.8355,  AlOH
     *   7.4400, 23.7389, -1.0179, -0.9947, -1.4353,  3.2530, -1.9224,  NH2
     *   3.0600, 11.4907, -0.4322, -0.6561, -0.5978,  2.4923, -2.4038,  AlH
     *   7.7600, 12.4438, -0.4756, -0.4909, -1.4623,  2.6823, -1.5396,  CN
     *  16.5382, 26.9571, -0.7464, -0.4921, -0.8506, -0.1365,  0.2358,  CO2
     *   5.3590, 12.3380, -0.4956, -0.2251, -0.1907, -0.2038,  0.2579,  SO
     *   6.8700, 11.9229, -1.4044,  0.7899, -0.7317, -0.0193, -0.4994,  TiO
     *   4.3693, 12.3190, -0.5050, -0.0290, -0.0266, -0.6002,  0.4572,  S2
c    *   2.4100, 12.1214,  0.9438,  2.2756, -0.1086,  4.1281, -1.9952,  FeH
c Dissociation energy from Dulick 2003
     *   1.5980, 12.1214,  0.9438,  2.2756, -0.1086,  4.1281, -1.9952,  FeH
     *  12.1388, 36.6661, -1.4062, -0.9258, -1.6969,  0.6005,  1.2302,  NH3
     *  13.2363, 25.1318, -0.5532, -0.0850, -0.9817,  0.6676,  0.3054,  HCN
     *  11.8560, 24.6414, -0.9415, -0.1856, -0.2948, -0.1630,  0.5836,  HCO
     *   5.1156, 12.8758, -0.4856, -0.5054, -0.0776, -0.0713,  0.2369,  O2
     *   7.9400, 23.8609, -1.0762, -0.4928, -0.4092,  0.0031,  0.3761,  CH2
     *   5.8690, 12.2896, -0.9180, -0.6238,  0.1243, -0.3525,  0.4767,  HF
c     *   0.0000, 18.8343, 12.4131, 11.9991,  6.8079,  8.4071,  2.6202,  H3+
     *   4.3730, 18.8343, 12.4131, 11.9991,  6.8079,  8.4071,  2.6202,  H3+
     *   1.7000, 10.1982, -0.9309,  1.8315, -5.6059,  6.9571, -3.5023,  CaH
     *  10.9653, 24.8807, -0.0033,  0.4796, -1.6979,  3.5631, -2.5414,  Al2O
     *   5.2700, 12.2132, -0.5246, -0.1918, -0.6810,  1.7287, -1.5839,  AlO
     *  12.6885, 36.6540, -1.3373, -1.0064, -0.5880, -0.2362,  0.8764,  CH3
     *   0.0000, 17.8513,-15.5361,-17.6144,-13.1604, -6.4819, -5.6361,  SiH2
     *   3.5300, 10.7940,  0.0122,  1.1189, -1.8758,  2.9976, -2.7758,  MgO
c     *   6.2100, 12.4672, -0.4452, -0.0100, -0.1868, -0.3860,  0.6230,  C2
     *   6.2970, 12.4672, -0.4452, -0.0100, -0.1868, -0.3860,  0.6230,  C2
c     *   6.3710, 12.4672, -0.4452, -0.0100, -0.1868, -0.3860,  0.6230,  C2
     *  13.2915, 25.9340, -1.4243,  1.6519, -0.7240, -0.7271,  0.7518,  TiO2
     *  12.9619, 25.9238, -1.2927,  1.3710, -2.4073,  2.2875, -0.5486,  VO2
     *   1.8800, 10.7184, -0.3642,  0.7843, -6.5309, 13.2912, -9.9502,  NaH
     *   5.1200, 11.8277, -0.3468, -1.0735,  1.8038, -1.7748,  0.4333,  AlCl
     *   6.8900, 12.2422, -0.4905, -0.4198,  0.0242,  0.3868, -0.5765,  AlF
     *   6.4100, 12.8108, -0.5811, -0.7895, -2.6766,  8.5158, -6.9993,  VO
     *   7.3550, 12.8487, -0.7627, -0.2538,  1.5240, -4.0119,  3.0234,  CS
     *   8.0735, 23.3256, -0.5884,  0.3637, -2.4401,  3.3936, -1.7121,  MgOH
     *  11.7451, 25.2051, -0.9105,  1.0031, -0.7207, -1.1064,  1.6239,  PO2
     *   8.7035, 23.1900, -1.0964,  2.5340, -5.9823,  5.3416, -1.1946,  CaOH
     *   6.4895, 23.0863, -1.3781,  0.2539, -0.6746, -1.2341,  1.5623/  PH2
      DATA C02/
     *  12.2087, 24.9752, -0.3204, -0.5640, -0.8997,  1.6927, -0.7771,  C2H
     *   6.9600, 12.5225, -1.2695,  1.7628, -2.0543, -1.2215,  2.3706,  ScO
     *  15.6364, 37.7022, -0.5885, -0.0823, -1.7283,  3.0502, -2.0176,  AlO2H
     *   3.8400, 11.9140, -0.5187, -0.1193, -0.3886,  1.1704, -1.2299,  AlS
     *   4.2000, 12.5326, -1.0657,  1.0360, -1.5641,  0.9560, -0.3218,  FeO
     *   4.4000, 11.0587, -1.3926,  1.4461, -2.1552,  3.3409, -3.1078,  CrO
     *  17.2173, 49.9426, -0.9720, -2.4957, -0.0017, -2.3299,  3.1042,  CH4
     *   4.8000, 11.9223, -0.6951,  0.1870, -0.7158,  0.4121,  0.0296,  NS
     *  11.1405, 25.9246, -0.5809,  0.0734, -0.3333,  0.1699,  0.0529,  SO2
     *   6.6880, 14.0972,  4.2904,  4.9608,  2.9390,  3.9789,  0.8908,  SiN
     *   4.7600, 19.9888, -6.7088, -4.3846, -2.8142, -2.3004, -0.3157,  OH-
     *   7.8500, 12.4674, -1.1280,  0.0368,  0.2221,  1.1043, -1.8804,  ZrO
     *  10.8500, 17.5169, 33.0097, 36.2110, 26.7396, 15.2392, 11.4130,  NO+
     *  14.4650, 25.6324, -1.5339,  1.1586, -0.9355,  1.6114, -1.2154,  ZrO2
     *   8.2800, 12.6246, -0.6966, -0.3874,  0.2531, -0.7582,  0.5307,  BO
     *  13.0355, 26.5610, -0.2891,  0.3006, -0.4009,  0.5864, -0.4006,  SiO2
     *  12.7425, 25.2283, -0.4780, -0.3611, -0.2189, -0.2108,  0.5883,  HBO
     *   4.6400, 11.8909, -0.8762,  0.1138,  0.0665, -0.5226,  0.3331,  SiC
     *  15.2000, 25.8617, -1.4050, -0.3896,  1.0805,  2.9269, -3.7531,  YO2
     *   4.7500, 11.6628, -1.4463,  1.3742, -0.8127, -0.4623,  0.2288,  TiS
     *  19.0991, 38.4541, -0.7808, -0.4220, -0.9239,  1.0793, -0.2304,  HBO2
     *  16.9704, 37.7481, -0.2529, -1.0622, -0.1485, -0.7058,  1.1910,  C2H2
     *  14.3762, 26.3815, -0.1712,  0.1197,  0.0059, -0.9891,  1.1946,  OCS
     *   0.0000,  2.5576, -0.5567, -4.5109, -4.3690, -0.1528, -3.1319,  ZrO+
     *   8.0150, 23.3420, -0.6139,  1.4091, -6.8466, 13.0407, -9.2977,  NaOH
     *   4.0900, 10.6268, -1.1367,  2.5278, -5.6022,  4.8741, -1.1616,  CaCl
     *  12.9003, 25.5751, -0.0730,  0.2808, -1.1757,  2.3733, -1.6726,  AlOF
     *   7.2900, 12.4422, -1.3547,  1.3087,  0.1688, -5.4106,  5.1158,  YO
     *   4.2300, 11.0864, -0.4463,  1.1926, -7.5820, 15.2552,-11.1116,  NaCl
     *  14.5371, 25.6134, -0.0508,  0.3710, -0.6246, -0.7682,  0.5868,  C2O
     *  11.4442, 24.7107, -0.5678, -0.0389,  1.0076, -4.6514,  4.3893,  CHP
     *   3.7900, 19.0227, -8.0668, -5.9821, -3.8685, -3.1838, -1.0364,  HS-
     *   0.7300, 19.7162, -5.0018, -2.7680, -1.2845, -0.9859, -0.3380,  H2-
     *   2.1200, 12.4717,  0.1601,  1.4596, -0.2012,  5.0788, -4.5487,  TiH
     *   9.7800, 35.8044, -1.3937, -0.2650, -0.6732, -2.5437,  2.9710,  PH3
     *   2.4000, 11.3146, -0.5595,  0.3619, -2.0065,  3.8766, -2.9900,  MgS
C 30-dec-2008 NP: added the dissociation energy from NIST
C
     *   0.0000,  4.5751,  3.4421,  0.7560, -1.7011,  1.4510, -1.3922,  TiO+
C    *  13.6890,  4.5751,  3.4421,  0.7560, -1.7011,  1.4510, -1.3922,  TiO+
     *  21.1510, 31.0805, 10.7070, 12.8687, 10.5799,  6.4414,  3.6171,  LaO2
     *   3.2100, 12.1817, -0.7102, -0.2403,  1.1042, -1.3644,  0.3198,  Si2
     *  13.2716, 48.6914, -1.0602, -1.2802, -0.8603,  0.1159, -0.0701,  SiH4
     *   8.2349, 24.0157, -0.6514, -0.6064, -0.6542,  0.9096, -0.5839,  BH2
     *  10.9011, 25.1839, -0.1060,  0.2530, -1.1850,  2.3355, -1.6111,  AlOCl
     *   8.2300, 12.1920,  0.1751, -0.7678, -1.3836,  1.7704, -0.0265,  LaO
     *  14.0629, 25.1475, -0.2270,  0.7024, -0.8499,  0.4583,  0.1889,  C2N
     *  20.0747, 38.6719, -0.2664,  0.2782, -1.2642,  1.6020, -0.5248,  AlBO2
     *   4.3400, 10.9561, -0.8720,  3.4218,-12.2306, 18.7863,-11.1011,  KCl
     *   3.2300, 19.3359, -5.7570, -3.5853, -1.3882, -2.3313, -0.4930,  SiH-
     *   5.4800, 11.0459, -0.8574,  2.3137, -4.6777,  4.4532, -1.1716,  CaF
     *  17.8875, 47.4921, -1.1390,  2.7534, -7.2248,  6.3242, -1.1381,  CaO2H2
     *   8.1892, 23.3129, -1.0581,  3.5131,-11.3115, 16.9078, -9.8867/  KOH
      DATA C03/
     *  10.3100, 21.7682, -5.8992, -3.8627, -4.0284,  1.2924, -2.5856,  CN-
     *  16.1405, 37.9519, -0.0230,  0.6639, -2.4910,  5.5385, -4.2945,  Al2O2
     *   9.0621, 23.3478, -2.1422,  1.7058, -1.6807, 10.3429,-14.0183,  BaOH
     *   8.6837, 23.1042, -1.2656,  3.2436, -7.2017,  6.5067, -1.7129,  SrOH
     *  13.9839, 25.6721, -0.0784,  0.0544, -0.2755,  0.6140, -0.3673,  BO2
     *   5.5700, 12.0158, -0.5187, -0.1216,  0.6738, -0.6377,  0.1588,  SiF
C
C 30-dec-2008 NP: added dissociation energy as dissociation energy of CH
C                 (3.465eV) + electron affinity of CH (1.238eV from NIST)
     *   0.0000, 16.4621,-13.8562,-13.1896, -9.2577, -6.3354, -2.5704,  CH-
C    *   4.7030, 16.4621,-13.8562,-13.1896, -9.2577, -6.3354, -2.5704,  CH-
     *  13.8610, 26.3081, -1.3134,  0.1185, -0.0461, -0.4056,  0.8088,  C3
     *   8.4800, 21.1413, -5.8697, -3.3745, -2.7491, -1.8902, -0.2441,  C2-
     *  17.1545, 48.1845, -0.5683,  0.1125, -3.0973,  4.3727, -2.1978,  MgO2H2
     *   9.3961, 23.7967, -0.6500,  0.2061, -1.9381,  2.1259, -0.6451,  BeOH
     *  10.4305, 24.8357, -0.4930, -0.4550,  0.8862, -2.7257,  2.4025,  HBS
     *  13.1966, 25.7392,  0.0961, -0.7979, -0.1515,  4.2750, -4.6336,  SiC2
     *  17.4231, 48.8561, -0.4831,  0.9575, -1.9798, -0.0476,  1.2346,  FeO2H2
     *  10.0930, 25.0689, -1.5784,  2.2605, -3.1152,  3.7375, -2.5596,  CrO2
     *  20.0817, 49.3051, -0.2203,  0.6123, -1.9159,  3.0362, -0.6588,  BeH2O2
     *  11.4541, 36.8342, -1.3068, -1.2283, -0.7130, -0.1039,  0.8121,  BH3
     *  12.5346, 24.2744, -0.4230,  2.1003, -7.6565, 14.5171,-10.4377,  NaCN
     *   6.5483, 23.5736, -0.7830, -0.0881, -2.2398,  2.7050, -1.5244,  BeH2
     *  10.1248, 24.8268, -0.3784,  0.5561, -0.7324,  1.7508, -1.6977,  Si2N
     *   9.3132, 22.5681, -0.7730,  3.2979, -6.3686,  5.5210, -0.9987,  CaCl2
     *  18.8913, 37.0212, -0.3881,  1.7934, -7.5472, 14.9782,-11.0505,  NaBO2
     *   0.0000, 19.8338,-46.6804,-50.9308,-35.9059,-13.5611,-23.8103,  C3H
     *  15.5315, 26.0301, -0.1824,  0.0109, -0.3944,  0.5184, -0.0882,  OBF
     *  11.9993, 26.2368, -0.1708,  0.2491,  0.4220, -2.2962,  2.2409,  CS2
     *   8.9381, 23.5703, -0.6263,  1.0060, -4.3983,  7.4665, -4.8955,  LiOH
     *   1.5500, 11.3681, -0.1946, -0.0669, -2.3347,  5.3477, -4.0343,  Al2
     *   4.8400, 11.3090, -0.5602,  0.5886, -3.9705,  7.3873, -5.2571,  LiCl
     *  11.3225, 25.4462, -1.0487,  1.8142, -1.5110,  0.4282, -0.0240,  TiOCl
     *  23.3326, 62.7915, -1.3095, -1.6903, -0.9624, -1.6171,  2.5521,  C2H4
     *   7.4689, 23.8059, -0.5629,  0.0019, -0.3896, -0.7781,  0.3890,  CHCl
     *   6.6900, 14.8883,  5.3193,  8.9551,  3.7271,  5.1452,  1.0391,  TiCl
     *  19.2284, 37.1933,  0.1308, -0.0614, -0.9981,  2.9770, -2.1833,  AlOF2
     *  18.9713, 36.8674, -0.8338,  3.8816,-11.3916, 16.8414, -9.6911,  KBO2
     *  11.2271, 25.9412,  0.1074, -0.8813, -0.2594,  4.4112, -4.4861,  Si2C
     *   9.2183, 24.5270, -0.6453, -1.0757, -0.7155,  2.2944, -1.4513,  CHF
     *   0.0000, 11.8175,-29.4442,-30.6402,-22.9279,-13.1209, -8.8023,  BO-
     *  10.9760, 27.6834,  5.5082,  6.6402,  5.5692,  2.7324,  1.9375,  AlO2
     *  18.0802, 47.0050, -2.3587,  2.3466, -2.2753,  8.4432,-11.3032,  BaO2H2
     *  12.8526, 25.8889, -1.0260,  1.8361, -1.5017,  0.3478,  0.0486,  OTiF
     *   6.5000, 20.6745, -7.9942, -5.7057, -2.6759, -6.1649,  1.2656,  CS-
     *  21.5636, 39.0495, -0.1190,  0.7088, -1.5184,  0.4914,  0.9277,  C2N2
     *  17.5958, 46.9386, -1.3295,  3.5725, -8.4710,  7.5694, -1.8456,  SrO2H2
     *  12.2076, 25.3442, -0.0379, -0.1189, -0.8276,  1.3188, -0.6986,  ClCN
     *  10.6135, 23.6489, -0.5207,  0.0519, -0.6538,  1.9149, -1.5058,  AlClF
     *  12.5010, 24.1386, -0.8692,  4.1888,-11.7377, 17.1662, -9.8522,  KCN
     *   8.8688, 23.5425, -0.5528,  0.0031, -0.7346,  2.3344, -1.9878,  AlCl2
     *   9.6070, 22.2204, -2.5275,  2.8555, -1.4987,  7.7865,-11.3039,  BaCl2
     *  12.3143, 24.3964, -0.4940,  0.0699, -0.5475,  1.6261, -1.2695,  AlF2
     *   8.1536, 22.9187, -0.1815,  0.6847, -2.4792,  4.3296, -2.7691/  MgCl2
      DATA C04/
     *   0.0000, 17.5598,-16.6727,-14.0707,-13.0780, -5.4193, -4.7856,  FeO-
     *  20.4537, 49.9913, -0.5362, -0.7176, -1.2169,  1.1206, -0.3773,  BO2H2
     *  14.1133, 48.5194, -0.8436, -1.0629, -0.7362,  0.3080, -0.3403,  SiH3Cl
     *   8.3239, 23.6272, -0.2108,  1.1105, -2.1105,  1.5380, -0.1684,  FeCl2
     *   7.3840, 24.8600, -0.1499, -0.1631,  0.1378,  1.6604, -1.9986,  Si3
     *  16.1268, 48.9782, -0.8260, -1.0380, -0.6452, -0.1029,  0.1199,  SiH3F
     *  16.2992, 49.7196, -1.2716, -1.4752, -1.1626,  0.6516, -0.0837,  CH3Cl
     *   9.1791, 22.1133, -1.4891,  4.1050, -7.6534,  6.6694, -1.5355,  SrCl2
     *  11.6845, 23.2600, -1.2039,  3.3661, -6.2828,  5.1661, -0.6547,  CaF2
     *  13.7563, 25.2856, -0.4137,  1.0746, -1.1248,  0.2935,  0.3807,  TiF2
     *  19.4163, 36.9346, -0.3977,  1.3814, -4.7577,  8.2956, -5.5779,  LiBO2
     *   9.5422, 23.6489, -0.6541,  0.7042, -2.5258,  4.5411, -3.0359,  MgClF
     *  19.3953, 37.4967, -0.4103,  0.6249, -2.5737,  3.7334, -2.0769,  BeBO2
     *  16.1988, 37.8077, -0.3545, -0.2428, -0.1731, -1.4896,  1.9844,  C2HCl
     *   9.9277, 24.6274, -0.5062,  0.9860, -1.3100,  0.8075, -0.0931,  TiCl2
     *  19.7168, 40.3256, -0.2533,  0.3731, -0.5863, -0.6939,  0.9337,  C4
     *  30.6562, 75.8041, -1.6269, -1.1205, -1.8109,  2.1354, -0.8357,  H3BO3
     *  10.7510, 23.8686, -0.6130,  0.7434, -2.6657,  5.0507, -3.5509,  MgF2
     *   0.0000, 13.8534,-28.5088,-27.6557,-25.0420, -4.2145,-21.0916,  BaClF
     *  13.3200, 24.6323, -0.2099,  0.5174, -1.9085,  2.9836, -1.7351,  BeF2
     *  16.6788, 38.1093, -0.3632, -0.2642, -0.4287, -0.5573,  0.9863,  C2HF
     *   9.6498, 23.7877, -0.2606,  0.4816, -1.7048,  2.1226, -0.8176,  BeCl2
     *  15.7352, 37.1910, -1.0480,  1.8371, -1.1420, -0.7526,  1.2880,  TiOCl2
     *  10.7683, 24.3508, -0.5859,  0.0972, -0.3635,  0.9082, -0.3338,  ZrCl2
     *  11.9101, 22.9073, -2.4413,  2.9420, -1.3655,  7.3312,-10.8692,  BaF2
     *  12.4073, 25.2586, -0.5256,  0.7548, -2.0655,  2.2598, -0.9944,  BeC2
     *   9.9676, 24.0020, -0.4765,  1.0925, -3.6131,  4.2582, -1.8225,  Be2O
     *  11.3542, 22.8132, -1.4157,  4.1790, -7.3508,  5.5696, -0.4507,  SrF2
     *  13.7587, 24.7160, -1.0103,  0.2376, -0.4664, -0.9114,  6.9672,  ZrF2
     *  13.0910, 27.6502,  6.5468,  8.2502,  7.3334,  4.1191,  1.2402,  FeF2
     *  12.5389, 37.9053, -1.3490,  3.1985, -1.1165, -6.7253,  7.3584,  P4
     *  19.0240, 49.7099, -0.5565, -0.7375, -0.2251, -1.1324,  1.2457,  SiH2F2
     *   3.2806, 41.7329, 32.0127, 34.5233, 27.1981, 13.3168, 13.4808,  H3O+
     *  27.0859, 54.0398,  0.0077,  0.4169, -0.9261, -0.3135,  0.6322,  C5
     *  19.7864, 37.9176, -0.7063,  1.7895, -1.5401,  0.9448, -0.6313,  TiF3
     *  14.3199, 37.3165, -0.8450,  1.6603, -1.6009,  0.8934, -0.5070,  TiCl3
     *  15.5540, 36.5254, -0.7361,  0.8503, -0.3688,  0.0324,  0.0881,  ZrCl3
     *  10.6603, 34.6664, -0.4567,  3.2641,-13.6211, 27.6173,-20.7914,  Na2Cl2
     *  18.1954, 60.7438, -0.7643,  2.2577,-14.4187, 28.3225,-20.4866,  (NaOH)2
     *  28.8149, 64.3940, -0.2174,  1.3367, -6.6368,  8.6309, -4.6284,  Be3O3
     *  10.8345, 33.9871, -1.3140,  7.4840,-21.9583, 33.6428,-20.3143,  K2Cl2
     *  18.3196, 60.4179, -1.6298,  6.4524,-22.9230, 33.8810,-20.0092,  (KOH)2
     *  20.4364, 49.7173, -0.6667,  0.8064, -0.1308, -0.4433,  0.8970,  ZrCl4
     *  27.1266, 62.7471, -0.3813,  3.6624,-15.0927, 27.0694,-18.7738,  (NaCN)2
     *  27.0557, 51.2712, -0.5271,  0.8930, -0.5666,  1.5292, -1.3568,  ZrF4
     *  20.3442, 61.3686, -0.8410,  1.3617, -9.5297, 16.1158,-11.1739,  (LiOH)2
     *   1.9300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CrH
     *   1.0499,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Li2
     *   2.8020,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  B2
     *   1.6060,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/  F2
      DATA C05/
     *   0.7368,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Na2
     *   0.0790,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Mg2
     *   5.0310,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  P2
     *   2.4740,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Cl2
     *   0.5520,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  K2
     *   2.0430,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Cu2
     *   3.9600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  As2
     *   3.3870,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Se2
     *   3.0880,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Sb2
     *   2.6330,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Te2
     *   1.5395,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  I2
     *   0.4167,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Cs2
     *   2.4456,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  He2+
     *   6.2020,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  C2+
     *   8.7076,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  N2+
     *   6.3670,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  O2+
     *   1.2600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  Ne2+
     *   4.9500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  P2+
     *   5.1430,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  S2+
     *   2.4286,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  LiH
     *   1.9730,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BeH
     *   3.5390,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BH
     *   3.0400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PH
     *   1.7708,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  KH
     *   2.6020,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  MnH
     *   2.4980,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CoH
     *   2.4510,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NiH
     *   2.6020,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CuH
     *   0.8500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ZnH
     *   2.8190,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GaH
     *   2.6890,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GeH
     *   2.8020,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AsH
     *   3.2200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SeH
     *   3.7560,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HBr
     *   1.7480,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  RbH
     *   1.6600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SrH
     *   2.0600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AgH
     *   0.6770,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CdH
     *   2.4810,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  InH
     *   2.6900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SnH
     *   2.4460,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SbH
     *   2.7670,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TeH
     *   3.0529,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HI
     *   1.7790,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CsH
     *   1.9500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BaH
     *   1.8600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  YbH
     *   3.3870,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PtH
     *   3.3610,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AuH
     *   0.3744,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HgH
     *   1.9870,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/  TlH
      DATA C06/
     *   1.5900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PbH
     *   2.9000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BiH
     *   1.8450,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HeH+
     *   3.1440,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BeH+
     *   4.0849,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CH+
     *   4.4770,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NH+
     *   5.0182,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  OH+
     *   3.4230,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HF+
     *   2.0800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NeH+
     *   1.9390,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  MgH+
     *   1.6310,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AlH+
     *   3.2440,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SiH+
     *   3.3790,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PH+
     *   3.5690,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SH+
     *   4.6569,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HCl+
     *   2.2000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ZnH+
     *   3.8920,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HBr+
     *   1.8220,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CdH+
     *   2.1080,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HgH+
     *   5.7110,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CF
     *   5.2800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CP
     *   4.0770,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CCl
     *   6.0800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CSe
     *   3.2570,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CBr
     *   5.9720,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  RhC
     *   6.5010,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  IrC
     *   6.2840,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PtC
     *   5.3950,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CN+
     *   8.3654,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CO+
     *   3.8770,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BN
     *   3.3000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NF
     *   3.7780,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AlN
     *   6.3600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PN
     *   3.4220,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NCl
     *   4.9000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TiN
     *   5.0310,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AsN
     *   3.8600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SeN
     *   5.8200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ZrN
     *   5.3000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NS+
     *   3.4910,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  LiO
     *   4.4900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BeO
     *   2.2420,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  FO
     *   2.7580,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NaO
     *   6.0670,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PO
     *   2.7337,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ClO
     *   2.7760,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  KO
     *   3.9860,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CaO
     *   3.7100,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  MnO
     *   3.7600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NiO
     *   2.9400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/  CuO
      DATA C07/
     *   3.8400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GaO
     *   6.8040,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GeO
     *   4.9740,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AsO
     *   4.4150,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SeO
     *   2.4289,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BrO
     *   2.8230,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  RbO
     *   4.3800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SrO
     *   7.4900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NbO
     *   3.5500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  InO
     *   5.4300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SnO
     *   4.4600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SbO
     *   3.8600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TeO
     *   2.4500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  IO
     *   5.7900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BaO
     *   7.1560,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TbO
     *   6.8950,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  LuO
     *   8.2600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HfO
     *   8.6560,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TaO
     *   7.4200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  WO
     *   4.0200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PtO
     *   3.8400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PbO
     *   3.4600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BiO
     *   9.0510,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ThO
     *   3.3400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BO+
     *   4.9100,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SiO+
     *   8.2400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PO+
     *   5.3950,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SO+
     *   5.0910,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AsO+
     *   7.8490,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TaO+
     *   5.9500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  LiF
     *   5.9020,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BeF
     *   7.5500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BF
     *   4.9090,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NaF
     *   4.7600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  MgF
     *   4.5620,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PF
     *   3.5220,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SF
     *   5.0310,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  KF
     *   6.1710,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ScF
     *   4.5750,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  MnF
     *   4.4000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  NiF
     *   4.4200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CuF
     *   3.7300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ZnF
     *   6.0200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GaF
     *   5.3800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GeF
     *   4.2000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AsF
     *   3.4700,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SeF
     *   2.8600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BrF
     *   5.0800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  RbF
     *   5.5400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SrF
     *   7.0600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/  YF
      DATA C08/
     *   3.5300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AgF
     *   3.1200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CdF
     *   5.3100,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  InF
     *   4.8920,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SnF
     *   4.5000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SbF
     *   2.7800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  IF
     *   5.3210,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CsF
     *   5.9800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BaF
     *   6.7900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  LaF
     *   5.5500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HoF
     *   5.3300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  YbF
     *   4.1600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  LuF
     *   1.8000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HgF
     *   4.5100,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TlF
     *   3.6400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PbF
     *   0.8650,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  LiNa
     *   4.4500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AsP
     *   3.6600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SbP
     *   3.2350,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BeS
     *   5.7110,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BS
     *   4.2890,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PS
     *   3.4300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CaS
     *   4.9200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ScS
     *   3.3900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CrS
     *   2.8100,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CuS
     *   5.4900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GeS
     *   3.8900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AsS
     *   3.8640,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SeS
     *   3.4700,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SrS
     *   5.4400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  YS
     *   4.8000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SnS
     *   3.4300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TeS
     *   4.3000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BaS
     *   5.9020,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  LaS
     *   4.0900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PbS
     *   3.2300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BiS
     *   3.8420,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BeCl
     *   5.3340,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BCl
     *   3.2000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  MgCl
     *   4.2810,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SiCl
     *   3.1920,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PCl
     *   3.3900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ScCl
     *   3.4700,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  MnCl
     *   3.3800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  FeCl
     *   3.8800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CuCl
     *   2.3330,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ZnCl
     *   4.7600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GaCl
     *   4.0100,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GeCl
     *   4.6000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AsCl
     *   3.3000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/  SeCl
      DATA C09/
     *   2.2346,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BrCl
     *   4.3930,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  RbCl
     *   4.2000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SrCl
     *   5.3800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  YCl
     *   3.2200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AgCl
     *   2.1200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CdCl
     *   4.4000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  InCl
     *   3.5900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SnCl
     *   3.6900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SbCl
     *   2.1514,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  ICl
     *   4.5800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  CsCl
     *   4.5500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  BaCl
     *   3.8400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  YbCl
     *   2.8700,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AuCl
     *   0.9150,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  HgCl
     *   3.8300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  TlCl
     *   3.0800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  PbCl
     *   3.2600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  AlSe
     *   5.5400,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SiSe
     *   4.9830,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GeSe
     *   3.8900,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  KBr
     *   3.9770,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  SiTe
     *   4.0720,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  GeTe
     *   3.3000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/  KI
C
C Coefficients for constructing partition functions (and then equilibrium
C   constants, perhaps). For diatomic molecules other than H2 and CO, the
C   data are from Sauval & Tatum (1984, ApJS, 56, 193). For H2 and CO, the
C   data are from Irwin (1987, A&A, 182, 348). For polyatomic molecules,
C   the coefficients are from Irwin (1988, A&AS, 74,145).
C Coefficients used to construct the partition function, as follows:
C
C     log10(Q) = Sum{i=0,9}{PCOEF(i+1)*log10(THETA)**i}
C                                                           Ioniz. pot.
      DATA P01/
     *   1.69179,      -1.72270,       0.798033,     -0.157089,         H2
     *  -0.535313,      1.75818,      -2.63895,       1.35708,          H2
     *   0.0,           0.0,                                 15.42593,  H2
     *   3.615300,     -1.773848,      0.3516181,     0.08620792,       CO
     *   0.2911791,    -1.141469,      2.513133,     -2.886502,         CO
     *   1.238932,      0.0,                                 14.01400,  CO
     *   4.344711818,  -3.6343233,     1.415963,      0.01594,          H2O
     *   0.56542,      -1.2583,        0.53796,       3*0.0, 12.62100,  H2O
     *   3.0929, -1.6778,  0.6743, -0.1874,  0.0000,  5*0.0, 13.01700,  OH
     *   3.2643, -1.7303,  0.4192,  0.0000,  0.0000,  5*0.0, 15.58100,  N2
     *   4.2275, -1.9144,  0.7201, -1.3099,  1.1657,  5*0.0, 11.49000,  SiO
     *  1.0, 9*0.,                                           10.42200,  HS
     *   5.117210341,  -3.94844146,    1.23193,       0.076156,         H2S
     *   0.42163,      -0.453534,      0.0,           3*0.0, 10.45700,  H2S
     *   3.0735, -1.8501,  0.9607, -0.3935,  0.0000,  5*0.0, 13.49000,  NH
     *   3.6908, -1.9801,  0.7704, -0.2247,  0.0000,  5*0.0,  7.91000,  SiH
     *   3.3586, -2.0656,  0.9624, -0.2239,  0.0000,  5*0.0, 10.64000,  CH
     *   2.5410, -2.4336,  1.4979,  0.0192, -0.7483,  5*0.0, -1.00000,  H2+
     *   4.3073, -1.8255,  0.3765,  0.0000,  0.0000,  5*0.0,  9.26420,  NO
     *   3.6704, -2.2682,  0.9354, -0.2597,  0.0000,  5*0.0,  7.20000,  MgH
     *   2.8005, -1.7476,  0.5310,  0.0000,  0.0000,  5*0.0, 12.74400,  HCl
     *   4.8026, -1.9753,  0.2600,  0.0000,  0.0000,  5*0.0, 10.53000,  SiS
     *   6.103792598,  -4.3938712,     0.662588,      0.3751,           AlOH
     *   0.38386,      -0.2147,        0.0,           3*0.0, -1.00000,  AlOH
     *   4.819621858,  -3.84200734,    1.5386462,     0.784399,         NH2
     *  -2.34404,       2.50803,      -1.13304,       3*0.0, 11.14000,  NH2
     *   3.3209, -2.5909,  1.7415, -0.7636,  0.0000,  5*0.0,  5.50000,  AlH
     *   4.0078, -2.1514,  0.9226, -0.1671,  0.0000,  5*0.0, 13.59800,  CN
     *   6.01081285,   -4.438833,      0.840462,      0.2945,           CO2
     *   0.3694,       -0.273,         0.0,           3*0.0, 13.77700,  CO2
     *   4.7963, -2.1308,  0.5224,  0.0000,  0.0000,  5*0.0, 10.29400,  SO
C The line with 5.7765 is from Alard and Hauschildt who artificially increased
C TiO parition function by a factor of 3. Also change in ionization energy
C according to the latest NIST data.
C    *   5.7765, -2.3739,  0.8940, -0.3641,  0.0000,  5*0.0,  6.40000,  TiO
     *   5.3051, -2.3739,  0.8940, -0.3641,  0.0000,  5*0.0,  6.81900,  TiO
     *   5.0796, -2.1967,  0.4101,  0.0000,  0.0000,  5*0.0,  9.35600,  S2
     *   4.6265980,    -2.5625800,     0.38885943,    0.40219820,       FeH
     *  -0.21386399,    0.027845045,   0.0,           3*0.0,  7.37000,  FeH
     *   5.884176216,  -5.8364867,     1.608417,      1.50876,          NH3
     *  -0.59607,      -0.58961,       0.2459,        3*0.0, -1.00000,  NH3
     *   5.434042379,  -4.2409874,     0.988745,      0.49464,          HCN
     *   0.03719,      -0.22924,       0.0,           3*0.0, 13.60000,  HCN
     *   6.298781639,  -3.85672804,    0.8551678,     0.321901,         HCO
     *   0.020274,      0.15254,      -0.25298,       3*0.0,  8.12000,  HCO
     *   4.0636, -2.0779,  0.7660, -0.2111,  0.0000,  5*0.0, 12.06970,  O2
     *  1.0, 9*0.,                                           10.39600,  CH2
     *   2.4164, -1.6132,  0.6357, -0.1767,  0.0000,  5*0.0, 16.03000,  HF
     *  1.0, 9*0.,                                           -1.00000,  H3+
     *   3.8411, -2.3891,  1.3578, -0.6893,  0.0000,  5*0.0,  5.86000,  CaH
     *  1.0, 9*0.,                                           -1.00000,  Al2O
     *   4.9191, -2.6291,  0.5831,  0.3163,  0.0000,  5*0.0,  9.46000,  AlO
     *  1.0, 9*0.,                                            9.84000,  CH3
     *  1.0, 9*0.,                                            8.80000,  SiH2
     *   5.3182, -2.6502, -0.2781, -0.7823,  1.3107,  5*0.0,  8.76000,  MgO
     *   4.3091, -2.2406,  0.4865, -0.2049,  0.0000,  5*0.0, 11.40000,  C2
     *  1.0, 9*0.,                                            9.50000,  TiO2
     *   8.457240767,  -4.1987868,     0.334575,      0.20744,          VO2
     *   0.18226,      -0.053465,      0.0,           3*0.0, -1.00000,  VO2
     *   3.5453, -2.3457,  0.8557, -0.1685,  0.0000,  5*0.0,  4.70000,  NaH
     *   5.1115, -2.2303,  0.8001, -0.5192,  0.0000,  5*0.0,  9.40000,  AlCl
     *   4.5405, -2.1033,  0.6208, -0.2930,  0.0000,  5*0.0, -1.00000,  AlF
     *   5.0687, -2.2186,  0.9545, -0.4592,  0.0000,  5*0.0,  7.23860,  VO
     *   4.1646, -1.9348,  0.8034, -1.3669,  1.1561,  5*0.0, 11.33000,  CS
     *   6.8401894714, -4.338616427,   0.71600166,    0.128126,         MgOH
     *   0.5978087,    -0.8658369,     0.385049,      3*0.0,  7.50000,  MgOH
     *  1.0, 9*0.,                                           11.90000,  PO2
     *   7.1623971155, -4.471282563,   1.1221899,    -0.558812,         CaOH
     *   0.2294,        1.78658,      -2.95118,       1.41591,          CaOH
     *   2*0.0,                                               5.80000,  CaOH
     *  1.0, 9*0.,                                            9.82400/  PH2
      DATA P02/
     *  1.0, 9*0.,                                           11.61000,  C2H
     *   4.8065, -2.2129,  0.9991, -0.5414,  0.0000,  5*0.0, -1.00000,  ScO
     *  1.0, 9*0.,                                           -1.00000,  AlO2H
     *   5.2461, -2.1319,  0.5340, -0.2309,  0.0000,  5*0.0, -1.00000,  AlS
     *   5.5642, -2.1947,  0.5065,  0.0000,  0.0000,  5*0.0,  8.90000,  FeO
     *   5.5270, -2.1311,  0.6523, -0.2533,  0.0000,  5*0.0,  7.85000,  CrO
     *  1.0, 9*0.,                                           12.61000,  CH4
     *   4.8052, -1.9619,  0.3140,  0.0000,  0.0000,  5*0.0,  8.87000,  NS
     *  1.0, 9*0.,                                           12.34900,  SO2
     *   4.6570, -2.3587,  0.8819, -0.1642,  0.0000,  5*0.0, -1.00000,  SiN
     *  1.0, 9*0.,                                           -1.00000,  OH-
     *   5.3279, -2.4694,  0.2164, -0.2313,  0.0000,  5*0.0,  6.00000,  ZrO
     *   3.5649, -1.7328,  0.4241,  0.0000,  0.0000,  5*0.0, -1.00000,  NO+
     *   8.72011985,   -4.247295,      0.2758,        0.20738,          ZrO2
     *   0.09406,       0.0,           0.0,           3*0.0, -1.00000,  ZrO2
     *   3.9953, -1.8665,  0.5965, -0.1617,  0.0000,  5*0.0, 13.30000,  BO
     *  1.0, 9*0.,                                           -1.00000,  SiO2
     *  1.0, 9*0.,                                           -1.00000,  HBO
     *   5.1477, -1.8671,  0.2404,  0.0000,  0.0000,  5*0.0,  9.20000,  SiC
     *  1.0, 9*0.,                                           -1.00000,  YO2
     *   5.8948, -2.2183,  0.5928, -0.3106,  0.0000,  5*0.0,  7.10000,  TiS
     *  1.0, 9*0.,                                           -1.00000,  HBO2
     *   7.1220464309, -6.966653604,   1.9668235,     0.362597,         C2H2
     *   0.608996,     -0.920435,      0.271892,      3*0.0, 11.40000,  C2H2
     *  1.0, 9*0.,                                           11.18500,  OCS
     *  1.0, 9*0.,                                           -1.00000,  ZrO+
     *  1.0, 9*0.,                                           -1.00000,  NaOH
     *   5.7494, -2.3340,  0.8685, -0.5306,  0.0000,  5*0.0,  5.86000,  CaCl
     *  1.0, 9*0.,                                           -1.00000,  AlOF
     *   4.9515, -2.0866,  0.6565, -0.3082,  0.0000,  5*0.0,  6.00000,  YO
     *   5.3364, -2.2844,  0.2820,  0.1185,  0.0000,  5*0.0, -1.00000,  NaCl
     *  1.0, 9*0.,                                           -1.00000,  C2O
     *  1.0, 9*0.,                                           10.79000,  CHP
     *  1.0, 9*0.,                                           -1.00000,  HS-
     *  1.0, 9*0.,                                           -1.00000,  H2-
     *  1.0, 9*0.,                                            6.00000,  TiH
     *  1.0, 9*0.,                                            9.86900,  PH3
     *   5.0367, -2.1625,  0.4859, -0.1780,  0.0000,  5*0.0, -1.00000,  MgS
     *  1.0, 9*0.,                                           -1.00000,  TiO+
     *  1.0, 9*0.,                                           -1.00000,  LaO2
     *   5.2617, -2.1485,  0.5647, -0.2985,  0.0000,  5*0.0, -1.00000,  Si2
     *  1.0, 9*0.,                                           -1.00000,  SiH4
     *  1.0, 9*0.,                                            9.80000,  BH2
     *  1.0, 9*0.,                                           -1.00000,  AlOCl
     *   5.1147, -2.5016,  1.0445, -0.3135,  0.0000,  5*0.0,  4.95000,  LaO
     *  1.0, 9*0.,                                           12.00000,  C2N
     *  1.0, 9*0.,                                           -1.00000,  AlBO2
     *   5.6860, -2.3016,  0.2086,  0.1763,  0.0000,  5*0.0, -1.00000,  KCl
     *  1.0, 9*0.,                                           -1.00000,  SiH-
     *   5.2010, -2.2653,  0.8941, -0.5384,  0.0000,  5*0.0, -1.00000,  CaF
     *  1.0, 9*0.,                                           -1.00000,  CaO2H2
     *  1.0, 9*0.,                                            7.50000/  KOH
      DATA P03/
     *  1.0, 9*0.,                                           -1.00000,  CN-
     *  1.0, 9*0.,                                           -1.00000,  Al2O2
     *  1.0, 9*0.,                                           -1.00000,  BaOH
     *  1.0, 9*0.,                                           -1.00000,  SrOH
     *  1.0, 9*0.,                                           -1.00000,  BO2
     *   5.0871, -2.0375,  0.4478, -0.1243,  0.0000,  5*0.0,  7.54000,  SiF
     *  1.0, 9*0.,                                           -1.00000,  CH-
     *   6.618407932,  -3.576399,      0.883642,      0.087548,         C3
     *   0.04817,      -0.16471,       0.0,           3*0.0, -1.00000,  C3
     *  1.0, 9*0.,                                           -1.00000,  C2-
     *  1.0, 9*0.,                                           -1.00000,  MgO2H2
     *  1.0, 9*0.,                                           -1.00000,  BeOH
     *  1.0, 9*0.,                                           -1.00000,  HBS
     *   7.54651307623,-5.075563869,   1.82960795,    0.0983258,        SiC2
     *  -6.335157,     14.33103,     -13.01689,       4.428233,         SiC2
     *   2*0.0,                                              10.20000,  SiC2
     *  1.0, 9*0.,                                           -1.00000,  FeO2H2
     *  1.0, 9*0.,                                           -1.00000,  CrO2
     *  1.0, 9*0.,                                           -1.00000,  BeH2O2
     *  1.0, 9*0.,                                           -1.00000,  BH3
     *  1.0, 9*0.,                                           -1.00000,  NaCN
     *  1.0, 9*0.,                                           -1.00000,  BeH2
     *  1.0, 9*0.,                                           -1.00000,  Si2N
     *  1.0, 9*0.,                                           -1.00000,  CaCl2
     *  1.0, 9*0.,                                           -1.00000,  NaBO2
     *  1.0, 9*0.,                                           -1.00000,  C3H
     *  1.0, 9*0.,                                           -1.00000,  OBF
     *  1.0, 9*0.,                                           10.07300,  CS2
     *  1.0, 9*0.,                                           -1.00000,  LiOH
     *   5.5538, -2.3365,  0.5754, -0.2119,  0.0000,  5*0.0,  5.40000,  Al2
     *   4.5605, -2.2216,  0.5760, -0.1706,  0.0000,  5*0.0,  9.57000,  LiCl
     *  1.0, 9*0.,                                           -1.00000,  TiOCl
     *  1.0, 9*0.,                                           -1.00000,  C2H4
     *  1.0, 9*0.,                                           -1.00000,  CHCl
     *  1.0, 9*0.,                                           -1.00000,  TiCl
     *  1.0, 9*0.,                                           -1.00000,  AlOF2
     *  1.0, 9*0.,                                           -1.00000,  KBO2
     *  1.0, 9*0.,                                           -1.00000,  Si2C
     *  1.0, 9*0.,                                           10.06000,  CHF
     *  1.0, 9*0.,                                           -1.00000,  BO-
     *  1.0, 9*0.,                                           -1.00000,  AlO2
     *  1.0, 9*0.,                                           -1.00000,  BaO2H2
     *  1.0, 9*0.,                                           -1.00000,  OTiF
     *  1.0, 9*0.,                                           -1.00000,  CS-
     *  1.0, 9*0.,                                           -1.00000,  C2N2
     *  1.0, 9*0.,                                           -1.00000,  SrO2H2
     *  1.0, 9*0.,                                           12.36000,  ClCN
     *  1.0, 9*0.,                                           -1.00000,  AlClF
     *  1.0, 9*0.,                                           -1.00000,  KCN
     *  1.0, 9*0.,                                           -1.00000,  AlCl2
     *  1.0, 9*0.,                                           -1.00000,  BaCl2
     *  1.0, 9*0.,                                           -1.00000,  AlF2
     *  1.0, 9*0.,                                           -1.00000/  MgCl2
      DATA P04/
     *  1.0, 9*0.,                                           -1.00000,  FeO-
     *  1.0, 9*0.,                                           -1.00000,  BO2H2
     *  1.0, 9*0.,                                           -1.00000,  SiH3Cl
     *  1.0, 9*0.,                                           -1.00000,  FeCl2
     *  1.0, 9*0.,                                           -1.00000,  Si3
     *  1.0, 9*0.,                                           -1.00000,  SiH3F
     *  1.0, 9*0.,                                           -1.00000,  CH3Cl
     *  1.0, 9*0.,                                           -1.00000,  SrCl2
     *  1.0, 9*0.,                                           -1.00000,  CaF2
     *  1.0, 9*0.,                                           -1.00000,  TiF2
     *  1.0, 9*0.,                                           -1.00000,  LiBO2
     *  1.0, 9*0.,                                           -1.00000,  MgClF
     *  1.0, 9*0.,                                           -1.00000,  BeBO2
     *  1.0, 9*0.,                                           -1.00000,  C2HCl
     *  1.0, 9*0.,                                           -1.00000,  TiCl2
     *  1.0, 9*0.,                                           -1.00000,  C4
     *  1.0, 9*0.,                                           -1.00000,  H3BO3
     *  1.0, 9*0.,                                           -1.00000,  MgF2
     *  1.0, 9*0.,                                           -1.00000,  BaClF
     *  1.0, 9*0.,                                           -1.00000,  BeF2
     *  1.0, 9*0.,                                           -1.00000,  C2HF
     *  1.0, 9*0.,                                           -1.00000,  BeCl2
     *  1.0, 9*0.,                                           -1.00000,  TiOCl2
     *  1.0, 9*0.,                                           -1.00000,  ZrCl2
     *  1.0, 9*0.,                                           -1.00000,  BaF2
     *  1.0, 9*0.,                                           -1.00000,  BeC2
     *  1.0, 9*0.,                                           -1.00000,  Be2O
     *  1.0, 9*0.,                                           -1.00000,  SrF2
     *  1.0, 9*0.,                                           -1.00000,  ZrF2
     *  1.0, 9*0.,                                           -1.00000,  FeF2
     *  1.0, 9*0.,                                           -1.00000,  P4
     *  1.0, 9*0.,                                           -1.00000,  SiH2F2
     *  1.0, 9*0.,                                           -1.00000,  H3O+
     *  1.0, 9*0.,                                           -1.00000,  C5
     *  1.0, 9*0.,                                           -1.00000,  TiF3
     *  1.0, 9*0.,                                           -1.00000,  TiCl3
     *  1.0, 9*0.,                                           -1.00000,  ZrCl3
     *  1.0, 9*0.,                                           -1.00000,  Na2Cl2
     *  1.0, 9*0.,                                           -1.00000,  Na2O2H2
     *  1.0, 9*0.,                                           -1.00000,  Be3O3
     *  1.0, 9*0.,                                           -1.00000,  K2Cl2
     *  1.0, 9*0.,                                           -1.00000,  K2O2H2
     *  1.0, 9*0.,                                           -1.00000,  ZrCl4
     *  1.0, 9*0.,                                           -1.00000,  Na2C2N2
     *  1.0, 9*0.,                                           -1.00000,  ZrF4
     *  1.0, 9*0.,                                           -1.00000,  Li2O2H2
     *  1.0, 9*0.,                                            7.33000,  CrH
     *  1.0, 9*0.,                                           -1.00000,  Li2
     *  1.0, 9*0.,                                           -1.00000,  B2
     *  1.0, 9*0.,                                           -1.00000/  F2
      DATA P05/
     *  1.0, 9*0.,                                           -1.00000,  Na2
     *  1.0, 9*0.,                                           -1.00000,  Mg2
     *  1.0, 9*0.,                                           -1.00000,  P2
     *  1.0, 9*0.,                                           -1.00000,  Cl2
     *  1.0, 9*0.,                                           -1.00000,  K2
     *  1.0, 9*0.,                                           -1.00000,  Cu2
     *  1.0, 9*0.,                                           -1.00000,  As2
     *  1.0, 9*0.,                                           -1.00000,  Se2
     *  1.0, 9*0.,                                           -1.00000,  Sb2
     *  1.0, 9*0.,                                           -1.00000,  Te2
     *  1.0, 9*0.,                                           -1.00000,  I2
     *  1.0, 9*0.,                                           -1.00000,  Cs2
     *  1.0, 9*0.,                                           -1.00000,  He2+
     *  1.0, 9*0.,                                           -1.00000,  C2+
     *  1.0, 9*0.,                                           -1.00000,  N2+
     *  1.0, 9*0.,                                           -1.00000,  O2+
     *  1.0, 9*0.,                                           -1.00000,  Ne2+
     *  1.0, 9*0.,                                           -1.00000,  P2+
     *  1.0, 9*0.,                                           -1.00000,  S2+
     *  1.0, 9*0.,                                           -1.00000,  LiH
     *  1.0, 9*0.,                                           -1.00000,  BeH
     *  1.0, 9*0.,                                           -1.00000,  BH
     *  1.0, 9*0.,                                           -1.00000,  PH
     *  1.0, 9*0.,                                           -1.00000,  KH
     *  1.0, 9*0.,                                           -1.00000,  MnH
     *  1.0, 9*0.,                                           -1.00000,  CoH
     *  1.0, 9*0.,                                           -1.00000,  NiH
     *  1.0, 9*0.,                                           -1.00000,  CuH
     *  1.0, 9*0.,                                           -1.00000,  ZnH
     *  1.0, 9*0.,                                           -1.00000,  GaH
     *  1.0, 9*0.,                                           -1.00000,  GeH
     *  1.0, 9*0.,                                           -1.00000,  AsH
     *  1.0, 9*0.,                                           -1.00000,  SeH
     *  1.0, 9*0.,                                           -1.00000,  HBr
     *  1.0, 9*0.,                                           -1.00000,  RbH
     *  1.0, 9*0.,                                           -1.00000,  SrH
     *  1.0, 9*0.,                                           -1.00000,  AgH
     *  1.0, 9*0.,                                           -1.00000,  CdH
     *  1.0, 9*0.,                                           -1.00000,  InH
     *  1.0, 9*0.,                                           -1.00000,  SnH
     *  1.0, 9*0.,                                           -1.00000,  SbH
     *  1.0, 9*0.,                                           -1.00000,  TeH
     *  1.0, 9*0.,                                           -1.00000,  HI
     *  1.0, 9*0.,                                           -1.00000,  CsH
     *  1.0, 9*0.,                                           -1.00000,  BaH
     *  1.0, 9*0.,                                           -1.00000,  YbH
     *  1.0, 9*0.,                                           -1.00000,  PtH
     *  1.0, 9*0.,                                           -1.00000,  AuH
     *  1.0, 9*0.,                                           -1.00000,  HgH
     *  1.0, 9*0.,                                           -1.00000/  TlH
      DATA P06/
     *  1.0, 9*0.,                                           -1.00000,  PbH
     *  1.0, 9*0.,                                           -1.00000,  BiH
     *  1.0, 9*0.,                                           -1.00000,  HeH+
     *  1.0, 9*0.,                                           -1.00000,  BeH+
     *  1.0, 9*0.,                                           -1.00000,  CH+
     *  1.0, 9*0.,                                           -1.00000,  NH+
     *  1.0, 9*0.,                                           -1.00000,  OH+
     *  1.0, 9*0.,                                           -1.00000,  HF+
     *  1.0, 9*0.,                                           -1.00000,  NeH+
     *  1.0, 9*0.,                                           -1.00000,  MgH+
     *  1.0, 9*0.,                                           -1.00000,  AlH+
     *  1.0, 9*0.,                                           -1.00000,  SiH+
     *  1.0, 9*0.,                                           -1.00000,  PH+
     *  1.0, 9*0.,                                           -1.00000,  SH+
     *  1.0, 9*0.,                                           -1.00000,  HCl+
     *  1.0, 9*0.,                                           -1.00000,  ZnH+
     *  1.0, 9*0.,                                           -1.00000,  HBr+
     *  1.0, 9*0.,                                           -1.00000,  CdH+
     *  1.0, 9*0.,                                           -1.00000,  HgH+
     *  1.0, 9*0.,                                           -1.00000,  CF
     *  1.0, 9*0.,                                           -1.00000,  CP
     *  1.0, 9*0.,                                           -1.00000,  CCl
     *  1.0, 9*0.,                                           -1.00000,  CSe
     *  1.0, 9*0.,                                           -1.00000,  CBr
     *  1.0, 9*0.,                                           -1.00000,  RhC
     *  1.0, 9*0.,                                           -1.00000,  IrC
     *  1.0, 9*0.,                                           -1.00000,  PtC
     *  1.0, 9*0.,                                           -1.00000,  CN+
     *  1.0, 9*0.,                                           -1.00000,  CO+
     *  1.0, 9*0.,                                           -1.00000,  BN
     *  1.0, 9*0.,                                           -1.00000,  NF
     *  1.0, 9*0.,                                           -1.00000,  AlN
     *  1.0, 9*0.,                                           -1.00000,  PN
     *  1.0, 9*0.,                                           -1.00000,  NCl
     *  1.0, 9*0.,                                           -1.00000,  TiN
     *  1.0, 9*0.,                                           -1.00000,  AsN
     *  1.0, 9*0.,                                           -1.00000,  SeN
     *  1.0, 9*0.,                                           -1.00000,  ZrN
     *  1.0, 9*0.,                                           -1.00000,  NS+
     *  1.0, 9*0.,                                           -1.00000,  LiO
     *  1.0, 9*0.,                                           -1.00000,  BeO
     *  1.0, 9*0.,                                           -1.00000,  FO
     *  1.0, 9*0.,                                           -1.00000,  NaO
     *  1.0, 9*0.,                                           -1.00000,  PO
     *  1.0, 9*0.,                                           -1.00000,  ClO
     *  1.0, 9*0.,                                           -1.00000,  KO
     *  1.0, 9*0.,                                           -1.00000,  CaO
     *  1.0, 9*0.,                                           -1.00000,  MnO
     *  1.0, 9*0.,                                           -1.00000,  NiO
     *  1.0, 9*0.,                                           -1.00000/  CuO
      DATA P07/
     *  1.0, 9*0.,                                           -1.00000,  GaO
     *  1.0, 9*0.,                                           -1.00000,  GeO
     *  1.0, 9*0.,                                           -1.00000,  AsO
     *  1.0, 9*0.,                                           -1.00000,  SeO
     *  1.0, 9*0.,                                           -1.00000,  BrO
     *  1.0, 9*0.,                                           -1.00000,  RbO
     *  1.0, 9*0.,                                           -1.00000,  SrO
     *  1.0, 9*0.,                                           -1.00000,  NbO
     *  1.0, 9*0.,                                           -1.00000,  InO
     *  1.0, 9*0.,                                           -1.00000,  SnO
     *  1.0, 9*0.,                                           -1.00000,  SbO
     *  1.0, 9*0.,                                           -1.00000,  TeO
     *  1.0, 9*0.,                                           -1.00000,  IO
     *  1.0, 9*0.,                                           -1.00000,  BaO
     *  1.0, 9*0.,                                           -1.00000,  TbO
     *  1.0, 9*0.,                                           -1.00000,  LuO
     *  1.0, 9*0.,                                           -1.00000,  HfO
     *  1.0, 9*0.,                                           -1.00000,  TaO
     *  1.0, 9*0.,                                           -1.00000,  WO
     *  1.0, 9*0.,                                           -1.00000,  PtO
     *  1.0, 9*0.,                                           -1.00000,  PbO
     *  1.0, 9*0.,                                           -1.00000,  BiO
     *  1.0, 9*0.,                                           -1.00000,  ThO
     *  1.0, 9*0.,                                           -1.00000,  BO+
     *  1.0, 9*0.,                                           -1.00000,  SiO+
     *  1.0, 9*0.,                                           -1.00000,  PO+
     *  1.0, 9*0.,                                           -1.00000,  SO+
     *  1.0, 9*0.,                                           -1.00000,  AsO+
     *  1.0, 9*0.,                                           -1.00000,  TaO+
     *  1.0, 9*0.,                                           -1.00000,  LiF
     *  1.0, 9*0.,                                           -1.00000,  BeF
     *  1.0, 9*0.,                                           -1.00000,  BF
     *  1.0, 9*0.,                                           -1.00000,  NaF
     *  1.0, 9*0.,                                           -1.00000,  MgF
     *  1.0, 9*0.,                                           -1.00000,  PF
     *  1.0, 9*0.,                                           -1.00000,  SF
     *  1.0, 9*0.,                                           -1.00000,  KF
     *  1.0, 9*0.,                                           -1.00000,  ScF
     *  1.0, 9*0.,                                           -1.00000,  MnF
     *  1.0, 9*0.,                                           -1.00000,  NiF
     *  1.0, 9*0.,                                           -1.00000,  CuF
     *  1.0, 9*0.,                                           -1.00000,  ZnF
     *  1.0, 9*0.,                                           -1.00000,  GaF
     *  1.0, 9*0.,                                           -1.00000,  GeF
     *  1.0, 9*0.,                                           -1.00000,  AsF
     *  1.0, 9*0.,                                           -1.00000,  SeF
     *  1.0, 9*0.,                                           -1.00000,  BrF
     *  1.0, 9*0.,                                           -1.00000,  RbF
     *  1.0, 9*0.,                                           -1.00000,  SrF
     *  1.0, 9*0.,                                           -1.00000/  YF
      DATA P08/
     *  1.0, 9*0.,                                           -1.00000,  AgF
     *  1.0, 9*0.,                                           -1.00000,  CdF
     *  1.0, 9*0.,                                           -1.00000,  InF
     *  1.0, 9*0.,                                           -1.00000,  SnF
     *  1.0, 9*0.,                                           -1.00000,  SbF
     *  1.0, 9*0.,                                           -1.00000,  IF
     *  1.0, 9*0.,                                           -1.00000,  CsF
     *  1.0, 9*0.,                                           -1.00000,  BaF
     *  1.0, 9*0.,                                           -1.00000,  LaF
     *  1.0, 9*0.,                                           -1.00000,  HoF
     *  1.0, 9*0.,                                           -1.00000,  YbF
     *  1.0, 9*0.,                                           -1.00000,  LuF
     *  1.0, 9*0.,                                           -1.00000,  HgF
     *  1.0, 9*0.,                                           -1.00000,  TlF
     *  1.0, 9*0.,                                           -1.00000,  PbF
     *  1.0, 9*0.,                                           -1.00000,  LiNa
     *  1.0, 9*0.,                                           -1.00000,  AsP
     *  1.0, 9*0.,                                           -1.00000,  SbP
     *  1.0, 9*0.,                                           -1.00000,  BeS
     *  1.0, 9*0.,                                           -1.00000,  BS
     *  1.0, 9*0.,                                           -1.00000,  PS
     *  1.0, 9*0.,                                           -1.00000,  CaS
     *  1.0, 9*0.,                                           -1.00000,  ScS
     *  1.0, 9*0.,                                           -1.00000,  CrS
     *  1.0, 9*0.,                                           -1.00000,  CuS
     *  1.0, 9*0.,                                           -1.00000,  GeS
     *  1.0, 9*0.,                                           -1.00000,  AsS
     *  1.0, 9*0.,                                           -1.00000,  SeS
     *  1.0, 9*0.,                                           -1.00000,  SrS
     *  1.0, 9*0.,                                           -1.00000,  YS
     *  1.0, 9*0.,                                           -1.00000,  SnS
     *  1.0, 9*0.,                                           -1.00000,  TeS
     *  1.0, 9*0.,                                           -1.00000,  BaS
     *  1.0, 9*0.,                                           -1.00000,  LaS
     *  1.0, 9*0.,                                           -1.00000,  PbS
     *  1.0, 9*0.,                                           -1.00000,  BiS
     *  1.0, 9*0.,                                           -1.00000,  BeCl
     *  1.0, 9*0.,                                           -1.00000,  BCl
     *  1.0, 9*0.,                                           -1.00000,  MgCl
     *  1.0, 9*0.,                                           -1.00000,  SiCl
     *  1.0, 9*0.,                                           -1.00000,  PCl
     *  1.0, 9*0.,                                           -1.00000,  ScCl
     *  1.0, 9*0.,                                           -1.00000,  MnCl
     *  1.0, 9*0.,                                           -1.00000,  FeCl
     *  1.0, 9*0.,                                           -1.00000,  CuCl
     *  1.0, 9*0.,                                           -1.00000,  ZnCl
     *  1.0, 9*0.,                                           -1.00000,  GaCl
     *  1.0, 9*0.,                                           -1.00000,  GeCl
     *  1.0, 9*0.,                                           -1.00000,  AsCl
     *  1.0, 9*0.,                                           -1.00000/  SeCl
      DATA P09/
     *  1.0, 9*0.,                                           -1.00000,  BrCl
     *  1.0, 9*0.,                                           -1.00000,  RbCl
     *  1.0, 9*0.,                                           -1.00000,  SrCl
     *  1.0, 9*0.,                                           -1.00000,  YCl
     *  1.0, 9*0.,                                           -1.00000,  AgCl
     *  1.0, 9*0.,                                           -1.00000,  CdCl
     *  1.0, 9*0.,                                           -1.00000,  InCl
     *  1.0, 9*0.,                                           -1.00000,  SnCl
     *  1.0, 9*0.,                                           -1.00000,  SbCl
     *  1.0, 9*0.,                                           -1.00000,  ICl
     *  1.0, 9*0.,                                           -1.00000,  CsCl
     *  1.0, 9*0.,                                           -1.00000,  BaCl
     *  1.0, 9*0.,                                           -1.00000,  YbCl
     *  1.0, 9*0.,                                           -1.00000,  AuCl
     *  1.0, 9*0.,                                           -1.00000,  HgCl
     *  1.0, 9*0.,                                           -1.00000,  TlCl
     *  1.0, 9*0.,                                           -1.00000,  PbCl
     *  1.0, 9*0.,                                           -1.00000,  AlSe
     *  1.0, 9*0.,                                           -1.00000,  SiSe
     *  1.0, 9*0.,                                           -1.00000,  GeSe
     *  1.0, 9*0.,                                           -1.00000,  KBr
     *  1.0, 9*0.,                                           -1.00000,  SiTe
     *  1.0, 9*0.,                                           -1.00000,  GeTe
     *  1.0, 9*0.,                                           -1.00000/  KI
C
C
C Try to find the input speicies name (SPNAME) in the list (SPLIST) of
C species for which we have equilibrium constant coefficients. Note that
C the index is stored in a new variable J, rather than using the loop
C variable I, because some optimizers don't save the loop variable after
C normal termination of the loop.
C
      DO 1 I=1,MSPEC
      J=I
      IF(SPLIST(J).EQ.SPNAME) GO TO 2
   1  CONTINUE
C
C Fall through to here, if requested molecule was not in SPLIST.
C Print a warning, but return anyway.
C
      WRITE(*,*) 'MOLCON: Don''t have the equilibrium constant for ',
     *           'molecule: "', SPNAME, '"'
      EQK =1.D20
      PART=1.D0
      RETURN
C
C Calculate independent variable for polynomial expansions.
C Note that the polynomial expansions in Sauval & Tatum (1984) and Irwin
C (1987,1988) are in terms of log10(5040/T), not log10(5039.7475/T), but
C the more accurate value of 5039.7475 should be used in converting the
C partition function into an equilibrium constant.
C
   2  TLIM=MAX(1250.,T)
      TH=5040.D0/TLIM
      LOGTH=LOG10(TH)
C
C Construct equilibrium constant from polynomial coefficients and
C dissociation constant. A "+1" term at the end would convert from
C pascals (i.e. N/m/m as in Sauval) to dynes/cm/cm.
C
c      if (t.lt.1600) logth=log10(5040.0/1600.0)
c      if (t.gt.7730) logth=log10(5040.0/7730.0)
      EQK=COEF(2,J)+LOGTH*(COEF(3,J)+LOGTH*(COEF(4,J)+
     &              LOGTH*(COEF(5,J)+LOGTH*(COEF(6,J)+
     &              LOGTH*(COEF(7,J))))))
     &             -TH*COEF(1,J)
C    &             +1.0D0
      EQK =10.D0**EQK
C
C Just for the reference, the relation between partition functions
C and equilibrium constant:
C
C            P(A)*P(B)*...      N(A)*N(B)*...
C K(AB...) = ------------- = kT-------------- =
C              P(AB...)           N(AB...)
C
C             2*pi*kT 3/2    M(A)*M(B)*... 3/2   Q(A)*Q(B)*...
C       = kT*(-------)    * (-------------)    * ------------- * exp(-D(AB)/kT)
C               h^2            M(AB...)           Q(AB...)
C
C where, K - equilibrium constant, Q - partition functions, M - masses
C        P - partial pressures, N - number densities, T - temperature,
C        D - complete dissociation energy, h - plank constant. Remember
C        to use masses in grams (1 amu = 1.660540E-24 g) and energy in
C        ergs (1 eV = 1.60219E-12 ergs). Also, k = 1.38065E-16 erg/K,
C        h = 6.626076E-27 erg s, and pi = 3.1415926536.
C
C Construct partition function from polynomial coefficients.
C
      PART=PCOEF(NPCOEF-1,J)
      DO 3 I=NPCOEF-2,1,-1
      PART=LOGTH*PART+PCOEF(I,J)
    3 CONTINUE
C
C Copy ionization potential
C
      PION=PCOEF(NPCOEF,J)
C
C Calculate equilibrium constant (EQK) from partition function, dissociation
C constant, and other information passed into subroutine. The constants used
C are:  79.733501 = 1.5*log10(2*pi/h/h)  [in cgs units] and
C      -15.859914 = alog10(k)            [in cgs units].
C       5039.7475 = alog10(e)*k*(eV/erg)
C
c      EQK_ST=(NTOT-1)*(79.733501D0+2.5D0*(LOG10(TLIM)-15.859914D0))+
c     &       1.5D0*RATIOM+QPRD-PART-COEF(1,J)*5039.7475D0/TLIM
C
C Convert equilibrium constant and partition function from logarithms.
C
c      EQK_ST=10.D0**EQK_ST
      PART=10.D0**PART
C
C Check if there is relevant data in Paul Barklem's tables
C
      CALL KP_Q_SPLN(SPNAME,T,Qm_spln,Kp_spln,BARKLEM)
      IF(BARKLEM) THEN
c        EQK =Kp_spln-COEF(1,J)*5039.7475D0/TLIM
        EQK =Kp_spln-COEF(1,J)*5040.D0/T
        EQK =10.D0**EQK
        PART=10.D0**Qm_spln
      ENDIF
      if(spname.eq.'H3O+') then
        EQK_ST=(NTOT-1)*(79.733501D0+2.5D0*(LOG10(T)-15.859914D0))+
     &         1.5D0*RATIOM+QPRD-PART-COEF(1,J)*5039.7475D0/T
        EQK=10.D0**EQK_ST
      endif
c      write(*,'(''cMOLCON:'',F10.1,A9,5G13.6)') T,SPNAME,EQK,
c     &                             PART,BARKLEM
c      if(spname.eq.'NO') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'C3') write(*,'(a,f10.2,1p6e16.8)')
c     &   spname,t , eqk, eqk_st, part, TH, LOGTH, TLIM
c      if(spname.eq.'H3O+') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'SiS') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'NO') write(*,'(a,f10.2,1p3e16.8,L)')
c     &   spname,t , eqk, eqk_st, part,barklem
c      if(spname.eq.'CH') write(*,'(a,f10.2,1p5e16.8,L)')
c     &   spname,t , eqk, eqk_st, part,COEF(1,J),Kp_spln,barklem
c      if(spname.eq.'CH-') write(*,'(a,f10.2,1p5e16.8,L)')
c     &   spname,t , eqk, eqk_st, part,COEF(1,J),Kp_spln,barklem
c      if(spname.eq.'CH-') write(*,'(a,f10.2,1p3e14.6,i3,1p2e14.6,L)')
c     &   spname,t , eqk, eqk_st, part,NTOT,QPRD,RATIOM,BARKLEM
c      if(spname.eq.'H2') write(*,'(a,f10.2,1p3e14.6,i3,1p2e14.6)')
c     &   spname,t , eqk, eqk_st, part,NTOT,Kp_spln,COEF(1,J)*5040.D0/T
c
c Don't use EQK_ST based on partition function - use direct fit to EQK.
c
c      EQK=EQK_ST
C
C Done.
C
      RETURN
      END
C---------------------- Start of Barklem subroutines ------------------------
C----------------------- End of Berklem subroutines ------------------------
      SUBROUTINE SPL_INIT(X,Y,Y2,U,N)
C
C  Computes second derivative approximations for cubic spline interpolation
C
      IMPLICIT NONE
      INTEGER N
      REAL*8 X(N),Y(N),Y2(N),U(N)
      INTEGER I
      REAL*8 SIG,P,YY1,YY2,YY3
C
C  Natural lower boundary condition
C
      Y2(1)=0.D0
      U(1)=0.D0
      DO 1 I=2,N-1
      SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
      P=SIG*Y2(I-1)+2.D0
      Y2(I)=(SIG-1.D0)/P
      YY1=Y(I-1)
      YY2=Y(I  )
      YY3=Y(I+1)
      U(I)=(6.D0*((YY3-YY2)/(X(I+1)-X(I))-(YY2-YY1)/
     /     (X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
  1   CONTINUE
C
C  Natural upper boundary condition
C
      Y2(N)=0.D0
      DO 2 I=N-1,1,-1
      Y2(I)=Y2(I)*Y2(I+1)+U(I)
  2   CONTINUE
C
      RETURN
      END

      REAL*8 FUNCTION SPL_INTERP(KLO,KHI,XA,YA,Y2A,N,X)
C
C  Performs cubic spline interpolation
C
      IMPLICIT NONE
      INTEGER KLO,KHI,N
      REAL*8 XA(N),YA(N),Y2A(N),X
      REAL*8 A,B,H,Y1,Y2
C
      H=XA(KHI)-XA(KLO)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y1=YA(KLO)
      Y2=YA(KHI)
      SPL_INTERP=A*Y1+B*Y2+((A*A-1.D0)*A*Y2A(KLO)+
     +                      (B*B-1.D0)*B*Y2A(KHI))*(H*H)/6.D0
C
      RETURN
      END

      SUBROUTINE XSAHA(IEL,TT,XNELEC,XNATOM,MAXION,POTI,FRCT,MODE)
C
C     MODE=1 returns ionization fractions/partition functions
C     MODE=2 returns ionization fractions
C     MODE=3 returns partition functions
C     MODE=4 returns total number of electrons produced
C     MODE=5 returns in MAXION(!) the number of ionization stages
C            available in XSAHA
C
C     ALL OF THE ABOVE IS FOR ALL IONIZATION STAGES UP TO MAXION
C
C  Parameters:
C     IEL    - (input) element atomic number (Hydrogen: 1)
C     TT     - (input) temperature (Kelvins)
C     XNELEC - (input) electron number density (cm^-3)
C     XNATOM - (input) particle number density (excluding electrons) (cm^-3)
C     MAXION - (input/output) size of the output arrays
C     POTI   - (output array of MAXION) ionization potential (eV)
C     FRCT   - (output array of MAXION) results according to MODE
C     MODE   - (input) see above
C
      INTEGER ELESIZ,IONSIZ,IEL
      PARAMETER (ELESIZ=100,IONSIZ=6)
      DOUBLE PRECISION FFF(IONSIZ),FEXARG,FRCT(MAXION),CF
      REAL IP(IONSIZ),PART(IONSIZ),POTLO(IONSIZ),SCALE(4),
     *     POTI(MAXION),TT
      INTEGER LOCZ(ELESIZ+1)
      LOGICAL FIRST

      INTEGER SIZ_H ,SIZ_He,SIZ_Li,SIZ_Be,SIZ_B ,SIZ_C ,SIZ_N ,SIZ_O ,
     1        SIZ_F ,SIZ_Ne,SIZ_Na,SIZ_Mg,SIZ_Al,SIZ_Si,SIZ_P ,SIZ_S ,
     2        SIZ_Cl,SIZ_Ar,SIZ_K ,SIZ_Ca,SIZ_Sc,SIZ_Ti,SIZ_V ,SIZ_Cr,
     3        SIZ_Mn,SIZ_Fe,SIZ_Co,SIZ_Ni,SIZ_Cu,SIZ_Zn,SIZ_Ga,SIZ_Ge,
     4        SIZ_As,SIZ_Se,SIZ_Br,SIZ_Kr,SIZ_Rb,SIZ_Sr,SIZ_Y ,SIZ_Zr,
     5        SIZ_Nb,SIZ_Mo,SIZ_Tc,SIZ_Ru,SIZ_Rh,SIZ_Pd,SIZ_Ag,SIZ_Cd,
     6        SIZ_In,SIZ_Sn,SIZ_Sb,SIZ_Te,SIZ_I ,SIZ_Xe,SIZ_Cs,SIZ_Ba,
     7        SIZ_La,SIZ_Ce,SIZ_Pr,SIZ_Nd,SIZ_Pm,SIZ_Sm,SIZ_Eu,SIZ_Gd,
     8        SIZ_Tb,SIZ_Dy,SIZ_Ho,SIZ_Er,SIZ_Tm,SIZ_Yb,SIZ_Lu,SIZ_Hf,
     9        SIZ_Ta,SIZ_W ,SIZ_Re,SIZ_Os,SIZ_Ir,SIZ_Pt,SIZ_Au,SIZ_Hg,
     A        SIZ_Tl,SIZ_Pb,SIZ_Bi,SIZ_Po,SIZ_At,SIZ_Rn,SIZ_Fr,SIZ_Ra,
     B        SIZ_Ac,SIZ_Th,SIZ_Pa,SIZ_U ,SIZ_Np,SIZ_Pu,SIZ_Am,SIZ_Cm,
     C        SIZ_Bk,SIZ_Cf,SIZ_Es
      INTEGER OFF_H ,OFF_He,OFF_Li,OFF_Be,OFF_B ,OFF_C ,OFF_N ,OFF_O ,
     1        OFF_F ,OFF_Ne,OFF_Na,OFF_Mg,OFF_Al,OFF_Si,OFF_P ,OFF_S ,
     2        OFF_Cl,OFF_Ar,OFF_K ,OFF_Ca,OFF_Sc,OFF_Ti,OFF_V ,OFF_Cr,
     3        OFF_Mn,OFF_Fe,OFF_Co,OFF_Ni,OFF_Cu,OFF_Zn,OFF_Ga,OFF_Ge,
     4        OFF_As,OFF_Se,OFF_Br,OFF_Kr,OFF_Rb,OFF_Sr,OFF_Y ,OFF_Zr,
     5        OFF_Nb,OFF_Mo,OFF_Tc,OFF_Ru,OFF_Rh,OFF_Pd,OFF_Ag,OFF_Cd,
     6        OFF_In,OFF_Sn,OFF_Sb,OFF_Te,OFF_I ,OFF_Xe,OFF_Cs,OFF_Ba,
     7        OFF_La,OFF_Ce,OFF_Pr,OFF_Nd,OFF_Pm,OFF_Sm,OFF_Eu,OFF_Gd,
     8        OFF_Tb,OFF_Dy,OFF_Ho,OFF_Er,OFF_Tm,OFF_Yb,OFF_Lu,OFF_Hf,
     9        OFF_Ta,OFF_W ,OFF_Re,OFF_Os,OFF_Ir,OFF_Pt,OFF_Au,OFF_Hg,
     A        OFF_Tl,OFF_Pb,OFF_Bi,OFF_Po,OFF_At,OFF_Rn,OFF_Fr,OFF_Ra,
     B        OFF_Ac,OFF_Th,OFF_Pa,OFF_U ,OFF_Np,OFF_Pu,OFF_Am,OFF_Cm,
     C        OFF_Bk,OFF_Cf,OFF_Es
C
C In order to add data for another ionization stage to a particular element
C one has to do two things: increase the value of SIZ_<elname> and add the
C data line(s) in the DATA NNN_<elname>
C
      PARAMETER (SIZ_H = 2, OFF_H = 1)
      INTEGER NNN_H (8*SIZ_H )
      PARAMETER (SIZ_He= 3, OFF_He=OFF_H +SIZ_H  )
      INTEGER NNN_He(8*SIZ_He)
      PARAMETER (SIZ_Li= 4, OFF_Li=OFF_He+SIZ_He)
      INTEGER NNN_Li(8*SIZ_Li)
      PARAMETER (SIZ_Be= 4, OFF_Be=OFF_Li+SIZ_Li)
      INTEGER NNN_Be(8*SIZ_Be)
      PARAMETER (SIZ_B = 4, OFF_B =OFF_Be+SIZ_Be)
      INTEGER NNN_B (8*SIZ_B )
      PARAMETER (SIZ_C = 6, OFF_C =OFF_B +SIZ_B )
      INTEGER NNN_C (8*SIZ_C )
      PARAMETER (SIZ_N = 6, OFF_N =OFF_C +SIZ_C )
      INTEGER NNN_N (8*SIZ_N )
      PARAMETER (SIZ_O = 6, OFF_O =OFF_N +SIZ_N )
      INTEGER NNN_O (8*SIZ_O )
      PARAMETER (SIZ_F = 6, OFF_F =OFF_O +SIZ_O )
      INTEGER NNN_F (8*SIZ_F )
      PARAMETER (SIZ_Ne= 6, OFF_Ne=OFF_F +SIZ_F )
      INTEGER NNN_Ne(8*SIZ_Ne)
      PARAMETER (SIZ_Na= 6, OFF_Na=OFF_Ne+SIZ_Ne)
      INTEGER NNN_Na(8*SIZ_Na)
      PARAMETER (SIZ_Mg= 6, OFF_Mg=OFF_Na+SIZ_Na)
      INTEGER NNN_Mg(8*SIZ_Mg)
      PARAMETER (SIZ_Al= 6, OFF_Al=OFF_Mg+SIZ_Mg)
      INTEGER NNN_Al(8*SIZ_Al)
      PARAMETER (SIZ_Si= 6, OFF_Si=OFF_Al+SIZ_Al)
      INTEGER NNN_Si(8*SIZ_Si)
      PARAMETER (SIZ_P = 6, OFF_P =OFF_Si+SIZ_Si)
      INTEGER NNN_P (8*SIZ_P )
      PARAMETER (SIZ_S = 6, OFF_S =OFF_P +SIZ_P )
      INTEGER NNN_S (8*SIZ_S )
      PARAMETER (SIZ_Cl= 5, OFF_Cl=OFF_S +SIZ_S )
      INTEGER NNN_Cl(8*SIZ_Cl)
      PARAMETER (SIZ_Ar= 5, OFF_Ar=OFF_Cl+SIZ_Cl)
      INTEGER NNN_Ar(8*SIZ_Ar)
      PARAMETER (SIZ_K = 5, OFF_K =OFF_Ar+SIZ_Ar)
      INTEGER NNN_K (8*SIZ_K )
      PARAMETER (SIZ_Ca= 5, OFF_Ca=OFF_K +SIZ_K )
      INTEGER NNN_Ca(8*SIZ_Ca)
      PARAMETER (SIZ_Sc= 5, OFF_Sc=OFF_Ca+SIZ_Ca)
      INTEGER NNN_Sc(8*SIZ_Sc)
      PARAMETER (SIZ_Ti= 5, OFF_Ti=OFF_Sc+SIZ_Sc)
      INTEGER NNN_Ti(8*SIZ_Ti)
      PARAMETER (SIZ_V = 5, OFF_V =OFF_Ti+SIZ_Ti)
      INTEGER NNN_V (8*SIZ_V )
      PARAMETER (SIZ_Cr= 5, OFF_Cr=OFF_V +SIZ_V )
      INTEGER NNN_Cr(8*SIZ_Cr)
      PARAMETER (SIZ_Mn= 5, OFF_Mn=OFF_Cr+SIZ_Cr)
      INTEGER NNN_Mn(8*SIZ_Mn)
      PARAMETER (SIZ_Fe= 5, OFF_Fe=OFF_Mn+SIZ_Mn)
      INTEGER NNN_Fe(8*SIZ_Fe)
      PARAMETER (SIZ_Co= 5, OFF_Co=OFF_Fe+SIZ_Fe)
      INTEGER NNN_Co(8*SIZ_Co)
      PARAMETER (SIZ_Ni= 5, OFF_Ni=OFF_Co+SIZ_Co)
      INTEGER NNN_Ni(8*SIZ_Ni)
      PARAMETER (SIZ_Cu= 3, OFF_Cu=OFF_Ni+SIZ_Ni)
      INTEGER NNN_Cu(8*SIZ_Cu)
      PARAMETER (SIZ_Zn= 3, OFF_Zn=OFF_Cu+SIZ_Cu)
      INTEGER NNN_Zn(8*SIZ_Zn)
      PARAMETER (SIZ_Ga= 3, OFF_Ga=OFF_Zn+SIZ_Zn)
      INTEGER NNN_Ga(8*SIZ_Ga)
      PARAMETER (SIZ_Ge= 3, OFF_Ge=OFF_Ga+SIZ_Ga)
      INTEGER NNN_Ge(8*SIZ_Ge)
      PARAMETER (SIZ_As= 3, OFF_As=OFF_Ge+SIZ_Ge)
      INTEGER NNN_As(8*SIZ_As)
      PARAMETER (SIZ_Se= 3, OFF_Se=OFF_As+SIZ_As)
      INTEGER NNN_Se(8*SIZ_Se)
      PARAMETER (SIZ_Br= 3, OFF_Br=OFF_Se+SIZ_Se)
      INTEGER NNN_Br(8*SIZ_Br)
      PARAMETER (SIZ_Kr= 3, OFF_Kr=OFF_Br+SIZ_Br)
      INTEGER NNN_Kr(8*SIZ_Kr)
      PARAMETER (SIZ_Rb= 3, OFF_Rb=OFF_Kr+SIZ_Kr)
      INTEGER NNN_Rb(8*SIZ_Rb)
      PARAMETER (SIZ_Sr= 3, OFF_Sr=OFF_Rb+SIZ_Rb)
      INTEGER NNN_Sr(8*SIZ_Sr)
      PARAMETER (SIZ_Y = 3, OFF_Y =OFF_Sr+SIZ_Sr)
      INTEGER NNN_Y (8*SIZ_Y )
      PARAMETER (SIZ_Zr= 3, OFF_Zr=OFF_Y +SIZ_Y )
      INTEGER NNN_Zr(8*SIZ_Zr)
      PARAMETER (SIZ_Nb= 3, OFF_Nb=OFF_Zr+SIZ_Zr)
      INTEGER NNN_Nb(8*SIZ_Nb)
      PARAMETER (SIZ_Mo= 3, OFF_Mo=OFF_Nb+SIZ_Nb)
      INTEGER NNN_Mo(8*SIZ_Mo)
      PARAMETER (SIZ_Tc= 3, OFF_Tc=OFF_Mo+SIZ_Mo)
      INTEGER NNN_Tc(8*SIZ_Tc)
      PARAMETER (SIZ_Ru= 3, OFF_Ru=OFF_Tc+SIZ_Tc)
      INTEGER NNN_Ru(8*SIZ_Ru)
      PARAMETER (SIZ_Rh= 3, OFF_Rh=OFF_Ru+SIZ_Ru)
      INTEGER NNN_Rh(8*SIZ_Rh)
      PARAMETER (SIZ_Pd= 3, OFF_Pd=OFF_Rh+SIZ_Rh)
      INTEGER NNN_Pd(8*SIZ_Pd)
      PARAMETER (SIZ_Ag= 3, OFF_Ag=OFF_Pd+SIZ_Pd)
      INTEGER NNN_Ag(8*SIZ_Ag)
      PARAMETER (SIZ_Cd= 3, OFF_Cd=OFF_Ag+SIZ_Ag)
      INTEGER NNN_Cd(8*SIZ_Cd)
      PARAMETER (SIZ_In= 3, OFF_In=OFF_Cd+SIZ_Cd)
      INTEGER NNN_In(8*SIZ_In)
      PARAMETER (SIZ_Sn= 3, OFF_Sn=OFF_In+SIZ_In)
      INTEGER NNN_Sn(8*SIZ_Sn)
      PARAMETER (SIZ_Sb= 3, OFF_Sb=OFF_Sn+SIZ_Sn)
      INTEGER NNN_Sb(8*SIZ_Sb)
      PARAMETER (SIZ_Te= 3, OFF_Te=OFF_Sb+SIZ_Sb)
      INTEGER NNN_Te(8*SIZ_Te)
      PARAMETER (SIZ_I = 3, OFF_I =OFF_Te+SIZ_Te)
      INTEGER NNN_I (8*SIZ_I )
      PARAMETER (SIZ_Xe= 3, OFF_Xe=OFF_I +SIZ_I )
      INTEGER NNN_Xe(8*SIZ_Xe)
      PARAMETER (SIZ_Cs= 3, OFF_Cs=OFF_Xe+SIZ_Xe)
      INTEGER NNN_Cs(8*SIZ_Cs)
      PARAMETER (SIZ_Ba= 3, OFF_Ba=OFF_Cs+SIZ_Cs)
      INTEGER NNN_Ba(8*SIZ_Ba)
      PARAMETER (SIZ_La= 3, OFF_La=OFF_Ba+SIZ_Ba)
      INTEGER NNN_La(8*SIZ_La)
      PARAMETER (SIZ_Ce= 4, OFF_Ce=OFF_La+SIZ_La)
      INTEGER NNN_Ce(8*SIZ_Ce)
      PARAMETER (SIZ_Pr= 4, OFF_Pr=OFF_Ce+SIZ_Ce)
      INTEGER NNN_Pr(8*SIZ_Pr)
      PARAMETER (SIZ_Nd= 4, OFF_Nd=OFF_Pr+SIZ_Pr)
      INTEGER NNN_Nd(8*SIZ_Nd)
      PARAMETER (SIZ_Pm= 3, OFF_Pm=OFF_Nd+SIZ_Nd)
      INTEGER NNN_Pm(8*SIZ_Pm)
      PARAMETER (SIZ_Sm= 3, OFF_Sm=OFF_Pm+SIZ_Pm)
      INTEGER NNN_Sm(8*SIZ_Sm)
      PARAMETER (SIZ_Eu= 4, OFF_Eu=OFF_Sm+SIZ_Sm)
      INTEGER NNN_Eu(8*SIZ_Eu)
      PARAMETER (SIZ_Gd= 3, OFF_Gd=OFF_Eu+SIZ_Eu)
      INTEGER NNN_Gd(8*SIZ_Gd)
      PARAMETER (SIZ_Tb= 3, OFF_Tb=OFF_Gd+SIZ_Gd)
      INTEGER NNN_Tb(8*SIZ_Tb)
      PARAMETER (SIZ_Dy= 3, OFF_Dy=OFF_Tb+SIZ_Tb)
      INTEGER NNN_Dy(8*SIZ_Dy)
      PARAMETER (SIZ_Ho= 3, OFF_Ho=OFF_Dy+SIZ_Dy)
      INTEGER NNN_Ho(8*SIZ_Ho)
      PARAMETER (SIZ_Er= 3, OFF_Er=OFF_Ho+SIZ_Ho)
      INTEGER NNN_Er(8*SIZ_Er)
      PARAMETER (SIZ_Tm= 3, OFF_Tm=OFF_Er+SIZ_Er)
      INTEGER NNN_Tm(8*SIZ_Tm)
      PARAMETER (SIZ_Yb= 3, OFF_Yb=OFF_Tm+SIZ_Tm)
      INTEGER NNN_Yb(8*SIZ_Yb)
      PARAMETER (SIZ_Lu= 3, OFF_Lu=OFF_Yb+SIZ_Yb)
      INTEGER NNN_Lu(8*SIZ_Lu)
      PARAMETER (SIZ_Hf= 3, OFF_Hf=OFF_Lu+SIZ_Lu)
      INTEGER NNN_Hf(8*SIZ_Hf)
      PARAMETER (SIZ_Ta= 3, OFF_Ta=OFF_Hf+SIZ_Hf)
      INTEGER NNN_Ta(8*SIZ_Ta)
      PARAMETER (SIZ_W = 3, OFF_W =OFF_Ta+SIZ_Ta)
      INTEGER NNN_W (8*SIZ_W )
      PARAMETER (SIZ_Re= 3, OFF_Re=OFF_W +SIZ_W )
      INTEGER NNN_Re(8*SIZ_Re)
      PARAMETER (SIZ_Os= 3, OFF_Os=OFF_Re+SIZ_Re)
      INTEGER NNN_Os(8*SIZ_Os)
      PARAMETER (SIZ_Ir= 3, OFF_Ir=OFF_Os+SIZ_Os)
      INTEGER NNN_Ir(8*SIZ_Ir)
      PARAMETER (SIZ_Pt= 3, OFF_Pt=OFF_Ir+SIZ_Ir)
      INTEGER NNN_Pt(8*SIZ_Pt)
      PARAMETER (SIZ_Au= 3, OFF_Au=OFF_Pt+SIZ_Pt)
      INTEGER NNN_Au(8*SIZ_Au)
      PARAMETER (SIZ_Hg= 3, OFF_Hg=OFF_Au+SIZ_Au)
      INTEGER NNN_Hg(8*SIZ_Hg)
      PARAMETER (SIZ_Tl= 3, OFF_Tl=OFF_Hg+SIZ_Hg)
      INTEGER NNN_Tl(8*SIZ_Tl)
      PARAMETER (SIZ_Pb= 3, OFF_Pb=OFF_Tl+SIZ_Tl)
      INTEGER NNN_Pb(8*SIZ_Pb)
      PARAMETER (SIZ_Bi= 3, OFF_Bi=OFF_Pb+SIZ_Pb)
      INTEGER NNN_Bi(8*SIZ_Bi)
      PARAMETER (SIZ_Po= 3, OFF_Po=OFF_Bi+SIZ_Bi)
      INTEGER NNN_Po(8*SIZ_Po)
      PARAMETER (SIZ_At= 3, OFF_At=OFF_Po+SIZ_Po)
      INTEGER NNN_At(8*SIZ_At)
      PARAMETER (SIZ_Rn= 3, OFF_Rn=OFF_At+SIZ_At)
      INTEGER NNN_Rn(8*SIZ_Rn)
      PARAMETER (SIZ_Fr= 3, OFF_Fr=OFF_Rn+SIZ_Rn)
      INTEGER NNN_Fr(8*SIZ_Fr)
      PARAMETER (SIZ_Ra= 3, OFF_Ra=OFF_Fr+SIZ_Fr)
      INTEGER NNN_Ra(8*SIZ_Ra)
      PARAMETER (SIZ_Ac= 3, OFF_Ac=OFF_Ra+SIZ_Ra)
      INTEGER NNN_Ac(8*SIZ_Ac)
      PARAMETER (SIZ_Th= 3, OFF_Th=OFF_Ac+SIZ_Ac)
      INTEGER NNN_Th(8*SIZ_Th)
      PARAMETER (SIZ_Pa= 3, OFF_Pa=OFF_Th+SIZ_Th)
      INTEGER NNN_Pa(8*SIZ_Pa)
      PARAMETER (SIZ_U = 3, OFF_U =OFF_Pa+SIZ_Pa)
      INTEGER NNN_U (8*SIZ_U )
      PARAMETER (SIZ_Np= 3, OFF_Np=OFF_U +SIZ_U )
      INTEGER NNN_Np(8*SIZ_Np)
      PARAMETER (SIZ_Pu= 3, OFF_Pu=OFF_Np+SIZ_Np)
      INTEGER NNN_Pu(8*SIZ_Pu)
      PARAMETER (SIZ_Am= 3, OFF_Am=OFF_Pu+SIZ_Pu)
      INTEGER NNN_Am(8*SIZ_Am)
      PARAMETER (SIZ_Cm= 3, OFF_Cm=OFF_Am+SIZ_Am)
      INTEGER NNN_Cm(8*SIZ_Cm)
      PARAMETER (SIZ_Bk= 3, OFF_Bk=OFF_Cm+SIZ_Cm)
      INTEGER NNN_Bk(8*SIZ_Bk)
      PARAMETER (SIZ_Cf= 3, OFF_Cf=OFF_Bk+SIZ_Bk)
      INTEGER NNN_Cf(8*SIZ_Cf)
      PARAMETER (SIZ_Es= 3, OFF_Es=OFF_Cf+SIZ_Cf)
      INTEGER NNN_Es(8*SIZ_Es)

      PARAMETER (NTABLE=OFF_Es+SIZ_Es-1)
      INTEGER NNNPFN(8,NTABLE)

      EQUIVALENCE (NNNPFN(1,OFF_H ),NNN_H (1))
      EQUIVALENCE (NNNPFN(1,OFF_He),NNN_He(1))
      EQUIVALENCE (NNNPFN(1,OFF_Li),NNN_Li(1))
      EQUIVALENCE (NNNPFN(1,OFF_Be),NNN_Be(1))
      EQUIVALENCE (NNNPFN(1,OFF_B ),NNN_B (1))
      EQUIVALENCE (NNNPFN(1,OFF_C ),NNN_C (1))
      EQUIVALENCE (NNNPFN(1,OFF_N ),NNN_N (1))
      EQUIVALENCE (NNNPFN(1,OFF_O ),NNN_O (1))
      EQUIVALENCE (NNNPFN(1,OFF_F ),NNN_F (1))
      EQUIVALENCE (NNNPFN(1,OFF_Ne),NNN_Ne(1))
      EQUIVALENCE (NNNPFN(1,OFF_Na),NNN_Na(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mg),NNN_Mg(1))
      EQUIVALENCE (NNNPFN(1,OFF_Al),NNN_Al(1))
      EQUIVALENCE (NNNPFN(1,OFF_Si),NNN_Si(1))
      EQUIVALENCE (NNNPFN(1,OFF_P ),NNN_P (1))
      EQUIVALENCE (NNNPFN(1,OFF_S ),NNN_S (1))
      EQUIVALENCE (NNNPFN(1,OFF_Cl),NNN_Cl(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ar),NNN_Ar(1))
      EQUIVALENCE (NNNPFN(1,OFF_K ),NNN_K (1))
      EQUIVALENCE (NNNPFN(1,OFF_Ca),NNN_Ca(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sc),NNN_Sc(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ti),NNN_Ti(1))
      EQUIVALENCE (NNNPFN(1,OFF_V ),NNN_V (1))
      EQUIVALENCE (NNNPFN(1,OFF_Cr),NNN_Cr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mn),NNN_Mn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Fe),NNN_Fe(1))
      EQUIVALENCE (NNNPFN(1,OFF_Co),NNN_Co(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ni),NNN_Ni(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cu),NNN_Cu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Zn),NNN_Zn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ga),NNN_Ga(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ge),NNN_Ge(1))
      EQUIVALENCE (NNNPFN(1,OFF_As),NNN_As(1))
      EQUIVALENCE (NNNPFN(1,OFF_Se),NNN_Se(1))
      EQUIVALENCE (NNNPFN(1,OFF_Br),NNN_Br(1))
      EQUIVALENCE (NNNPFN(1,OFF_Kr),NNN_Kr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rb),NNN_Rb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sr),NNN_Sr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Y ),NNN_Y (1))
      EQUIVALENCE (NNNPFN(1,OFF_Zr),NNN_Zr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Nb),NNN_Nb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mo),NNN_Mo(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tc),NNN_Tc(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ru),NNN_Ru(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rh),NNN_Rh(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pd),NNN_Pd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ag),NNN_Ag(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cd),NNN_Cd(1))
      EQUIVALENCE (NNNPFN(1,OFF_In),NNN_In(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sn),NNN_Sn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sb),NNN_Sb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Te),NNN_Te(1))
      EQUIVALENCE (NNNPFN(1,OFF_I ),NNN_I (1))
      EQUIVALENCE (NNNPFN(1,OFF_Xe),NNN_Xe(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cs),NNN_Cs(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ba),NNN_Ba(1))
      EQUIVALENCE (NNNPFN(1,OFF_La),NNN_La(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ce),NNN_Ce(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pr),NNN_Pr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Nd),NNN_Nd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pm),NNN_Pm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sm),NNN_Sm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Eu),NNN_Eu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Gd),NNN_Gd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tb),NNN_Tb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Dy),NNN_Dy(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ho),NNN_Ho(1))
      EQUIVALENCE (NNNPFN(1,OFF_Er),NNN_Er(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tm),NNN_Tm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Yb),NNN_Yb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Lu),NNN_Lu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Hf),NNN_Hf(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ta),NNN_Ta(1))
      EQUIVALENCE (NNNPFN(1,OFF_W ),NNN_W (1))
      EQUIVALENCE (NNNPFN(1,OFF_Re),NNN_Re(1))
      EQUIVALENCE (NNNPFN(1,OFF_Os),NNN_Os(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ir),NNN_Ir(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pt),NNN_Pt(1))
      EQUIVALENCE (NNNPFN(1,OFF_Au),NNN_Au(1))
      EQUIVALENCE (NNNPFN(1,OFF_Hg),NNN_Hg(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tl),NNN_Tl(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pb),NNN_Pb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Bi),NNN_Bi(1))
      EQUIVALENCE (NNNPFN(1,OFF_Po),NNN_Po(1))
      EQUIVALENCE (NNNPFN(1,OFF_At),NNN_At(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rn),NNN_Rn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Fr),NNN_Fr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ra),NNN_Ra(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ac),NNN_Ac(1))
      EQUIVALENCE (NNNPFN(1,OFF_Th),NNN_Th(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pa),NNN_Pa(1))
      EQUIVALENCE (NNNPFN(1,OFF_U ),NNN_U (1))
      EQUIVALENCE (NNNPFN(1,OFF_Np),NNN_Np(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pu),NNN_Pu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Am),NNN_Am(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cm),NNN_Cm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Bk),NNN_Bk(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cf),NNN_Cf(1))
      EQUIVALENCE (NNNPFN(1,OFF_Es),NNN_Es(1))
      SAVE NNNPFN,LOCZ,SCALE,FIRST,FFF
C      ( 1)( 2)  ( 3)( 4)  ( 5)( 6)  ( 7)( 8)  ( 9)(10)   ( IP )   G Ion REF
      DATA NNN_H/
     1 200020001,200020011,201620881,231228281,378953411, 1359502, 1,00,D+F H  1
     2 100010001,100010001,100010001,100010001,100010001, 1359500, 1,01/G   H  2
      DATA NNN_He/
     1 100010001,100010011,102111241,145022061,363059451, 2458104, 2,00,D+F He 1
     2 200020001,200020071,208524971,382669341,128222452, 5440302, 2,01,D+F He 2
     3 100010001,100010001,100010001,100010001,100010001, 5440300, 2,02/G   He 3
      DATA NNN_Li/
     1 200020011,201220481,212922881,258731081,394251691,  538901, 3,00,D+F Li 1
     2 100010001,100010201,126225521, 67216512,351165562, 7561907, 3,01,D+F Li 2
     3 200020001,200020211,227936571, 69610342,137217102,12241800, 3,02,D+F Li 3
     4 100010001,100010001,100010001,100010001,100010001,12241800, 3,03/G   Li 4
      DATA NNN_Be/
     1 100010051,104311441,131615641,190623681,298037691,  931900, 4,00,AEL Be 1
     2 200120231,211422771,249627631,309034911,398545051, 1820600, 4,01,AEL Be 2
     3 100010001,100010201,126225521, 67216512,351165562,15385000, 4,02,AEL Be 3
     4 200020001,200020011,201220661,223426161,332644691,21765700, 4,03/AEL Be 4
      DATA NNN_B/
     1 600060001,600560281,608761991,637466191,693973361,  829500, 5,00,AEL B  1
     2 100310831,132016901,214226411,315736741,419147071, 2514900, 5,01,AEL B  2
     3 200721061,233526401,297533311,369040481,440747651, 3792000, 5,02,AEL B  3
     4 100010001,100010001,100010001,100010001,100010001,25929800, 5,03/G   B  4
      DATA NNN_C/
     1 893292271, 96110042,105311262,126315202,196126432, 1125508, 6,00,D+F C  1
     2 595060251,620865751,713280191, 95712292,167623542, 2437501, 6,01,D+F C  2
     3 105513201,180324851,341851341, 88416332,296550722, 4787101, 6,02,D+F C  3
     4 204922771,262630421,350941931,494556971,644872001, 6447600, 6,03,D+F C  4
     5 100010001,100010001,100010001,100010001,100010001,39207700, 6,04,G   C  5
     6 200020001,200020001,200020001,200020001,200020001,48998100, 6,05/G   C  6
      DATA NNN_N/
     1 403141851,457051681,594071181, 92913362,203331152, 1452915, 7,00,D+F N  1
     2 919899541,107211512,124914302,182526232,403762662, 2959202, 7,01,D+F N  2
     3 596862721,684177081, 88110342,128317062,239334312, 4742501, 7,02,D+F N  3
     4 112816481,240733751,462068491,116419932,283736822, 7744900, 7,03,D+F N  4
     5 210124681,293634211,391145791,539862151,703178471, 9786200, 7,04,D+F N  5
     6 100010001,100010001,100010001,100010001,100010001,55205700, 7,05/G   N  6
      DATA NNN_O/
     1 874789691,924795711, 99410492,115213492,169022242, 1361307, 8,00,D+F O  1
     2 424151091,622874781, 91312832,221842502, 79914013, 3510711, 8,01,D+F O  2
     3  95610702,118113032,149619922,329761642,101914173, 5488500, 8,02,D+F O  3
     4 603567171,775391141,106612482,143716252,181420032, 7739300, 8,03,D+F O  4
     5 124420321,306943181,606281181,101712232,142916342,11387300, 8,04,D+F O  5
     6 215026541,323137551,421546491,508255151,594863811,13807900, 8,05/AEL O  6
      DATA NNN_F/
     1 575958511,589859231,595860671,636470031,815199581, 1741802, 9,00,D+F F  1
     2 900296401,102610802,113912542,152921152,318348952, 3498003, 9,01,D+F F  2
     3 469162651,791295541,121419552,402686872,154822203, 6264500, 9,02,D+F F  3
     4  99511422,129214572,170523002,320140922,498458762, 8713900, 9,03,D+F F  4
     5 615472711, 87710602,127215002,172919582,218624152,11421300, 9,04,D+F F  5
     6 135324181,377252001,661580261, 94410852,122613672,15711700, 9,05/AEL F  6
      DATA NNN_Ne/
     1 100010001,100010051,105313051,210239461, 74013022, 2155808,10,00,D+F Ne 1
     2 580158751,591759741,642687101,159332652, 64111533, 4106907,10,01,D+F Ne 2
     3  93510272,110411662,127116062,257647882, 75110223, 6350000,10,02,D+F Ne 3
     4 529774371, 94611322,135816202,188221442,240626682, 9701900,10,03,D+F Ne 4
     5 103312152,140616092,181320182,222224262,263128352,12630000,10,04,AEL Ne 5
     6 629178711, 98311802,136715512,173619202,210422892,15790900,10,05/AEL Ne 6
      DATA NNN_Na/
     1 200020001,200320211,207322131,253031421,417657451,  513802,11,00,D+F Na 1
     2 100010001,100010161,119621261, 50711872,246445382, 4728901,11,01,D+F Na 2
     3 580158751,591860351, 71813142,321968812,106014333, 7165000,11,02,D+F Na 3
     4  96910772,116012242,130714232,153916552,177118872, 9888000,11,03,D+F Na 4
     5 601386081,108812932,148916832,187820722,226624612,13836900,11,04,AEL Na 5
     6 105712442,144616652,189221182,234425702,279630222,17209000,11,05/AEL Na 6
      DATA NNN_Mg/
     1 100010011,101410621,118414581,204831781,509479731,  764404,12,00,D+F Mg 1
     2 200120051,202921001,226926901,368457091, 92814872, 1503101,12,01,D+F Mg 2
     3 100010001,100110611,177455431,176546012, 99718753, 8011905,12,02,D+F Mg 3
     4 579758751,591459501,600560591,611461681,622362781,10928900,12,03,AEL Mg 4
     5 100611232,120612752,134214102,147815462,161416822,14122900,12,04,AEL Mg 5
     6 674896701,121814462,167018942,211723412,256527892,18648900,12,05/AEL Mg 6
      DATA NNN_Al/
     1 558857701,583558761,593260591,635969541,796790971,  598400,13,00,D+F Al 1
     2 100310211,110313021,172828201, 55311252,215637942, 1882203,13,01,D+F Al 2
     3 200320201,208622331,250530971,410251081,611571211, 2844000,13,02,D+F Al 3
     4 100010001,100210881,207436531,523168101,838999681,11996000,13,03,D+F Al 4
     5 577758651,591259631,604461351,622563161,640764981,15377000,13,04,AEL Al 5
     6 103511582,124713242,140014772,155316292,170517812,19042000,13,05/AEL Al 6
      DATA NNN_Si/
     1 825189211, 95210052,106211532,134317202,237934082,  814913,14,00,D+F Si 1
     2 563057761,588160311,631768671,791097651,127817282, 1634000,14,01,D+F Si 2
     3 101110771,126716471,232438081, 71914052,262045302, 3346001,14,02,D+F Si 3
     4 200720521,217224081,284439171,551370951, 86810262, 4513000,14,03,D+F Si 4
     5 100010001,100210881,207436531,523168101,838999681,16672900,14,04,FAK Si 5
     6 575458521,591459851,610063201,672674071,843698661,20510900,14,05/AEL Si 6
      DATA NNN_P/
     1 402643441,496757481,658274401,833492941,103511532, 1048300,15,00,AEL P  1
     2 874497931,106011282,119812802,138415142,164717802, 1972000,15,01,AEL P  2
     3 564058061,604164611,709579551, 90410172,112912422, 3015500,15,02,AEL P  3
     4 100811411,149720221,280936121,441552181,602168241, 5135400,15,03,AEL P  4
     5 200420781,227025361,281430911,336936471,392542021, 6500700,15,04,AEL P  5
     6 100010001,100010001,100010001,100010001,100010001,22041300,15,05/G   P  6
      DATA NNN_S/
     1 822887891,930697831,102610932,121614492,185124742, 1035708,16,00,D+F S  1
     2 443056011,694982961, 96911522,144218572,227326892, 2339900,16,01,D+F S  2
     3  91610392,113512242,136416942,233429882,364242962, 3500000,16,02,D+F S  3
     4 560058861,633871081, 82410062,123314602,168619132, 4728900,16,03,D+F S  4
     5 104512901,177025421,375163021,122420462,286036742, 7250000,16,04,D+F S  5
     6 202321571,241428261,358355061, 78310152,124814802, 8802800,16,05/D+F S  6
      DATA NNN_Cl/
     1 538155931,571657911,598067191, 89013782,227737172, 1300916,17,00,D+F Cl 1
     2 873396771,104411072,118513532,175525872,406763932, 2379903,17,01,D+F Cl 2
     3 506569571, 87610522,134421682,439092662,182132573, 3990006,17,02,D+F Cl 3
     4  95110872,120013232,154921252,345149322,641378942, 5350000,17,03,D+F Cl 4
     5 558960371,677779341, 95311692,138816082,182720472, 6780000,17,04/D+F Cl 5
      DATA NNN_Ar/
     1 100010001,100010051,106913911,240147261, 90716112, 1575411,18,00,D+F Ar 1
     2 550256831,578158781,636585461,151530162, 58010303, 2762007,18,01,D+F Ar 2
     3  92110362,112412002,133216772,254443722, 76512833, 4090003,18,02,D+F Ar 3
     4 582082081,103112292,149920212,309750502,720793642, 5978900,18,03,D+F Ar 4
     5  97111072,123213982,172625622,463976582,106413633, 7500000,18,04/D+F Ar 5
      DATA NNN_K/
     1 200020011,200720361,211923291,280137141,525575741,  433803,19,00,D+F K  1
     2 100010001,100110341,135929551, 79119282,405274892, 3180905,19,01,D+F K  2
     3 554657081,581260301, 73012702,285363872,129023363, 4600005,19,02,D+F K  3
     4  96010862,118413212,180836632, 90321023,416863253, 6090000,19,03,D+F K  4
     5 657793361,119515082,195826322,352944302,533162332, 8259900,19,04/D+F K  5
      DATA NNN_Ca/
     1 100110061,104311741,145919971,294345051, 69010322,  611003,20,00,D+F Ca 1
     2 205822781,279234761,427553061,688994901,136319772, 1186701,20,01,D+F Ca 2
     3 100010001,100510821,168744821,130232522, 69012813, 5121003,20,02,D+F Ca 3
     4 555157161,585662471, 82816862, 42510013,168423663, 6700000,20,03,D+F Ca 4
     5  99411262,123814062,182930402,484766392, 84310223, 8438900,20,04/D+F Ca 5
      DATA NNN_Sc/
     1 924696691,105212282,151219062,240530032,368944512,  653900,21,00,AEL Sc 1
     2 190424662,297634542,391743752,482952832,573761912, 1280000,21,01,AEL Sc 2
     3 976799291,101110322,105810882,111911502,118112122, 2475000,21,02,AEL Sc 3
     4 100010001,100510821,168744821,130232522, 69012813, 7390000,21,03,FAK Sc 4
     5 555157161,585662471, 82816862, 42510013,168423663, 9200000,21,04/FAK Sc 5
      DATA NNN_Ti/
     1 181021172,260333222,430155582,710089242,110213293,  681900,22,00,D+F Ti 1
     2 474659872,721284672, 98211413,134515623,177919963, 1356900,22,01,D+F Ti 2
     3 228327012,308134272,381143862,534563472,734983512, 2747000,22,02,D+F Ti 3
     4 971498311, 99210032,102610572,108711172,114711782, 4324000,22,03,D+F Ti 4
     5 100010001,100510821,168744821,130232522, 69012813, 9980000,22,04/FAK Ti 5
      DATA NNN_V/
     1 272835172,425851532,632278322, 97212013,146817723,  674000,23,00,AEL V  1
     2 373954132,743597002,121414713,173920143,229225713, 1464900,23,01,AEL V  2
     3 323142642,519660272,679975352,824789522, 96610363, 2930900,23,02,AEL V  3
     4 248329302,324234952,373439752,421744582,469949412, 4800000,23,03,AEL V  4
     5 970698231,990699881,100710152,102410322,104010482, 6500000,23,04/AEL V  5
      DATA NNN_Cr/
     1 717277611, 92911652,152620872,295141952,550468122,  676400,24,00,D+F Cr 1
     2  71611552,205635512,558281952,115315823,205625293, 1649000,24,01,D+F Cr 2
     3 280639822,538369722, 87610823,129115003,170919183, 3095000,24,02,D+F Cr 3
     4 377150952,616070292,791788382, 97610683,116012523, 5000000,24,03,D+F Cr 4
     5 264730962,341436462,394042872,463549832,533056782, 7300000,24,04/D+F Cr 5
      DATA NNN_Mn/
     1 600060321,629270891, 86911302,151020222,267534752,  743100,25,00,AEL Mn 1
     2 739594821,139921212,309342852,567372412, 97112553, 1563600,25,01,AEL Mn 2
     3  98417472,265535782,454754842,641973532,828792212, 3369000,25,02,AEL Mn 3
     4 328847052,586668342,771785912, 94710343,112112093, 5300000,25,03,AEL Mn 4
     5 422055132,636770792,779285062,921999322,106411363, 7600000,25,04/AEL Mn 5
      DATA NNN_Fe/
C    1 197023222,274433302,416753952,723799822,139419053,  787038,26,00,D+F Fe 1
     1 197023222,274433302,416753952,723799822,139419053,  790024,26,00,D+F Fe 1! Ion. potential from NIST J. Sugar and C. Corliss, J. Phys. Chem. Ref. Data 14, 1-664 (1985).
     2 409453722,686687452,110213823,174322233,286437043, 1618792,26,01,D+F Fe 2! Kurucz
c    2 409453722,686687452,110213823,174322233,286437043, 1617902,26,01,D+F Fe 2
c    3 262136422,501167232, 87911303,138916483,190721673, 3064300,26,02,D+F Fe 3
     3 262136422,501167232, 87911303,138916483,190721673, 3065200,26,02,D+F Fe 3 ! Kurucz
     4  98723522,420363072, 87011423,145117913,215925463, 5700000,26,03,AEL Fe 4
     5 388854482,666275742,846693572,102511143,120312923, 7900000,26,04/D+F Fe 5
      DATA NNN_Co/
c    1 199427202,335740022,474957182,708090462,118315403,  786000,27,00,D+F Co 1
     1 199427202,335740022,474957182,708090462,118315403,  788100,27,00,D+F Co 1
     2 279739202,490858232,684582472,104713233,159818733, 1704900,27,01,D+F Co 2
     3 279836622,461857562,720693022,124915873,192522633, 3349000,27,02,D+F Co 3
     4 262136422,501167232, 87911303,138916483,190821673, 5300000,27,03,FAK Co 4
     5  98723522,420363072, 87011423,145117913,215925463, 8300000,27,04/FAK Co 5
      DATA NNN_Ni/
c    1 227027622,306233052,356839222,446052912,652382292,  763314,28,00,D+F Ni 1
     1 227027622,306233052,356839222,446052912,652382292,  763996,28,00,D+F Ni 1
     2 108416342,222428472,353944332,577378932,110314303, 1814900,28,01,D+F Ni 2
     3 198724282,293236452,468362702, 86511123,136016073, 3516000,28,02,D+F Ni 3
     4 279836622,461857562,720693022,124915873,192522633, 5600000,28,03,FAK Ni 4
     5 262136422,501167232, 87911303,138916483,190721673, 7900000,28,04/FAK Ni 5
      DATA NNN_Cu/
     1 201620781,231026761,314737361,450555381,692386911,  772301,29,00,D+F Cu 1
     2 109415761,247938311, 58910042,190937022, 68311693, 2028903,29,01,D+F Cu 2
     3 897195961,107212972,165021182,260230862,356940532, 3682900,29,02/D+F Cu 3
      DATA NNN_Zn/
     1 100010001,100410231,108712611,167124841,388460411,  939102,30,00,D+F Zn 1
     2 200020021,201620761,223726341,351352061, 80812472, 1796001,30,01,D+F Zn 2
     3 100610471,122617301,300566361,149924112,332342352, 3970000,30,02/D+F Zn 3
      DATA NNN_Ga/
     1 403245601,493151431,529654331,559358091,611065171,  600000,31,00,AEL Ga 1
     2  99710051,104511541,135016501,208226431,321837921, 2050900,31,01,AEL Ga 2
     3 199820071,204521391,229124761,266028451,302932131, 3070000,31,02/AEL Ga 3
      DATA NNN_Ge/
     1 502665261,755183501,901496201,102410942,117912812,  787900,32,00,AEL Ge 1
     2 422848161,512153401,557458941,636270361,794489061, 1593000,32,01,AEL Ge 2
     3 100010261,114613921,175221251,249828711,324436181, 3421000,32,02/AEL Ge 3
      DATA NNN_As/
     1 403143241,491856701,649173781,840396751,113013392,  981000,33,00,AEL As 1
     2 593676641,884697521,105911572,129515012,180322212, 1858700,33,01,AEL As 2
     3 484470541, 91510972,125614082,157017612,199722912, 2829900,33,02/AEL As 3
      DATA NNN_Se/
     1 630172361,799686381,919797221,102810942,117712832,  975000,34,00,AEL Se 1
     2 438055511,691582151, 94510732,121413672,152016732, 2150000,34,01,AEL Se 2
     3 651982921, 94610382,113212492,139515462,169718482, 3200000,34,02/AEL Se 3
      DATA NNN_Br/
     1 437347431,498951671,538559501, 74710812,169126672, 1183910,35,00,D+F Br 1
     2 705183611, 93510092,111614162,222932532,427652992, 2160000,35,01,D+F Br 2
     3 510869921, 87410312,123116552,236530712,377744832, 3590000,35,02/D+F Br 3
      DATA NNN_Kr/
     1 100010001,100010051,105012781,198535971, 65911422, 1399507,36,00,D+F Kr 1
     2 461049811,522254261,609088131,168935052, 68612253, 2455908,36,01,D+F Kr 2
     3 759990901,101911142,129017782,302856642, 99414333, 3690000,36,02/D+F Kr 3
      DATA NNN_Rb/
     1 200020011,200720361,211523021,269434141,459163351,  417502,37,00,D+F Rb 1
     2 100010001,100110321,129524961, 61014202,291753192, 2750004,37,01,D+F Rb 2
     3 473650891,533156051, 66810932,232950852, 99915303, 4000000,37,02/D+F Rb 3
      DATA NNN_Sr/
     1 100110041,104111741,146019721,281941411,607785251,  569202,38,00,D+F Sr 1
     2 202621931,255331271,384347931,624085761,122417632, 1102600,38,01,D+F Sr 2
     3 100010001,100110321,129524961, 61014202,291753192, 4300000,38,02/FAK Sr 3
      DATA NNN_Y/
c    1 791587851,100012192,155119942,254031782,389946932,  637900,39,00,AEL Y  1
     1 791587851,100012192,155119942,254031782,389946932,  621710,39,00,AEL Y  1 ! From Kurucz
     2 118217102,220827002,319036792,416646512,513256072, 1223000,39,01,AEL Y  2
     3  92510012,104710862,112311612,120212472,132814282, 2050000,39,02/AEL Y  3
      DATA NNN_Zr/
     1 141320802,291439702,531170262, 92712273,162521053,  663400,40,00,D+F Zr 1 ! Ion. potential from NIST P.A. Hackett, M.R. Humphries, S.A. Mitchell, and D.M. Rayner, J. Chem. Phys. 85, 3194-3197 (1986)
     2 354454352,724689652,107212643,148517093,193321573, 1312900,40,01,D+F Zr 2
     3 209727032,324537052,415446282,510255752,604965222, 2298000,40,02/D+F Zr 3
      DATA NNN_Nb/
     1 256636022,465759302,749693962,116514243,171520333,  687900,41,00,AEL Nb 1
c     1 256636022,465759302,749693962,116514243,171520333,  675890,41,00,AEL Nb 1 ! From Kurucz
     2 335157222, 84511463,147718363,221826083,299933893, 1431900,41,01,AEL Nb 2
     3 223725352,280830972,340937362,406844002,473150632, 2503900,41,02/AEL Nb 3
      DATA NNN_Mo/
c    1 703972941, 82610822,154822682,327244912,571469372,  709900,42,00,D+F Mo 1
     1 703972941, 82610822,154822682,327244912,571469372,  709250,42,00,D+F Mo 1 ! From Kurucz
     2  69113342,270146932, 71810043,131916543,200323603, 1614900,42,01,NPk Mo 2 ! PFs are calculated using energy levels from Nilsson & Pickering, 2003, Phys. Scr., 67, 223
     3 267645462,669890262,115514323,173620673,242528083, 2714900,42,02/AEL Mo 3
      DATA NNN_Tc/
     1  90113722,190525812,348647032,631684102,110714373,  728000,43,00,Pal Tc 1 ! PFs are taken from Palmeri et al. 2007, MNRAS, 374, 63
     2 132521482,335250142, 72110033,135517843,229929083, 1525900,43,01,Pal Tc 2 ! PFs are taken from Palmeri et al. 2007, MNRAS, 374, 63
     3  80117462,174618952,189518952,189518952,189518952, 3000000,43,02/Pal Tc 3 ! PFs are taken from Palmeri et al. 2007, MNRAS, 374, 63
      DATA NNN_Ru/
     1 176824122,318941082,515263202,761790472,106112303,  736400,44,00,AEL Ru 1
     2 221934642,501968372, 88911173,136316243,189221613, 1675900,44,01,AEL Ru 2
     3 210622722,241025422,267928262,297731272,327834282, 2846000,44,02/AEL Ru 3
      DATA NNN_Rh/
     1 148520202,255230902,364942462,489656082,638872352,  746000,45,00,AEL Rh 1
     2 153421292,288137912,484660322,720187062,101011483, 1807000,45,01,AEL Rh 2
     3 254537212,492362292,770592182,107312243,137615273, 3104900,45,02/AEL Rh 3
      DATA NNN_Pd/
     1 115919651,320746011,607576761, 95011642,141817172,  832900,46,00,AEL Pd 1
     2 755087211,105913442,173122222,282034722,412247732, 1941900,46,01,AEL Pd 2
     3 180223462,289735212,414247632,538460052,662672472, 3292000,46,02/AEL Pd 3
      DATA NNN_Ag/
     1 200020001,200220141,206422141,257633021,455164681,  757403,47,00,D+F Ag 1
     2 100810581,125817401,260641031, 66210072,135316982, 2148000,47,01,D+F Ag 2
     3 795887491, 97711762,156620252,248329422,340038582, 3481900,47,02/D+F Ag 3
      DATA NNN_Cd/
     1 100010001,100410241,109212891,176827421,444268771,  899003,48,00,D+F Cd 1
     2 200020021,201720921,233329881,451475371,127520782, 1690301,48,01,D+F Cd 2
     3 100310281,114815371,246138311,519265531,791492761, 3747000,48,02/D+F Cd 3
      DATA NNN_In/
     1 252431921,368440461,433746521,512259221,723389021,  578400,49,00,D+F In 1
     2 100110071,104611651,146118581,225426511,304734431, 1886000,49,01,D+F In 2
     3 200120111,205021611,243628031,317035371,390442701, 2802900,49,02/D+F In 3
      DATA NNN_Sn/
     1 232637101,488058571,669074381,816189091, 97210632,  734200,50,00,AEL Sn 1
     2 286335941,408144471,479351961,571862901,686274341, 1462700,50,01,AEL Sn 2
     3 100010251,114013811,175321601,256829751,338337901, 3049000,50,02/AEL Sn 3
      DATA NNN_Sb/
     1 404043481,494656811,646772781,813490751,101411372,  863900,51,00,AEL Sb 1
     2 303147981,618472951,827392621,103711702,131214532, 1650000,51,01,AEL Sb 2
     3 313037601,429347901,536260591,689477591,862494881, 2529900,51,02/AEL Sb 3
      DATA NNN_Te/
     1 526258801,657372351,784284071,897095741,102711082,  900900,52,00,AEL Te 1
     2 440855541,686481251, 93810792,125414792,176321132, 1860000,52,01,AEL Te 2
     3 349054751,699883081, 96611302,134216202,197724212, 2800000,52,02/AEL Te 3
      DATA NNN_I/
     1 405342041,438645621,475751071,587974491,102214572, 1045404,53,00,D+F I  1
     2 568567471,773485861, 94510362,112712182,130914002, 1909000,53,01,D+F I  2
     3 514269581, 86910562,130716652,215327742,351843662, 3200000,53,02/AEL I  3
      DATA NNN_Xe/
     1 100010001,100010091,109515351,291060661,119621482, 1212716,54,00,D+F Xe 1
     2 414844131,465649111,538464651, 87112232,158019362, 2120000,54,01,D+F Xe 2
     3 615475101,867797531,112213462,157618062,203622662, 3209900,54,02/D+F Xe 3
      DATA NNN_Cs/
     1 200020001,201020501,215623871,283536181,462756261,  389300,55,00,D+F Cs 1
     2 100010001,100310371,119016501,269146361, 77912412, 2510000,55,01,D+F Cs 2
     3 424445601,481750061,516953311,549356551,581759791, 3500000,55,02/D+F Cs 3
      DATA NNN_Ba/
     1 101210791,135119351,282340571,574580391,111015062,  521002,56,00,D+F Ba 1
     2 262638611,504160621,698579371, 91010692,129115952, 1000000,56,01,D+F Ba 2
     3 100010001,100310351,118416321,264945521, 76512182, 3700000,56,02/FAK Ba 3
      DATA NNN_La/
     1  71111992,172323592,312540402,510763182,765791012,  557700,57,00,AEL La 1
     2 204529582,383647882,582469262,807992692,104911723, 1106000,57,01,AEL La 2
     3  94712552,148416582,179819212,203621522,227424042, 1917700,57,02/AEL La 3
      DATA NNN_Ce/
     1 516771922,101415733,230431963,422563713,661579353,  553870,58,00,AEL Ce 1 ! PFs are taken from Palmeri et al. 2000, Phys. Scr., 61, 323
     2  71918863,305242193,538665523,771988853,100511224, 1085000,58,01,MZH Ce 2 ! PFs are taken from Palmeri et al. 2000, Phys. Scr., 61, 323
     3 506183092,108612923,146416133,174418603,196520603, 2020000,58,02,CCB Ce 3 ! PFs are taken from Cowley & Barisciano 1994, Obs., 114, 308
     4 118012722,134214202,152616852,191722342,264131332, 3690600,58,03/RW  Ce 4 ! PFs are calculated using energy levels from Reader & Wyart 2009, Phys. Rev. A, 80, 042517
      DATA NNN_Pr/
     1 146526632,508289352,142720943,287237333,465456163,  547300,59,00,Sne Pr 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2  53615083,324256453, 86012064,159720354,251930474, 1055000,59,01,ISA Pr 2 ! PFs are calculated using energy levels from Mashonkina et al. 2009, A&A, 495, 297
     3 421093902,165924663,331041793,507660143,700980743, 2162400,59,02,ISA Pr 3 ! PFs are calculated using energy levels from Mashonkina et al. 2009, A&A, 495, 297
     4 373649462,593368882,785988552, 98810923,119813043, 3900000,59,03/AEL Pr 4 ! PFs are calculated using NIST energy levels
      DATA NNN_Nd/
     1 145623072,410172132,120218793,276138313,505263693,  552500,60,00,Sne Nd 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2  47511303,223037433,559777223,100512564,151817894, 1073000,60,01,ISA Nd 2 ! PFs are calculated using energy levels from Mashonkina et al. 2005, A&A, 441, 309
     3 432699302,204835193,525971403, 90710984,128314614, 2218000,60,02,ISA Nd 3 ! PFs are calculated using energy levels from Ryabchikova et al. 2006, A&A, 456, 329
     4 104717683,241529543,339937663,407343323,455447453, 4042000,60,03/Wyt Nd 4 ! PFs are calculated using energy levels from Wyart et al. 2006, J. Phys. B39, L77
      DATA NNN_Pm/
     1 293029302,339657372, 97415223,219529733,383647633,  558200,61,00,Fiv Pm 1 ! PFs are taken from Fivet at al. 2007, MNRAS, 380, 771
     2  53611273,274552953, 86912833,176222974,288035004, 1090000,61,01,Fiv Pm 2 ! PFs are taken from Fivet at al. 2007, MNRAS, 380, 771
     3  49012373,262048233,482348233,519661563,709279783, 2230000,61,02/Fiv Pm 3 ! PFs are taken from Fivet at al. 2007, MNRAS, 380, 771
      DATA NNN_Sm/
     1  92915672,222431062,444763802, 89612173,159520253,  564370,62,00,AEL Sm 1
     2 315059662, 97114563,204627093,342541693,490556383, 1106900,62,01,AEL Sm 2
     3 269037812,520270372, 91111273,133915483,172719093, 2340000,62,02/AEL Sm 3
      DATA NNN_Eu/
     1 800080571,851699301,127617362,240433032,444958442,  567045,63,00,AEL Eu 1
     2 125416052,211828182,375549622,644381732,101112213, 1124100,63,01,AEL Eu 2
     3  82514782, 47913863,315459503, 98114674,204226924, 2492000,63,02,ISA Eu 3 ! PFs are calculated using energy levels from Wyart et al. 2008, A&A, 483, 339
     4 353543472,487852542,553557522,592460632,617962762, 4265000,63,03/AEL Eu 4 ! PFs are calculated using NIST energy levels
      DATA NNN_Gd/
     1 244232982,441460242, 82611223,149719523,247930643,  615000,64,00,Sne Gd 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 534793262,139219123,247730843,371043333,495055893, 1209000,64,01,AEL Gd 2
     3 364145232,514756362,604864112,673870372,732276072, 2063000,64,02/AEL Gd 3
      DATA NNN_Tb/
     1 546880382,113515623,209227313,347543173,524362333,  586390,65,00,Sne Tb 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2  56510823,163922043,279234353,417550623,615575303, 1151900,65,01,Sne Tb 2 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     3  53713323,276551143, 85012894,181224014,304037114, 2191000,65,02/ISA Tb 3 ! PFs are calculated using Wyart & Ryabtsev extended energy levels analysis (Ryabtsev, private communication)
      DATA NNN_Dy/
     1 175219662,262038952,604693902,142320733,288338103,  593890,66,00,Sne Dy 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 347359162,108619003,300742453,533359923,606555733, 1167000,66,01,Sne Dy 2 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     3 320279972,191238513, 66810374,148019834,253331184, 2280000,66,02/ISA Dy 3 ! PFs are calculated using Wyart & Ryabtsev extended energy levels analysis (Ryabtsev, private communication)
      DATA NNN_Ho/
     1 222635002,542276772,100312353,145716713,187020703,  602160,67,00,FAK Ho 1
     2 321455092,112322203,401966563,102014674,200226144, 1180000,67,01,Bor Ho 2 ! PFs are taken from Bord & Cowley 2002, Sol. Phys., 211, 3
     3 222635002,542276772,100312353,145716713,187020703, 2284000,67,02/AEL Ho 3
      DATA NNN_Er/
     1 131715322,213632462,504577482,115416533,226829683,  610780,68,00,Sne Er 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 282946962, 81713443,201827463,339638403,399938623, 1193000,68,01,Sne Er 2 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     3 801281851, 91511592,166126662,472591362,190642503, 2274000,68,02/Irw Er 3 ! PFs are calculated using polynomial approximation from Irwin 1981, ApJS, 45, 621
      DATA NNN_Tm/
     1 800381111, 87510702,147621462,310343462,585475982,  618436,69,00,AEL Tm 1
     2 156718872,279244452,678196342,128316243,197823443, 1205000,69,01,AEL Tm 2
     3  93517192,364666132,103414613,192624193,293334613, 2368000,69,02/AEL Tm 3
      DATA NNN_Yb/
     1 104410001,100011021,142920191,299545391, 68910342,  625394,70,00,Sne Yb 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 200120901,270345231, 81714042,223533112,461959862, 1218400,70,01,AEL Yb 2
     3 100312561,250851931, 91914182,198626022,323638692, 2505000,70,02/AEL Yb 3
      DATA NNN_Lu/
     1 514664441,759086851, 99211442,133315612,182721252,  542589,71,00,AEL Lu 1
     2 125924831,438667801, 98714112,199727872,380850742, 1389900,71,01,AEL Lu 2
C     2 112718911,335853801,742987841,895879721,626944081, 1389900,71,01,Sne Lu 2
     3 323948621,661297271,158626482,426865032, 93712843, 2095960,71,02/AEL Lu 3
      DATA NNN_Hf/
     1 659294081,128016962,222528952,372047062,585171462,  700000,72,00,AEL Hf 1
     2  99117882,274638812,520867322, 84410313,123314453, 1489900,72,01,AEL Hf 2
     3 187427702,343739872,448049452,539358282,625266642, 2329900,72,02/AEL Hf 3
      DATA NNN_Ta/
     1  65210892,171325762,373552252,705192012,116414343,  787900,73,00,AEL Ta 1
     2 192837842,600784802,111113823,165419233,218524383, 1620000,73,01,AEL Ta 2
     3  99117872,274638812,520867312, 84410313,123314453, 2400000,73,02/FAK Ta 3
      DATA NNN_W/
     1 398981651,130019172,273438022,516168382, 88411163,  797900,74,00,AEL W  1
     2 131429482,523279952,111414623,183422233,262130233, 1770000,74,01,AEL W  2
     3 192837842,600784792,111113823,165419233,218524383, 2500000,74,02/FAK W  3
      DATA NNN_Re/
     1 600963001, 75910412,150121572,301940972,539168952,  787000,75,00,AEL Re 1
     2  73710852,190731262,464964142, 83810503,127315053, 1660000,75,01,AEL Re 2
     3 131429482,523279952,111414623,183422233,262130233, 2600000,75,02/FAK Re 3
      DATA NNN_Os/
     1 110815502,216829732,398752322,672484682,104612673,  850000,76,00,AEL Os 1
     2 168225972,362046562,566766422,757484612, 93010103, 1700000,76,01,AEL Os 2
     3  73710852,190731262,464964142, 83810503,127315053, 2700000,76,02/FAK Os 3
      DATA NNN_Ir/
     1 128117692,236030402,381847322,582671422, 87110533,  896700,77,00,AEL Ir 1 ! IP=8.96702 eV according to NIST
     2 216133402,476163702,811599542,118413753,156417503, 1691000,77,01,VKM Ir 2 ! PFs are calculated from energy levels of van Kleef & Metsch 1978, Physica C95, 251; IP=16.91 eV from Carlson et al. 1970, Atomic Data and Nuclear Data Table, 2, 63
     3 168225972,362046562,566766422,757484612, 93010103, 2800000,77,02/FAK Ir 3
      DATA NNN_Pt/
     1 158918512,207523002,254328242,316335762,407246582,  900000,78,00,AEL Pt 1
     2  98115462,224930742,401150612,623475412, 89910583, 1855900,78,01,AEL Pt 2
     3 110815502,216829732,398752322,672484682,104612673, 2900000,78,02/FAK Pt 3
      DATA NNN_Au/
     1 203222611,265731251,364042301,494958601,702084731,  922000,79,00,AEL Au 1
     2 120521331,357753801, 75310062,130516572,206925452, 2050000,79,01,AEL Au 2
     3 651780821,108814772,195925252,316338622,460853882, 3000000,79,02/AEL Au 3
      DATA NNN_Hg/
     1 100010001,100110111,105211851,152122101,341552811, 1043002,80,00,D+F Hg 1
     2 200320211,210023021,268834231,480472341,111416912, 1875000,80,01,D+F Hg 2
     3 104012871,186129471,458664151, 82410072,119013732, 3420000,80,02/D+F Hg 3
      DATA NNN_Tl/
     1 200420711,222424271,265429161,325637371,442853911,  610500,81,00,AEL Tl 1
     2 100010021,101910801,121414641,189525811,358949721, 2041900,81,01,AEL Tl 2
     3 200020311,216624611,296337451,489064791, 85711212, 2979900,81,02/AEL Tl 3
      DATA NNN_Pb/
     1 103411711,147819101,244331781,434862751, 93113762,  741404,82,00,D+F Pb 1
     2 204122231,248227841,311535621,429153941,651976431, 1502800,82,01,D+F Pb 2
     3 100210131,106812201,154522671,381665951, 95512512, 3192900,82,02/D+F Pb 3
      DATA NNN_Bi/
     1 400140351,416944121,474851591,564362181,690477231,  728700,83,00,AEL Bi 1
     2 106814451,204427341,350744811,586879131,108314772, 1667900,83,01,AEL Bi 2
     3 205523051,264830231,345439921,469156001,675281671, 2555900,83,02/AEL Bi 3
      DATA NNN_Po/
     1 500950661,518153561,559058941,628968071,748483501,  843000,84,00,AEL Po 1
     2 443756241,696282451, 95411012,128615262,182922012, 1900000,84,01,FAK Po 2
     3 336953201,682481011, 93810882,127915272,184622442, 2700000,84,02/FAK Po 3
      DATA NNN_At/
     1 402841621,431544771,463148311,520059491,734896851,  930000,85,00,FAK At 1
     2 576168741,788387631, 96910642,116012552,135014462, 2000000,85,01,FAK At 2
     3 490265341,812797201,116614322,179622692,285035302, 2900000,85,02/FAK At 3
      DATA NNN_Rn/
     1 100010001,100010031,102311051,133018071,264539391, 1074500,86,00,AEL Rn 1
     2 402841621,431544771,463148311,520059491,734996851, 2000000,86,01,FAK Rn 2
     3 576168741,788387631, 96910642,116012552,135014462, 3000000,86,02/FAK Rn 3
      DATA NNN_Fr/
     1 200020011,201220591,218124481,296538611,488859141,  400000,87,00,FAK Fr 1
     2 100010001,100010031,102311051,133018071,264539401, 2200000,87,01,FAK Fr 2
     3 421645151,477449611,511852711,542455761,572958821, 3300000,87,02/FAK Fr 3
      DATA NNN_Ra/
     1 104110411,105712431,155420871,293741981,596683361,  527800,88,00,Qui Ra 1 ! PFs are taken from Quinet et al. 2007, A&A, 474, 307
     2 198321961,258631331,381946231,552565051,754486211, 1015000,88,01,Qui Ra 2 ! PFs are taken from Quinet et al. 2007, A&A, 474, 307
     3 100010001,100010031,102311051,133018071,264539391, 3400000,88,02/FAK Ra 3
      DATA NNN_Ac/
     1 441654441,664281721,101912862,163320772,263333182,  517000,89,00,Qui Ac 1 ! PFs are taken from Quinet et al. 2007, A&A, 474, 307
     2 195142621, 72610952,153420412,261732632,397747612, 1175000,89,01,Qui Ac 2 ! PFs are taken from Quinet et al. 2007, A&A, 474, 307
     3 723989131,103511752,130814352,155416652,177018682, 2000000,89,02/AEL Ac 3
      DATA NNN_Th/
     1  63810522,177929162,457168312, 97513353,175722323,  630670,90,00,Sne Th 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 167142052, 79912843,186125143,322539763,475155383, 1190000,90,01,BWt Th 2 ! PFs are calculated from 508 energy levels of Blaise & Wyart 1992, Energy Levels and Atomic Spectra of Actinides, Paris
     3 491281082,108913303,154717483,193921253,230924903, 1830000,90,02/BWt Th 3 ! PFs are calculated from 175 energy levels of Blaise & Wyart 1992, Energy Levels and Atomic Spectra of Actinides, Paris
      DATA NNN_Pa/
     1 347877992,129318323,240730533,380546863,570368573,  600000,91,00,AEL Pa 1
     2 347877992,129318323,240730533,380546863,570368573, 1200000,91,01,FAK Pa 2
     3 347777992,129318323,240730533,380546863,570368573, 2000000,91,02/FAK Pa 3
      DATA NNN_U/
     1 209530092,450866762, 96613623,186524763,318839893,  619400,92,00,AEL U  1
     2  51311613,230239873,615986563,112513714,158317444, 1060000,92,01,Sne U  2 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     3 211130612,456267402, 94912483,151817063,177417123, 2000000,92,02/Irw U  3 ! PFs are calculated using polynomial approximation from Irwin 1981, ApJS, 45, 621
      DATA NNN_Np/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,93,00,FAK Np 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,93,01,FAK Np 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,93,02/FAK Np 3
      DATA NNN_Pu/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,94,00,FAK Pu 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,94,01,FAK Pu 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,94,02/FAK Pu 3
      DATA NNN_Am/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,95,00,FAK Am 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,95,01,FAK Am 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,95,02/FAK Am 3
      DATA NNN_Cm/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,96,00,FAK Cm 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,96,01,FAK Cm 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,96,02/FAK Cm 3
      DATA NNN_Bk/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,97,00,FAK Bk 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,97,01,FAK Bk 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,97,02/FAK Bk 3
      DATA NNN_Cf/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,98,00,FAK Cf 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,98,01,FAK Cf 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,98,02/FAK Cf 3
      DATA NNN_Es/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,99,00,FAK Es 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,99,01,FAK Es 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,99,02/FAK Es 3
      DATA SCALE/0.001,0.01,0.1,1.0/,FIRST/.TRUE./
C
C  First time XSAHA is called find the starting locations for each element
C
      IF(FIRST) THEN
        FIRST=.FALSE.
        IZ=0
        DO 1 N=1,NTABLE
        IF(NNNPFN(7,N).NE.IZ.AND.IZ.LE.ELESIZ) THEN
          IZ=NNNPFN(7,N)
          LOCZ(IZ)=N
        ENDIF
   1    CONTINUE
        LOCZ(IZ+1)=NTABLE+1
      ENDIF
C
C  Find starting row in the partition table and the number of ionization
C  stages available for a given element IEL
C
      N=LOCZ(IEL)
      NIONS=LOCZ(IEL+1)-N
C
C  For MODE=5 return the number of ionizations available for IEL
C
      IF(MODE.EQ.5) THEN
        MAXION=NIONS
        RETURN
      ENDIF
C
C  Compute T and kT in eV
C
      TTKEV=8.6171E-5*TT
      TV=TTKEV
      TTK=1.38065E-16*TT
C
C  Lowering of the ionization potential in Volts for unit Zeff
C
      CHARGE=2.*XNELEC
      EXCESS=XNELEC-XNATOM
C
C  Special allowance for doubly ionized Helium
C
      IF(EXCESS.GT.0.) CHARGE=CHARGE-EXCESS+4.*(2.*EXCESS)
C
C  Original code:
C     DEBYE=SQRT(TTK/(2.8965E-18*CHARGE))
C     POTLOW=MIN(1.,1.44E-7/DEBYE)
C
C  Compute the inverse of Debye radius to avoid division by zero at low temperatures
C
      DEBYE=SQRT(2.8965E-18*CHARGE/TTK)
      POTLOW=MIN(1.,1.44E-7*DEBYE)
C
C  Solve the Saha equation
C
      NION2=NIONS
      N=N-1
      DO 2 IONN=1,NION2
      Z=IONN
      POTLO(IONN)=POTLOW*Z
C      write(*,*) IP(IONN)-POTLO(IONN)
      N=N+1
      NNN100=NNNPFN(6,N)/100
      IP(IONN)=FLOAT(NNN100)/1000.
      G=NNNPFN(6,N)-NNN100*100
      IF(N.EQ.1) THEN
        PART(1)=2.
c        IF(TT.LT.9000.) GO TO 2
        PART(1)=PART(1)+8.*EXP(-10.196/TV)+18.*EXP(-12.084/TV)+32.*
     *          EXP(-12.745/TV)+50.*EXP(-13.051/TV)+72.*EXP(-13.217/TV)
        D1=13.595/6.5/6.5/TV
        D2=POTLO(1)/TV
      ELSE
        T2000=IP(IONN)*2000./11.
        IT=MAX(1,MIN(9,INT(TT/T2000-.5)))
        DT=TT/T2000-FLOAT(IT)-.5
        PMIN=1.
        I=(IT+1)/2
        K1=NNNPFN(I,N)/100000
        K2=NNNPFN(I,N)-K1*100000
        K3=K2/10
        KSCALE=K2-K3*10
        IF(MOD(IT,2).EQ.0) THEN
          P1=K3*SCALE(KSCALE)
          K1=NNNPFN(I+1,N)/100000
          KSCALE=MOD(NNNPFN(I+1,N),10)
          P2=K1*SCALE(KSCALE)
        ELSE
          P1=K1*SCALE(KSCALE)
          P2=K3*SCALE(KSCALE)
          IF(DT.LT.0.AND.KSCALE.LE.1) KP1=P1
          IF(DT.LT.0.AND.KSCALE.LE.1.AND.KP1.EQ.INT(P2+.5)) PMIN=KP1
        END IF
        PART(IONN)=MAX(PMIN,P1+(P2-P1)*DT)
c        write(*,*) (NNNPFN(I,N),I=1,6),PART(IONN),IP(IONN),G,IONN
        IF(G.EQ.0.0.OR.POTLO(IONN).LT.0.1.OR.TT.LT.T2000*4.0) GO TO 2
        IF(TT.GT.(T2000*11.)) TV=(T2000*11.)*8.6171E-5
        D1=0.1/TV
      END IF
      D2=POTLO(IONN)/TV
      PART(IONN)=PART(IONN)+G*EXP(-IP(IONN)/TV)*
     *           (SQRT(13.595*Z*Z/TV/D2)**3*
     *           (1./3.+(1.-(.5+(1./18.+D2/120.)*D2)*D2)*D2)-
     -           SQRT(13.595*Z*Z/TV/D1)**3*
     *           (1./3.+(1.-(.5+(1./18.+D1/120.)*D1)*D1)*D1))
c      TV=TTKEV
   2  CONTINUE
C
      IF(MODE.NE.3) THEN
        CF=2.*2.4148D15*TT*SQRT(TT)/XNELEC
        FFF(1)=1.
        DO 3 IONN=2,NION2
C
C  IF is to avoid annoying floating point underflows
C
        FEXARG=(IP(IONN-1)-POTLO(IONN-1))/TV
c        write(*,*) IONN,NION2,PART(IONN)/PART(IONN-1),FEXARG
c        IF(FEXARG.GT.80.) THEN
c          FFF(IONN)=0.
c        ELSE
          FFF(IONN)=CF*PART(IONN)/PART(IONN-1)*EXP(-FEXARG)
c        END IF
   3    CONTINUE
        DO 4 IONN=NION2,2,-1
        FFF(1)=1.+FFF(IONN)*FFF(1)
   4    CONTINUE
        FFF(1)=1./FFF(1)
        DO 5 IONN=2,NION2
        FFF(IONN)=FFF(IONN-1)*FFF(IONN)
   5    CONTINUE
        DO 6 IONN=1,MAXION
        FRCT(IONN)=1.
   6    CONTINUE
      ELSE
        DO 7 IONN=1,MAXION
        FRCT(IONN)=0.
   7    CONTINUE
      END IF
C
C  Formulate the answer according to MODE
C
      NIONS=MIN(MAXION,NION2)
      IF(MODE.EQ.1) THEN
        FRCT(1)=FFF(1)/PART(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 8 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
          FRCT(IONN)=FFF(IONN)/PART(IONN)
   8      CONTINUE
        END IF
      ELSE IF(MODE.EQ.2) THEN
        FRCT(1)=FFF(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 9 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
          FRCT(IONN)=FFF(IONN)
   9      CONTINUE
        END IF
      ELSE IF(MODE.EQ.3) THEN
        FRCT(1)=PART(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 10 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
          FRCT(IONN)=PART(IONN)
  10      CONTINUE
        END IF
      ELSE IF(MODE.EQ.4) THEN
        FRCT(1)=0
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 11 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
          FRCT(1)=FRCT(1)+FFF(IONN)*(IONN-1)
  11      CONTINUE
        END IF
      END IF
C
      RETURN
      END
