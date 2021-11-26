// This just contains all the definitions necessary for the parallel version

// The SME library version (and compilation date)
#ifndef VERSION
#define VERSION "6.03, July 2019"
#endif

/* Datafile locations */
// DATA_DIR is defined in platform.h

#define DATAFILE_FE "Fe1_Bautista2017.dat.INTEL"
#define DATAFILE_NH "NH_Stancil2018.dat.INTEL"
#define DATAFILE_STEHLE "stehle_long.dat.INTEL"
#define DATAFILE_BPO "bpo_self.grid.INTEL"
#define DATAFILE_VCS "vcsbalmer.dat"


#ifndef STATE_DEFINE
#define STATE_DEFINE

/* Global static variables and arrays */
#define MAX_ELEM 100 // Maximum Number of elements
#define MOSIZE 288   // Maximum Number of layers in the Model Atmosphere
#define MUSIZE 77    // Maximum Number of mu angles
#define MAX_PATHLEN 512
#define MAX_OUT_LEN 511
#define SP_LEN 8

extern "C" typedef struct
{
  /* IMPORTANT NOTE

    The internal notation for the model mode is inconsistent with
    the krz convention (in the krz 0 is RHOX and 1 is TAU):

    MOTYPE==0   means depth scale is "Tau", plane-parralel
    MOTYPE==1   means depth scale is "Rhox", plane-parralel
    MOTYPE==3   means depth scale is "RhoX", spherical
    MOTYPE==-1  fake value used with the call to OPMTRX get just
                just the line opacities
  */
  short NRHOX;
  short NRHOX_allocated;
  short MOTYPE;
  double TEFF, GRAV, WLSTD, RADIUS;
  int NumberSpectralSegments, NLINES, NWAVE_C;
  double WFIRST, WLAST, VW_scale;
  int N_SPLIST, IXH1, IXH2, IXH2mol, IXH2pl, IXHMIN,
      IXHE1, IXHE2, IXHE3, IXC1, IXAL1, IXSI1, IXSI2, IXCA1,
      IXMG1, IXMG2, IXCA2, IXN1, IXFE1, IXO1, IXCH, IXNH, IXOH;
  int PATHLEN, change_byte_order;
  int allocated_NLTE_lines;
  // consistency flags
  short flagMODEL, flagWLRANGE, flagABUND, flagLINELIST,
      flagIONIZ, flagCONTIN, lineOPACITIES, flagH2broad,
      initNLTE;

  /* Global pointers for dynamically allocated arrays */

  // statically sized arrays
  short IFOP[20];
  float ABUND[MAX_ELEM];
  double RHOX[MOSIZE], T[MOSIZE], XNE[MOSIZE], XNA[MOSIZE],
      RHO[MOSIZE], VTURB[MOSIZE], RAD_ATMO[MOSIZE];
  double XNA_eos[MOSIZE], XNE_eos[MOSIZE], RHO_eos[MOSIZE];
  double AHYD[MOSIZE], AH2P[MOSIZE], AHMIN[MOSIZE], SIGH[MOSIZE],
      AHE1[MOSIZE], AHE2[MOSIZE], AHEMIN[MOSIZE],
      SIGHE[MOSIZE], ACOOL[MOSIZE], ALUKE[MOSIZE],
      AHOT[MOSIZE], SIGEL[MOSIZE], SIGH2[MOSIZE];
  double TKEV[MOSIZE], TK[MOSIZE], HKT[MOSIZE], TLOG[MOSIZE];
  double FREQ, FREQLG, EHVKT[MOSIZE], STIM[MOSIZE], BNU[MOSIZE];
  float H1FRACT[MOSIZE], HE1FRACT[MOSIZE], H2molFRACT[MOSIZE];
  double COPBLU[MOSIZE], COPRED[MOSIZE], COPSTD[MOSIZE];
  double *LINEOP[MOSIZE], *AVOIGT[MOSIZE], *VVOIGT[MOSIZE];
  double LTE_b[MOSIZE];
  char PATH[MAX_PATHLEN];
  int debug_print;

  // dynamic arrays
  double **ATOTAL;
  int *INDX_C;
  double *YABUND, *XMASS, *EXCUP, *ENU4, *ENL4;
  double **BNLTE_low, **BNLTE_upp;
  float **FRACT, **PARTITION_FUNCTIONS, *POTION, *MOLWEIGHT;
  short *MARK, *AUTOION, *IDHEL;
  int *ION, *ANSTEE;
  double *WLCENT, *EXCIT, *GF,
      *GAMRAD, *GAMQST, *GAMVW, *ALMAX,
      *Wlim_left, *Wlim_right;
  char *SPLIST, *spname;
  int *SPINDEX;
  /* Consistency flags */
  short *flagNLTE;
  char result[MAX_OUT_LEN + 1]; /* leave a space for a '\0' */
} GlobalState;
#endif

extern "C" const char * GetInterfaces();

// define global parameter access
GlobalState * _NewState();
const char * _FreeState(short clean_pointers, GlobalState *state);
GlobalState * _CopyState(short clean_pointers, GlobalState *state);

// define the external methods
const char * _SMELibraryVersion(); /* Return SME library version */
const char * _GetDataFiles();      /* Return the required data files */
const char * _GetLibraryPath(GlobalState *state);    /* Return the current data file directory */
const char * _SetLibraryPath(const char * path, int pathlen, GlobalState *state);    /* Set the data file directory */
const char * _InputWaveRange(double wfirst, double wlast, GlobalState *state);    /* Read in Wavelength range */
const char * _SetVWscale(double vw_scale, GlobalState *state);        /* Set van der Waals scaling factor */
const char * _SetH2broad(short flag, GlobalState *state);        /* Set flag for H2 molecule */
const char * _InputLineList(int nlines, int slen, const char * species, double * linelist, GlobalState *state);     /* Read in line list */
const char * _OutputLineList(int nlines, double * linelist, GlobalState *state);    /* Return line list */
const char * _UpdateLineList(short nlines, int slen, short * index, const char * species,  double * linelist, GlobalState *state);    /* Change line list parameters */
const char * _InputModel(short ndepth, double teff, double grav, double wlstd, const char * motype_str,int motype_slen,short * opflag, double * depth, double * temp, double * xne, double * xna, double * rho, double * vt, double radius, double * height, GlobalState *state);        /* Read in model atmosphere */
const char * _InputDepartureCoefficients(double * bmat, int lineindex, GlobalState *state);
const char * _GetDepartureCoefficients(double * bmat, int nrhox, int line, GlobalState *state);   /* Get NLTE b's for specific line */
const char * _GetNLTEflags(short * nlte_flags, int nlines, GlobalState *state);               /* Get line list NLTE flags */
const char * _ResetDepartureCoefficients(GlobalState *state); /* Reset LTE */
const char * _InputAbund(double * abund, int nelements, GlobalState *state);                 /* Read in abundances */
const char * _Opacity(short request_output, short nout, double * out1, double * out2, double * out3, GlobalState *state);                    /* Calculate opacities */
const char * _GetOpacity(short ifop, short length, double * result, const char * species, int slen, const char * type, int tlen, GlobalState *state);                 /* Returns specific cont. opacity */
const char * _Ionization(short ion, GlobalState *state);                 /* Perfrom EOS calculations */
const char * _GetDensity(short length, double * result, GlobalState *state);                 /* Returns density in g/cm^3 */
const char * _GetNatom(short length, double * result, GlobalState *state);                   /* Returns atomic number density */
const char * _GetNelec(short length, double * result, GlobalState *state);                   /* Returns electron number density */
const char * _Transf(short nmu, double * mu, double * cint_seg, double * cintr_seg, int nwmax, int nw, double * wint_seg, double * sint_seg, double accrt, double accwi, short keep_lineop, short long_continuum, GlobalState *state);                     /* Computes spectral synthesis */
const char * _CentralDepth(int nmu, double * mu, int nwsize, float * table, double accrt, GlobalState *state);               /* Computes line central depths */
const char * _GetLineOpacity(double wave, short nrhox, double * lop, double * cop, double * scr, double * tsf, double * csf, GlobalState *state);             /* Returns specific line opacity */
const char * _GetLineRange(double * linerange, int nlines, GlobalState *state);               /* Get validity range for every line */
