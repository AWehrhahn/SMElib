#ifndef SME_DLL
#ifdef BUILDING_SME_WIN_DLL
#define SME_DLL __declspec(dllexport)
#else
#define SME_DLL
#endif
#endif

#ifndef IDL_DEFINE
#define IDL_DEFINE
// Define IDL String
typedef int IDL_STRING_SLEN_T;
#define IDL_STRING_MAX_SLEN 2147483647

typedef struct
{                         /* Define string descriptor */
  IDL_STRING_SLEN_T slen; /* Length of string, 0 for null */
  short stype;            /* type of string, static or dynamic */
  char *s;                /* Addr of string */
} IDL_STRING;
#endif

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

// define global parameter access
extern "C" int SME_DLL Parallel_GetNLINES(int n, void *arg[], GlobalState *state);
extern "C" short SME_DLL Parallel_GetNRHOX(int n, void *arg[], GlobalState *state);
extern "C" char *SME_DLL Parallel_GetSPNAME(int n, void *arg[], GlobalState *state);

extern "C" GlobalState *Parallel_NewState();
extern "C" const char *Parallel_FreeState(short clean_pointers, GlobalState *state);
extern "C" GlobalState *Parallel_CopyState(short clean_pointers, GlobalState *state);

// define the external methods
extern "C" const char *SME_DLL Parallel_SMELibraryVersion(int n, void *arg[], GlobalState *state); /* Return SME library version */
extern "C" const char *SME_DLL Parallel_GetDataFiles(int n, void *arg[], GlobalState *state);      /* Return the required data files */
extern "C" const char *SME_DLL Parallel_GetLibraryPath(int n, void *arg[], GlobalState *state);    /* Return the current data file directory */
extern "C" const char *SME_DLL Parallel_SetLibraryPath(int n, void *arg[], GlobalState *state);    /* Set the data file directory */
extern "C" const char *SME_DLL Parallel_InputWaveRange(int n, void *arg[], GlobalState *state);    /* Read in Wavelength range */
extern "C" const char *SME_DLL Parallel_SetVWscale(int n, void *arg[], GlobalState *state);        /* Set van der Waals scaling factor */
extern "C" const char *SME_DLL Parallel_SetH2broad(int n, void *arg[], GlobalState *state);        /* Set flag for H2 molecule */
extern "C" const char *SME_DLL Parallel_ClearH2broad(int n, void *arg[], GlobalState *state);      /* Clear flag for H2 molecule */
extern "C" const char *SME_DLL Parallel_InputLineList(int n, void *arg[], GlobalState *state);     /* Read in line list */
extern "C" const char *SME_DLL Parallel_OutputLineList(int n, void *arg[], GlobalState *state);    /* Return line list */
extern "C" const char *SME_DLL Parallel_UpdateLineList(int n, void *arg[], GlobalState *state);    /* Change line list parameters */
extern "C" const char *SME_DLL Parallel_InputModel(int n, void *arg[], GlobalState *state);        /* Read in model atmosphere */
extern "C" const char *SME_DLL Parallel_InputDepartureCoefficients(int n, void *arg[], GlobalState *state);
extern "C" const char *SME_DLL Parallel_GetDepartureCoefficients(int n, void *arg[], GlobalState *state);   /* Get NLTE b's for specific line */
extern "C" const char *SME_DLL Parallel_GetNLTEflags(int n, void *arg[], GlobalState *state);               /* Get line list NLTE flags */
extern "C" const char *SME_DLL Parallel_ResetDepartureCoefficients(int n, void *arg[], GlobalState *state); /* Reset LTE */
extern "C" const char *SME_DLL Parallel_InputAbund(int n, void *arg[], GlobalState *state);                 /* Read in abundances */
extern "C" const char *SME_DLL Parallel_Opacity(int n, void *arg[], GlobalState *state);                    /* Calculate opacities */
extern "C" const char *SME_DLL Parallel_GetOpacity(int n, void *arg[], GlobalState *state);                 /* Returns specific cont. opacity */
extern "C" const char *SME_DLL Parallel_Ionization(int n, void *arg[], GlobalState *state);                 /* Perfrom EOS calculations */
extern "C" const char *SME_DLL Parallel_GetDensity(int n, void *arg[], GlobalState *state);                 /* Returns density in g/cm^3 */
extern "C" const char *SME_DLL Parallel_GetNatom(int n, void *arg[], GlobalState *state);                   /* Returns atomic number density */
extern "C" const char *SME_DLL Parallel_GetNelec(int n, void *arg[], GlobalState *state);                   /* Returns electron number density */
extern "C" const char *SME_DLL Parallel_Transf(int n, void *arg[], GlobalState *state);                     /* Computes spectral synthesis */
extern "C" const char *SME_DLL Parallel_CentralDepth(int n, void *arg[], GlobalState *state);               /* Computes line central depths */
extern "C" const char *SME_DLL Parallel_GetLineOpacity(int n, void *arg[], GlobalState *state);             /* Returns specific line opacity */
extern "C" const char *SME_DLL Parallel_GetLineRange(int n, void *arg[], GlobalState *state);               /* Get validity range for every line */
