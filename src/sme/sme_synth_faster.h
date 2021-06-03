
#ifdef BUILDING_SME_WIN_DLL
#define SME_DLL __declspec(dllexport)
#else
#define SME_DLL
#endif

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

// Define IDL String
typedef int IDL_STRING_SLEN_T;
#define IDL_STRING_MAX_SLEN 2147483647

typedef struct {		/* Define string descriptor */
  IDL_STRING_SLEN_T slen;	/* Length of string, 0 for null */
  short stype;			/* type of string, static or dynamic */
  char *s;			/* Addr of string */
} IDL_STRING;

// define the external methods
extern "C" char const *SME_DLL SMELibraryVersion(int n, void *arg[]); /* Return SME library version */
extern "C" char const *SME_DLL GetDataFiles(int n, void *arg[]);      /* Return the required data files */
extern "C" char const *SME_DLL GetLibraryPath(int n, void *arg[]);    /* Return the current data file directory */
extern "C" char const *SME_DLL SetLibraryPath(int n, void *arg[]);    /* Set the data file directory */
extern "C" char const *SME_DLL InputWaveRange(int n, void *arg[]);    /* Read in Wavelength range */
extern "C" char const *SME_DLL SetVWscale(int n, void *arg[]);        /* Set van der Waals scaling factor */
extern "C" char const *SME_DLL SetH2broad(int n, void *arg[]);        /* Set flag for H2 molecule */
extern "C" char const *SME_DLL ClearH2broad(int n, void *arg[]);      /* Clear flag for H2 molecule */
extern "C" char const *SME_DLL InputLineList(int n, void *arg[]);     /* Read in line list */
extern "C" char const *SME_DLL OutputLineList(int n, void *arg[]);    /* Return line list */
extern "C" char const *SME_DLL UpdateLineList(int n, void *arg[]);    /* Change line list parameters */
extern "C" char const *SME_DLL InputModel(int n, void *arg[]);        /* Read in model atmosphere */
extern "C" char const *SME_DLL InputDepartureCoefficients(int n, void *arg[]);
extern "C" char const *SME_DLL GetDepartureCoefficients(int n, void *arg[]);   /* Get NLTE b's for
                                                                                 specific line */
extern "C" char const *SME_DLL GetNLTEflags(int n, void *arg[]);               /* Get line list NLTE flags */
extern "C" char const *SME_DLL ResetDepartureCoefficients(int n, void *arg[]); /* Reset LTE */
extern "C" char const *SME_DLL InputAbund(int n, void *arg[]);                 /* Read in abundances */
extern "C" char const *SME_DLL Opacity(int n, void *arg[]);                    /* Calculate opacities */
extern "C" char const *SME_DLL GetOpacity(int n, void *arg[]);                 /* Returns specific cont. opacity */
extern "C" char const *SME_DLL Ionization(int n, void *arg[]);                 /* Perfrom EOS calculations */
extern "C" char const *SME_DLL GetDensity(int n, void *arg[]);                 /* Returns density in g/cm^3 */
extern "C" char const *SME_DLL GetNatom(int n, void *arg[]);                   /* Returns atomic number density */
extern "C" char const *SME_DLL GetNelec(int n, void *arg[]);                   /* Returns electron number density */
extern "C" char const *SME_DLL Transf(int n, void *arg[]);                     /* Computes spectral synthesis */
extern "C" char const *SME_DLL CentralDepth(int n, void *arg[]);               /* Computes line central depths */
extern "C" char const *SME_DLL GetLineOpacity(int n, void *arg[]);             /* Returns specific line opacity */
extern "C" char const *SME_DLL GetLineRange(int n, void *arg[]);               /* Get validity range for every line */
