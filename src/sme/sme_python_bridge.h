#ifndef IDL_DEFINE
#define IDL_DEFINE
// Define IDL String
typedef int IDL_STRING_SLEN_T;
#define IDL_STRING_MAX_SLEN 2147483647

typedef struct {		/* Define string descriptor */
  IDL_STRING_SLEN_T slen;	/* Length of string, 0 for null */
  short stype;			/* type of string, static or dynamic */
  char *s;			/* Addr of string */
} IDL_STRING;
#endif

// define the external methods
const char * Python_SMELibraryVersion(); /* Return SME library version */
const char * Python_GetDataFiles();      /* Return the required data files */
const char * Python_GetLibraryPath();    /* Return the current data file directory */
const char * Python_SetLibraryPath(IDL_STRING * path);    /* Set the data file directory */
const char * Python_InputWaveRange(double wmin, double wmax);    /* Read in Wavelength range */
const char * Python_SetVWscale(double gamma6);        /* Set van der Waals scaling factor */
const char * Python_SetH2broad();        /* Set flag for H2 molecule */
const char * Python_ClearH2broad();      /* Clear flag for H2 molecule */
const char * Python_InputLineList(int nlines, IDL_STRING * species, double * atomic);     /* Read in line list */
const char * Python_OutputLineList(int nlines, double * atomic);    /* Return line list */
const char * Python_UpdateLineList(short nlines, IDL_STRING * species, double * atomic, short * index);    /* Change line list parameters */
const char * Python_InputModel(short ndepth, double teff, double grav, double wlstd, IDL_STRING * motype, short * opflag, double * depth, double * temp, double * xne, double * xna, double * rho, double * vt, double radius, double * height);        /* Read in model atmosphere */
const char * Python_InputDepartureCoefficients(double * bmat, int lineindex);
const char * Python_GetDepartureCoefficients(double * bmat, int nrhox, int line);   /* Get NLTE b's for specific line */
const char * Python_GetNLTEflags(short * nlte_flags, int nlines);               /* Get line list NLTE flags */
const char * Python_ResetDepartureCoefficients(); /* Reset LTE */
const char * Python_InputAbund(double * abund);                 /* Read in abundances */
const char * Python_Opacity();                    /* Calculate opacities */
const char * Python_GetOpacity(short ifop, short length, double * result, IDL_STRING * species, IDL_STRING * key);                 /* Returns specific cont. opacity */
const char * Python_Ionization(short ion);                 /* Perfrom EOS calculations */
const char * Python_GetDensity(short length, double * result);                 /* Returns density in g/cm^3 */
const char * Python_GetNatom(short length, double * result);                   /* Returns atomic number density */
const char * Python_GetNelec(short length, double * result);                   /* Returns electron number density */
const char * Python_Transf(short nmu, double * mu, double * cint_seg, double * cintr_seg, int nwmax, int nw, double * wint_seg, double * sint_seg, double accrt, double accwi, short keep_lineop, short long_continuum);                     /* Computes spectral synthesis */
const char * Python_CentralDepth(int nmu, double * mu, int nwsize, float * table, double accrt);               /* Computes line central depths */
const char * Python_GetLineOpacity(double wave, short nrhox, double * lop, double * cop, double * scr, double * tsf, double * csf);             /* Returns specific line opacity */
const char * Python_GetLineRange(double * linerange, int nlines);               /* Get validity range for every line */
