#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "platform.h"
#include "sme_synth_parallel.h"

/* Constants */
#define PI 3.14159265358979e0
#define SQRTPI 1.7724538509e0
#define CLIGHT 2.99792458e18
#define CLIGHTcm 2.99792458e10

#define pow10(x) exp(2.30258509299405e0 * (x))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define round(x) (x >= 0) ? (int)(x + 0.5) : (int)(x - 0.5)

/* Useful data */

float AMASS[MAX_ELEM] = {0.,
                         1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999,
                         18.998, 20.179, 22.990, 24.305, 26.982, 28.086, 30.974, 32.060,
                         35.453, 39.948, 39.102, 40.080, 44.956, 47.900, 50.941, 51.996,
                         54.938, 55.847, 58.933, 58.710, 63.546, 65.370, 69.720, 72.590,
                         74.922, 78.960, 79.904, 83.800, 85.468, 87.620, 88.906, 91.220,
                         92.906, 95.940, 98.906, 101.070, 102.905, 106.400, 107.868, 112.400,
                         114.820, 118.690, 121.750, 127.600, 126.905, 131.300, 132.905, 137.340,
                         138.906, 140.120, 140.908, 144.240, 146.000, 150.400, 151.960, 157.250,
                         158.925, 162.500, 164.930, 167.260, 168.934, 170.040, 174.970, 178.490,
                         180.948, 183.850, 186.200, 190.200, 192.200, 195.090, 196.967, 200.590,
                         204.370, 207.190, 208.981, 210.000, 210.000, 222.000, 223.000, 226.025,
                         227.000, 232.038, 230.040, 238.029, 237.048, 242.000, 242.000, 245.000,
                         248.000, 252.000, 253.000};
char ELEMEN[MAX_ELEM][3] = {" ",
                            "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne",
                            "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca",
                            "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                            "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr",
                            "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
                            "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
                            "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                            "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg",
                            "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                            "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cs", "Es"};

/* Default OK response */
const char OK_response = '\0';

/*
FREE macro to avoid freeing empty pointers
The second version below can be used to trace any attempts to
to do such a terrible thing!
*/
#if DEBUG_MEMORY_CHECK
#define CALLOC(ptr, varlen, vartype)                                \
  if (ptr != NULL)                                                  \
  {                                                                 \
    printf("Attempt to re-allocate %s line #%d\n", #ptr, __LINE__); \
    exit(99);                                                       \
  }                                                                 \
  ptr = (vartype *)calloc(varlen, sizeof(vartype));

#define FREE(ptr)                                                                    \
  if (ptr != NULL)                                                                   \
  {                                                                                  \
    free((char *)ptr);                                                               \
    ptr = NULL;                                                                      \
  }                                                                                  \
  else                                                                               \
  {                                                                                  \
    printf("Attempt to free unallocated variable %s at line #%d\n", #ptr, __LINE__); \
    exit(98);                                                                        \
  }
#else
#define CALLOC(ptr, varlen, vartype) ptr = (vartype *)calloc(varlen, sizeof(vartype))

#define FREE(ptr)      \
  if (ptr != NULL)     \
  {                    \
    free((char *)ptr); \
    ptr = NULL;        \
  }
#endif
/* Modules */

void ALAM(double *, GlobalState *state);
void CONTOP(double, double *, GlobalState *state);
void HOP(double *, int, int, GlobalState *state);
void H2PLOP(double *, int, int, GlobalState *state);
void HMINOP(double *, int, int, GlobalState *state);
void HMINOP_old(double *, int, int, GlobalState *state);
void HRAYOP(double *, int, GlobalState *state);
void HE1OP(double *, int, int, GlobalState *state);
void HE1OP_new(double *, int, int, GlobalState *state);
void HE2OP(double *, int, int, GlobalState *state);
void HEMIOP(double *, int, GlobalState *state);
void HERAOP(double *, int, GlobalState *state);
void COOLOP(double *, GlobalState *state);
double C1OP(int, GlobalState *state);
double MG1OP(int, GlobalState *state);
double AL1OP(int, GlobalState *state);
double SI1OP(int, GlobalState *state);
double FE1OP(int, GlobalState *state);
double C1OP_new(int, GlobalState *state);
double MG1OP_new(int, GlobalState *state);
double N1OP(int, GlobalState *state);
double O1OP(int, GlobalState *state);
double MG2OP(int, GlobalState *state);
double SI2OP(int, GlobalState *state);
double CA2OP(int, GlobalState *state);

void LUKEOP(double *, GlobalState *state);
void HOTOP(double *, GlobalState *state);
void ELECOP(double *, GlobalState *state);
void H2RAOP(double *, int, GlobalState *state);
int RKINTS(double *, int, double, double, double *, double *, double *,
           int, int &, double *, short, GlobalState *state);
int RKINTS_sph(double rhox[][2 * MOSIZE], int, int NRHOXs[], double, double,
               double *, double *, double *, int, int &,
               double *, short, int grazing[], GlobalState *state);
double FCINTG(double, double, double *, GlobalState *state);
void TBINTG(int, double *, double *, double *, double *, GlobalState *state);
void TBINTG_sph(int, double *, double *, double *, double *, int, GlobalState *state);
void CENTERINTG(double *, int, int, double *, double *, GlobalState *state);
void LINEOPAC(int, GlobalState *state);
void OPMTRX(double, double *, double *, double *, double *, int, int, GlobalState *state);
void OPMTRXn(double, double *, double *, double *, GlobalState *state);
void OPMTRX1(int, double *, GlobalState *state);
void GAMHE(short, double, double, double, double &, double &, GlobalState *state);
double HFNM(int, int);
double VCSE1F(double);
double VACAIR(double);
double SOFBET(double, double, int, int);

/* EOS FORTRAN routines */

extern "C" void xsaha_(int &, float &, float &, float &, int &, float *,
                       double *, int &);
extern "C" int eqcount_(char[][3], char *, int *, int &, int &, int &, int, int);
extern "C" int eqlist_(float *, char[][3], char *, int *, int *, char *, int &,
                       int &, int &, int &, int, int, int);

extern "C" void eqstat_(int &, float &, float &, float &, float *, char[][3],
                        float *, int &, int *, char *, float *, float *, float *,
                        float *, int &, int &, float &, float &, float &, int &,
                        int, int);
extern "C" void eqpf_(float &, float &, float &, float *, char[][3],
                      float *, int &, char *, int &, float *, int, int);

/* H-lines FORTRAN routines */

extern "C" float hlinop_(double &, int &, int &, double &, float &, float &,
                         float &, float &, float &);
extern "C" void hlinprof_(double &, double &, float &, float &, int &, int &,
                          float &, float &, float &, float &, char *, int *,
                          int *);

/* Code */

char *ByteSwap(char *s, int n)
{
  char c;
  int i, j;

  for (i = 0, j = n - 1; i < n / 2; i++, j--)
  {
    c = s[i];
    s[i] = s[j];
    s[j] = c;
  }
  return s;
}

char *Terminator(char *s, int len)
{
  static char tmpstore[128];
  strncpy(tmpstore, s, min(len, 127));
  tmpstore[127] = '\0';
  return tmpstore;
}

char *strtrim(char *s)
{
  int i, j, l = strlen(s);
  for (i = 0; i < l; i++)
    if (!isspace(s[i]))
      break;
  for (j = l - 1; j >= i; j--)
    if (isspace(s[j]))
      s[j] = '\0';
  return s + i;
}

int compress(char *target, char *source)
{
  /*
  This funcion  copies string "source" to string "target" elliminating
  all white spaces (space, tab, NL). All other characters are moved to
  the left, so normally "target" has the same or smaller length than
  source.
  "compress" returns the length of the compressed string.

   Author: N.Piskunov

   LAST UPDATE: October 24, 1994
   C++ Version: October 25, 1994
  */
  int s = 0, t = 0;
  do
    if (!isspace(source[s]))
      target[t++] = source[s];
  while (source[s++] != '\0');
  return t - 1;
}

extern "C" GlobalState *SME_DLL NewState(int n, void *args[])
{
  GlobalState *state = (GlobalState *)malloc(sizeof(GlobalState));

  state->FREQ = 0;
  state->FREQLG = 0;
  state->NRHOX = -1;
  state->NRHOX_allocated = -1;
  state->MOTYPE = 0;
  state->TEFF = -1;
  state->GRAV = -1;
  state->WLSTD = -1;
  state->RADIUS = -1;
  state->NumberSpectralSegments = -1;
  state->NLINES = -1;
  state->NWAVE_C = -1;
  state->WFIRST = -1;
  state->WLAST = -1;
  state->VW_scale = -1;
  state->N_SPLIST = -1;
  state->IXH1 = state->IXH2 = state->IXH2mol = state->IXH2pl = state->IXHMIN =
      state->IXHE1 = state->IXHE2 = state->IXHE3 = state->IXC1 = state->IXAL1 = state->IXSI1 =
          state->IXSI2 = state->IXCA1 = state->IXMG1 = state->IXMG2 = state->IXCA2 = state->IXN1 =
              state->IXFE1 = state->IXO1 = state->IXCH = state->IXNH = state->IXOH = -1;
  state->PATHLEN = 0;
  state->change_byte_order = 0;
  state->allocated_NLTE_lines = -1;
  // consistency flags
  state->flagMODEL = state->flagWLRANGE = state->flagABUND = state->flagLINELIST =
      state->flagIONIZ = state->flagCONTIN = state->lineOPACITIES = state->flagH2broad =
          state->initNLTE = 0;

  /* Global pointers for dynamically allocated arrays */

  // statically sized arrays
  // short IFOP[20];
  // float ABUND[MAX_ELEM];
  // double RHOX[MOSIZE], T[MOSIZE], XNE[MOSIZE], XNA[MOSIZE],
  //     RHO[MOSIZE], VTURB[MOSIZE], RAD_ATMO[MOSIZE];
  // double XNA_eos[MOSIZE], XNE_eos[MOSIZE], RHO_eos[MOSIZE];
  // double AHYD[MOSIZE], AH2P[MOSIZE], AHMIN[MOSIZE], SIGH[MOSIZE],
  //     AHE1[MOSIZE], AHE2[MOSIZE], AHEMIN[MOSIZE],
  //     SIGHE[MOSIZE], ACOOL[MOSIZE], ALUKE[MOSIZE],
  //     AHOT[MOSIZE], SIGEL[MOSIZE], SIGH2[MOSIZE];
  // double TKEV[MOSIZE], TK[MOSIZE], HKT[MOSIZE], TLOG[MOSIZE];
  // double FREQ, FREQLG, EHVKT[MOSIZE], STIM[MOSIZE], BNU[MOSIZE];
  // float H1FRACT[MOSIZE], HE1FRACT[MOSIZE], H2molFRACT[MOSIZE];
  // double COPBLU[MOSIZE], COPRED[MOSIZE], COPSTD[MOSIZE];
  // double LTE_b[MOSIZE];
  // char PATH[MAX_PATHLEN];

  // double *LINEOP[MOSIZE], *AVOIGT[MOSIZE], *VVOIGT[MOSIZE];
  for (int i = 0; i < MOSIZE; i++)
  {
    state->LINEOP[i] = NULL;
    state->AVOIGT[i] = NULL;
    state->VVOIGT[i] = NULL;
  }
  state->debug_print = 0;

  // dynamic arrays
  state->ATOTAL = NULL;
  state->INDX_C = NULL;
  state->YABUND = NULL;
  state->XMASS = NULL;
  state->EXCUP = NULL;
  state->ENU4 = NULL;
  state->ENL4 = NULL;
  state->BNLTE_low = NULL;
  state->BNLTE_upp = NULL;
  state->FRACT = NULL;
  state->PARTITION_FUNCTIONS = NULL;
  state->POTION = NULL;
  state->MOLWEIGHT = NULL;
  state->MARK = NULL;
  state->AUTOION = NULL;
  state->IDHEL = NULL;
  state->ION = NULL;
  state->ANSTEE = NULL;
  state->WLCENT = NULL;
  state->EXCIT = NULL;
  state->GF = NULL;
  state->GAMRAD = NULL;
  state->GAMQST = NULL;
  state->GAMVW = NULL;
  state->ALMAX = NULL;
  state->Wlim_left = NULL;
  state->Wlim_right = NULL;
  state->SPLIST = NULL;
  state->spname = NULL;
  state->SPINDEX = NULL;
  /* Consistency flags */
  state->flagNLTE = NULL;
  return state;
}

extern "C" const char *SME_DLL FreeState(int n, void *args[], GlobalState *state)
{
  free(state);
  return OK_response;
}

extern "C" GlobalState *SME_DLL CopyState(int n, void *args[], GlobalState *state)
{
  // NOTE: This is a shallow copy
  GlobalState *new_state = NewState(0, NULL);
  new_state->FREQ = state->FREQ;
  new_state->FREQLG = state->FREQLG;
  new_state->NRHOX = state->NRHOX;
  new_state->NRHOX_allocated = state->NRHOX_allocated;
  new_state->MOTYPE = state->MOTYPE;
  new_state->TEFF = state->TEFF;
  new_state->GRAV = state->GRAV;
  new_state->WLSTD = state->WLSTD;
  new_state->RADIUS = state->RADIUS;
  new_state->NumberSpectralSegments = state->NumberSpectralSegments;
  new_state->NLINES = state->NLINES;
  new_state->NWAVE_C = state->NWAVE_C;
  new_state->WFIRST = state->WFIRST;
  new_state->WLAST = state->WLAST;
  new_state->VW_scale = state->VW_scale;
  new_state->N_SPLIST = state->N_SPLIST;
  new_state->IXH1 = state->IXH1;
  new_state->IXH2 = state->IXH2;
  new_state->IXH2mol = state->IXH2mol;
  new_state->IXH2pl = state->IXH2pl;
  new_state->IXHMIN = state->IXHMIN;
  new_state->IXHE1 = state->IXHE1;
  new_state->IXHE2 = state->IXHE2;
  new_state->IXHE3 = state->IXHE3;
  new_state->IXC1 = state->IXC1;
  new_state->IXAL1 = state->IXAL1;
  new_state->IXSI1 = state->IXSI1;
  new_state->IXSI2 = state->IXSI2;
  new_state->IXCA1 = state->IXCA1;
  new_state->IXMG1 = state->IXMG1;
  new_state->IXMG2 = state->IXMG2;
  new_state->IXCA2 = state->IXCA2;
  new_state->IXN1 = state->IXN1;
  new_state->IXFE1 = state->IXFE1;
  new_state->IXO1 = state->IXO1;
  new_state->IXCH = state->IXCH;
  new_state->IXNH = state->IXNH;
  new_state->IXOH = state->IXOH;
  new_state->PATHLEN = state->PATHLEN;
  new_state->change_byte_order = state->change_byte_order;
  new_state->allocated_NLTE_lines = state->allocated_NLTE_lines;
  // consistency flags
  new_state->flagMODEL = state->flagMODEL;
  new_state->flagWLRANGE = state->flagWLRANGE;
  new_state->flagABUND = state->flagABUND;
  new_state->flagLINELIST = state->flagLINELIST;
  new_state->flagIONIZ = state->flagIONIZ;
  new_state->flagCONTIN = state->flagCONTIN;
  new_state->lineOPACITIES = state->lineOPACITIES;
  new_state->flagH2broad = state->flagH2broad;
  new_state->initNLTE = state->initNLTE;

  /* Global pointers for dynamically allocated arrays */

  // statically sized arrays
  // These have been assigned new memory, so the contents are copied
  for (int i = 0; i < 20; i++)
  {
    new_state->IFOP[i] = state->IFOP[i];
  }
  for (int i = 0; i < MAX_ELEM; i++)
  {
    new_state->ABUND[i] = state->ABUND[i];
  }
  for (int i = 0; i < MOSIZE; i++)
  {
    new_state->RHOX[i] = state->RHOX[i];
    new_state->T[i] = state->T[i];
    new_state->XNE[i] = state->XNE[i];
    new_state->XNA[i] = state->XNA[i];
    new_state->RHO[i] = state->RHO[i];
    new_state->VTURB[i] = state->VTURB[i];
    new_state->RAD_ATMO[i] = state->RAD_ATMO[i];
    new_state->XNA_eos[i] = state->XNA_eos[i];
    new_state->XNE_eos[i] = state->XNE_eos[i];
    new_state->RHO_eos[i] = state->RHO_eos[i];
    new_state->AHYD[i] = state->AHYD[i];
    new_state->AH2P[i] = state->AH2P[i];
    new_state->AHMIN[i] = state->AHMIN[i];
    new_state->SIGH[i] = state->SIGH[i];
    new_state->AHE1[i] = state->AHE1[i];
    new_state->AHE2[i] = state->AHE2[i];
    new_state->AHEMIN[i] = state->AHEMIN[i];
    new_state->SIGHE[i] = state->SIGHE[i];
    new_state->ACOOL[i] = state->ACOOL[i];
    new_state->ALUKE[i] = state->ALUKE[i];
    new_state->AHOT[i] = state->AHOT[i];
    new_state->SIGEL[i] = state->SIGEL[i];
    new_state->SIGH2[i] = state->SIGH2[i];
    new_state->TKEV[i] = state->TKEV[i];
    new_state->TK[i] = state->TK[i];
    new_state->HKT[i] = state->HKT[i];
    new_state->TLOG[i] = state->TLOG[i];
    new_state->EHVKT[i] = state->EHVKT[i];
    new_state->STIM[i] = state->STIM[i];
    new_state->BNU[i] = state->BNU[i];
    new_state->H1FRACT[i] = state->H1FRACT[i];
    new_state->HE1FRACT[i] = state->HE1FRACT[i];
    new_state->H2molFRACT[i] = state->H2molFRACT[i];
    new_state->COPBLU[i] = state->COPBLU[i];
    new_state->COPRED[i] = state->COPRED[i];
    new_state->COPSTD[i] = state->COPSTD[i];
    new_state->LTE_b[i] = state->LTE_b[i];

    // Allocate space for the line opacities and Voigt parameters
    for (int i = 0; i < state->NRHOX; i++)
    {
      new_state->LINEOP[i] = state->LINEOP[i];
      new_state->AVOIGT[i] = state->AVOIGT[i];
      new_state->VVOIGT[i] = state->VVOIGT[i];
    }
  }

  strncpy(new_state->PATH, state->PATH, MAX_PATHLEN);

  new_state->debug_print = state->debug_print;

  // dynamic arrays
  // for those we don't know the size, we just copy the pointer;
  new_state->ATOTAL = state->ATOTAL;
  new_state->INDX_C = state->INDX_C;
  new_state->YABUND = state->YABUND;
  new_state->XMASS = state->XMASS;
  new_state->EXCUP = state->EXCUP;
  new_state->ENU4 = state->ENU4;
  new_state->ENL4 = state->ENL4;
  new_state->BNLTE_low = state->BNLTE_low;
  new_state->BNLTE_upp = state->BNLTE_upp;
  new_state->flagNLTE = state->flagNLTE;
  new_state->FRACT = state->FRACT;
  new_state->PARTITION_FUNCTIONS = state->PARTITION_FUNCTIONS;
  new_state->POTION = state->POTION;
  new_state->MOLWEIGHT = state->MOLWEIGHT;
  new_state->MARK = state->MARK;
  new_state->AUTOION = state->AUTOION;
  new_state->IDHEL = state->IDHEL;
  new_state->ION = state->ION;
  new_state->ANSTEE = state->ANSTEE;
  new_state->WLCENT = state->WLCENT;
  new_state->EXCIT = state->EXCIT;
  new_state->GF = state->GF;
  new_state->GAMRAD = state->GAMRAD;
  new_state->GAMQST = state->GAMQST;
  new_state->GAMVW = state->GAMVW;
  new_state->ALMAX = state->ALMAX;
  new_state->Wlim_left = state->Wlim_left;
  new_state->Wlim_right = state->Wlim_right;
  new_state->SPLIST = state->SPLIST;
  new_state->spname = state->spname;
  new_state->SPINDEX = state->SPINDEX;
  return new_state;
}

extern "C" int SME_DLL GetNLINES(int n, void *args[], GlobalState *state)
{
  return state->NLINES;
}

extern "C" short SME_DLL GetNRHOX(int n, void *args[], GlobalState *state)
{
  return state->NRHOX;
}

extern "C" char *SME_DLL GetSPNAME(int n, void *args[], GlobalState *state)
{
  return state->spname;
}

extern "C" char const *SME_DLL SMELibraryVersion(int n, void *arg[], GlobalState *state) /* Return SME library version */
{
  sprintf(state->result, "SME Library version: %s, %s", VERSION, PLATFORM);
  return state->result;
}

extern "C" char const *SME_DLL GetDataFiles(int n, void *arg[], GlobalState *state) /* Return SME library version */
{
  sprintf(state->result, "%s;%s;%s;%s;%s", DATAFILE_FE, DATAFILE_NH, DATAFILE_STEHLE, DATAFILE_VCS, DATAFILE_BPO);
  return state->result;
}

extern "C" char const *SME_DLL GetLibraryPath(int n, void *arg[], GlobalState *state)
{
  sprintf(state->result, "%s", state->PATH);
  return state->result;
}

/*
  Set SME library datafile location
  If smelib was installed using make install the default location should point to the data files already
*/
extern "C" char const *SME_DLL SetLibraryPath(int n, void *arg[], GlobalState *state)
{
  state->PATHLEN = 0;
  if (n == 1)
  {
    state->PATHLEN = (*(IDL_STRING *)arg[0]).slen;
    strncpy(state->PATH, (*(IDL_STRING *)arg[0]).s, state->PATHLEN); /* Copy path to the Hydrogen line data files */
    state->PATH[state->PATHLEN] = '\0';
    state->change_byte_order = 1;
    state->change_byte_order = (*((char *)(&state->change_byte_order))) ? 0 : 1; /* Check if big-endian than need to change byte order */
    return &OK_response;
  }
  strcpy(state->result, "No path was specified");
  return state->result;
}

extern "C" char const *SME_DLL InputWaveRange(int n, void *arg[], GlobalState *state) /* Read in Wavelength range */
{
  int i;

  if (n < 2)
  {
    strcpy(state->result, "Only one argument found");
    return state->result;
  }
  if (state->flagWLRANGE)
  {
    if (fabs(state->WFIRST - *(double *)arg[0]) < 1.e-3 &&
        fabs(state->WLAST - *(double *)arg[1]) < 1.e-3)
      return &OK_response;
  }
  state->WFIRST = *(double *)arg[0];
  state->WLAST = *(double *)arg[1];
  if (state->WFIRST >= state->WLAST || state->WFIRST <= 0.0 || state->WLAST <= 0.)
  {
    state->flagWLRANGE = 0;
    strcpy(state->result, "Wrong wavelength range");
    return state->result;
  }
  else
  {
    state->flagWLRANGE = 1;
    state->flagCONTIN = 0;
    return &OK_response;
  }
}

extern "C" char const *SME_DLL SetVWscale(int n, void *arg[], GlobalState *state) /* Set van der Waals scaling factor */
{
  if (n < 1)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  state->VW_scale = *(double *)arg[0];
  state->VW_scale = fabs(state->VW_scale);
  return &OK_response;
}

extern "C" char const *SME_DLL SetH2broad(int n, void *arg[], GlobalState *state) /* Set flag for H2 molecule */
{
  state->flagH2broad = 1;
  return &OK_response;
}

extern "C" char const *SME_DLL ClearH2broad(int n, void *arg[], GlobalState *state) /* Clear flag for H2 molecule */
{
  state->flagH2broad = 0;
  return &OK_response;
}

extern "C" char const *SME_DLL InputLineList(int n, void *arg[], GlobalState *state) /* Read in line list */
{
  short l;
  int LINE, i;
  IDL_STRING *a0;
  double GFLOG, GRLG10, GSLG10, GWLG10,
      *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8;
  /*
   state->NLINES - NUMBERS OF SPECTRAL LINES;
   For each line:
   state->ION    - IONIZATION STAGE (1 - neutral, 2 - single ion, etc.)
   state->WLCENT - UNSHIFTED CENTRAL WAVELENGTH (Angstroems);
   state->EXCIT  - LOW LEVEL EXCITATION POTENTIAL IN eV;
   GFLOG  - log(state->GF);
   state->GAMRAD - RADIATION DAMPING (C1);
   state->GAMQST - QUADRATIC STARK DUMPING (C4);
   state->GAMVW  - VAN DER WAALS DUMPING (C6);
  */
  if (n < 2)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  if (state->flagLINELIST)
  {
    if (state->spname != NULL)
      FREE(state->spname);
    if (state->SPINDEX != NULL)
      FREE(state->SPINDEX);
    if (state->ION != NULL)
      FREE(state->ION);
    if (state->MARK != NULL)
      FREE(state->MARK);
    if (state->AUTOION != NULL)
      FREE(state->AUTOION);
    if (state->WLCENT != NULL)
      FREE(state->WLCENT);
    if (state->EXCIT != NULL)
      FREE(state->EXCIT);
    if (state->GF != NULL)
      FREE(state->GF);
    if (state->GAMRAD != NULL)
      FREE(state->GAMRAD);
    if (state->GAMQST != NULL)
      FREE(state->GAMQST);
    if (state->GAMVW != NULL)
      FREE(state->GAMVW);
    if (state->ANSTEE != NULL)
      FREE(state->ANSTEE);
    if (state->IDHEL != NULL)
      FREE(state->IDHEL);
    if (state->ALMAX != NULL)
      FREE(state->ALMAX);
    if (state->Wlim_left != NULL)
      FREE(state->Wlim_left);
    if (state->Wlim_right != NULL)
      FREE(state->Wlim_right);
    state->flagLINELIST = 0;
  }

  if (state->lineOPACITIES)
  {
    for (i = 0; i < state->NRHOX; i++)
    {
      if (state->LINEOP[i] != NULL)
        FREE(state->LINEOP[i]);
      if (state->AVOIGT[i] != NULL)
        FREE(state->AVOIGT[i]);
      if (state->VVOIGT[i] != NULL)
        FREE(state->VVOIGT[i]);
    }
    state->lineOPACITIES = 0;
  }

  state->NLINES = *(int *)arg[0];
  if (state->NLINES < 1)
  {
    state->flagLINELIST = 0;
    strcpy(state->result, "No line list");
    return state->result;
  }

  a3 = (double *)arg[2]; /* Setup pointers to line parameters */
  a3 += 2 * state->NLINES;
  for (LINE = 0; LINE < state->NLINES - 1; LINE++)
  {
    if (a3[LINE] > a3[LINE + 1]) /* Check that central wavelength are monotoneously increasing */
    {
      state->flagLINELIST = 0;
      strcpy(state->result, "Line list is not sorted in wavelength ascending order");
      return state->result;
    }
  }

  CALLOC(state->spname, state->NLINES * 8, char);
  CALLOC(state->SPINDEX, state->NLINES, int);
  CALLOC(state->ION, state->NLINES, int);
  CALLOC(state->MARK, state->NLINES, short);
  CALLOC(state->AUTOION, state->NLINES, short);
  CALLOC(state->WLCENT, state->NLINES, double);
  CALLOC(state->EXCIT, state->NLINES, double);
  CALLOC(state->GF, state->NLINES, double);
  CALLOC(state->GAMRAD, state->NLINES, double);
  CALLOC(state->GAMQST, state->NLINES, double);
  CALLOC(state->GAMVW, state->NLINES, double);
  CALLOC(state->ANSTEE, state->NLINES, int);
  CALLOC(state->IDHEL, state->NLINES, short);
  CALLOC(state->ALMAX, state->NLINES, double);
  CALLOC(state->Wlim_left, state->NLINES, double);
  CALLOC(state->Wlim_right, state->NLINES, double);

  if (state->Wlim_right == NULL)
  {
    if (state->spname != NULL)
    {
      FREE(state->spname);
    }
    if (state->SPINDEX != NULL)
      FREE(state->SPINDEX);
    if (state->ION != NULL)
      FREE(state->ION);
    if (state->MARK != NULL)
      FREE(state->MARK);
    if (state->AUTOION != NULL)
      FREE(state->AUTOION);
    if (state->WLCENT != NULL)
      FREE(state->WLCENT);
    if (state->EXCIT != NULL)
      FREE(state->EXCIT);
    if (state->GF != NULL)
      FREE(state->GF);
    if (state->GAMRAD != NULL)
      FREE(state->GAMRAD);
    if (state->GAMQST != NULL)
      FREE(state->GAMQST);
    if (state->GAMVW != NULL)
      FREE(state->GAMVW);
    if (state->ANSTEE != NULL)
      FREE(state->ANSTEE);
    if (state->IDHEL != NULL)
      FREE(state->IDHEL);
    if (state->ALMAX != NULL)
      FREE(state->ALMAX);
    if (state->Wlim_left != NULL)
      FREE(state->Wlim_left);
    if (state->Wlim_right != NULL)
      FREE(state->Wlim_right);
    state->flagLINELIST = 0;
    strcpy(state->result, "Not enough memory");
    return state->result;
  }

  a0 = (IDL_STRING *)arg[1]; /* Pointer to the list of species    */
  a1 = (double *)arg[2];     /* Setup pointers to line parameters */
  a2 = a1 + state->NLINES;
  a3 = a2 + state->NLINES;
  a4 = a3 + state->NLINES;
  a5 = a4 + state->NLINES;
  a6 = a5 + state->NLINES;
  a7 = a6 + state->NLINES;
  a8 = a7 + state->NLINES;

  state->VW_scale = 1;
  for (LINE = 0; LINE < state->NLINES; LINE++)
  {
    /* state->spname will be passed to FORTRAN, so no trailing zero's, fixed length
   padded with spaces instead */
    memcpy(state->spname + 8 * LINE, a0[LINE].s, a0[LINE].slen);
    if (a0[LINE].slen < 8)
      for (l = a0[LINE].slen; l < 8; l++)
        state->spname[8 * LINE + l] = ' ';
    //    state->ION[LINE]   =(int)a2[LINE]; /* Ionization            */
    for (l = 0; l < a0[LINE].slen; l++)
      if (*(a0[LINE].s + l) == ' ')
        break;
    state->ION[LINE] = (l == a0[LINE].slen) ? 1 : atoi(a0[LINE].s + l + 1);
    state->WLCENT[LINE] = a3[LINE];                      /* Central wavelength    */
    state->EXCIT[LINE] = a4[LINE];                       /* Excitation            */
    GFLOG = a5[LINE];                                    /* Oscillator strength   */
    state->GAMRAD[LINE] = a6[LINE];                      /* Radiative damping     */
    state->GAMQST[LINE] = a7[LINE];                      /* Stark damping         */
    state->GAMVW[LINE] = a8[LINE];                       /* Van der Waals damping */
    state->MARK[LINE] = -1;                              /* Initialize line flag  */
    state->Wlim_left[LINE] = state->WLCENT[LINE] - 150.; /* Initialize line contribution limits */
    state->Wlim_right[LINE] = state->WLCENT[LINE] + 150.;

    if (state->EXCIT[LINE] > 100.)
      state->EXCIT[LINE] = state->EXCIT[LINE] / 8065.544;
    if (state->GAMRAD[LINE] < 20. && state->GAMRAD[LINE] > 0.)
      state->GAMRAD[LINE] = pow10(state->GAMRAD[LINE]);
    GRLG10 = 0.;
    if (state->GAMRAD[LINE] > 0.)
      GRLG10 = log10(state->GAMRAD[LINE]);
    if (strncmp(state->spname + 8 * LINE, "H 1", 3)) /* Non-Hydrogen line */
    {
      if (state->GAMQST[LINE] < 0.)
        state->GAMQST[LINE] = pow10(state->GAMQST[LINE]);
      GSLG10 = 0.;
      if (state->GAMQST[LINE] > 0.)
        GSLG10 = log10(state->GAMQST[LINE]);
      if (state->GAMVW[LINE] < 0.)
      {
        state->GAMVW[LINE] = pow10(state->GAMVW[LINE]);
        GWLG10 = 0.;
        if (state->GAMVW[LINE] > 0.)
          GWLG10 = log10(state->GAMVW[LINE]);
        state->ANSTEE[LINE] = 0;
      }
      else if (state->GAMVW[LINE] > 10.)
      {
        GWLG10 = 0.;
        state->ANSTEE[LINE] = 1;
      }
    }
    else /* For hydrogen lines state->GAMQST & state->GAMVW have special meaning */
    {
      int nLO, nUP;
      nLO = GSLG10 = state->GAMQST[LINE];
      nUP = GWLG10 = state->GAMVW[LINE];
      if (nUP <= nLO || nLO <= 0) // Incorrect Hydrogen line format. Ignore it.
      {
        printf("SME will not compute H I line at %g A because energy level numbers are incorrect:\n",
               state->WLCENT[LINE]);
        printf("n_lower=%d, n_upper=%d\n", nLO, nUP);
        state->MARK[LINE] = 2;
      }
    }

    state->GF[LINE] = pow10(GFLOG);
  }
  state->flagLINELIST = 1;
  return &OK_response;
}

extern "C" char const *SME_DLL OutputLineList(int n, void *arg[], GlobalState *state) /* Return line list */
{
  int LINE, Nlines;
  double *a1;
  /*
   state->NLINES - NUMBERS OF SPECTRAL LINES;
   For each line:
   state->GAMRAD - RADIATION DAMPING (C1);
   state->GAMQST - QUADRATIC STARK DUMPING (C4);
   state->GAMVW  - VAN DER WAALS DUMPING (C6);
*/

  if (n < 2)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  if (!state->flagLINELIST)
  {
    strcpy(state->result, "No line list");
    return state->result;
  }
  Nlines = *(int *)arg[0];
  if (state->NLINES < 1)
  {
    state->flagLINELIST = 0;
    strcpy(state->result, "No line list");
    return state->result;
  }
  a1 = (double *)arg[1];

  for (LINE = 0; LINE < min(Nlines, state->NLINES); LINE++)
  {
    a1[6 * LINE] = state->WLCENT[LINE];
    a1[6 * LINE + 1] = state->GF[LINE];
    a1[6 * LINE + 2] = state->EXCIT[LINE];
    a1[6 * LINE + 3] = (state->GAMRAD[LINE] > 0.) ? log10(state->GAMRAD[LINE]) : 0.; /* Radiative damping     */
    if (strncmp(state->spname + 8 * LINE, "H ", 2))                                  /* Non-Hydrogen line     */
    {
      a1[6 * LINE + 4] = (state->GAMQST[LINE] > 0.) ? log10(state->GAMQST[LINE]) : 0.; /* Stark damping         */
      a1[6 * LINE + 5] = (state->GAMVW[LINE] > 0. &&
                          state->GAMVW[LINE] < 5.)
                             ? log10(state->GAMVW[LINE])
                             : state->GAMVW[LINE]; /* Van der Waals damping */
    }
    else /* Hydrogen line         */
    {
      a1[6 * LINE + 4] = state->GAMQST[LINE]; /* Stark damping         */
      a1[6 * LINE + 5] = state->GAMVW[LINE];  /* Van der Waals damping */
    }
  }
  return &OK_response;
}

extern "C" char const *SME_DLL UpdateLineList(int n, void *arg[], GlobalState *state) /* Change line list parameters */
{
  static char ERRMES[60];
  char tmpname[8];
  short LINE, NUPDTE, *INDEX;
  double GFLOG, GRLG10, GSLG10, GWLG10,
      *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8;
  IDL_STRING *a0;
  /*
   NUPDTE - NUMBERS OF SPECTRAL LINES;
   INDEX  - ARRAY OF INDICES IN EXISTING LINE LIST;
   For each line:
   state->ION    - IONIZATION STAGE (1 - neutral)
   state->WLCENT - UNSHIFTED CENTRAL WAVELENGTH (ANGSTREMS);
   state->EXCIT  - LOW LEVEL EXCITATION POTENTIAL IN EV;
   GFLOG  - log(state->GF);
   state->GAMRAD - RADIATION DAMPING (C1);
   state->GAMQST - QUADRATIC STARK DUMPING (C4);
   state->GAMVW  - VAN DER WAALS DUMPING (C6).
*/

  if (n < 4)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  if (!state->flagLINELIST)
  {
    strcpy(state->result, "Line list was not set. Cannot update.");
    return state->result;
  }
  NUPDTE = *(short *)arg[0];
  if (NUPDTE < 1)
    return &OK_response;

  a0 = (IDL_STRING *)arg[1]; /* Setup pointers for species        */
  a1 = (double *)arg[2];     /* Setup pointers to line parameters */
  a2 = a1 + NUPDTE;
  a3 = a2 + NUPDTE;
  a4 = a3 + NUPDTE;
  a5 = a4 + NUPDTE;
  a6 = a5 + NUPDTE;
  a7 = a6 + NUPDTE;
  a8 = a7 + NUPDTE;
  INDEX = (short *)arg[3];
  for (LINE = 0; LINE < NUPDTE; LINE++)
  {
    double WW, EXC;
    short i, l;

    i = INDEX[LINE];
    if (i < 0 || i >= state->NLINES)
    {
      strcpy(state->result, "Replacement index is out of range");
      return state->result;
    }

    /* state->spname will be passed to FORTRAN, so no trailing
   zero's, fixed length padded with spaces instead */

    memcpy(tmpname, a0[LINE].s, a0[LINE].slen);
    if (a0[LINE].slen < 8)
      for (l = a0[LINE].slen; l < 8; l++)
        tmpname[l] = ' ';
    WW = a3[LINE]; /* Wavelength */
    EXC = a4[LINE];
    if (EXC > 100.)
      EXC /= 8065.544; /* Excitation */

    /* Make sure we are talking about the same line.
   Check species name and excitation potential */

    if (strncmp(state->spname + 8 * i, tmpname, 8) || fabs(EXC - state->EXCIT[i]) > 0.005)
    {
      sprintf(ERRMES, "Attempt to replace line %d with another line", i);
      printf("Subst: %10.4f, '%s', %f, %f\n", WW, tmpname, EXC, a5[LINE]);
      printf("Orig:  %10.4f, '%4s', %f, %f\n", state->WLCENT[i], state->spname + 8 * i, state->EXCIT[i],
             log10(state->GF[i]));
      return ERRMES;
    }

    state->WLCENT[i] = WW;
    GFLOG = a5[LINE];
    state->GAMRAD[i] = a6[LINE];
    state->GAMQST[i] = a7[LINE];
    state->GAMVW[i] = a8[LINE];
    if (state->GAMRAD[i] < 20. && state->GAMRAD[i] > 0.)
      state->GAMRAD[i] = pow10(state->GAMRAD[i]);
    GRLG10 = 0.;
    if (state->GAMRAD[i] > 0.)
      GRLG10 = log10(state->GAMRAD[i]);
    if (strncmp(state->spname + 8 * i, "H ", 2)) /* Non-Hydrogen line */
    {
      if (state->GAMQST[i] < 0.)
        state->GAMQST[i] = pow10(state->GAMQST[i]);
      GSLG10 = 0.;
      if (state->GAMQST[i] > 0.)
        GSLG10 = log10(state->GAMQST[i]);
      if (state->GAMVW[i] < 0.)
        state->GAMVW[i] = pow10(state->GAMVW[i]);
      GWLG10 = 0.;
      if (state->GAMVW[i] > 0.)
        GWLG10 = log10(state->GAMVW[i]);
    }
    else /* For hydrogen lines this parameters have special meaning */
    {
      GSLG10 = state->GAMQST[i];
      GWLG10 = state->GAMVW[i];
    }
    state->GF[i] = pow10(GFLOG);
    state->MARK[i] = -1;                                     /* Mark line for is unknown in terms of opacity contribution */
    state->Wlim_left[i] = max(state->WLCENT[i] - 1000., 0.); /* Initialize line contribution limits */
    state->Wlim_right[i] = min(state->WLCENT[i] + 1000., 20000000.);
  }
  return &OK_response;
}

extern "C" char const *SME_DLL InputModel(int n, void *arg[], GlobalState *state) /* Read in model atmosphere */
{
  int IM, im, i, arg_offset;
  short *ifop, l;
  char motype[5];
  IDL_STRING *s;
  double TAU, DTAU1, DTAU2;
  double *a1, *a2, *a3, *a4, *a5, *a6, *a7;
  int L;

  if (n < 12)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }

  // Free invalidated arrays
  if (state->lineOPACITIES)
  {
    for (L = 0; L < state->NRHOX; L++)
    {
      FREE(state->LINEOP[L]);
      FREE(state->AVOIGT[L]);
      FREE(state->VVOIGT[L]);
    }
  }

  state->flagMODEL = 0;
  state->flagCONTIN = 0;
  state->lineOPACITIES = 0;

  state->NRHOX = *(short *)arg[0];
  if (state->NRHOX > MOSIZE)
  {
    sprintf(state->result, "SME library supports atmospheric model with maximum %d depth layers", MOSIZE);
    return state->result;
  }

  state->TEFF = *(double *)arg[1];
  state->GRAV = *(double *)arg[2];
  state->WLSTD = *(double *)arg[3];

  s = (IDL_STRING *)arg[4];
  l = min(4, s->slen);
  strncpy(motype, s->s, l);
  motype[l] = 0;
  for (i = 0; i < strlen(motype); i++)
    motype[i] = toupper(motype[i]);

  // Adding provision for spherical models
  if (!strncmp(motype, "TAU", 3))
  {
    state->MOTYPE = 0;
    arg_offset = 0;
    state->RADIUS = -1.;
  }
  else if (!strncmp(motype, "RHOX", 4))
  {
    state->MOTYPE = 1;
    arg_offset = 0;
    state->RADIUS = -1.;
  }
  else if (!strncmp(motype, "SPH", 3))
  {
    state->MOTYPE = 3;
    arg_offset = 1;
    state->RADIUS = *(double *)arg[5];
  }

  ifop = (short *)arg[5 + arg_offset];
  for (i = 0; i < 20; i++)
    state->IFOP[i] = ifop[i];

  // Allocate space for the line opacities and Voigt parameters
  if (!state->lineOPACITIES)
  {
    for (L = 0; L < state->NRHOX; L++)
    {
      CALLOC(state->LINEOP[L], state->NLINES, double);
      CALLOC(state->AVOIGT[L], state->NLINES, double);
      CALLOC(state->VVOIGT[L], state->NLINES, double);
    }
    state->lineOPACITIES = 1;
  }

  a1 = (double *)arg[6 + arg_offset];
  a2 = (double *)arg[7 + arg_offset];
  a3 = (double *)arg[8 + arg_offset];
  a4 = (double *)arg[9 + arg_offset];
  a5 = (double *)arg[10 + arg_offset];
  a6 = (double *)arg[11 + arg_offset];
  if (state->MOTYPE == 3)
    a7 = (double *)arg[12 + arg_offset];

  for (IM = im = 0; IM < state->NRHOX; im++, IM++) /* Copy model on the original grid */
  {                                                /* Intermediate points are found   */
    state->RHOX[IM] = a1[im];                      /* by iterpolation                 */
    state->T[IM] = a2[im];
    state->XNE[IM] = a3[im];
    state->XNA[IM] = a4[im];
    state->RHO[IM] = a5[im];
    state->VTURB[IM] = a6[im];
    if (state->MOTYPE == 3)
      state->RAD_ATMO[IM] = a7[im];
  }

  for (IM = 0; IM < state->NRHOX; IM++)
  {
    state->TKEV[IM] = 8.6171e-5 * state->T[IM];  // Temperature in eV
    state->TK[IM] = 1.38054e-16 * state->T[IM];  // Temperature times Boltzmann factor kT
                                                 // NP changed the value of the Planck constant from 6.6256e-27 in the line below 22-Jan-2018
    state->HKT[IM] = 6.6261e-27 / state->TK[IM]; // Plank constant divided by kT h/kT (h is in erg*s)
    state->TLOG[IM] = log(state->T[IM]);
  }
  state->flagMODEL = 1;
  return &OK_response;
}

extern "C" char const *SME_DLL InputDepartureCoefficients(int n, void *arg[], GlobalState *state)
{
  /* Reads in NLTE b's  for one transition at a time. The calling sequence
   requires a pointer to a double array of the size 2*state->NRHOX and an integer
   with the transition number. The logic of handling NLTE is the following:

   1) The first call is detected using a global static flag initNBLTE.
      At this moment we set the "default" departure coefficients state->LTE_b to 1,
      allocate the the vector of pointer the size of the line list and set them
      all to default and allocate the vector of flags state->flagNLTE all set to 0 (false)
   2) The initialization flag (state->initNLTE) is set to true
   3) The state->BNLTE_low and state->BNLTE_upp corresponding to the specified line are allocated
      state->NRHOX memory and the input array is copied there. The corresponding state->flagNLTE
      is set to 1 (true)
   4) Subsequent calls to the routine may allocate memory to other pointers or reset
      already existing once. In this case memory is reallocated to avoid leaks if
      state->NRHOX changes
   5) There no need to reset NLTE system in a given run, only in the end of calculations
  */
  int im, line;
  double *b;

  if (n < 2) // We assume that the caller will provide 2*state->NRHOX element array, so
             // be careful on the IDL side. The other argument is the line number.
  {
    strcpy(state->result, "No arguments found");
    return state->result;
  }
  if (!state->flagMODEL)
  {
    strcpy(state->result, "Model atmosphere must be set before departure coefficients");
    return state->result;
  }
  if (!state->flagLINELIST)
  {
    strcpy(state->result, "Line list must be set before departure coefficients");
    return state->result;
  }

  if (!state->initNLTE) // Initialize the departure arrays for the first time
  {
    for (im = 0; im < MOSIZE; im++)
      state->LTE_b[im] = 1.; // Initialize the default LTE b's

    CALLOC(state->BNLTE_low, state->NLINES, double *);
    CALLOC(state->BNLTE_upp, state->NLINES, double *);
    CALLOC(state->flagNLTE, state->NLINES, short);
    for (line = 0; line < state->NLINES; line++) // Set all lines to LTE first
    {
      state->BNLTE_low[line] = state->LTE_b;
      state->BNLTE_upp[line] = state->LTE_b;
      state->flagNLTE[line] = 0;
    }
    state->allocated_NLTE_lines = state->NLINES;
    state->initNLTE = 1;
  } // End of initialization

  b = (double *)arg[0];
  line = *(int *)arg[1];

  if (line < 0 || line >= state->allocated_NLTE_lines)
  {
    strcpy(state->result, "Attempt to set departure coefficients for non-existing transition");
    return state->result;
  }

  if (state->flagNLTE[line])
  {
    FREE(state->BNLTE_low[line]);
    FREE(state->BNLTE_upp[line]);
  }

  CALLOC(state->BNLTE_low[line], state->NRHOX, double); // Allocate departure coefficient arrays
  CALLOC(state->BNLTE_upp[line], state->NRHOX, double);

  for (im = 0; im < state->NRHOX; im++) // Copy departure coefficients
  {
    state->BNLTE_low[line][im] = *b++;
    state->BNLTE_upp[line][im] = *b++;
  }
  state->flagNLTE[line] = 1;

  return &OK_response;
}

extern "C" char const *SME_DLL GetDepartureCoefficients(int n, void *arg[], GlobalState *state) /* Get NLTE b's for specific line */
{
  int im;
  int nrhox, line;
  double *b;

  if (n < 3) // Check if arguments are present
  {
    strcpy(state->result, "Requires an array pointer, its length and line number");
    return state->result;
  }

  if (!state->initNLTE)
  {
    strcpy(state->result, "NLTE mode was not initialized. No departure coefficients available.");
    return state->result;
  }

  line = *(int *)arg[2];
  if (line < 0 || line >= state->NLINES)
  {
    strcpy(state->result, "Attempt to set departure coefficients for non-existing transition");
    return state->result;
  }

  b = (double *)arg[0];
  nrhox = *(int *)arg[1];

  if (state->flagNLTE[line])
  {
    for (im = 0; im < min(nrhox, state->NRHOX); im++)
    {
      *b++ = state->BNLTE_low[line][im];
      *b++ = state->BNLTE_upp[line][im];
    }
  }
  else
  {
    for (im = 0; im < min(nrhox, state->NRHOX); im++)
    {
      *b++ = 1.e0;
      *b++ = 1.e0;
    }
  }

  return &OK_response;
}

extern "C" char const *SME_DLL GetNLTEflags(int n, void *arg[], GlobalState *state) /* Get NLTE flag for every line */
{
  int nlines, line;
  short *b;

  if (n < 2) // Check if arguments are present
  {
    strcpy(state->result, "GetNLTELines: Requires an array pointer and its length");
    return state->result;
  }

  b = (short *)arg[0];
  nlines = *(int *)arg[1];

  if (!state->initNLTE)
  {
    for (line = 0; line < min(nlines, state->NLINES); line++)
    {
      b[line] = 0;
    }
    return &OK_response;
    ;
  }

  for (line = 0; line < min(nlines, state->NLINES); line++)
  {
    b[line] = state->flagNLTE[line];
  }

  return &OK_response;
}

extern "C" char const *SME_DLL ResetDepartureCoefficients(int n, void *arg[], GlobalState *state) /* Reset LTE */
{
  int line;

  if (!state->initNLTE)
    return &OK_response;

  for (line = 0; line < state->allocated_NLTE_lines; line++)
  {
    if (state->flagNLTE[line])
    {
      FREE(state->BNLTE_low[line]);
      FREE(state->BNLTE_upp[line]);
    }
  }
  FREE(state->flagNLTE);
  FREE(state->BNLTE_low);
  FREE(state->BNLTE_upp);
  state->allocated_NLTE_lines = 0;
  state->initNLTE = 0;

  return &OK_response;
}

extern "C" char const *SME_DLL InputAbund(int n, void *arg[], GlobalState *state) /* Read in abundances */
{
  int i;
  double *a;

  if (n < 1)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  a = (double *)arg[0];
  state->ABUND[0] = 1;
  for (i = 1; i < MAX_ELEM; i++)
  {
    state->ABUND[i] = (a[i - 1] >= 0.) ? a[i - 1] : pow10(a[i - 1]);
  }
  state->flagABUND = 1;
  state->flagCONTIN = 0;
  return &OK_response;
}

extern "C" char const *SME_DLL Opacity(int n, void *arg[], GlobalState *state) /* Calculate opacities */
{
  short i, nrhox;
  double *a1, *a2, *a3;

  if (n > 0)
  {
    if ((state->MOTYPE != 0 && n < 3) ||
        (state->MOTYPE == 0 && n < 4))
    {
      strcpy(state->result, "Opacity: Not enough arguments");
      return state->result;
    }
  }
  if (!state->flagMODEL)
  {
    strcpy(state->result, "Model atmosphere not set");
    return state->result;
  }
  if (!state->flagWLRANGE)
  {
    strcpy(state->result, "Wavelength interval was not specified");
    return state->result;
  }
  if (!state->flagABUND)
  {
    strcpy(state->result, "Abundances were not set");
    return state->result;
  }
  if (!state->flagIONIZ)
  {
    strcpy(state->result, "Molecular-ionization equilibrium was not computed");
    return state->result;
  }
  state->flagCONTIN = 0;

  // Continuous opacity at the red edge

  CONTOP(state->WLAST, state->COPRED, state);

  if (state->MOTYPE == 0)
    CONTOP(state->WLSTD, state->COPSTD, state); // Compute special opacity vector

  // Continuous opacity at the blue edge

  CONTOP(state->WFIRST, state->COPBLU, state);

  if (n >= 3)
  {
    i = *(short *)arg[0]; /* Length of IDL arrays */
    nrhox = min(state->NRHOX, i);
    a1 = (double *)arg[1];
    a2 = (double *)arg[2];
    if (state->MOTYPE == 0)
      a3 = (double *)arg[3];
    for (i = 0; i < nrhox; i++)
    {
      a1[i] = state->COPBLU[i];
      a2[i] = state->COPRED[i];
      if (n >= 4 && state->MOTYPE == 0)
        a3[i] = state->COPSTD[i];
    }
  }

  state->flagCONTIN = 1;
  return &OK_response;
}

void CONTOP(double WLCONT, double *opacity, GlobalState *state)
{
  /*  This subroutine computes the continuous opacity vector for one
    or two wavelengths.

     AUTHOR: N.Piskunov

     LAST UPDATE: January 12, 1992

    IF state->MOTYPE!= 0 - Kurucz type model with state->RHOX as depth scale
             == 0 - Depth parameter is TAUSTD

    WLCONT    - continuum wavelength
    opacity   - depth array of continuous opacity
*/
  double FREQ15;
  int j;

  state->FREQ = 2.997925e18 / WLCONT;
  state->FREQLG = log(state->FREQ);
  for (j = 0; j < state->NRHOX; j++)
  {
    state->EHVKT[j] = exp(-state->FREQ * state->HKT[j]);
    FREQ15 = state->FREQ * 1.e-15;
    state->STIM[j] = 1. - state->EHVKT[j];
    state->BNU[j] = 1.47439e-2 * FREQ15 * FREQ15 * FREQ15 * state->EHVKT[j] / state->STIM[j];
  }
  ALAM(opacity, state);
}

void ALAM(double *opacity, GlobalState *state)
{
  /*  THIS SUBROUTINE COMPUTES CONTINUOUS OPACITY USING
    KURUCZ's ATLAS-9 SUBROUTINES.
  */
  int J;

  /*  CLEAR OPACITY ACCUMULATORS */

  for (J = 0; J < state->NRHOX; J++)
  {
    state->AHYD[J] = 0;
    state->AH2P[J] = 0;
    state->AHMIN[J] = 0;
    state->SIGH[J] = 0;
    state->AHE1[J] = 0;
    state->AHE2[J] = 0;
    state->AHEMIN[J] = 0;
    state->SIGHE[J] = 0;
    state->ACOOL[J] = 0;
    state->ALUKE[J] = 0;
    state->AHOT[J] = 0;
    state->SIGEL[J] = 0;
    state->SIGH2[J] = 0;
  }

  if (state->IFOP[0] == 1)
    HOP(state->AHYD, state->IXH1, state->IXH2, state);
  if (state->IFOP[1] == 1)
    H2PLOP(state->AH2P, state->IXH1, state->IXH2, state);
  if (state->IFOP[2] == 1)
    HMINOP(state->AHMIN, state->IXH1, state->IXHMIN, state);
  if (state->IFOP[3] == 1)
    HRAYOP(state->SIGH, state->IXH1, state);
  if (state->IFOP[4] == 1)
    HE1OP_new(state->AHE1, state->IXHE1, state->IXHE2, state);
  if (state->IFOP[5] == 1)
    HE2OP(state->AHE2, state->IXHE2, state->IXHE3, state);
  if (state->IFOP[6] == 1)
    HEMIOP(state->AHEMIN, state->IXHE1, state);
  if (state->IFOP[7] == 1)
    HERAOP(state->SIGHE, state->IXHE1, state);
  if (state->IFOP[8] == 1)
    COOLOP(state->ACOOL, state);
  if (state->IFOP[9] == 1)
    LUKEOP(state->ALUKE, state);
  if (state->IFOP[10] == 1)
    HOTOP(state->AHOT, state);
  if (state->IFOP[11] == 1)
    ELECOP(state->SIGEL, state);
  if (state->IFOP[12] == 1)
    H2RAOP(state->SIGH2, state->IXH2mol, state);

  /*  CALCULATE THE TOTAL CONTINUOUS OPACITY */

  for (J = 0; J < state->NRHOX; J++)
  {
    opacity[J] = state->AHYD[J] + state->AH2P[J] + state->AHMIN[J] + state->SIGH[J] + state->AHE1[J] + state->AHE2[J] +
                 state->AHEMIN[J] + state->SIGHE[J] + state->ACOOL[J] + state->ALUKE[J] + state->AHOT[J] + state->SIGEL[J] +
                 state->SIGH2[J];
  }
  return;
}

double SEATON(double FREQ, double FREQ0, double XSECT, double POWER, double A)
{
  return XSECT * (A + (1. - A) * (FREQ0 / FREQ)) *
         pow(sqrt(FREQ0 / FREQ), floor(2. * POWER + 0.01));
}

double COULBF1S(double FREQ, double Z)
{
  static int kw = 72, mion = 1006;
  static double GAUNT1S[151] =
      {
          0.7973, 0.8094, 0.8212, 0.8328, 0.8439, 0.8548, 0.8653, 0.8754, 0.8852,
          0.8946, 0.9035, 0.9120, 0.9201, 0.9278, 0.9351, 0.9420, 0.9484, 0.9544,
          0.9601, 0.9653, 0.9702, 0.9745, 0.9785, 0.9820, 0.9852, 0.9879, 0.9903,
          0.9922, 0.9938, 0.9949, 0.9957, 0.9960, 0.9960, 0.9957, 0.9949, 0.9938,
          0.9923, 0.9905, 0.9884, 0.9859, 0.9832, 0.9801, 0.9767, 0.9730, 0.9688,
          0.9645, 0.9598, 0.9550, 0.9499, 0.9445, 0.9389, 0.9330, 0.9269, 0.9206,
          0.9140, 0.9071, 0.9001, 0.8930, 0.8856, 0.8781, 0.8705, 0.8627, 0.8546,
          0.8464, 0.8381, 0.8298, 0.8213, 0.8128, 0.8042, 0.7954, 0.7866, 0.7777,
          0.7685, 0.7593, 0.7502, 0.7410, 0.7318, 0.7226, 0.7134, 0.7042, 0.6951,
          0.6859, 0.6767, 0.6675, 0.6584, 0.6492, 0.6401, 0.6310, 0.6219, 0.6129,
          0.6039, 0.5948, 0.5859, 0.5769, 0.5680, 0.5590, 0.5502, 0.5413, 0.5324,
          0.5236, 0.5148, 0.5063, 0.4979, 0.4896, 0.4814, 0.4733, 0.4652, 0.4572,
          0.4493, 0.4415, 0.4337, 0.4261, 0.4185, 0.4110, 0.4035, 0.3962, 0.3889,
          0.3818, 0.3749, 0.3680, 0.3611, 0.3544, 0.3478, 0.3413, 0.3348, 0.3285,
          0.3222, 0.3160, 0.3099, 0.3039, 0.2980, 0.2923, 0.2866, 0.2810, 0.2755,
          0.2701, 0.2648, 0.2595, 0.2544, 0.2493, 0.2443, 0.2394, 0.2345, 0.2298,
          0.2251, 0.2205, 0.2160, 0.2115, 0.2072, 0.2029, 0.1987};
  double coulbf1s, elog;
  int I;

  coulbf1s = 0.;
  if (FREQ / (Z * Z) < 3.28805e15)
    return 0.;
  elog = log10(FREQ / (Z * Z) / 3.28805e15);
  I = (int)(elog / 0.02);
  I = max(min(I, 149), 0);
  coulbf1s = GAUNT1S[I] + (GAUNT1S[I + 1] - GAUNT1S[I]) / 0.02 * (elog - I * 0.02);
  return coulbf1s;
}

void LINTER(double XOLD[], double YOLD[], int NOLD,
            double XNEW[], double YNEW[], int NNEW)
{ // Assuning sorted in XOLD ind XNEW ascending order
  int IOLD, INEW;

  IOLD = 1;
  for (INEW = 0; INEW < NNEW; INEW++)
  {
    while (XNEW[INEW] >= XOLD[IOLD])
    {
      if (IOLD == NOLD - 1)
        break;
      IOLD++;
    }
    YNEW[INEW] = YOLD[IOLD - 1] + (YOLD[IOLD] - YOLD[IOLD - 1]) /
                                      (XOLD[IOLD] - XOLD[IOLD - 1]) *
                                      (XNEW[INEW] - XOLD[IOLD - 1]);
  }
  return;
}

int MAP1(double XOLD[], double FOLD[], int NOLD,
         double XNEW[], double FNEW[], int NNEW)
{
  int L, L1, L2, LL, K;
  double A, B, C, D, CBAC, CFOR, BBAC, BFOR, ABAC, AFOR, WT;

  L = 1;
  LL = -1;
  CFOR = BFOR = AFOR = 0.;
  for (K = 0; K < NNEW; K++)
  {
    while (L < NOLD)
    {
      if (XNEW[K] < XOLD[L])
        break;
      L++;
    }
    if (L == LL)
    {
      FNEW[K] = A + (B + C * XNEW[K]) * XNEW[K];
      continue;
    }
    if (L == NOLD)
    {
      L = min(NOLD - 1, L);
      C = 0.;
      B = (FOLD[L] - FOLD[L - 1]) / (XOLD[L] - XOLD[L - 1]);
      A = FOLD[L] - XOLD[L] * B;
      LL = L;
      FNEW[K] = A + (B + C * XNEW[K]) * XNEW[K];
      continue;
    }
    if (L > 2)
    {
      L1 = L - 1;
      if (L <= LL + 1 && (L != 2 || L != 3))
      {
        CBAC = CFOR;
        BBAC = BFOR;
        ABAC = AFOR;
      }
      else
      {
        L2 = L - 2;
        D = (FOLD[L1] - FOLD[L2]) / (XOLD[L1] - XOLD[L2]);
        CBAC = FOLD[L] / ((XOLD[L] - XOLD[L1]) * (XOLD[L] - XOLD[L2])) +
               (FOLD[L2] / (XOLD[L] - XOLD[L2]) - FOLD[L1] / (XOLD[L] - XOLD[L1])) /
                   (XOLD[L1] - XOLD[L2]);
        BBAC = D - (XOLD[L1] + XOLD[L2]) * CBAC;
        ABAC = FOLD[L2] - XOLD[L2] * D + XOLD[L1] * XOLD[L2] * CBAC;
      }
      if (L == NOLD)
      {
        C = CBAC;
        B = BBAC;
        A = ABAC;
        LL = L;
        FNEW[K] = A + (B + C * XNEW[K]) * XNEW[K];
        continue;
      }
      D = (FOLD[L] - FOLD[L1]) / (XOLD[L] - XOLD[L1]);
      CFOR = FOLD[L + 1] / ((XOLD[L + 1] - XOLD[L]) * (XOLD[L + 1] - XOLD[L1])) +
             (FOLD[L1] / (XOLD[L + 1] - XOLD[L1]) - FOLD[L] / (XOLD[L + 1] - XOLD[L])) /
                 (XOLD[L] - XOLD[L1]);
      BFOR = D - (XOLD[L] + XOLD[L1]) * CFOR;
      AFOR = FOLD[L1] - XOLD[L1] * D + XOLD[L] * XOLD[L1] * CFOR;
      WT = 0.;
      if (fabs(CFOR) != 0.)
        WT = fabs(CFOR) / (fabs(CFOR) + fabs(CBAC));
      A = AFOR + WT * (ABAC - AFOR);
      B = BFOR + WT * (BBAC - BFOR);
      C = CFOR + WT * (CBAC - CFOR);
      LL = L;
      FNEW[K] = A + (B + C * XNEW[K]) * XNEW[K];
    }
    else
    {
      L = min(NOLD - 1, L);
      C = 0.;
      B = (FOLD[L] - FOLD[L - 1]) / (XOLD[L] - XOLD[L - 1]);
      A = FOLD[L] - XOLD[L] * B;
      LL = L;
      FNEW[K] = A + (B + C * XNEW[K]) * XNEW[K];
    }
  }
  return LL - 1;
}

double XKARZAS(double FREQ, double ZEFF2, int N, int L)
{
  //     Karzas, W.J. and Latter, R. 1961, ApJS 6, 167-212.
  static float XN[15][29] =
      {{-30.274422, -29.048572, -28.181067, -26.962272, -25.437868, // X1
        -24.444170, -23.404269, -22.248421, -21.454163, -20.858944,
        -20.390346, -19.694283, -19.200905, -18.835387, -18.556686,
        -18.339364, -18.168213, -18.030238, -17.826632, -17.633456,
        -17.461067, -17.322353, -17.245241, -17.223162, -17.211266,
        -17.204840, -17.202587, -17.200999, -17.199715},
       {-31.779474, -30.553459, -29.685827, -28.466543, -26.940432, // X2
        -25.943993, -24.898608, -23.729491, -22.917021, -22.298979,
        -21.803393, -21.042629, -20.473370, -20.025469, -19.660029,
        -19.355246, -19.098003, -18.876442, -18.517855, -18.127425,
        -17.714170, -17.308930, -17.038908, -16.953361, -16.905447,
        -16.879127, -16.869826, -16.863085, -16.857754},
       {-32.659912, -31.433874, -30.566210, -29.346836, -27.820290, // X3
        -26.823453, -25.777089, -24.605440, -23.789519, -23.167057,
        -22.666147, -21.891933, -21.306393, -20.839041, -20.451712,
        -20.122889, -19.840361, -19.591597, -19.176587, -18.699419,
        -18.149566, -17.533628, -17.049033, -16.875774, -16.773227,
        -16.714935, -16.693926, -16.678663, -16.666369},
       {-33.284599, -32.058554, -31.190879, -29.971473, -28.444826, // X4
        -27.447836, -26.401066, -25.228582, -24.411413, -23.787317,
        -23.284581, -22.505775, -21.914353, -21.439606, -21.044235,
        -20.705972, -20.413135, -20.153596, -19.714525, -19.197426,
        -18.576241, -17.824248, -17.155428, -16.887819, -16.719154,
        -16.619216, -16.582315, -16.555295, -16.533276},
       {-33.769146, -32.543097, -31.675417, -30.455996, -28.929303, // X5
        -27.932243, -26.885239, -25.712408, -24.894628, -24.269941,
        -23.766226, -22.985245, -22.390846, -21.912586, -21.513577,
        -21.170761, -20.873304, -20.608270, -20.156957, -19.619181,
        -18.958075, -18.121143, -17.308727, -16.951892, -16.712503,
        -16.563827, -16.507488, -16.465627, -16.431184},
       {-34.165051, -32.939000, -32.071317, -30.851888, -29.325169, // X6
        -28.328071, -27.280986, -26.107892, -25.289843, -24.664705,
        -24.160564, -23.378190, -22.782394, -22.302428, -21.901012,
        -21.555896, -21.255472, -20.987585, -20.529803, -19.979782,
        -19.295022, -18.402541, -17.482757, -17.047424, -16.737838,
        -16.536084, -16.457331, -16.397931, -16.348398},
       {-34.499784, -33.273731, -32.406047, -31.186614, -29.659879, // X7
        -28.662758, -27.615624, -26.442410, -25.624138, -24.998790,
        -24.494343, -23.711394, -23.114332, -22.633333, -22.230699,
        -21.884181, -21.582185, -21.312152, -20.849982, -20.292819,
        -19.593097, -18.663739, -17.663648, -17.161477, -16.785637,
        -16.528798, -16.425342, -16.345983, -16.278790},
       {-34.789743, -33.563690, -32.696004, -31.476568, -29.949823, // X8
        -28.952576, -27.905521, -26.732230, -25.913849, -25.288312,
        -24.783697, -24.000359, -23.402741, -22.921064, -22.517235,
        -22.169801, -21.866776, -21.595595, -21.130798, -20.568503,
        -19.858590, -18.903358, -17.843146, -17.285660, -16.849210,
        -16.537235, -16.407454, -16.306014, -16.218699},
       {-35.045505, -33.819451, -32.951765, -31.732326, -30.205575, // X9
        -29.208318, -28.161241, -26.987832, -26.169441, -25.543807,
        -25.039029, -24.255440, -23.657439, -23.175297, -22.770919,
        -22.422852, -22.118723, -21.846749, -21.380133, -20.814545,
        -20.097359, -19.123314, -18.017622, -17.414518, -16.923750,
        -16.558183, -16.401026, -16.275647, -16.165911},
       {-35.274293, -34.048238, -33.180551, -31.961111, -30.434355, // X10
        -29.437090, -28.389998, -27.216550, -26.398051, -25.772354,
        -25.267495, -24.483312, -23.885464, -23.402587, -22.997820,
        -22.649302, -22.344664, -22.072514, -21.604193, -21.035827,
        -20.313639, -19.326284, -18.184568, -17.544349, -17.005732,
        -16.588554, -16.403642, -16.253350, -16.118795},
       {-35.481256, -34.255201, -33.387514, -32.168073, -30.641313, // X11
        -29.644043, -28.596939, -27.423463, -26.604924, -25.979176,
        -25.474255, -24.689915, -24.091864, -23.608739, -23.203681,
        -22.854826, -22.549810, -22.276842, -21.807547, -21.237407,
        -20.511071, -19.513620, -18.342986, -17.672186, -17.092253,
        -16.625647, -16.412652, -16.237373, -16.076228},
       {-35.670198, -34.444144, -33.576456, -32.357014, -30.830251, // X12
        -29.832977, -28.785864, -27.612367, -26.793798, -26.168012,
        -25.663043, -24.878583, -24.280378, -23.797065, -23.391784,
        -23.042673, -22.737368, -22.464078, -21.994040, -21.422148,
        -20.692935, -19.687256, -18.494545, -17.795069, -17.182159,
        -16.669643, -16.429381, -16.227310, -16.037494},
       {-35.844009, -34.617954, -33.750266, -32.530823, -31.004058, // X13
        -30.006781, -28.959661, -27.786148, -26.967555, -26.341739,
        -25.836687, -25.051753, -24.453445, -23.969994, -23.564544,
        -23.215236, -22.909707, -22.636559, -22.165546, -21.592592,
        -20.861125, -19.849269, -18.640363, -17.921966, -17.273191,
        -16.719020, -16.451969, -16.222218, -16.001878},
       {-36.004932, -34.778877, -33.911189, -32.691746, -31.164979, // X14
        -30.167699, -29.120574, -27.947047, -27.128436, -26.502596,
        -25.997515, -25.212506, -24.614103, -24.130536, -23.724949,
        -23.375482, -23.069774, -22.796032, -22.324557, -21.750758,
        -21.017491, -20.000677, -18.777116, -18.041065, -17.364348,
        -16.772813, -16.479089, -16.221551, -15.968930},
       {-36.154748, -34.928693, -34.061005, -32.841561, -31.314793, // X15
        -30.317511, -29.270382, -28.096844, -27.278218, -26.652358,
        -26.147254, -25.362186, -24.763705, -24.280044, -23.874346,
        -23.524751, -23.218899, -22.944996, -22.473148, -21.898667,
        -21.163944, -20.143099, -18.907170, -18.155759, -17.454858,
        -16.827663, -16.509932, -16.224591, -15.938340}};
  static float FREQN[15][29] = {
      {19.516982, 19.164810, 18.915052, 18.563043, 18.120083, // FREQ1
       17.828904, 17.521260, 17.174377, 16.931912, 16.747387,
       16.600083, 16.377277, 16.215909, 16.094200, 15.999955,
       15.925518, 15.866216, 15.817969, 15.745954, 15.676626,
       15.613849, 15.562692, 15.533972, 15.525713, 15.521260,
       15.518864, 15.518023, 15.517421, 15.516939},
      {19.516949, 19.164737, 18.914922, 18.562750, 18.119270, // FREQ2
       17.827313, 17.518023, 17.167149, 16.919200, 16.727792,
       16.572317, 16.329852, 16.145327, 15.998094, 15.876964,
       15.775097, 15.688665, 15.613849, 15.492095, 15.358548,
       15.215909, 15.074566, 14.979337, 14.948961, 14.931912,
       14.922531, 14.919200, 14.916804, 14.914879},
      {19.516943, 19.164723, 18.914898, 18.562696, 18.119119, // FREQ3
       17.827018, 17.517421, 17.165797, 16.916804, 16.724064,
       16.566974, 16.320472, 16.130898, 15.977703, 15.849803,
       15.740463, 15.646019, 15.562696, 15.423010, 15.261631,
       15.074579, 14.863704, 14.696235, 14.635934, 14.600123,
       14.579728, 14.572359, 14.567017, 14.562696},
      {19.516941, 19.164719, 18.914889, 18.562677, 18.119066, // FREQ4
       17.826915, 17.517210, 17.165323, 16.915963, 16.722752,
       16.565089, 16.317140, 16.125732, 15.970333, 15.839881,
       15.727658, 15.630046, 15.543267, 15.395977, 15.221861,
       15.011789, 14.756488, 14.527662, 14.435545, 14.377277,
       14.342650, 14.329852, 14.320471, 14.312819},
      {19.516940, 19.164717, 18.914886, 18.562668, 18.119042, // FREQ5
       17.826867, 17.517112, 17.165103, 16.915573, 16.722143,
       16.564213, 16.315589, 16.123320, 15.966880, 15.835211,
       15.721601, 15.622449, 15.533972, 15.382871, 15.202143,
       14.979337, 14.696203, 14.420029, 14.298047, 14.215909,
       14.164752, 14.145327, 14.130897, 14.118999},
      {19.516940, 19.164715, 18.914883, 18.562663, 18.119029, // FREQ6
       17.826841, 17.517059, 17.164984, 16.915361, 16.721812,
       16.563737, 16.314744, 16.122004, 15.964992, 15.832652,
       15.718275, 15.618265, 15.528838, 15.375583, 15.191044,
       14.960636, 14.659571, 14.348026, 14.199875, 14.094175,
       14.025088, 13.998063, 13.977668, 13.960636},
      {19.516939, 19.164715, 18.914882, 18.562661, 18.119021, // FREQ7
       17.826825, 17.517027, 17.164912, 16.915233, 16.721612,
       16.563450, 16.314234, 16.121209, 15.963850, 15.831103,
       15.716257, 15.615723, 15.525712, 15.371128, 15.184212,
       14.948958, 14.635891, 14.298034, 14.127792, 13.999929,
       13.912303, 13.876929, 13.849764, 13.826742},
      {19.516939, 19.164714, 18.914881, 18.562659, 18.119016, // FREQ8
       17.826815, 17.517006, 17.164865, 16.915150, 16.721482,
       16.563263, 16.313903, 16.120692, 15.963107, 15.830094,
       15.714942, 15.614066, 15.523672, 15.368212, 15.179720,
       14.941207, 14.619801, 14.262209, 14.073663, 13.925602,
       13.819464, 13.775217, 13.740590, 13.710759},
      {19.516939, 19.164714, 18.914881, 18.562657, 18.119012, // FREQ9
       17.826808, 17.516992, 17.164833, 16.915093, 16.721394,
       16.563135, 16.313676, 16.120337, 15.962597, 15.829401,
       15.714039, 15.612925, 15.522267, 15.366202, 15.176613,
       14.935812, 14.608414, 14.235819, 14.032225, 13.866132,
       13.741981, 13.688539, 13.645876, 13.608454},
      {19.516939, 19.164714, 18.914880, 18.562657, 18.119009, // FREQ10
       17.826803, 17.516982, 17.164810, 16.915052, 16.721330,
       16.563043, 16.313513, 16.120083, 15.962231, 15.828904,
       15.713391, 15.612108, 15.521260, 15.364758, 15.174377,
       14.931912, 14.600083, 14.215909, 13.999955, 13.817969,
       13.676626, 13.613849, 13.562692, 13.516939},
      {19.516939, 19.164713, 18.914880, 18.562656, 18.119008, // FREQ11
       17.826799, 17.516974, 17.164793, 16.915022, 16.721283,
       16.562976, 16.313392, 16.119895, 15.961961, 15.828537,
       15.712911, 15.611502, 15.520513, 15.363687, 15.172715,
       14.929003, 14.593814, 14.200566, 13.974434, 13.778545,
       13.621032, 13.548931, 13.488931, 13.434153},
      {19.516939, 19.164713, 18.914880, 18.562655, 18.119006, // FREQ12
       17.826796, 17.516969, 17.164780, 16.914999, 16.721247,
       16.562924, 16.313301, 16.119752, 15.961755, 15.828257,
       15.712546, 15.611041, 15.519944, 15.362870, 15.171447,
       14.926778, 14.588984, 14.188523, 13.953966, 13.745966,
       13.573403, 13.492115, 13.423028, 13.358576},
      {19.516939, 19.164713, 18.914880, 18.562655, 18.119005, // FREQ13
       17.826794, 17.516964, 17.164770, 16.914981, 16.721219,
       16.562884, 16.313230, 16.119641, 15.961595, 15.828039,
       15.712262, 15.610681, 15.519501, 15.362233, 15.170457,
       14.925038, 14.585188, 14.178914, 13.937343, 13.718804,
       13.532347, 13.442104, 13.363780, 13.289052},
      {19.516939, 19.164713, 18.914879, 18.562655, 18.119004, // FREQ14
       17.826792, 17.516961, 17.164762, 16.914967, 16.721197,
       16.562852, 16.313173, 16.119552, 15.961468, 15.827866,
       15.712036, 15.610396, 15.519149, 15.361728, 15.169670,
       14.923652, 14.582152, 14.171135, 13.923684, 13.695974,
       13.496762, 13.397869, 13.310243, 13.224682},
      {19.516939, 19.164713, 18.914879, 18.562654, 18.119003, // FREQ15
       17.826791, 17.516958, 17.164756, 16.914956, 16.721179,
       16.562826, 16.313127, 16.119481, 15.961365, 15.827726,
       15.711854, 15.610166, 15.518864, 15.361319, 15.169034,
       14.922532, 14.579688, 14.164756, 13.912343, 13.676639,
       13.465764, 13.358576, 13.261657, 13.164756}};
  static float XL[6][6][29] = {
      {{-30.274422, -29.048572, -28.181067, -26.962272, -25.437868, // X1s
        -24.444170, -23.404269, -22.248421, -21.454163, -20.858944,
        -20.390346, -19.694283, -19.200905, -18.835387, -18.556686,
        -18.339364, -18.168213, -18.030238, -17.826632, -17.633456,
        -17.461067, -17.322353, -17.245241, -17.223162, -17.211266,
        -17.204840, -17.202587, -17.200999, -17.199715},
       {-31.177414, -29.951530, -29.083850, -27.864712, -26.339031, // X2s
        -25.343652, -24.299685, -23.134693, -22.327692, -21.716473,
        -21.228927, -20.487480, -19.941059, -19.517455, -19.178033,
        -18.899376, -18.668043, -18.471683, -18.160149, -17.830286,
        -17.492277, -17.172499, -16.965517, -16.901255, -16.865263,
        -16.845632, -16.838714, -16.833696, -16.829681},
       {-31.705705, -30.479739, -29.612265, -28.392746, -26.866974, // X3s
        -25.871133, -24.826672, -23.659806, -22.850344, -22.235989,
        -21.744734, -20.993964, -20.435725, -19.998364, -19.643303,
        -19.347420, -19.097776, -18.881962, -18.529746, -18.137370,
        -17.701228, -17.231454, -16.873769, -16.748412, -16.674666,
        -16.633129, -16.617776, -16.606984, -16.598091},
       {-32.080641, -30.854674, -29.986801, -28.767697, -27.241693, // X4s
        -26.245685, -25.200974, -24.033538, -23.223063, -22.607845,
        -22.115266, -21.360872, -20.798453, -20.355878, -19.995174,
        -19.692644, -19.435600, -19.211713, -18.841933, -18.420428,
        -17.932110, -17.363567, -16.873130, -16.680219, -16.559751,
        -16.488746, -16.462241, -16.443053, -16.427763},
       {-32.371142, -31.145245, -30.277611, -29.058332, -27.532386, // X5s
        -26.536299, -25.491539, -24.323724, -23.512880, -22.897091,
        -22.403960, -21.648140, -21.083702, -20.638728, -20.275002,
        -19.969127, -19.708598, -19.480857, -19.102318, -18.665521,
        -18.148008, -17.516456, -16.921283, -16.663742, -16.492247,
        -16.386117, -16.345903, -16.316173, -16.291778},
       {-32.608820, -31.382756, -30.515126, -29.295866, -27.769793, // X6s
        -26.773814, -25.728819, -24.560932, -23.750086, -23.133811,
        -22.640288, -21.883631, -21.318035, -20.871913, -20.506426,
        -20.198858, -19.936428, -19.706400, -19.322760, -18.877373,
        -18.342274, -17.669792, -16.995256, -16.680122, -16.457336,
        -16.312694, -16.256489, -16.214113, -16.178612}},
      {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {-35.779538, -34.184208, -33.083933, -31.512708, -29.543604, // X2p
        -28.256123, -26.903279, -25.387738, -24.333408, -23.531477,
        -22.889415, -21.907557, -21.178842, -20.610306, -20.152156,
        -19.774043, -19.458248, -19.189136, -18.759267, -18.299831,
        -17.823327, -17.365980, -17.066362, -16.972218, -16.919695,
        -16.890892, -16.880696, -16.873357, -16.867478},
       {-36.234105, -34.655854, -33.538432, -31.967064, -29.997698, // X3p
        -28.709867, -27.356451, -25.839127, -24.782259, -23.977343,
        -23.331485, -22.340276, -21.599900, -21.017917, -20.544424,
        -20.149344, -19.815760, -19.527654, -19.058410, -18.538322,
        -17.967020, -17.364676, -16.918642, -16.765111, -16.675798,
        -16.625318, -16.607492, -16.594210, -16.583614},
       {-36.585694, -35.007703, -33.890016, -32.318668, -30.349350, // X4p
        -29.061334, -27.707618, -26.189677, -25.132040, -24.325956,
        -23.678826, -22.684226, -21.939671, -21.352566, -20.873369,
        -20.471723, -20.130813, -19.835172, -19.348733, -18.800381,
        -18.178384, -17.480038, -16.904760, -16.685329, -16.550262,
        -16.471169, -16.442151, -16.420831, -16.403759},
       {-36.866137, -35.287883, -34.170413, -32.599199, -30.629663, // X5p
        -29.341564, -27.987755, -26.469536, -25.411517, -24.604882,
        -23.957191, -22.961135, -22.214481, -21.625034, -21.142933,
        -20.738297, -20.393941, -20.094254, -19.599261, -19.036165,
        -18.385686, -17.626125, -16.948476, -16.665818, -16.480643,
        -16.367024, -16.324502, -16.292865, -16.266917},
       {-37.098169, -35.519950, -34.402525, -32.831070, -30.861699, // X6p
        -29.573885, -28.219694, -26.701459, -25.643044, -24.836230,
        -24.188105, -23.191275, -22.443490, -21.852666, -21.369042,
        -20.962634, -20.616374, -20.314553, -19.814673, -19.242970,
        -18.575541, -17.775947, -17.020568, -16.681448, -16.445735,
        -16.294606, -16.235710, -16.191866, -16.154983}},
      {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {-41.364414, -39.434006, -38.066663, -36.143204, -33.730242, // X3d
        -32.150245, -30.487089, -28.617809, -27.311427, -26.313205,
        -25.509946, -24.270587, -23.339149, -22.602299, -21.924436,
        -21.493723, -21.063954, -20.691590, -20.080654, -19.397357,
        -18.637161, -17.823176, -17.209853, -16.996234, -16.871214,
        -16.800539, -16.775144, -16.756765, -16.741919},
       {-41.585694, -39.655304, -38.288039, -36.364454, -33.951410, // X4d
        -32.371226, -30.707789, -28.837992, -27.530994, -26.531796,
        -25.727043, -24.484484, -23.549206, -22.807462, -22.198909,
        -21.686891, -21.250382, -20.870478, -20.243060, -19.532238,
        -18.722925, -17.815346, -17.075994, -16.798160, -16.628568,
        -16.529588, -16.493472, -16.467238, -16.445815},
       {-41.816885, -39.886598, -38.519116, -36.595706, -34.182651, // X5d
        -32.602365, -30.938792, -29.068803, -27.761491, -26.761551,
        -25.956256, -24.712472, -23.775049, -23.031086, -22.420027,
        -21.905038, -21.464940, -21.081321, -20.445565, -19.720393,
        -18.883701, -17.916497, -17.077571, -16.738117, -16.519620,
        -16.387033, -16.337715, -16.301341, -16.271391},
       {-42.024362, -40.094064, -38.726686, -36.803137, -34.390124, // X6d
        -32.809866, -31.146180, -29.276029, -27.968300, -26.968324,
        -26.162701, -24.918051, -23.979662, -23.234506, -22.621799,
        -22.105162, -21.663212, -21.277514, -20.637026, -19.903484,
        -19.050185, -18.044511, -17.129904, -16.735338, -16.467566,
        -16.298269, -16.232977, -16.184230, -16.143922}},
      {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {-47.062815, -44.780358, -43.163100, -40.887314, -38.030685, // X4f
        -36.158301, -34.185235, -31.963719, -30.407089, -29.214529,
        -28.252197, -26.761810, -25.634821, -24.737662, -23.998757,
        -23.374580, -22.839980, -22.373323, -21.598611, -20.713453,
        -19.693804, -18.530997, -17.563112, -17.193424, -16.965517,
        -16.832288, -16.783370, -16.747717, -16.718672},
       {-47.128880, -44.846322, -43.229046, -40.953347, -38.096716, // X5f
        -36.224291, -34.250943, -32.029199, -30.472360, -29.279276,
        -28.316408, -26.824527, -25.695751, -24.796176, -24.054627,
        -23.427631, -22.889877, -22.419401, -21.636478, -20.737351,
        -19.690904, -18.469715, -17.404053, -16.973748, -16.697901,
        -16.531879, -16.469784, -16.423961, -16.386588},
       {-47.267412, -44.984913, -43.367636, -41.091842, -38.235239, // X6f
        -36.362731, -34.389528, -32.167518, -30.610443, -29.417223,
        -28.453971, -26.961283, -25.831491, -24.930907, -24.187725,
        -23.559075, -23.019383, -22.547066, -21.759545, -20.852145,
        -19.789541, -18.530522, -17.390884, -16.906727, -16.582667,
        -16.380139, -16.302886, -16.245236, -16.197380}},
      {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {-52.894711, -50.260082, -48.392958, -45.765034, -42.464679, // X5g
        -40.300146, -38.017153, -35.443424, -33.636754, -32.250427,
        -31.129593, -29.389103, -28.068001, -27.012118, -26.138711,
        -25.398332, -24.761042, -24.202462, -23.268415, -22.188504,
        -20.919298, -19.415147, -18.073478, -17.521544, -17.163795,
        -16.946562, -16.865194, -16.805098, -16.755865},
       {-52.845039, -50.210247, -48.343069, -45.715131, -42.414728, // X6g
        -40.250164, -37.967149, -35.393156, -33.586496, -32.199833,
        -31.078643, -29.337458, -27.969702, -26.958401, -26.083595,
        -25.341555, -24.702345, -24.141808, -23.203287, -22.115356,
        -20.830007, -19.288694, -17.874057, -17.268729, -16.863465,
        -16.610369, -16.513883, -16.442010, -16.382570}},
      {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {-58.850334, -55.863542, -53.746437, -50.766409, -47.022317, // X6h
        -44.565391, -41.972509, -39.046704, -36.990356, -35.410261,
        -34.131188, -32.140740, -30.626018, -29.411767, -28.404701,
        -27.548439, -26.808936, -26.159088, -25.067378, -23.795088,
        -22.279431, -20.436907, -18.711058, -17.957760, -17.446882,
        -17.124901, -17.001376, -16.909196, -16.832806}}};
  static float EKARZAS[29] = {10000., 4444., 2500., 1111., 400., 204.1, 100., 44.44,
                              25., 16., 11.11, 6.25, 4., 2.778, 2.041, 1.562, 1.235, 1., 0.6944, 0.4444,
                              0.25, 0.1111, 0.04, 0.02041, 0.01, 0.004444, 0.0025, 0.001111, 0.};
  double FREQLG, X, FREQN15[29];
  int I;

  FREQLG = log10(FREQ / ZEFF2);
  if (N <= 15)
  {
    if (L >= N || N > 6)
    {
      if (FREQLG < FREQN[N - 1][28])
        return 0.;
      for (I = 2; I < 30; I++)
      {
        if (FREQLG > FREQN[N - 1][I - 1])
          break;
      }
      X = (FREQLG - FREQN[N - 1][I - 1]) / (FREQN[N - 1][I - 2] - FREQN[N - 1][I - 1]) *
              (XN[N - 1][I - 2] - XN[N - 1][I - 1]) +
          XN[N - 1][I - 1];
      return exp(X * 2.30258509299405e0) / ZEFF2;
    }
    if (FREQLG < FREQN[N - 1][28])
      return 0.;

    for (I = 2; I < 30; I++)
    {
      if (FREQLG > FREQN[N - 1][I - 1])
        break;
    }
    X = (FREQLG - FREQN[N - 1][I - 1]) / (FREQN[N - 1][I - 2] - FREQN[N - 1][I - 1]) *
            (XL[L][N - 1][I - 2] - XL[L][N - 1][I - 1]) +
        XL[L][N - 1][I - 1];
    return exp(X * 2.30258509299405e0) / ZEFF2;
  }

  FREQN15[28] = log10(109677.576 * 2.99792458e10 / (N * N));
  if (FREQLG < FREQN15[28])
    return 0.;
  for (I = 2; I < 29; I++)
  {
    FREQN15[I - 1] = log10((EKARZAS[I - 1] + 1. / (N * N)) * 109677.576 * 2.99792458e10);
    if (FREQLG > FREQN15[I - 1])
      break;
  }

  X = (FREQLG - FREQN15[I - 1]) / (FREQN15[I - 2] - FREQN15[I - 1]) *
          (XN[14][I - 2] - XN[14][I - 1]) +
      XN[14][I - 1];
  return exp(X * 2.30258509299405e0) / ZEFF2;
}

double COULX(int N, double freq, double Z, GlobalState *state)
{
  static double A[6] = {0.9916, 1.105, 1.101, 1.101, 1.102, 1.0986},
                B[6] = {2.719e3, -2.375e4, -9.863e3, -5.765e3, -3.909e3, -2.704e3},
                C[6] = {-2.268e10, 4.077e8, 1.035e8, 4.593e7, 2.371e7, 1.229e7};
  double CLX, FREQ1;
  int n;

  n = (N + 1) * (N + 1);
  if (freq < Z * Z * 3.28805e15 / n)
    return 0.;

  FREQ1 = freq * 1.e-10;
  CLX = 0.2815 / FREQ1 / FREQ1 / FREQ1 / n / n / (N + 1) * Z * Z * Z * Z;
  if (N >= 6)
    return CLX;
  if (N == 0)
  {
    CLX *= COULBF1S(state->FREQ, Z);
    return CLX;
  }
  CLX *= (A[N] + (B[N] + C[N] * (Z * Z / FREQ1)) * (Z * Z / FREQ1));
  return CLX;
}

double COULFF(int J, int NZ, GlobalState *state)
{
  static double Z4LOG[6] = {0., 1.20412, 1.90849, 2.40824, 2.79588, 3.11261},
                A[12][11] = {
                    {5.53, 5.49, 5.46, 5.43, 5.40, 5.25, 5.00, 4.69, 4.48, 4.16, 3.85},
                    {4.91, 4.87, 4.84, 4.80, 4.77, 4.63, 4.40, 4.13, 3.87, 3.52, 3.27},
                    {4.29, 4.25, 4.22, 4.18, 4.15, 4.02, 3.80, 3.57, 3.27, 2.98, 2.70},
                    {3.64, 3.61, 3.59, 3.56, 3.54, 3.41, 3.22, 2.97, 2.70, 2.45, 2.20},
                    {3.00, 2.98, 2.97, 2.95, 2.94, 2.81, 2.65, 2.44, 2.21, 2.01, 1.81},
                    {2.41, 2.41, 2.41, 2.41, 2.41, 2.32, 2.19, 2.02, 1.84, 1.67, 1.50},
                    {1.87, 1.89, 1.91, 1.93, 1.95, 1.90, 1.80, 1.68, 1.52, 1.41, 1.30},
                    {1.33, 1.39, 1.44, 1.49, 1.55, 1.56, 1.51, 1.42, 1.33, 1.25, 1.17},
                    {0.90, 0.95, 1.00, 1.08, 1.17, 1.30, 1.32, 1.30, 1.20, 1.15, 1.11},
                    {0.55, 0.58, 0.62, 0.70, 0.85, 1.01, 1.15, 1.18, 1.15, 1.11, 1.08},
                    {0.33, 0.36, 0.39, 0.46, 0.59, 0.76, 0.97, 1.09, 1.13, 1.10, 1.08},
                    {0.19, 0.21, 0.24, 0.28, 0.38, 0.53, 0.76, 0.96, 1.08, 1.09, 1.09}};
  double GAMLOG, HVKTLG, P, Q, CLFF;
  int IGAM, IHVKT;

  GAMLOG = 10.39638 - state->TLOG[J] / 1.15129 + Z4LOG[NZ - 1];
  IGAM = min((int)(GAMLOG + 7.), 10);
  if (IGAM < 1)
    IGAM = 1;

  HVKTLG = (state->FREQLG - state->TLOG[J]) / 1.15129 - 20.63764;
  IHVKT = min((int)(HVKTLG + 9.), 11);
  if (IHVKT < 1)
    IHVKT = 1;
  P = GAMLOG - (IGAM - 7);
  Q = HVKTLG - (IHVKT - 9);
  CLFF = (1. - P) * ((1. - Q) * A[IHVKT - 1][IGAM - 1] + Q * A[IHVKT][IGAM - 1]) +
         P * ((1. - Q) * A[IHVKT - 1][IGAM] + Q * A[IHVKT][IGAM]);
  return CLFF;
}

void HOP(double *ahyd, int iH1, int iH2, GlobalState *state) /* REQUIRES FUNCTIONS COULX AND COULFF */
{
  double BOLT[MOSIZE][8], EXLIM[MOSIZE], FREET[MOSIZE], BOLTEX[MOSIZE];
  double CONT[8], H, CFREE, XR, EX, C, nH1;
  int J, N;

  for (J = 0; J < state->NRHOX; J++)
  {
    nH1 = state->FRACT[J][iH1];
    for (N = 0; N < 8; N++)
      BOLT[J][N] = exp(-13.595 * (1. - 1. / (N + 1) / (N + 1)) / state->TKEV[J]) *
                   2. * (N + 1) * (N + 1) * nH1 / state->RHO[J];
    FREET[J] = state->XNE[J] * state->FRACT[J][iH2] / (sqrt(state->T[J]) * state->RHO[J]);
    XR = nH1 / 13.595 * state->TKEV[J] / state->RHO[J];
    BOLTEX[J] = exp(-13.427 / state->TKEV[J]) * XR;
    EXLIM[J] = exp(-13.595 / state->TKEV[J]) * XR;
  }
  for (N = 0; N < 8; N++)
    CONT[N] = COULX(N, state->FREQ, 1., state);
  CFREE = 3.6919e8 / (state->FREQ * state->FREQ);
  C = ((2.815e29 / state->FREQ) / state->FREQ) / state->FREQ;
  for (J = 0; J < state->NRHOX; J++)
  {
    EX = BOLTEX[J];
    if (state->FREQ < 4.05933e13)
      EX = EXLIM[J] / state->EHVKT[J];
    H = (CONT[6] * BOLT[J][6] + CONT[7] * BOLT[J][7] + (EX - EXLIM[J]) * C +
         COULFF(J, 1, state) * FREET[J] / state->FREQ * CFREE) *
        state->STIM[J];
    for (N = 0; N < 6; N++)
      H += CONT[N] * BOLT[J][N] * (1. - state->EHVKT[J]);
    ahyd[J] = H;
  }
  return;
}

void HRAYOP(double *sigh, int iH1, GlobalState *state)
{
  double WAVE, WW, SIG, nH1;
  int J;

  WAVE = CLIGHT / min(state->FREQ, 2.463e15); // Wavelength in Angstroems
  WW = WAVE * WAVE;
  SIG = (5.799e-13 + 1.422e-6 / WW + 2.784 / (WW * WW)) / (WW * WW);
  for (J = 0; J < state->NRHOX; J++)
    sigh[J] = SIG * state->FRACT[J][iH1] * 2. / state->RHO[J];
  return;
}

void H2PLOP(double *ah2p, int iH1, int iH2, GlobalState *state)
{
  double FR, ES, FREQ15, nH1;
  int J;

  if (state->FREQ > 3.28805e15)
    return;
  FR = -3.0233e3 + (3.7797e2 + (-1.82496e1 + (3.9207e-1 - 3.1672e-3 * state->FREQLG) *
                                                 state->FREQLG) *
                                   state->FREQLG) *
                       state->FREQLG;
  FREQ15 = state->FREQ * 1.e-15;
  ES = -7.342e-3 + (-2.409 + (1.028 + (-0.4230 + (0.1224 - 0.01351 * FREQ15) *
                                                     FREQ15) *
                                          FREQ15) *
                                 FREQ15) *
                       FREQ15;
  for (J = 0; J < state->NRHOX; J++)
  {
    //  ah2p[J]=exp(-ES/state->TKEV[J]+FR)*2.*state->FRACT[J][iH1]*state->FRACT[J][iH2]/state->RHO[J]*state->STIM[J];
    nH1 = state->FRACT[J][iH1] * 2;
    ah2p[J] = exp(-ES / state->TKEV[J] + FR) * nH1 * state->FRACT[J][iH2] / state->RHO[J] * state->STIM[J];
    //    printf("%d %10.5g %10.5g\n",J,ah2p[J],state->STIM[J]);
  }
  return;
}

void HMINOP_old(double *ahmin, int iH1, int iHmin, GlobalState *state)
{
  double HMINBF, HMINFF, H, FREQ1, B, C, HMINFR, nH1;
  int J;

  FREQ1 = state->FREQ * 1.e-10;
  B = (1.3727e-15 + 4.3748 / state->FREQ) / FREQ1;
  C = -2.5993e-7 / (FREQ1 * FREQ1);
  if (state->FREQ <= 1.8259e14)
    HMINBF = 0.;
  else if (state->FREQ >= 2.111e14)
    HMINBF = 6.801e-10 + (5.358e-3 + (1.481e3 + (-5.519e7 +
                                                 4.808e11 / FREQ1) /
                                                    FREQ1) /
                                         FREQ1) /
                             FREQ1;
  else
    HMINBF = 3.695e-6 + (-1.251e-1 + 1.052e3 / FREQ1) / FREQ1;
  for (J = 0; J < state->NRHOX; J++)
  {
    //    HMINFF=(B+C/state->T[J])*state->FRACT[J][iH1]*state->XNE[J]*2.e-20;
    nH1 = state->FRACT[J][iH1] * 2;
    HMINFF = (B + C / state->T[J]) * nH1 * state->XNE[J] * 1.e-20;
    if (state->T[J] > 7730.)
      HMINFR = exp(0.7552 / state->TKEV[J]) / (2. * 2.4148E15 * state->T[J] * sqrt(state->T[J])) * state->FRACT[J][iH1] * state->XNE[J];
    // Bug fixed 2007-12-15: Partition function of H- is 1 and not 2 as we used
    // before:
    else
      HMINFR = state->FRACT[J][iHmin];
    //    printf("state->T: %10.1f Kurucz: %11.6e EOS: %11.6e\n",state->T[J],
    //                          exp(0.7552/state->TKEV[J])/(2.*2.4148E15*state->T[J]*
    //                          sqrt(state->T[J]))*state->FRACT[J][iH1]*state->XNE[J],state->FRACT[J][iHmin]);
    H = HMINBF * (1. - state->EHVKT[J]) * HMINFR * 1.e-10;
    ahmin[J] = (H + HMINFF) / state->RHO[J];
  }
  return;
}

void HMINOP(double *ahmin, int iH1, int iHmin, GlobalState *state)
{
  //From Mathisen (1984), after Wishart (1979) and  Broad & Reinhardt (1976)
  static double WBF[85] = {18.00, 19.60, 21.40, 23.60, 26.40, 29.80, 34.30,
                           40.40, 49.10, 62.60, 111.30, 112.10, 112.67, 112.95, 113.05,
                           113.10, 113.20, 113.23, 113.50, 114.40, 121.00, 139.00, 164.00,
                           175.00, 200.00, 225.00, 250.00, 275.00, 300.00, 325.00, 350.00,
                           375.00, 400.00, 425.00, 450.00, 475.00, 500.00, 525.00, 550.00,
                           575.00, 600.00, 625.00, 650.00, 675.00, 700.00, 725.00, 750.00,
                           775.00, 800.00, 825.00, 850.00, 875.00, 900.00, 925.00, 950.00,
                           975.00, 1000.00, 1025.00, 1050.00, 1075.00, 1100.00, 1125.00, 1150.00,
                           1175.00, 1200.00, 1225.00, 1250.00, 1275.00, 1300.00, 1325.00, 1350.00,
                           1375.00, 1400.00, 1425.00, 1450.00, 1475.00, 1500.00, 1525.00, 1550.00,
                           1575.00, 1600.00, 1610.00, 1620.00, 1630.00, 1643.91};
  static double BF[85] = {0.067, 0.088, 0.117, 0.155, 0.206, 0.283, 0.414,
                          0.703, 1.24, 2.33, 11.60, 13.90, 24.30, 66.70, 95.00,
                          56.60, 20.00, 14.60, 8.50, 7.10, 5.43, 5.91, 7.29,
                          7.918, 9.453, 11.08, 12.75, 14.46, 16.19, 17.92, 19.65,
                          21.35, 23.02, 24.65, 26.24, 27.77, 29.23, 30.62, 31.94,
                          33.17, 34.32, 35.37, 36.32, 37.17, 37.91, 38.54, 39.07,
                          39.48, 39.77, 39.95, 40.01, 39.95, 39.77, 39.48, 39.06,
                          38.53, 37.89, 37.13, 36.25, 35.28, 34.19, 33.01, 31.72,
                          30.34, 28.87, 27.33, 25.71, 24.02, 22.26, 20.46, 18.62,
                          16.74, 14.85, 12.95, 11.07, 9.211, 7.407, 5.677, 4.052,
                          2.575, 1.302, 0.8697, 0.4974, 0.1989, 0.};
  // Bell and Berrington J.Phys.B,vol. 20, 801-806,1987.
  static double WAVEK[22] = {.50, .40, .35, .30, .25, .20, .18, .16, .14, .12, .10, .09,
                             .08, .07, .06, .05, .04, .03, .02, .01, .008, .006};
  static double THETAFF[11] = {
      0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6};
  static double FF[22][11] = {
      //  FFBEG=
      {.0178, .0222, .0308, .0402, .0498, .0596, .0695, .0795, .0896, .131, .172}, //  1823
      {.0228, .0280, .0388, .0499, .0614, .0732, .0851, .0972, .110, .160, .211},  //  2278
      {.0277, .0342, .0476, .0615, .0760, .0908, .105, .121, .136, .199, .262},    //  2604
      {.0364, .0447, .0616, .0789, .0966, .114, .132, .150, .169, .243, .318},     //  3038
      {.0520, .0633, .0859, .108, .131, .154, .178, .201, .225, .321, .418},       //  3645
      {.0791, .0959, .129, .161, .194, .227, .260, .293, .327, .463, .602},        //  4557
      {.0965, .117, .157, .195, .234, .272, .311, .351, .390, .549, .711},         //  5063
      {.121, .146, .195, .241, .288, .334, .381, .428, .475, .667, .861},          //  5696
      {.154, .188, .249, .309, .367, .424, .482, .539, .597, .830, 1.07},          //  6510
      {.208, .250, .332, .409, .484, .557, .630, .702, .774, 1.06, 1.36},          //  7595
      {.293, .354, .468, .576, .677, .777, .874, .969, 1.06, 1.45, 1.83},          //  9113
                                                                                   // FFEND=
      {.358, .432, .572, .702, .825, .943, 1.06, 1.17, 1.28, 1.73, 2.17},          // 10126
      {.448, .539, .711, .871, 1.02, 1.16, 1.29, 1.43, 1.57, 2.09, 2.60},          // 11392
      {.579, .699, .924, 1.13, 1.33, 1.51, 1.69, 1.86, 2.02, 2.67, 3.31},          // 13019
      {.781, .940, 1.24, 1.52, 1.78, 2.02, 2.26, 2.48, 2.69, 3.52, 4.31},          // 15189
      {1.11, 1.34, 1.77, 2.17, 2.53, 2.87, 3.20, 3.51, 3.80, 4.92, 5.97},          // 18227
      {1.73, 2.08, 2.74, 3.37, 3.90, 4.50, 5.01, 5.50, 5.95, 7.59, 9.06},          // 22784
      {3.04, 3.65, 4.80, 5.86, 6.86, 7.79, 8.67, 9.50, 10.3, 13.2, 15.6},          // 30378
      {6.79, 8.16, 10.7, 13.1, 15.3, 17.4, 19.4, 21.2, 23.0, 29.5, 35.0},          // 45567
      {27.0, 32.4, 42.6, 51.9, 60.7, 68.9, 76.8, 84.2, 91.4, 117., 140.},          // 91134
      {42.3, 50.6, 66.4, 80.8, 94.5, 107., 120., 131., 142., 183., 219.},          //113918
      {75.1, 90.0, 118., 144., 168., 191., 212., 234., 253., 325., 388.}};         //151890

  double WFFLOG[22], FFLOG[11][22], FFTT[11], THETA[MOSIZE], FFTHETA[MOSIZE];
  double WAVE[1], WAVELOG[1], XHMIN[MOSIZE], FFTLOG[1], H, HMINBF[1], HMINFF;
  int J, IWAVE, ITHETA, MAXWAVE;

  for (IWAVE = 0; IWAVE < 22; IWAVE++)
  {
    // 91.134 number taken from Bell and Berrington
    WFFLOG[IWAVE] = log(91.134e0 / WAVEK[IWAVE]);
    for (ITHETA = 0; ITHETA < 11; ITHETA++)
      FFLOG[ITHETA][IWAVE] = log(FF[IWAVE][ITHETA] * 1.e-26);
  }

  for (J = 0; J < state->NRHOX; J++)
  {
    THETA[J] = 5040. / state->T[J];
    // .754209 Hotop and Lineberger J. Phys. Chem. Ref. Data Vol. 14, 731-752, 1985
    XHMIN[J] = exp(0.754209 / state->TKEV[J]) / (2. * 2.4148e15 * state->T[J] * sqrt(state->T[J])) * state->FRACT[J][iH1] * state->XNE[J];
  }
  WAVE[0] = CLIGHT / state->FREQ * 0.1; // Wavelength in nanometers
  WAVELOG[0] = log(WAVE[0]);
  for (ITHETA = 0; ITHETA < 11; ITHETA++)
  {
    LINTER(WFFLOG, FFLOG[ITHETA], 22, WAVELOG, FFTLOG, 1);
    FFTT[ITHETA] = exp(FFTLOG[0]) / THETAFF[ITHETA] * 5040. * 1.380658e-16;
  }

  HMINBF[0] = 0.;
  if (state->FREQ > 1.82365E14)
    MAXWAVE = MAP1(WBF, BF, 85, WAVE, HMINBF, 1);
  for (J = 0; J < state->NRHOX; J++)
  {
    LINTER(THETAFF, FFTT, 11, THETA + J, FFTHETA + J, 1);
    HMINFF = FFTHETA[J] * state->FRACT[J][iH1] * 2. * state->XNE[J] / state->RHO[J];
    //    H=HMINBF[0]*1.e-18*(1.-state->EHVKT[J])*XHMIN[J]/state->RHO[J];
    H = HMINBF[0] * 1.e-18 * (1. - state->EHVKT[J]) * state->FRACT[J][iHmin] * state->PARTITION_FUNCTIONS[J][iHmin] / state->RHO[J];
    ahmin[J] = H + HMINFF;
  }
  return;
}

void HE1OP(double *ahe1, int iHe1, int iHe2, GlobalState *state) /* REQUIRES FUNCTION COULFF. Needs update!!! */
{
  double BOLT[MOSIZE][10], EXLIM[MOSIZE], BOLTEX[MOSIZE], FREET[MOSIZE], TRANS[10];
  double FREQ3, CFREE, C, HE1, EX, XRLOG;
  static double G[10] = {1., 3., 1., 9., 3., 3., 1., 9., 20., 3.},
                HEFREQ[10] = {5.9452090E15, 1.1528440E15, 0.9603331E15, 0.8761076E15,
                              0.8147104E15, 0.4519048E15, 0.4030971E15, 0.3821191E15,
                              0.3660215E15, 0.3627891E15},
                CHI[10] = {0., 19.819, 20.615, 20.964, 21.217, 22.718, 22.920, 23.006,
                           23.073, 23.086};
  int J, N, NMIN, IMIN;

  for (J = 0; J < state->NRHOX; J++)
  {
    for (N = 0; N < 10; N++)
    {
      BOLT[J][N] = exp(-CHI[N] / state->TKEV[J] + log(state->FRACT[J][iHe1]) - log(state->RHO[J])) * G[N];
    }
    FREET[J] = state->XNE[J] * 1.e-10 * state->FRACT[J][iHe2] * 1.e-10 / state->RHO[J] / sqrt(state->T[J]) * 1.e-10;
    /*    XRLOG=log(state->FRACT[J][iHe1]*(4/2/13.595)*state->TKEV[J]/state->RHO[J]); */
    XRLOG = log(state->FRACT[J][iHe1] * (2. / 13.595) * state->TKEV[J] / state->RHO[J]);
    BOLTEX[J] = exp(-23.730 / state->TKEV[J] + XRLOG);
    EXLIM[J] = exp(-24.587 / state->TKEV[J] + XRLOG);
  }
  FREQ3 = state->FREQ * 1.e-10;
  FREQ3 = FREQ3 * FREQ3 * FREQ3;
  CFREE = 3.6919e8 / FREQ3;
  C = 2.815e-1 / FREQ3;
  for (NMIN = 0; NMIN < 10; NMIN++)
  {
    TRANS[NMIN] = 0;
    IMIN = NMIN;
    if (HEFREQ[NMIN] <= state->FREQ)
      break;
  }

  switch (IMIN)
  {
  case 0:
    TRANS[0] = exp(33.32e0 - 2.e0 * state->FREQLG);
  case 1:
    TRANS[1] = exp(-390.026e0 + (21.035e0 - 0.318e0 * state->FREQLG) * state->FREQLG);
  case 2:
    TRANS[2] = exp(26.83e0 - 1.91e0 * state->FREQLG);
  case 3:
    TRANS[3] = exp(61.21e0 - 2.9e0 * state->FREQLG);
  case 4:
    TRANS[4] = exp(81.35e0 - 3.5e0 * state->FREQLG);
  case 5:
    TRANS[5] = exp(12.69e0 - 1.54e0 * state->FREQLG);
  case 6:
    TRANS[6] = exp(23.85e0 - 1.86e0 * state->FREQLG);
  case 7:
    TRANS[7] = exp(49.30e0 - 2.60e0 * state->FREQLG);
  case 8:
    TRANS[8] = exp(85.20e0 - 3.69e0 * state->FREQLG);
  case 9:
    TRANS[9] = exp(58.81e0 - 2.89e0 * state->FREQLG);
  default:
    break;
  }

  for (J = 0; J < state->NRHOX; J++)
  {
    EX = BOLTEX[J];
    if (state->FREQ < 2.055e14)
      EX = EXLIM[J] / state->EHVKT[J];
    HE1 = (EX - EXLIM[J]) * C;
    for (N = 0; N < 10; N++)
      HE1 += TRANS[N] * BOLT[J][N];
    ahe1[J] = (HE1 + COULFF(J, 1, state) * FREET[J] * CFREE) * state->STIM[J];
  }
  return;
}

double CROSSHE(double FREQ)
{
  // Marr, G.V. and West, J.B. Atomic Data and Nuclear Data Tables,
  //   vol 18, 497-508, 1976.
  static double X505[92] = {7.58, 7.46, 7.33, 7.19, 7.06, 6.94, 6.81,
                            6.68, 6.55, 6.43, 6.30, 6.18, 6.05, 5.93, 5.81, 5.69, 5.57,
                            5.45, 5.33, 5.21, 5.10, 4.98, 4.87, 4.76, 4.64, 4.53, 4.42,
                            4.31, 4.20, 4.09, 4.00, 3.88, 3.78, 3.68, 3.57, 3.47, 3.37,
                            3.27, 3.18, 3.08, 2.98, 2.89, 2.80, 2.70, 2.61, 2.52, 2.44,
                            2.35, 2.26, 2.18, 2.10, 2.02, 1.94, 1.86, 1.78, 1.70, 1.63,
                            1.55, 1.48, 1.41, 1.34, 1.28, 1.21, 1.14, 1.08, 1.02, .961,
                            .903, .847, .792, .738, .687, .637, .588, .542, .497, .454,
                            .412, .373, .335, .299, .265, .233, .202, .174, .147, .123,
                            .100, .0795, .0609, .0443, .0315},
                X50[16] = {.0315, .0282, .0250, .0220, .0193, .0168, .0145, .0124,
                           .0105, .00885, .00736, .00604, .00489, .00389, .00303, .00231},
                X20[11] = {.00231, .00199, .00171, .00145, .00122, .00101, .000832,
                           .000673, .000535, .000417, .000318},
                X10[21] = {.000318, .000274, .000235, .000200, .000168, .000139, .000115,
                           .000093, .000074, .000057, .000044, .000032, .000023, .000016, .000010,
                           .000006, .000003, .000001, .0000006, .0000003, 0.};
  double WAVE;
  int i;

  if (FREQ < 5.945209e15)
    return 0.;
  WAVE = CLIGHT / FREQ;
  if (WAVE > 50.)
  {
    i = 93. - (WAVE - 50.) / 5.;
    i = min(92, max(2, i));
    return ((WAVE - (92 - i) * 5 - 50) / 5. * (X505[i - 2] - X505[i - 1]) + X505[i - 1]) * 1.e-18;
  }
  if (WAVE > 20.)
  {
    i = 17. - (WAVE - 20.) / 2.;
    i = min(16, max(2, i));
    return ((WAVE - (16 - i) * 2 - 20) / 2. * (X50[i - 2] - X50[i - 1]) + X50[i - 1]) * 1.e-18;
  }
  if (WAVE > 10.)
  {
    i = 12. - (WAVE - 10.) / 1.;
    i = min(11, max(2, i));
    return ((WAVE - (11 - i) * 1 - 10) / 1. * (X20[i - 2] - X20[i - 1]) + X20[i - 1]) * 1.e-18;
  }
  i = 22. - WAVE / 0.5;
  i = min(21, max(2, i));
  return ((WAVE - (21 - i) * 0.5) / 0.5 * (X10[i - 2] - X10[i - 1]) + X10[i - 1]) * 1.e-18;
}

double HE111S(double FREQ)
{
  // Following Mathisen
  static double W[64] = {
      504.3, 501.5, 498.7, 493.3, 488.1, 480.3, 477.8, 454.0, 443.0,
      395.0, 356.4, 348.2, 324.6, 302.0, 298.1, 275.6, 260.6, 256.2,
      239.4, 224.6, 220., 215, 210., 205., 200., 195., 190.,
      185., 180., 175., 170., 165., 160., 155., 150., 145.,
      135., 130., 125., 120., 115., 110., 105., 100., 95.,
      90., 85., 80., 75., 70., 65., 60., 55., 50.,
      45., 40., 35., 30., 25., 20., 15., 10., 5., 0.},
                X[64] = {7.346, 7.317, 7.259, 7.143, 7.030, 6.857, 6.800, 6.284, 6.041, 4.977, 4.138, 3.961, 3.474, 3.025, 2.945, 2.522, 2.259, 2.179, 1.901, 1.684, 1.61, 1.53, 1.45, 1.38, 1.30, 1.22, 1.14, 1.08, 1.02, 0.961, 0.903, 0.847, 0.792, 0.738, 0.687, 0.637, 0.542, 0.497, 0.454, 0.412, 0.373, 0.335, 0.299, 0.265, 0.233, 0.202, 0.174, 0.147, 0.124, 0.103, 0.0840, 0.0676, 0.0535, 0.0414, .0311, .0266, .0158, .0104, .00637, .00349, .00161, .00054, .000083, 0.};
  double WAVE;
  int i;

  if (FREQ < 5.945209e15)
    return 0.;
  WAVE = CLIGHT / FREQ;
  for (i = 1; i < 64; i++)
    if (WAVE > W[i])
      break;
  return ((WAVE - W[i]) / (W[i - 1] - W[i]) * (X[i - 1] - X[i]) + X[i]) * 1.e-18;
}

double HE12s1S(double FREQ)
{
  static double FREQ1S[16] = {
      15.947182, 15.913654, 15.877320, 15.837666, 15.794025,
      15.745503, 15.690869, 15.628361, 15.555317, 15.467455,
      15.357189, 15.289399, 15.251073, 15.209035, 15.162487,
      14.982421},
                X1S[16] = {-19.635557, -19.159345, -18.958474, -18.809535, -18.676481, -18.546006, -18.410962, -18.264821, -18.100205, -17.909165, -17.684370, -17.557867, -17.490360, -17.417876, -17.349386, -17.084441};
  double FREQLG10, WAVENO, EK, EPS, X;
  int i;

  if (FREQ < 32033.214e0 * CLIGHTcm)
    return 0;

  if (FREQ > 2.4 * 109722.267e0 * CLIGHTcm)
  {
    WAVENO = FREQ / CLIGHTcm;
    EK = (WAVENO - 32033.214e0) / 109722.267e0;
    EPS = 2. * (EK - 2.612316e0) / 0.00322e0;
    return 0.008175e0 * pow(484940. / WAVENO, 2.71) * 8.067e-18 *
           (EPS + 76.21) * (EPS + 76.21) / (1. + EPS * EPS);
  }

  FREQLG10 = log10(FREQ);
  for (i = 1; i < 16; i++)
    if (FREQLG10 > FREQ1S[i])
      break;
  X = (FREQLG10 - FREQ1S[i]) / (FREQ1S[i - 1] - FREQ1S[i]) *
          (X1S[i - 1] - X1S[i]) +
      X1S[i];
  return pow10(X);
}

double HE12s3S(double FREQ)
{
  static double FREQ3S[16] = {
      15.956523, 15.923736, 15.888271, 15.849649, 15.807255,
      15.760271, 15.707580, 15.647601, 15.577992, 15.495055,
      15.392451, 15.330345, 15.295609, 15.257851, 15.216496,
      15.061770},
                X3S[16] = {-18.426022, -18.610700, -18.593051, -18.543304, -18.465513, -18.378707, -18.278574, -18.164329, -18.033346, -17.882435, -17.705542, -17.605584, -17.553459, -17.500667, -17.451318, -17.266686};
  double FREQLG10, WAVENO, EK, EPS, X;
  int i;

  if (FREQ < 38454.691 * CLIGHTcm)
    return 0.;

  if (FREQ > 2.4 * 109722.267 * CLIGHTcm)
  {
    WAVENO = FREQ / CLIGHTcm;
    EK = (WAVENO - 38454.691e0) / 109722.267e0;
    EPS = 2. * (EK - 2.47898e0) / 0.000780e0;
    return 0.01521e0 * pow(470310.e0 / WAVENO, 3.12) *
           8.067e-18 * (EPS - 122.4e0) * (EPS - 122.4e0) / (1. + EPS * EPS);
  }

  FREQLG10 = log10(FREQ);
  for (i = 1; i < 16; i++)
    if (FREQLG10 > FREQ3S[i])
      break;
  X = (FREQLG10 - FREQ3S[i]) / (FREQ3S[i - 1] - FREQ3S[i]) *
          (X3S[i - 1] - X3S[i]) +
      X3S[i];
  return pow10(X);
}

double HE12p1P(double FREQ)
{
  static double FREQ1P[16] = {
      15.939981, 15.905870, 15.868850, 15.828377, 15.783742,
      15.733988, 15.677787, 15.613218, 15.537343, 15.445346,
      15.328474, 15.255641, 15.214064, 15.168081, 15.116647,
      14.911002},
                X1P[16] = {-18.798876, -19.685922, -20.011664, -20.143030, -20.091354, -19.908333, -19.656788, -19.367745, -19.043016, -18.674484, -18.240861, -17.989700, -17.852015, -17.702677, -17.525347, -16.816344};
  double FREQLG10, WAVENO, X, EK, EPS1S, EPS1D;
  int i;

  if (FREQ < 27175.76 * CLIGHTcm)
    return 0;

  if (FREQ > 2.4 * 109722.267 * CLIGHTcm)
  {
    WAVENO = FREQ / CLIGHTcm;
    EK = (WAVENO - 27175.76e0) / 109722.267e0;
    EPS1S = 2. * (EK - 2.446534e0) / 0.01037e0;
    EPS1D = 2. * (EK - 2.59427e0) / 0.00538e0;
    return 0.9487e-3 * pow(466750. / WAVENO, 3.69) * 8.067e-18 *
           ((EPS1S - 29.30) * (EPS1S - 29.30) / (1. + EPS1S * EPS1S) +
            (EPS1D + 172.4) * (EPS1D + 172.4) / (1. + EPS1D * EPS1D));
  }

  FREQLG10 = log10(FREQ);
  for (i = 1; i < 16; i++)
    if (FREQLG10 > FREQ1P[i])
      break;
  X = (FREQLG10 - FREQ1P[i]) / (FREQ1P[i - 1] - FREQ1P[i]) *
          (X1P[i - 1] - X1P[i]) +
      X1P[i];
  return pow10(X);
}

double HE12p3P(double FREQ)
{
  static double FREQ3P[16] = {
      15.943031, 15.909169, 15.872441, 15.832318, 15.788107,
      15.738880, 15.683351, 15.619667, 15.545012, 15.454805,
      15.340813, 15.270195, 15.230054, 15.185821, 15.136567,
      14.942557},
                X3P[16] = {-19.791021, -19.697886, -19.591421, -19.471855, -19.337053, -19.183958, -19.009750, -18.807990, -18.570571, -18.288361, -17.943476, -17.738737, -17.624154, -17.497163, -17.403183, -17.032999};
  double FREQLG10, X;
  int i;

  if (FREQ < 29223.753 * CLIGHTcm)
    return 0.;
  FREQLG10 = log10(FREQ);
  for (i = 1; i < 16; i++)
    if (FREQLG10 > FREQ3P[i])
      break;
  X = (FREQLG10 - FREQ3P[i]) / (FREQ3P[i - 1] - FREQ3P[i]) *
          (X3P[i - 1] - X3P[i]) +
      X3P[i];
  return pow10(X);
}

void HE1OP_new(double *ahe1, int iHe1, int iHe2, GlobalState *state)
{
  static double G[10] = {1., 3., 1., 9., 3., 3., 1., 9., 20., 3.},
                HEFREQ[10] = {5.945209e15, 1.152844e15, .9603331e15,
                              .8761076e15, .8147104e15, .4519048e15, .4030971e15,
                              .3821191e15, .3660215e15, .3627891E15},
                CHI[10] = {0., 19.819, 20.615, 20.964, 21.217,
                           22.718, 22.920, 23.006, 23.073, 23.086};
  double BOLT[10][MOSIZE], EXLIM[MOSIZE], TRANS[10], TRANS1S[10],
      TRANSN[27], BOLTN[27][MOSIZE], BOLTEX[MOSIZE],
      FREET[MOSIZE];
  double RYD, XR, XRLOG, FREQ3, FREQHE, ELIM, ZEFF2, CFREE, C, HE1, EX;
  int J, N, IMIN, NMIN;

  RYD = 109722.273 * CLIGHTcm;
  for (J = 0; J < state->NRHOX; J++)
  {
    for (N = 0; N < 10; N++)
      BOLT[N][J] = exp(-CHI[N] / state->TKEV[J]) * G[N] * state->FRACT[J][iHe1] / state->RHO[J];
    for (N = 3; N < 27; N++)
      BOLTN[N][J] = exp(-24.587 * (1. - 1. / (N * N)) / state->TKEV[J]) * 4. * N * N * state->FRACT[J][iHe1] / state->RHO[J];
    //  FREET[J]=state->XNE[J]*XNF(J,4)/state->RHO(J)/SQRT(state->T(J))
    FREET[J] = state->XNE[J] * 1.e-10 * state->FRACT[J][iHe2] * state->PARTITION_FUNCTIONS[J][iHe2] *
               1.e-10 / state->RHO[J] / sqrt(state->T[J]) * 1.e-10;
    //  XR=XNFP(J,3)*(4./2./13.595)*state->TKEV(J)/state->RHO(J)
    XRLOG = log(state->FRACT[J][iHe1] * (2. / 13.595) * state->TKEV[J] / state->RHO[J]);
    BOLTEX[J] = exp(-23.730 / state->TKEV[J] + XRLOG);
    EXLIM[J] = exp(-24.587 / state->TKEV[J] + XRLOG);
    //    ahe1[J]=0.1;
  }
  FREQ3 = state->FREQ * 1.e-10;
  FREQ3 = FREQ3 * FREQ3 * FREQ3;
  CFREE = 3.6919e8 / FREQ3;
  C = 2.815e-1 / FREQ3;

  for (NMIN = 0; NMIN < 10; NMIN++)
  {
    TRANS[NMIN] = 0;
    IMIN = NMIN + 1;
    if (HEFREQ[NMIN] <= state->FREQ)
      break;
    IMIN = 0;
  }
  switch (IMIN)
  {
  case 0:
  {
    for (J = 0; J < state->NRHOX; J++)
    {
      EX = (state->FREQ < 2.055e14) ? EXLIM[J] / state->EHVKT[J] : BOLTEX[J];
      HE1 = (EX - EXLIM[J]) * C;
      ahe1[J] = (HE1 + COULFF(J, 1, state) * FREET[J] * CFREE) * state->STIM[J];
    }
    return;
  }
  case 1:
    TRANS[0] = CROSSHE(state->FREQ);
  case 2:
    TRANS[1] = HE12s3S(state->FREQ);
  case 3:
    TRANS[2] = HE12s1S(state->FREQ);
  case 4:
    TRANS[3] = HE12p3P(state->FREQ);
  case 5:
    TRANS[4] = HE12p1P(state->FREQ);
  case 6:
    TRANS[5] = XKARZAS(state->FREQ, 1.236439e0, 3, 0); // 1s3s 3S
  case 7:
    TRANS[6] = XKARZAS(state->FREQ, 1.102898e0, 3, 0); // 1s3s 1S
  case 8:
    TRANS[7] = XKARZAS(state->FREQ, 1.045499e0, 3, 1); // 1s3p 3P
  case 9:
    TRANS[8] = XKARZAS(state->FREQ, 1.001427e0, 3, 2); // 1s3d 3D+1D
  case 10:
    TRANS[9] = XKARZAS(state->FREQ, 0.9926e0, 3, 1); // 1s3p 1P
  default:
    break;
  }
  // HeII n=2
  ELIM = 527490.06e0;
  FREQHE = (ELIM - 171135.00e0) * CLIGHTcm;
  if (state->FREQ >= FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[4] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
    FREQHE = (ELIM - 169087.e0) * CLIGHTcm;
  }
  if (state->FREQ >= FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[3] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
    FREQHE = (ELIM - 166277.546e0) * CLIGHTcm;
  }
  if (state->FREQ >= FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[2] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
    FREQHE = (ELIM - 159856.069e0) * CLIGHTcm;
  }
  if (state->FREQ < FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[1] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
  }

  // HeII n=3
  ELIM = 588451.59e0;
  FREQHE = (ELIM - 186209.471e0) * CLIGHTcm;
  if (state->FREQ >= FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[9] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
    FREQHE = (ELIM - 186101.e0) * CLIGHTcm;
  }
  if (state->FREQ >= FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[8] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
    FREQHE = (ELIM - 185564.e0) * CLIGHTcm;
  }
  if (state->FREQ >= FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[7] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
    FREQHE = (ELIM - 184864.e0) * CLIGHTcm;
  }
  if (state->FREQ >= FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[6] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
    FREQHE = (ELIM - 183236.e0) * CLIGHTcm;
  }
  if (state->FREQ >= FREQHE)
  {
    ZEFF2 = FREQHE / RYD;
    TRANS[5] += XKARZAS(state->FREQ, ZEFF2, 1, 0);
    if (state->FREQ >= 1.25408e16)
    {
      for (N = 4; N < 28; N++)
      {
        ZEFF2 = 4.e0 - 3.e0 / (N * N);
        TRANSN[N - 1] = XKARZAS(state->FREQ, ZEFF2, 1, 0);
      }
    }
  }
  //  printf("IMIN=%d, state->FREQ=%g\n",IMIN,state->FREQ);
  //  return;
  for (J = 0; J < state->NRHOX; J++)
  {
    EX = (state->FREQ < 2.055e14) ? EXLIM[J] / state->EHVKT[J] : BOLTEX[J];
    HE1 = (EX - EXLIM[J]) * C;
    for (N = IMIN - 1; N < 10; N++)
      HE1 += TRANS[N] * BOLT[N][J];
    if (state->FREQ >= 1.25408e16)
    {
      for (N = 3; N < 27; N++)
        HE1 += TRANSN[N] * BOLTN[N][J];
    }
    ahe1[J] = (HE1 + COULFF(J, 1, state) * FREET[J] * CFREE) * state->STIM[J];
  }
}

void HE2OP(double *ahe2, int iHe2, int iHe3, GlobalState *state) /*  REQUIRES FUNCTIONS COULX AND COULFF */
{
  /* FREQUENCIES ARE 4X HYDROGEN, CHI ARE FOR state->ION POT=54.403 */
  double HE2, C, CFREE, EX, FREQ3, BLTARG, BLTLOG, EXLLOG,
      XRLOG;
  double CONT[9], BOLT[MOSIZE][9], EXLIM[MOSIZE], FREET[MOSIZE], BOLTEX[MOSIZE];
  int J, N;

  for (J = 0; J < state->NRHOX; J++)
  {
    for (N = 0; N < 9; N++)
    {
      BLTARG = (54.403 - 54.403 / (N + 1) / (N + 1)) / state->TKEV[J] + log(state->RHO[J]);
      BOLT[J][N] = (state->FRACT[J][iHe2] == 0.0 || BLTARG > 80.) ? 0. : exp(-BLTARG) * 2. * (N + 1) * (N + 1) * state->FRACT[J][iHe2];
    }
    FREET[J] = state->XNE[J] * state->FRACT[J][iHe3] / sqrt(state->T[J]) / state->RHO[J];
    /*    XRLOG=log(state->TKEV[J]*(2/2/13.595)/state->RHO[J]); */
    XRLOG = log(state->TKEV[J] / 13.595 / state->RHO[J]);
    BLTLOG = 53.859 / state->TKEV[J] - XRLOG;
    BOLTEX[J] = (state->FRACT[J][iHe2] == 0.0 || BLTLOG > 80.) ? 0. : state->FRACT[J][iHe2] * exp(-BLTLOG);
    EXLLOG = 54.403 / state->TKEV[J] - XRLOG;
    EXLIM[J] = (state->FRACT[J][iHe2] == 0.0 || EXLLOG > 80.) ? 0. : state->FRACT[J][iHe2] * exp(-EXLLOG);
  }
  //  for(N=0; N<9; N++) CONT[N]=COULX(N, state->FREQ, 2.);
  for (N = 0; N < 9; N++)
    CONT[N] = XKARZAS(state->FREQ, 4.e0, N + 1, N + 1);
  FREQ3 = (state->FREQ * 1.e-05);
  FREQ3 = FREQ3 * FREQ3 * FREQ3;
  CFREE = 3.6919e-07 / FREQ3 * 4.;
  C = 2.815e14 * 2. * 2. / FREQ3;
  for (J = 0; J < state->NRHOX; J++)
  {
    EX = BOLTEX[J];
    if (state->FREQ < 1.31522e14)
      EX = EXLIM[J] / state->EHVKT[J];
    HE2 = (EX - EXLIM[J]) * C;
    for (N = 0; N < 9; N++)
      HE2 = HE2 + CONT[N] * BOLT[J][N];
    HE2 = (HE2 + COULFF(J, 2, state) * CFREE * FREET[J]) * state->STIM[J];
    ahe2[J] = (HE2 < 1.e-30) ? 0. : HE2;
  }
  return;
}

void HEMIOP(double *ahemin, int iHe1, GlobalState *state)
{
  double A, B, C;
  int J;

  A = 3.397e-26 + (-5.216e-11 + 7.039e05 / state->FREQ) / state->FREQ;
  B = -4.116e-22 + (1.067e-06 + 8.135e09 / state->FREQ) / state->FREQ;
  C = 5.081e-17 + (-8.724e-03 - 5.659e12 / state->FREQ) / state->FREQ;
  for (J = 0; J < state->NRHOX; J++)
    ahemin[J] = (A * state->T[J] + B + C / state->T[J]) * state->XNE[J] * state->FRACT[J][iHe1] / state->RHO[J] * 1.E-20;
  return;
}

void HERAOP(double *sighe, int iHe1, GlobalState *state)
{
  double WAVE, WW, SIG, S1;
  int J;

  WAVE = 2.997925e3 / min(state->FREQ * 1.e-15, 5.15); // wavelength in Angstroems
  WW = WAVE * WAVE;
  S1 = 1. + (2.44e5 + 5.94e10 / (WW - 2.90e5)) / WW;
  SIG = 5.484e-14 / WW / WW * S1 * S1;
  for (J = 0; J < state->NRHOX; J++)
    sighe[J] = SIG * state->FRACT[J][iHe1] / state->RHO[J];
  return;
}

double C1OP(int J, GlobalState *state) /* CROSS-SECTION */
{
  double C1240, C1444, X1240, X1444, X1100;

  C1240 = 5. * exp(-1.264 / state->TKEV[J]);
  C1444 = exp(-2.683 / state->TKEV[J]);
  X1444 = 0.;
  X1240 = 0.;
  X1100 = 0.;
  if (state->FREQ >= 2.7254e15)
    X1100 = SEATON(state->FREQ, 2.7254e15, 1.219e-17, 2.0, 3.317);
  if (state->FREQ >= 2.4196e15)
    X1240 = SEATON(state->FREQ, 2.4196e15, 1.030e-17, 1.5, 2.789);
  if (state->FREQ >= 2.0761e15)
    X1444 = SEATON(state->FREQ, 2.0761e15, 9.590e-18, 1.5, 3.501);
  return X1100 * 9. + X1240 * C1240 + X1444 * C1444;
}

double C1OP_new(int J, GlobalState *state) /* Cross-section                                */
{                                          /* This routine is based on R.L. Kurucz Atlas12 */
  static double ELEV[25] = {79314.86, 78731.27, 78529.62, 78309.76, 78226.35,
                            77679.82, 73975.91, 72610.72, 71374.90, 70743.95,
                            69722.00, 68856.33, 61981.82, 60373.00, 21648.01,
                            10192.63, 43.42, 16.42, 0.00, 119878.00,
                            105798.70, 97878.00, 75254.93, 64088.85, 33735.20},
                GLEV[25] = {9., 3., 7., 15., 21., 5., 1., 5., 9., 3., 15., 3., 3., 9., 1., 5., 5.,
                            3., 1., 3., 3., 5., 12., 15., 5.},
                RYD = 109732.298;
  double BOLT[25], X[25], Z, FREQ3, Z2FREQ, ZEFF2, ELIM, HCKT, WAVENO;
  double A, B, EPS, XS0, XS1, XD0, XD1, XD2, GFACTOR, H;
  int i, DEGEN;

  HCKT = state->HKT[J] * CLIGHTcm;
  for (i = 0; i < 25; i++)
  {
    BOLT[i] = GLEV[i] * exp(-ELEV[i] * HCKT);
    X[i] = 0.;
  }
  WAVENO = state->FREQ / CLIGHTcm;
  Z = 1.;
  FREQ3 = 2.815e29 / state->FREQ / state->FREQ / state->FREQ * Z * Z * Z * Z;
  Z2FREQ = 1.e20 * state->FREQ / (Z * Z);
  //  ELIM=90820.42
  //  C II 2P average
  ELIM = 90862.70;
  while (1)
  {
    //  2s2 2p3d 3P
    //  ELEV=79314.86
    if (WAVENO < ELIM - ELEV[0])
      break;
    //  GLEV=9.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[0]);
    X[0] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    //  2s2 2p3d 1P
    //  ELEV=78731.27
    if (WAVENO < ELIM - ELEV[1])
      break;

    //  GLEV=3.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[1]);
    X[1] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    //  2s2 2p3d 1F
    //  ELEV=78529.62
    if (WAVENO < ELIM - ELEV[2])
      break;

    //  GLEV=7.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[2]);
    X[2] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    //  2s2 2p3d 3D
    //  ELEV=78309.76
    if (WAVENO < ELIM - ELEV[3])
      break;

    //  GLEV=15.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[3]);
    X[3] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    //  2s2 2p3d 3F
    //  ELEV=78226.35
    if (WAVENO < ELIM - ELEV[4])
      break;

    //  GLEV=21.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[4]);
    X[4] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    //  2s2 2p3d 1D
    //  ELEV=77679.82
    if (WAVENO < ELIM - ELEV[5])
      break;

    //  GLEV=5.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[5]);
    X[5] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    //  2s2 2p3p 1S
    //  ELEV=73975.91
    if (WAVENO < ELIM - ELEV[6])
      break;

    //  GLEV=1.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[6]);
    X[6] = XKARZAS(state->FREQ, ZEFF2, 3, 1);
    //  2s2 2p3p 1D
    //  ELEV=72610.72
    if (WAVENO < ELIM - ELEV[7])
      break;

    //  GLEV=5.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[7]);
    X[7] = XKARZAS(state->FREQ, ZEFF2, 3, 1);
    //  2s2 2p3p 3P
    //  ELEV=71374.90
    if (WAVENO < ELIM - ELEV[8])
      break;

    //  GLEV=9.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[8]);
    X[8] = XKARZAS(state->FREQ, ZEFF2, 3, 1);
    //  2s2 2p3p 3S
    //  ELEV=70743.95
    if (WAVENO < ELIM - ELEV[9])
      break;

    //  GLEV=3.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[9]);
    X[9] = XKARZAS(state->FREQ, ZEFF2, 3, 1);
    //  2s2 2p3p 3D
    //  ELEV=69722.00
    if (WAVENO < ELIM - ELEV[10])
      break;

    //  GLEV=15.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[10]);
    X[10] = XKARZAS(state->FREQ, ZEFF2, 3, 1);
    //  2s2 2p3p 1P
    //  ELEV=68856.33
    if (WAVENO < ELIM - ELEV[11])
      break;

    //  GLEV=3.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[11]);
    X[11] = XKARZAS(state->FREQ, ZEFF2, 3, 1);
    //  2s2 2p3s 1P
    //  ELEV=61981.82
    if (WAVENO < ELIM - ELEV[12])
      break;

    //  GLEV=3.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[12]);
    X[12] = XKARZAS(state->FREQ, ZEFF2, 3, 0);
    //  2s2 2p3s 3P
    //  ELEV=60373.00
    if (WAVENO < ELIM - ELEV[13])
      break;

    //  GLEV=9.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[13]);
    X[13] = XKARZAS(state->FREQ, ZEFF2, 3, 0);
    break;
  }

  //  C II 2s2 2p 2P1/2
  ELIM = 90820.42;
  while (1)
  {
    //  2s2 2p2 1S
    //  ELEV=21648.01
    if (WAVENO < ELIM - ELEV[14])
      break;

    //  GLEV=1.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //    XS0=10.^(-16.80-(WAVENO-69172.400)/3.00/RYD)
    XS0 = pow10(-16.80 - (WAVENO - ELIM + ELEV[14]) / 3.00 / RYD);
    EPS = (WAVENO - 97700.) * 2. / 2743.;
    A = 68.e-18;
    B = 118.e-18;
    // Fit to Burke, P.G. and Taylor, K.state->T. 1979, J. Phys. B, 12, 2971-2984.
    XS1 = (A * EPS + B) / (EPS * EPS + 1.);
    X[14] = (XS0 + XS1) / 3.;
    //  2s2 2p2 1D
    //  ELEV=10192.63
    if (WAVENO < ELIM - ELEV[15])
      break;

    //  GLEV=5.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //     XD0=10.^(-16.80-(WAVENO-80627.760)/3.00/RYD)
    XD0 = pow10(-16.80 - (WAVENO - ELIM + ELEV[15]) / 3.00 / RYD);
    // Fit to Burke, P.G. and Taylor, K.state->T. 1979, J. Phys. B, 12, 2971-2984.
    EPS = (WAVENO - 93917.) * 2. / 9230.;
    A = 22.e-18;
    B = 26.e-18;
    XD1 = (A * EPS + B) / (EPS * EPS + 1.);
    // Fit to Burke, P.G. and Taylor, K.state->T. 1979, J. Phys. B, 12, 2971-2984.
    EPS = (WAVENO - 111130.) * 2. / 2743.;
    A = -10.5e-18;
    B = 46.e-18;
    XD2 = (A * EPS + B) / (EPS * EPS + 1.);
    X[15] = (XD0 + XD1 + XD2) * 1. / 3.;
    //  2s2 2p2 3P2
    //  ELEV=43.42
    if (WAVENO < ELIM - ELEV[16])
      break;

    //  GLEV=5.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //     X(16)=10.^(-16.80-(WAVENO-90777.000)/3.00/RYD)*1./3.
    X[16] = pow10(-16.80 - (WAVENO - ELIM + ELEV[16]) / 3.00 / RYD) / 3.;
    //  2s2 2p2 3P1
    //  ELEV=16.42
    if (WAVENO < ELIM - ELEV[17])
      break;

    //  GLEV=3.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //     X(17)=10.^(-16.80-(WAVENO-90777.000)/3.00/RYD)*1./3.
    X[17] = pow10(-16.80 - (WAVENO - ELIM + ELEV[17]) / 3.00 / RYD) / 3.;
    //  2s2 2p2 3P0
    //  ELEV=0.
    if (WAVENO < ELIM - ELEV[18])
      break;

    //  GLEV=1.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //     X(18)=10.^(-16.80-(WAVENO-90777.000)/3.00/RYD)*1./3.
    X[18] = pow10(-16.80 - (WAVENO - ELIM + ELEV[18]) / 3.00 / RYD) / 3.;
    break;
  }

  //  C II 2s2 2p 2P3/2
  ELIM = 90820.42 + 63.42;
  while (1)
  {
    //  2s2 2p2 1S
    //  ELEV=21648.01
    if (WAVENO < ELIM - ELEV[14])
      break;

    //  GLEV=1.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //     XS0=10.^(-16.80-(WAVENO-69172.400)/3.00/RYD)
    XS0 = pow10(-16.80 - (WAVENO - ELIM + ELEV[14]) / 3.00 / RYD);
    EPS = (WAVENO - 97700.) * 2. / 2743.;
    A = 68.e-18;
    B = 118.e-18;
    // Fit to Burke, P.G. and Taylor, K.state->T. 1979, J. Phys. B, 12, 2971-2984.
    XS1 = (A * EPS + B) / (EPS * EPS + 1.);
    X[14] += (XS0 + XS1) * 2. / 3.;
    //  2s2 2p2 1D
    //  ELEV=10192.63
    if (WAVENO < ELIM - ELEV[15])
      break;

    //  GLEV=5.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //     XD0=10.^(-16.80-(WAVENO-80627.760)/3.00/RYD)
    XD0 = pow10(-16.80 - (WAVENO - ELIM + ELEV[15]) / 3.00 / RYD);
    // Fit to Burke, P.G. and Taylor, K.state->T. 1979, J. Phys. B, 12, 2971-2984.
    EPS = (WAVENO - 93917.) * 2. / 9230.;
    A = 22.e-18;
    B = 26.e-18;
    XD1 = (A * EPS + B) / (EPS * EPS + 1.);
    // Fit to Burke, P.G. and Taylor, K.state->T. 1979, J. Phys. B, 12, 2971-2984.
    EPS = (WAVENO - 111130.) * 2. / 2743.;
    A = -10.5e-18;
    B = 46.e-18;
    XD2 = (A * EPS + B) / (EPS * EPS + 1.);
    X[15] += (XD0 + XD1 + XD2) * 2. / 3.;
    //  2s2 2p2 3P2
    //  ELEV=43.42
    if (WAVENO < ELIM - ELEV[16])

      //  GLEV=5.
      // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
      //     X(16)=10.^(-16.80-(WAVENO-90777.000)/3.00/RYD)*2./3.
      X[16] += pow10(-16.80 - (WAVENO - ELIM + ELEV[16]) / 3.00 / RYD) * 2. / 3.;
    //  2s2 2p2 3P1
    //  ELEV=16.42
    if (WAVENO < ELIM - ELEV[17])
      break;

    //  GLEV=3.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //     X(17)=10.^(-16.80-(WAVENO-90777.000)/3.00/RYD)*2./3.
    X[17] += pow10(-16.80 - (WAVENO - ELIM + ELEV[17]) / 3.00 / RYD) * 2. / 3.;
    //  2s2 2p2 3P0
    //  ELEV=0.
    if (WAVENO < ELIM - ELEV[18])
      break;

    //  GLEV=1.
    // Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
    //     X(18)=10.^(-16.80-(WAVENO-90777.000)/3.00/RYD)/3.
    //    X[18]+=pow10(-16.80-(WAVENO-ELIM+ELEV[18])/3.00/RYD)*2./3.;
    // Corrected to match the reference above
    X[18] += pow10(-16.80 - (WAVENO - ELIM + ELEV[18]) / 3.00 / RYD) * 2. / 3.;
    break;
  }

  //  C II 2s 2p2 4P1/2
  ELIM = 90820.42 + 43003.3;
  while (1)
  {
    //  2s2p3 1P
    //  ELEV=119878.
    if (WAVENO < ELIM - ELEV[19])
      break;

    //  GLEV=3.
    DEGEN = 3;
    ZEFF2 = 4. / RYD * (ELIM - ELEV[19]);
    X[19] = XKARZAS(state->FREQ, ZEFF2, 2, 1) * DEGEN;
    //  2s2p3 3S
    //  ELEV=105798.7
    if (WAVENO < ELIM - ELEV[20])
      break;

    //  GLEV=3.
    DEGEN = 3;
    ZEFF2 = 4. / RYD * (ELIM - ELEV[20]);
    X[20] = XKARZAS(state->FREQ, ZEFF2, 2, 1) * DEGEN;
    //  2s2p3 1D
    //  ELEV=97878.
    if (WAVENO < ELIM - ELEV[21])
      break;

    //  GLEV=5.
    DEGEN = 3;
    ZEFF2 = 4. / RYD * (ELIM - ELEV[21]);
    X[21] = XKARZAS(state->FREQ, ZEFF2, 2, 1) * DEGEN;
    //  2s2p3 3P
    //  ELEV=75254.93
    if (WAVENO < ELIM - ELEV[22])
      break;

    //  GLEV=12.
    DEGEN = 3;
    ZEFF2 = 4. / RYD * (ELIM - ELEV[22]);
    X[22] = XKARZAS(state->FREQ, ZEFF2, 2, 1) * DEGEN;
    //  2s2p3 3D
    //  ELEV=64088.85
    if (WAVENO < ELIM - ELEV[23])
      break;

    //  GLEV=15.
    DEGEN = 3;
    ZEFF2 = 4. / RYD * (ELIM - ELEV[23]);
    X[23] = XKARZAS(state->FREQ, ZEFF2, 2, 1) * DEGEN;
    //  2s2p3 5S
    //  ELEV=33735.20
    if (WAVENO < ELIM - ELEV[24])
      break;

    //  GLEV=5.
    DEGEN = 3;
    ZEFF2 = 4. / RYD * (ELIM - ELEV[24]);
    X[24] = XKARZAS(state->FREQ, ZEFF2, 2, 1) * DEGEN;
    break;
  }

  ELIM = 90820.42e0;
  GFACTOR = 6.;
  //  N=4 TO INFINITY
  H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
      (exp(-max(ELIM - RYD * Z * Z / 16., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
  //  printf("%d %g %g %g %g %g\n", J, H, ELIM, WAVENO, ELIM-WAVENO, HCKT);
  //  C II 2s 2p2 4P1/2
  //  ELIM=90820.42+43003.3
  for (i = 0; i < 25; i++)
    H += X[i] * BOLT[i];
  return H;
}

double MG1OP(int J, GlobalState *state) // CROSS-SECTION TIMES THE PARTITION FUNCTION
{
  static double PEACH[15][7] =
      {
          // TEMP: 4000     5000     6000     7000     8000     9000    10000     WAVE(A)
          {-42.474, -42.350, -42.109, -41.795, -41.467, -41.159, -40.883},  // 1500
          {-41.808, -41.735, -41.582, -41.363, -41.115, -40.866, -40.631},  // 1550
          {-41.273, -41.223, -41.114, -40.951, -40.755, -40.549, -40.347},  // 1621
          {-45.583, -44.008, -42.957, -42.205, -41.639, -41.198, -40.841},  // 1622
          {-44.324, -42.747, -41.694, -40.939, -40.370, -39.925, -39.566},  // 2513
          {-50.969, -48.388, -46.630, -45.344, -44.355, -43.568, -42.924},  // 2514
          {-50.633, -48.026, -46.220, -44.859, -43.803, -42.957, -42.264},  // 3756
          {-53.028, -49.643, -47.367, -45.729, -44.491, -43.520, -42.736},  // 3757
          {-51.785, -48.352, -46.050, -44.393, -43.140, -42.157, -41.363},  // 6549
          {-52.285, -48.797, -46.453, -44.765, -43.486, -42.480, -41.668},  // 6550
          {-52.028, -48.540, -46.196, -44.507, -43.227, -42.222, -41.408},  // 7234
          {-52.384, -48.876, -46.513, -44.806, -43.509, -42.488, -41.660},  // 7235
          {-52.363, -48.856, -46.493, -44.786, -43.489, -42.467, -41.639},  // 7291
          {-54.704, -50.772, -48.107, -46.176, -44.707, -43.549, -42.611},  // 7292
          {-54.359, -50.349, -47.643, -45.685, -44.198, -43.027, -42.418}}; // 9000
  static double FREQMG[7] = {1.9341452e15, 1.8488510e15, 1.1925797e15,
                             7.9804046e14, 4.5772110e14, 4.1440977e14,
                             4.1113514e14};
  static double FLOG[9] = {35.32123, 35.19844, 35.15334, 34.71490, 34.31318,
                           33.75728, 33.65788, 33.64994, 33.43947};
  static double TLG[7] = {8.29405, 8.51719, 8.69951, 8.85367,
                          8.98720, 9.10498, 9.21034};
  double XWL1, XWL2, D, D1, DT;
  int N, NT;

  NT = min(6, (int)floor(state->T[J] / 1000.) - 3);
  if (NT < 1)
    NT = 1;
  DT = (state->TLOG[J] - TLG[NT - 1]) / (TLG[NT] - TLG[NT - 1]);
  for (N = 0; N < 7; N++)
    if (state->FREQ > FREQMG[N])
      break;
  D = (state->FREQLG - FLOG[N]) / (FLOG[N + 1] - FLOG[N]);
  if (N > 1)
    N = 2 * N - 1;
  D1 = 1.0 - D;
  XWL1 = PEACH[N + 1][NT - 1] * D + PEACH[N][NT - 1] * D1;
  XWL2 = PEACH[N + 1][NT] * D + PEACH[N][NT] * D1;
  return exp(XWL1 * (1.0 - DT) + XWL2 * DT);
}

double MG1OP_new(int J, GlobalState *state) /* Cross-section                                */
{                                           /* This routine is based on R.L. Kurucz Atlas12 */
  static double ELEV[15] = {54676.710, 54676.438, 54192.284, 53134.642, 49346.729,
                            47957.034, 47847.797, 46403.065, 43503.333, 41197.043,
                            35051.264, 21919.178, 21870.464, 21850.405, 0.};
  static double GLEV[15] = {21., 7., 15., 5., 3., 15., 9., 5., 1., 3., 3., 5., 3., 1., 1.};
  static double RYD = 109732.298e0, ELIM = 61671.02e0, Z = 1., GFACTOR = 2.;
  double BOLT[15], X[15], FREQ3, WAVENO, H, HCKT, ZEFF2;
  int i;

  HCKT = state->HKT[J] * CLIGHTcm;
  for (i = 0; i < 15; i++)
  {
    BOLT[i] = GLEV[i] * exp(-ELEV[i] * HCKT);
    X[i] = 0.;
  }
  FREQ3 = 2.815e29 / state->FREQ / state->FREQ / state->FREQ * Z * Z * Z * Z;
  WAVENO = state->FREQ / CLIGHTcm;

  //  3s4f 3F
  //     ELEV=54676.710
  if (WAVENO < ELIM - ELEV[0])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    // Commented out because all X are zero.
    //    for(i=0; i<15; i++) H+=X[i]*BOLT[i];
    return H;
  }
  //     GLEV=21.
  ZEFF2 = 16. / RYD * (ELIM - ELEV[0]);
  X[0] = XKARZAS(state->FREQ, ZEFF2, 4, 3);
  //  3s4f 1F
  //     ELEV=54676.438
  if (WAVENO < ELIM - ELEV[1])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 1; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=7.
  ZEFF2 = 16. / RYD * (ELIM - ELEV[1]);
  X[1] = XKARZAS(state->FREQ, ZEFF2, 4, 3);
  //  3s4d 3D
  //     ELEV=54192.284
  if (WAVENO < ELIM - ELEV[2])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 2; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=15.
  ZEFF2 = 16. / RYD * (ELIM - ELEV[2]);
  X[2] = XKARZAS(state->FREQ, ZEFF2, 4, 2);
  //  3s4d 1D
  //     ELEV=53134.642
  if (WAVENO < ELIM - ELEV[3])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 3; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=5.
  ZEFF2 = 16. / RYD * (ELIM - ELEV[3]);
  X[3] = XKARZAS(state->FREQ, ZEFF2, 4, 2);
  //  3s4p 1P
  //     ELEV=49346.729
  if (WAVENO < ELIM - ELEV[4])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 4; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=3.
  ZEFF2 = 16. / RYD * (ELIM - ELEV[4]);
  X[4] = XKARZAS(state->FREQ, ZEFF2, 4, 1);
  //  3s3d 3D
  //     ELEV=47957.034
  if (WAVENO < ELIM - ELEV[5])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 5; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=15.
  X[5] = 25.e-18 * pow(13713.986e0 / WAVENO, 2.7);
  //  3s4p 3P
  //     ELEV=47847.797
  if (WAVENO < ELIM - ELEV[6])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 6; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=9.
  X[6] = 33.8e-18 * pow((13823.223e0 / WAVENO), 2.8);
  //  3s3d 1D
  //     ELEV=46403.065
  if (WAVENO < ELIM - ELEV[7])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 7; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=5.
  X[7] = 45.e-18 * pow((15267.955e0 / WAVENO), 2.7);
  //  3s4s 1S
  //     ELEV=43503.333
  if (WAVENO < ELIM - ELEV[8])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 8; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=1.
  X[8] = 0.43e-18 * pow((18167.687e0 / WAVENO), 2.6);
  //  3s4s 3S
  //     ELEV=41197.043
  if (WAVENO < ELIM - ELEV[9])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 9; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=3.
  X[9] = 2.1e-18 * pow((20473.617e0 / WAVENO), 2.6);
  //  2s3p 1P
  //     ELEV=35051.264
  if (WAVENO < ELIM - ELEV[10])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 10; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=3.
  X[10] = 16.e-18 * pow((26619.756e0 / WAVENO), 2.1) -
          7.8e-18 * pow((26619.756e0 / WAVENO), 9.5);
  //  3s3p 3P
  //     ELEV=21911.178
  if (WAVENO < ELIM - ELEV[11])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 11; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=5.
  ZEFF2 = 9. / RYD * (ELIM - ELEV[11]);
  X[11] = 20.e-18 * pow(39759.842e0 / WAVENO, 2.7);
  X[11] = max(X[11], 40.e-18 * pow(39759.842e0 / WAVENO, 14.));
  //  3s3p 3P
  //     ELEV=21870.464
  if (WAVENO < ELIM - ELEV[12])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 12; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=3.
  ZEFF2 = 9. / RYD * (ELIM - ELEV[12]);
  X[12] = 20.e-18 * pow((39759.842 / WAVENO), 2.7);
  X[12] = max(X[12], 40.e-18 * pow((39759.842e0 / WAVENO), 14.));
  //  3s3p 3P0
  //     ELEV=21850.405
  if (WAVENO < ELIM - ELEV[13])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 13; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=1.
  ZEFF2 = 9. / RYD * (ELIM - ELEV[13]);
  X[13] = 20.e-18 * pow((39759.842e0 / WAVENO), 2.7);
  X[13] = max(X[13], 40.e-18 * pow((39759.842e0 / WAVENO), 14.));
  //  3s2 1S
  //     ELEV=0.
  if (WAVENO < ELIM - ELEV[14])
  {
    H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
        (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
    for (i = 0; i < 14; i++)
      H += X[i] * BOLT[i];
    return H;
  }
  //     GLEV=1.
  X[14] = 1.1e-18 * pow((ELIM - ELEV[14]) / WAVENO, 10.);
  H = FREQ3 * GFACTOR * 2. / 2. / (RYD * Z * Z * HCKT) *
      (exp(-max(ELIM - RYD * Z * Z / 25., ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
  for (i = 0; i < 15; i++)
    H += X[i] * BOLT[i];
  return H;
}

double AL1OP(int J, GlobalState *state)
{
  return (state->FREQ >= 1.443e15) ? 2.1e-17 * pow(1.443e15 / state->FREQ, 3.) * 6 : 0.;
}

double AL1OP_new(int J, GlobalState *state) /* Cross-section                                */
{                                           /* This routine is based on R.L. Kurucz Atlas12 */
  double ELIM, WAVENO, F1, F2, al1op;

  WAVENO = state->FREQ / CLIGHTcm;
  ELIM = 48278.37e0;

  if (WAVENO < (ELIM - 112.061e0))
  {
    al1op = 0.;
  }
  else if (WAVENO >= (ELIM - 112.061e0) && WAVENO < ELIM)
  {
    //     3s2 3p 2P3/2
    //    al1op=6.5e-17*((ELIM-112.061e0)/WAVENO)^5*4.
    F1 = (ELIM - 112.061e0) / WAVENO;
    F1 = F1 * F1 * F1 * F1 * F1 * 4.;
    al1op = 6.5e-17 * F1;
  }
  else
  {
    //     3s2 3p 2P1/2
    //    al1op=6.5e-17*((ELIM-112.061e0)/WAVENO)^5*4.+
    //          6.5E-17*(ELIM/WAVENO)^5*2.;
    F1 = (ELIM - 112.061e0) / WAVENO;
    F1 = F1 * F1 * F1 * F1 * F1 * 4.;
    F2 = ELIM / WAVENO;
    F2 = F2 * F2 * F2 * F2 * F2 * 2.;
    al1op = 6.5e-17 * (F1 + F2);
  }
  return al1op;
}

double SI1OP(int J, GlobalState *state) /* Cross-section                                */
{
  static double PEACH[19][9] =
      /* TEMP:4000   5000   6000   7000   8000   9000  10000  11000  12000   WAVE(A)*/
      {{38.136, 38.138, 38.140, 38.141, 38.143, 38.144, 38.144, 38.145, 38.145},  /* 1200 */
       {37.834, 37.839, 37.843, 37.847, 37.850, 37.853, 37.855, 37.857, 37.858},  /* 1400 */
       {37.898, 37.898, 37.897, 37.897, 37.897, 37.896, 37.895, 37.895, 37.894},  /* 1519 */
       {40.737, 40.319, 40.047, 39.855, 39.714, 39.604, 39.517, 39.445, 39.385},  /* 1520 */
       {40.581, 40.164, 39.893, 39.702, 39.561, 39.452, 39.366, 39.295, 39.235},  /* 1676 */
       {45.521, 44.456, 43.753, 43.254, 42.878, 42.580, 42.332, 42.119, 41.930},  /* 1677 */
       {45.520, 44.455, 43.752, 43.251, 42.871, 42.569, 42.315, 42.094, 41.896},  /* 1978 */
       {55.068, 51.783, 49.553, 47.942, 46.723, 45.768, 44.997, 44.360, 43.823},  /* 1979 */
       {53.868, 50.369, 48.031, 46.355, 45.092, 44.104, 43.308, 42.652, 42.100},  /* 5379 */
       {54.133, 50.597, 48.233, 46.539, 45.261, 44.262, 43.456, 42.790, 42.230},  /* 5380 */
       {54.051, 50.514, 48.150, 46.454, 45.176, 44.175, 43.368, 42.702, 42.141},  /* 5624 */
       {54.442, 50.854, 48.455, 46.733, 45.433, 44.415, 43.592, 42.912, 42.340},  /* 5625 */
       {54.320, 50.722, 48.313, 46.583, 45.277, 44.251, 43.423, 42.738, 42.160},  /* 6260 */
       {55.691, 51.965, 49.444, 47.615, 46.221, 45.119, 44.223, 43.478, 42.848},  /* 6261 */
       {55.661, 51.933, 49.412, 47.582, 46.188, 45.085, 44.189, 43.445, 42.813},  /* 6349 */
       {55.973, 52.193, 49.630, 47.769, 46.349, 45.226, 44.314, 43.555, 42.913},  /* 6350 */
       {55.922, 52.141, 49.577, 47.715, 46.295, 45.172, 44.259, 43.500, 42.858},  /* 6491 */
       {56.828, 52.821, 50.110, 48.146, 46.654, 45.477, 44.522, 43.730, 43.061},  /* 6492 */
       {56.657, 52.653, 49.944, 47.983, 46.491, 45.315, 44.360, 43.569, 42.901}}; /*6900 */
                                                                                  /*     3P,1D,1S,1D,3D,3F,1D,3P */
  static double FREQSI[9] = {2.1413750e15, 1.97231650e15, 1.7879689e15,
                             1.5152920e15, 0.55723927e15, 5.3295914e14,
                             4.7886458e14, 4.72164220e14, 4.6185133e14};
  static double FLOG[11] = {35.45438, 35.30022, 35.21799, 35.11986, 34.95438,
                            33.95402, 33.90947, 33.80244, 33.78835, 33.76626,
                            33.70518};
  static double TLG[9] = {8.29405, 8.51719, 8.69951, 8.85367, 8.98720,
                          9.10498, 9.21034, 9.30565, 9.39266};
  double D, DT, DD, XWL1, XWL2;
  int NT, N;

  NT = min(8, (int)floor(state->T[J] / 1000.) - 3);
  if (NT < 1)
    NT = 1;
  DT = (state->TLOG[J] - TLG[NT - 1]) / (TLG[NT] - TLG[NT - 1]);
  for (N = 0; N < 9; N++)
    if (state->FREQ > FREQSI[N])
      break;
  D = (state->FREQLG - FLOG[N]) / (FLOG[N + 1] - FLOG[N]);
  if (N > 1)
    N = 2 * N - 1;
  DD = 1. - D;
  XWL1 = PEACH[N + 1][NT - 1] * D + PEACH[N][NT - 1] * DD;
  XWL2 = PEACH[N + 1][NT] * D + PEACH[N][NT] * DD;
  return exp(-(XWL1 * (1. - DT) + XWL2 * DT)) * 9.;
}

double SI1OP_new(int J, GlobalState *state) /* Cross-section                                */
{                                           /* This routine is based on R.L. Kurucz Atlas12 */
  static double ELEV[33] = {
      59962.284, 59100., 59077.112, 58893.40, 58801.529,
      58777., 57488.974, 56503.346, 54225.621, 53387.34,
      53362.24, 51612.012, 50533.424, 50189.389, 49965.894,
      49399.670, 49128.131, 48161.459, 47351.554, 47284.061,
      40991.884, 39859.920, 15394.370, 6298.850, 223.157,
      77.115, 0.000, 94000., 79664.0, 72000.,
      56698.738, 45303.310, 33326.053};
  static double GLEV[33] = {
      9., 56., 15., 7., 3., 28., 21., 5., 15., 3., 7., 1., 9., 5., 21.,
      3., 9., 15., 5., 3., 3., 9., 1., 5., 5., 3., 1., 3., 3., 5., 12., 15., 5.};
  double BOLT[33], X[33], HCKT, FREQ3, WAVENO, ELIM, RYD, ZEFF2, EPS, RESON1,
      DEGEN, GFACTOR, aSi1op;
  int I;

  HCKT = state->HKT[J] * CLIGHTcm;
  FREQ3 = 2.815E29 / state->FREQ / state->FREQ / state->FREQ;
  WAVENO = state->FREQ / CLIGHTcm;
  RYD = 109732.298e0;

  for (I = 0; I < 33; I++)
  {
    BOLT[I] = GLEV[I] * exp(-ELEV[I] * HCKT);
    X[I] = 0.;
  }

  while (1)
  {
    //Si II 3s2 3p 2P average
    ELIM = 65939.18e0;

    // 3s2 3p4d 3P
    // ELEV=59962.284
    if (WAVENO < ELIM - ELEV[0])
      break;

    // GLEV=9.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[0]);
    X[0] = XKARZAS(state->FREQ, ZEFF2, 4, 2);
    // 3s2 3p4f (2P3/2)4f
    // ELEV=59100.
    if (WAVENO < ELIM - ELEV[1])
      break;

    // GLEV=56.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[1]);
    X[1] = XKARZAS(state->FREQ, ZEFF2, 4, 3);
    // 3s2 3p4d 3D
    // ELEV=59077.112
    if (WAVENO < ELIM - ELEV[2])
      break;

    // GLEV=15.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[2]);
    X[2] = XKARZAS(state->FREQ, ZEFF2, 4, 2);
    // 3s2 3p4d 1F
    // ELEV=58893.40
    if (WAVENO < ELIM - ELEV[3])
      break;

    // GLEV=7.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[3]);
    X[3] = XKARZAS(state->FREQ, ZEFF2, 4, 2);
    // 3s2 3p4d 1P
    // ELEV=58801.529
    if (WAVENO < ELIM - ELEV[4])
      break;

    // GLEV=3.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[4]);
    X[4] = XKARZAS(state->FREQ, ZEFF2, 4, 2);
    // 3s2 3p4f (2P1/2)4f
    // ELEV=58777.
    if (WAVENO < ELIM - ELEV[5])
      break;

    // GLEV=28.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[5]);
    X[5] = XKARZAS(state->FREQ, ZEFF2, 4, 3);
    // 3s2 3p4d 3F
    // ELEV=57488.974
    if (WAVENO < ELIM - ELEV[6])
      break;

    // GLEV=21.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[6]);
    X[6] = XKARZAS(state->FREQ, ZEFF2, 4, 2);
    // 3s2 3p4d 1D
    // ELEV=56503.346
    if (WAVENO < ELIM - ELEV[7])
      break;

    // GLEV=5.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[7]);
    X[7] = XKARZAS(state->FREQ, ZEFF2, 4, 2);
    // 3s2 3p3d 3D
    // ELEV=54225.621
    if (WAVENO < ELIM - ELEV[8])
      break;

    // GLEV=15.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[8]);
    X[8] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    // 3s2 3p3d 1P
    // ELEV=53387.34
    if (WAVENO < ELIM - ELEV[9])
      break;

    // GLEV=3.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[9]);
    X[9] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    // 3s2 3p3d 1F
    // ELEV=53362.24
    if (WAVENO < ELIM - ELEV[10])
      break;

    // GLEV=7.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[10]);
    X[10] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    // 3s2 3p4p 1S
    // ELEV=51612.012
    if (WAVENO < ELIM - ELEV[11])
      break;

    // GLEV=1.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[11]);
    X[11] = XKARZAS(state->FREQ, ZEFF2, 4, 1);
    // 3s2 3p3d 3P
    // ELEV=50533.424
    if (WAVENO < ELIM - ELEV[12])
      break;

    // GLEV=9.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[12]);
    X[12] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    // 3s2 3p4p 1D
    // ELEV=50189.389
    if (WAVENO < ELIM - ELEV[13])
      break;

    // GLEV=5.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[13]);
    X[13] = XKARZAS(state->FREQ, ZEFF2, 4, 1);
    // 3s2 3p3d 3F
    // ELEV=49965.894
    if (WAVENO < ELIM - ELEV[14])
      break;

    // GLEV=21.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[14]);
    X[14] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    // 3s2 3p4p 3S
    // ELEV=49399.670
    if (WAVENO < ELIM - ELEV[15])
      break;

    // GLEV=3.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[15]);
    X[15] = XKARZAS(state->FREQ, ZEFF2, 4, 1);
    // 3s2 3p4p 3P
    // ELEV=49128.131
    if (WAVENO < ELIM - ELEV[16])
      break;

    // GLEV=9.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[16]);
    X[16] = XKARZAS(state->FREQ, ZEFF2, 4, 1);
    // 3s2 3p4p 3D
    // ELEV=48161.459
    if (WAVENO < ELIM - ELEV[17])
      break;

    // GLEV=15.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[17]);
    X[17] = XKARZAS(state->FREQ, ZEFF2, 4, 1);
    // 3s2 3p3d 1D
    // ELEV=47351.554
    if (WAVENO < ELIM - ELEV[18])
      break;

    // GLEV=5.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[18]);
    X[18] = XKARZAS(state->FREQ, ZEFF2, 3, 2);
    // 2s2 3p4p 1P
    // ELEV=47284.061
    if (WAVENO < ELIM - ELEV[19])
      break;

    // GLEV=3.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[19]);
    X[19] = XKARZAS(state->FREQ, ZEFF2, 4, 1);
    // 3s2 3p4s 1P
    // ELEV=40991.884
    if (WAVENO < ELIM - ELEV[20])
      break;

    // GLEV=3.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[20]);
    X[20] = XKARZAS(state->FREQ, ZEFF2, 4, 0);
    // 3s2 3p4s 3P
    // ELEV=39859.920
    if (WAVENO < ELIM - ELEV[21])
      break;

    // GLEV=9.
    ZEFF2 = 16. / RYD * (ELIM - ELEV[21]);
    X[21] = XKARZAS(state->FREQ, ZEFF2, 4, 0);
    break;
  }

  //Si II 3s2 3p 2P1/2
  ELIM = 65747.55e0;

  while (1)
  {
    // 3s2 3p2 1S
    // ELEV=15394.370
    if (WAVENO < ELIM - ELEV[22])
      break;

    // GLEV=1.
    EPS = (WAVENO - 70000.e0) * 2.e0 / 6500.e0;
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    RESON1 = (97.e-18 * EPS + 94.e-18) / (EPS * EPS + 1.);
    X[22] = (37.e-18 * pow(50353.180e0 / WAVENO, 2.40) + RESON1) / 3.;
    // 3s2 3p2 1D
    // ELEV=6298.850
    if (WAVENO < ELIM - ELEV[23])
      break;

    // GLEV=5.
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    EPS = (WAVENO - 78600.) * 2. / 13000.;
    RESON1 = (-10.e-18 * EPS + 77.e-18) / (EPS * EPS + 1.);
    X[23] = (24.5e-18 * pow(59448.70e0 / WAVENO, 1.85) + RESON1) / 3.;
    // 3s2 3p2 3P2
    // ELEV=223.157
    if (WAVENO < ELIM - ELEV[24])
      break;

    // GLEV=5.
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    if (WAVENO <= 74000.e0)
      X[24] = 72.e-18 * pow(65524.393e0 / WAVENO, 1.90) / 3.;
    else
      X[24] = 93.e-18 * pow(65524.393e0 / WAVENO, 4.00) / 3.;
    // 3s2 3p2 3P1
    // ELEV=77.115
    if (WAVENO < ELIM - ELEV[25])
      break;

    // GLEV=3.
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    if (WAVENO <= 74000.e0)
      X[25] = 72.e-18 * pow(65524.393e0 / WAVENO, 1.90) * 2. / 3.;
    else
      X[25] = 93.e-18 * pow(65524.393e0 / WAVENO, 4.00) * 2. / 3.;
    // 3s2 3p2 3P0
    // ELEV=0.00
    if (WAVENO < ELIM - ELEV[26])
      break;

    // GLEV=1.
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    if (WAVENO <= 74000.e0)
      X[26] = 72.e-18 * pow(65524.393e0 / WAVENO, 1.90) / 3.;
    else
      X[26] = 93.e-18 * pow(65524.393e0 / WAVENO, 4.00) / 3.;
    break;
  }

  // Si II 3s2 3p 2P3/2
  ELIM = 65747.55e0 + 287.45e0;

  while (1)
  {
    // 3s2 3p2 1S
    // ELEV=15394.370
    if (WAVENO < ELIM - ELEV[22])
      break;
    // GLEV=1.
    EPS = (WAVENO - 70000.e0) * 2. / 6500.e0;
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    RESON1 = (97.e-18 * EPS + 94.e-18) / (EPS * EPS + 1.);
    X[22] += (37.e-18 * pow(50353.180e0 / WAVENO, 2.40) + RESON1) * 2. / 3.;
    // 3s2 3p2 1D
    // ELEV=6298.850
    if (WAVENO < ELIM - ELEV[23])
      break;

    // GLEV=5.
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    EPS = (WAVENO - 78600.e0) * 2. / 13000.e0;
    RESON1 = (-10.e-18 * EPS + 77.e-18) / (EPS * EPS + 1.);
    X[23] += (24.5e-18 * pow(59448.700e0 / WAVENO, 1.85) + RESON1) * 2. / 3.;
    // 3s2 3p2 3P2
    // ELEV=223.157
    if (WAVENO < ELIM - ELEV[24])
      break;

    // GLEV=5.
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    if (WAVENO <= 74000.e0)
      X[24] += 72.e-18 * pow(65524.393e0 / WAVENO, 1.90) * 2. / 3.;
    else
      X[24] += 93.e-18 * pow(65524.393e0 / WAVENO, 4.00) * 2. / 3.;
    // 3s2 3p2 3P1
    // ELEV=77.115
    if (WAVENO < ELIM - ELEV[25])
      break;

    // GLEV=3.
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    if (WAVENO <= 74000.e0)
      X[25] += 72.e-18 * pow(65524.393e0 / WAVENO, 1.90) * 2. / 3.;
    else
      X[25] += 93.e-18 * pow(65524.393e0 / WAVENO, 4.00) * 2. / 3.;
    // 3s2 3p2 3P0
    // ELEV=0.00
    if (WAVENO < ELIM - ELEV[26])
      break;

    // GLEV=1.
    // fits to Nahar, S.N. and Pradhan, A.K. J.Phys.B 26, 1109-1127, 1993.
    if (WAVENO <= 74000.e0)
      X[26] += 72.e-18 * pow(65524.393e0 / WAVENO, 1.90) * 2. / 3.;
    else
      X[26] += 93.e-18 * pow(65524.393e0 / WAVENO, 4.00) * 2. / 3.;
    break;
  }

  //Si II 3s 3p2 4P1/2
  ELIM = 65747.5e0 + 42824.35e0;

  while (1)
  {
    // 3s3p3 1P
    // ELEV=94000.
    if (WAVENO < ELIM - ELEV[27])
      break;

    // GLEV=3.
    DEGEN = 3.;
    ZEFF2 = 9. / RYD * (ELIM - ELEV[27]);
    X[27] = XKARZAS(state->FREQ, ZEFF2, 3, 1) * DEGEN;
    // 3s3p3 3S
    // guess
    // ELEV=79664.0
    if (WAVENO < ELIM - ELEV[28])
      break;

    // GLEV=3.
    DEGEN = 3.;
    ZEFF2 = 9. / RYD * (ELIM - ELEV[28]);
    X[28] = XKARZAS(state->FREQ, ZEFF2, 3, 1) * DEGEN;
    // 3s3p3 1D
    // guess
    // ELEV=72000.
    if (WAVENO < ELIM - ELEV[29])
      break;

    // GLEV=5.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[29]);
    X[29] = XKARZAS(state->FREQ, ZEFF2, 3, 1) * DEGEN;
    // 3s3p3 3P
    // ELEV=56698.738
    if (WAVENO < ELIM - ELEV[30])
      break;

    // GLEV=12.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[30]);
    X[30] = XKARZAS(state->FREQ, ZEFF2, 3, 1) * DEGEN;
    // 2s2p3 3D
    // ELEV=45303.310
    if (WAVENO < ELIM - ELEV[31])
      break;

    // GLEV=15.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[31]);
    X[31] = XKARZAS(state->FREQ, ZEFF2, 3, 1) * DEGEN;
    // 2s2p3 5S
    // ELEV=33326.053
    if (WAVENO < ELIM - ELEV[32])
      break;

    // GLEV=5.
    ZEFF2 = 9. / RYD * (ELIM - ELEV[32]);
    X[32] = XKARZAS(state->FREQ, ZEFF2, 3, 1) * DEGEN;
    break;
  }

  ELIM = 65747.55e0;
  GFACTOR = 6.;

  // N=5 TO INFINITY
  aSi1op = FREQ3 * GFACTOR * 2. / 2. / (RYD * HCKT) *
           (exp(-max(ELIM - RYD, ELIM - WAVENO) * HCKT) - exp(-ELIM * HCKT));
  for (I = 0; I < 33; I++)
    aSi1op += X[I] * BOLT[I];
  return aSi1op;
}

double FE1OP(int J, GlobalState *state)
{
  /* 
  Cross-section time partition functions 
  This routine is based on R.L. Kurucz Atlas12 
  */
  static double G[48] = {25., 35., 21., 15., 9., 35., 33., 21., 27., 49., 9., 21.,
                         27., 9., 9., 25., 33., 15., 35., 3., 5., 11., 15., 13.,
                         15., 9., 21., 15., 21., 25., 35., 9., 5., 45., 27., 21.,
                         15., 21., 15., 25., 21., 35., 5., 15., 45., 35., 55., 25.};
  static double E[48] = {500., 7500., 12500., 17500., 19000., 19500., 19500.,
                         21000., 22000., 23000., 23000., 24000., 24000., 24500.,
                         24500., 26000., 26500., 26500., 27000., 27500., 28500.,
                         29000., 29500., 29500., 29500., 30000., 31500., 31500.,
                         33500., 33500., 34000., 34500., 34500., 35000., 35500.,
                         37000., 37000., 37000., 38500., 40000., 40000., 41000.,
                         41000., 43000., 43000., 43000., 43000., 44000.};
  static double WNO[48] = {63500., 58500., 53500., 59500., 45000., 44500., 44500.,
                           43000., 58000., 41000., 54000., 40000., 40000., 57500.,
                           55500., 38000., 57500., 57500., 37000., 54500., 53500.,
                           55000., 34500., 34500., 34500., 34000., 32500., 32500.,
                           32500., 32500., 32000., 29500., 29500., 31000., 30500.,
                           29000., 27000., 54000., 27500., 24000., 47000., 23000.,
                           44000., 42000., 42000., 21000., 42000., 42000.};
  double BOLT, XSECT, WAVENO, FE1OPACITY, XXX;
  int I;

  WAVENO = state->FREQ / CLIGHTcm;
  if (WAVENO < 21000.)
    return 0.;
  FE1OPACITY = 0.;
  for (I = 0; I < 48; I++)
  {
    BOLT = G[I] * exp(-E[I] * CLIGHTcm * state->HKT[J]);
    if (WNO[I] < WAVENO)
    {
      XXX = ((WNO[I] + 3000. - WAVENO) / WNO[I] / .1);
      XSECT = 3.e-18 / (1. + XXX * XXX * XXX * XXX);
    }
    else
      XSECT = 0.;
    FE1OPACITY += XSECT * BOLT;
  }
  return FE1OPACITY;
}

double FE1OP_new(int J, GlobalState *state)
{
  /*
  Cross-sections of Fe 1 photoionization time
  This routine is based on data provided by Bautista
  described in Bautista et al. 2017, A&A 606, 127
  */
  static double WN0 = 10000.000, WNSTEP = 20.000;
  static int n_WN = 12001, n_Ebin = 78, first = 1;
  static double Ebin[78], GCROSS[2401][78];
  double WAVENO, BOLT, FACTOR, kT_eV, fe1op;
  int i_wn, i_en, i;

  if (first)
  {
    char path[512];
    int headlen;
    char head[2048];
    float delta;
    FILE *fe1op_data;

    strncpy(path, state->PATH, state->PATHLEN + 1);
    strcat(path, DATAFILE_FE);
    fe1op_data = fopen(path, "rb");

    i = fread(&headlen, sizeof(int), 1, fe1op_data);
    if (state->change_byte_order)
      headlen = *(int *)ByteSwap((char *)&headlen, 4);
    i = fread(head, 1, headlen, fe1op_data);

    i = fread(&delta, sizeof(float), 1, fe1op_data);
    if (state->change_byte_order)
      delta = *(float *)ByteSwap((char *)&delta, 4);

    i = fread(&n_Ebin, sizeof(int), 1, fe1op_data);
    if (state->change_byte_order)
      n_Ebin = *(int *)ByteSwap((char *)&n_Ebin, 4);
    i = fread(Ebin, sizeof(double), n_Ebin, fe1op_data);
    if (state->change_byte_order)
    {
      for (i_en = 0; i_en < n_Ebin; i_en++)
        Ebin[i_en] = *(double *)ByteSwap((char *)(Ebin + i_en), 8);
    }

    i = fread(&n_WN, sizeof(int), 1, fe1op_data);
    if (state->change_byte_order)
      n_WN = *(int *)ByteSwap((char *)&n_WN, 4);

    i = fread(&WN0, sizeof(double), 1, fe1op_data);
    if (state->change_byte_order)
      WN0 = *(double *)ByteSwap((char *)&WN0, 8);

    i = fread(&WNSTEP, sizeof(double), 1, fe1op_data);
    if (state->change_byte_order)
      WNSTEP = *(double *)ByteSwap((char *)&WNSTEP, 8);

    i = fread(GCROSS, sizeof(double), n_Ebin * n_WN, fe1op_data);
    if (state->change_byte_order)
    {
      for (i_en = 0; i_en < n_Ebin; i_en++)
        for (i_wn = 0; i_wn < n_WN; i_wn++)
          GCROSS[i_en][i_wn] = *(double *)ByteSwap((char *)(GCROSS + i_wn * 78 + i_en), 8);
    }
    fclose(fe1op_data);
    first = 0;
  }

  WAVENO = state->FREQ / CLIGHTcm;
  kT_eV = state->TK[J] / 1.602176565e-12; // Changing kT from erg/K to eV/K
  if (WAVENO < WN0 || WAVENO > WN0 + WNSTEP * (n_WN - 1))
    return 0.;
  i_wn = (WAVENO - WN0) / WNSTEP;
  FACTOR = (WAVENO - WN0 - i_wn * WNSTEP) / WNSTEP;
  fe1op = 0.e0;
  for (i_en = 0; i_en < n_Ebin; i_en++)
  {
    BOLT = exp(-Ebin[i_en] / kT_eV);
    fe1op += ((GCROSS[i_wn + 1][i_en] - GCROSS[i_wn][i_en]) * FACTOR + GCROSS[i_wn][i_en]) * BOLT;
  }
  return fe1op; ///state->PARTITION_FUNCTIONS[J][state->IXFE1];
}

double CHOP(int J, GlobalState *state) /* Cross-section for CH molecule */
{
  static double CROSSCH[105][15] =
      {{-38.000, -38.000, -38.000, -38.000, -38.000, -38.000, -38.000,            // 0.1
        -38.000, -38.000, -38.000, -38.000, -38.000, -38.000, -38.000, -38.000},  // 0.1
       {-32.727, -31.151, -30.133, -29.432, -28.925, -28.547, -28.257,            // 0.2
        -28.030, -27.848, -27.701, -27.580, -27.479, -27.395, -27.322, -27.261},  // 0.2
       {-31.588, -30.011, -28.993, -28.290, -27.784, -27.405, -27.115,            // 0.3
        -26.887, -26.705, -26.558, -26.437, -26.336, -26.251, -26.179, -26.117},  // 0.3
       {-30.407, -28.830, -27.811, -27.108, -26.601, -26.223, -25.932,            // 0.4
        -25.705, -25.523, -25.376, -25.255, -25.154, -25.069, -24.997, -24.935},  // 0.4
       {-29.513, -27.937, -26.920, -26.218, -25.712, -25.334, -25.043,            // 0.5
        -24.816, -24.635, -24.487, -24.366, -24.266, -24.181, -24.109, -24.047},  // 0.5
       {-28.910, -27.341, -26.327, -25.628, -25.123, -24.746, -24.457,            // 0.6
        -24.230, -24.049, -23.902, -23.782, -23.681, -23.597, -23.525, -23.464},  // 0.6
       {-28.517, -26.961, -25.955, -25.261, -24.760, -24.385, -24.098,            // 0.7
        -23.873, -23.694, -23.548, -23.429, -23.329, -23.245, -23.174, -23.113},  // 0.7
       {-28.213, -26.675, -25.680, -24.993, -24.497, -24.127, -23.843,            // 0.8
        -23.620, -23.443, -23.299, -23.181, -23.082, -22.999, -22.929, -22.869},  // 0.8
       {-27.942, -26.427, -25.446, -24.769, -24.280, -23.915, -23.635,            // 0.9
        -23.416, -23.241, -23.100, -22.983, -22.887, -22.805, -22.736, -22.677},  // 0.9
       {-27.706, -26.210, -25.241, -24.572, -24.088, -23.728, -23.451,            // 1.0
        -23.235, -23.063, -22.923, -22.808, -22.713, -22.633, -22.565, -22.507},  // 1.0
       {-27.475, -26.000, -25.043, -24.382, -23.905, -23.548, -23.275,            // 1.1
        -23.062, -22.891, -22.753, -22.640, -22.546, -22.467, -22.400, -22.343},  // 1.1
       {-27.221, -25.783, -24.844, -24.193, -23.723, -23.372, -23.102,            // 1.2
        -22.892, -22.724, -22.588, -22.476, -22.384, -22.306, -22.240, -22.184},  // 1.2
       {-26.863, -25.506, -24.607, -23.979, -23.523, -23.182, -22.919,            // 1.3
        -22.714, -22.550, -22.417, -22.309, -22.218, -22.142, -22.078, -22.023},  // 1.3
       {-26.685, -25.347, -24.457, -23.835, -23.382, -23.044, -22.784,            // 1.4
        -22.580, -22.418, -22.286, -22.178, -22.089, -22.014, -21.950, -21.896},  // 1.4
       {-26.085, -24.903, -24.105, -23.538, -23.120, -22.805, -22.561,            // 1.5
        -22.370, -22.217, -22.093, -21.991, -21.906, -21.835, -21.775, -21.723},  // 1.5
       {-25.902, -24.727, -23.936, -23.376, -22.964, -22.654, -22.415,            // 1.6
        -22.227, -22.076, -21.955, -21.855, -21.772, -21.702, -21.644, -21.593},  // 1.6
       {-25.215, -24.196, -23.510, -23.019, -22.655, -22.378, -22.163,            // 1.7
        -21.992, -21.855, -21.744, -21.653, -21.577, -21.513, -21.459, -21.412},  // 1.7
       {-24.914, -23.937, -23.284, -22.820, -22.475, -22.212, -22.007,            // 1.8
        -21.845, -21.715, -21.609, -21.522, -21.449, -21.388, -21.336, -21.292},  // 1.8
       {-24.519, -23.637, -23.039, -22.606, -22.281, -22.030, -21.834,            // 1.9
        -21.678, -21.552, -21.450, -21.365, -21.295, -21.236, -21.185, -21.142},  // 1.9
       {-24.086, -23.222, -22.650, -22.246, -21.948, -21.722, -21.546,            // 2.0
        -21.407, -21.296, -21.205, -21.131, -21.070, -21.018, -20.974, -20.937},  // 2.0
       {-23.850, -23.018, -22.472, -22.088, -21.805, -21.590, -21.422,            // 2.1
        -21.289, -21.182, -21.095, -21.024, -20.964, -20.914, -20.872, -20.835},  // 2.1
       {-23.136, -22.445, -21.994, -21.676, -21.440, -21.259, -21.117,            // 2.2
        -21.004, -20.912, -20.837, -20.775, -20.723, -20.679, -20.642, -20.611},  // 2.2
       {-23.199, -22.433, -21.927, -21.573, -21.314, -21.119, -20.969,            // 2.3
        -20.851, -20.758, -20.682, -20.621, -20.571, -20.529, -20.493, -20.463},  // 2.3
       {-22.696, -22.020, -21.585, -21.286, -21.071, -20.912, -20.791,            // 2.4
        -20.697, -20.622, -20.563, -20.514, -20.475, -20.442, -20.414, -20.391},  // 2.4
       {-22.119, -21.557, -21.194, -20.943, -20.761, -20.624, -20.518,            // 2.5
        -20.434, -20.367, -20.313, -20.268, -20.231, -20.201, -20.175, -20.153},  // 2.5
       {-21.855, -21.300, -20.931, -20.673, -20.485, -20.344, -20.235,            // 2.6
        -20.151, -20.084, -20.031, -19.988, -19.953, -19.924, -19.900, -19.880},  // 2.6
       {-21.126, -20.673, -20.382, -20.184, -20.044, -19.943, -19.868,            // 2.7
        -19.811, -19.769, -19.736, -19.710, -19.690, -19.674, -19.662, -19.652},  // 2.7
       {-20.502, -20.150, -19.922, -19.766, -19.657, -19.578, -19.520,            // 2.8
        -19.478, -19.446, -19.422, -19.404, -19.390, -19.379, -19.371, -19.365},  // 2.8
       {-20.030, -19.724, -19.530, -19.399, -19.309, -19.245, -19.199,            // 2.9
        -19.166, -19.142, -19.125, -19.112, -19.103, -19.096, -19.091, -19.088},  // 2.9
       {-19.640, -19.364, -19.189, -19.074, -18.996, -18.943, -18.906,            // 3.0
        -18.881, -18.863, -18.852, -18.844, -18.839, -18.837, -18.836, -18.836},  // 3.0
       {-19.333, -19.092, -18.939, -18.838, -18.770, -18.725, -18.695,            // 3.1
        -18.675, -18.662, -18.655, -18.651, -18.649, -18.649, -18.651, -18.653},  // 3.1
       {-19.070, -18.880, -18.756, -18.674, -18.621, -18.585, -18.562,            // 3.2
        -18.548, -18.540, -18.536, -18.536, -18.537, -18.539, -18.542, -18.546},  // 3.2
       {-18.851, -18.708, -18.617, -18.558, -18.521, -18.498, -18.484,            // 3.3
        -18.477, -18.475, -18.476, -18.478, -18.482, -18.487, -18.493, -18.498},  // 3.3
       {-18.709, -18.599, -18.533, -18.494, -18.471, -18.459, -18.454,            // 3.4
        -18.454, -18.457, -18.462, -18.469, -18.476, -18.483, -18.490, -18.498},  // 3.4
       {-18.656, -18.572, -18.524, -18.497, -18.485, -18.480, -18.482,            // 3.5
        -18.486, -18.493, -18.501, -18.510, -18.519, -18.527, -18.536, -18.544},  // 3.5
       {-18.670, -18.613, -18.582, -18.566, -18.561, -18.562, -18.568,            // 3.6
        -18.575, -18.583, -18.592, -18.601, -18.610, -18.619, -18.627, -18.635},  // 3.6
       {-18.728, -18.700, -18.687, -18.683, -18.685, -18.691, -18.698,            // 3.7
        -18.706, -18.715, -18.723, -18.731, -18.739, -18.745, -18.752, -18.758},  // 3.7
       {-18.839, -18.835, -18.836, -18.842, -18.849, -18.857, -18.865,            // 3.8
        -18.872, -18.878, -18.883, -18.888, -18.892, -18.895, -18.898, -18.900},  // 3.8
       {-19.034, -19.041, -19.049, -19.057, -19.064, -19.069, -19.071,            // 3.9
        -19.071, -19.070, -19.068, -19.065, -19.061, -19.058, -19.054, -19.051},  // 3.9
       {-19.372, -19.378, -19.382, -19.380, -19.372, -19.359, -19.341,            // 4.0
        -19.321, -19.300, -19.280, -19.261, -19.243, -19.227, -19.212, -19.199},  // 4.0
       {-19.780, -19.777, -19.763, -19.732, -19.686, -19.631, -19.573,            // 4.1
        -19.517, -19.465, -19.419, -19.379, -19.344, -19.314, -19.288, -19.265},  // 4.1
       {-20.151, -20.133, -20.087, -20.009, -19.911, -19.810, -19.715,            // 4.2
        -19.631, -19.559, -19.497, -19.446, -19.402, -19.365, -19.333, -19.306},  // 4.2
       {-20.525, -20.454, -20.312, -20.138, -19.970, -19.825, -19.705,            // 4.3
        -19.607, -19.528, -19.464, -19.411, -19.367, -19.330, -19.300, -19.274},  // 4.3
       {-20.869, -20.655, -20.366, -20.104, -19.894, -19.731, -19.604,            // 4.4
        -19.505, -19.426, -19.363, -19.312, -19.271, -19.236, -19.208, -19.184},  // 4.4
       {-21.179, -20.768, -20.380, -20.081, -19.856, -19.686, -19.556,            // 4.5
        -19.454, -19.375, -19.311, -19.260, -19.218, -19.184, -19.155, -19.131},  // 4.5
       {-21.167, -20.601, -20.206, -19.925, -19.719, -19.565, -19.447,            // 4.6
        -19.355, -19.283, -19.226, -19.180, -19.143, -19.112, -19.087, -19.066},  // 4.6
       {-20.918, -20.348, -19.976, -19.720, -19.536, -19.401, -19.299,            // 4.7
        -19.220, -19.159, -19.112, -19.073, -19.043, -19.018, -18.998, -18.981},  // 4.7
       {-20.753, -20.204, -19.847, -19.602, -19.427, -19.299, -19.203,            // 4.8
        -19.129, -19.072, -19.028, -18.993, -18.965, -18.942, -18.924, -18.909},  // 4.8
       {-20.456, -19.987, -19.677, -19.460, -19.302, -19.186, -19.098,            // 4.9
        -19.030, -18.978, -18.937, -18.904, -18.878, -18.857, -18.841, -18.827},  // 4.9
       {-20.154, -19.734, -19.461, -19.272, -19.136, -19.035, -18.960,            // 5.0
        -18.902, -18.858, -18.824, -18.797, -18.775, -18.759, -18.745, -18.735},  // 5.0
       {-19.941, -19.544, -19.288, -19.114, -18.992, -18.903, -18.837,            // 5.1
        -18.788, -18.751, -18.723, -18.701, -18.684, -18.671, -18.661, -18.654},  // 5.1
       {-19.657, -19.321, -19.104, -18.956, -18.853, -18.779, -18.724,            // 5.2
        -18.684, -18.655, -18.632, -18.615, -18.602, -18.592, -18.585, -18.579},  // 5.2
       {-19.388, -19.109, -18.930, -18.810, -18.725, -18.664, -18.620,            // 5.3
        -18.586, -18.562, -18.543, -18.529, -18.518, -18.510, -18.503, -18.498},  // 5.3
       {-19.201, -18.953, -18.794, -18.686, -18.611, -18.556, -18.515,            // 5.4
        -18.485, -18.462, -18.446, -18.433, -18.423, -18.416, -18.410, -18.406},  // 5.4
       {-18.923, -18.719, -18.588, -18.500, -18.439, -18.396, -18.365,            // 5.5
        -18.344, -18.328, -18.318, -18.311, -18.307, -18.304, -18.303, -18.302},  // 5.5
       {-18.614, -18.458, -18.361, -18.298, -18.258, -18.232, -18.216,            // 5.6
        -18.206, -18.202, -18.201, -18.202, -18.205, -18.208, -18.213, -18.218},  // 5.6
       {-18.419, -18.295, -18.222, -18.178, -18.153, -18.139, -18.132,            // 5.7
        -18.131, -18.133, -18.138, -18.143, -18.150, -18.157, -18.164, -18.172},  // 5.7
       {-18.296, -18.201, -18.148, -18.118, -18.101, -18.094, -18.091,            // 5.8
        -18.093, -18.096, -18.101, -18.107, -18.113, -18.120, -18.126, -18.132},  // 5.8
       {-18.021, -17.992, -17.977, -17.970, -17.967, -17.968, -17.970,            // 5.9
        -17.974, -17.978, -17.983, -17.989, -17.994, -18.000, -18.005, -18.011},  // 5.9
       {-17.694, -17.686, -17.686, -17.691, -17.698, -17.708, -17.718,            // 6.0
        -17.729, -17.740, -17.750, -17.761, -17.771, -17.781, -17.790, -17.798},  // 6.0
       {-17.374, -17.384, -17.400, -17.420, -17.440, -17.462, -17.483,            // 6.1
        -17.503, -17.523, -17.541, -17.558, -17.575, -17.590, -17.603, -17.616},  // 6.1
       {-17.169, -17.199, -17.230, -17.262, -17.293, -17.323, -17.351,            // 6.2
        -17.378, -17.404, -17.427, -17.449, -17.469, -17.488, -17.505, -17.520},  // 6.2
       {-17.151, -17.184, -17.217, -17.250, -17.282, -17.313, -17.342,            // 6.3
        -17.369, -17.395, -17.418, -17.440, -17.461, -17.480, -17.497, -17.513},  // 6.3
       {-17.230, -17.260, -17.290, -17.320, -17.348, -17.375, -17.401,            // 6.4
        -17.425, -17.448, -17.469, -17.489, -17.508, -17.525, -17.541, -17.556},  // 6.4
       {-17.379, -17.403, -17.425, -17.446, -17.467, -17.486, -17.505,            // 6.5
        -17.524, -17.541, -17.558, -17.574, -17.588, -17.602, -17.615, -17.627},  // 6.5
       {-17.596, -17.604, -17.609, -17.612, -17.616, -17.622, -17.628,            // 6.6
        -17.636, -17.644, -17.652, -17.661, -17.670, -17.679, -17.687, -17.695},  // 6.6
       {-17.846, -17.823, -17.795, -17.770, -17.750, -17.735, -17.725,            // 6.7
        -17.719, -17.716, -17.715, -17.716, -17.719, -17.722, -17.726, -17.730},  // 6.7
       {-18.089, -18.015, -17.942, -17.882, -17.836, -17.802, -17.777,            // 6.8
        -17.760, -17.748, -17.740, -17.736, -17.734, -17.733, -17.734, -17.736},  // 6.8
       {-18.299, -18.156, -18.038, -17.947, -17.881, -17.833, -17.798,            // 6.9
        -17.774, -17.757, -17.745, -17.738, -17.733, -17.730, -17.729, -17.729},  // 6.9
       {-18.441, -18.243, -18.096, -17.991, -17.915, -17.860, -17.821,            // 7.0
        -17.792, -17.772, -17.757, -17.746, -17.738, -17.733, -17.730, -17.728},  // 7.0
       {-18.474, -18.262, -18.111, -18.004, -17.926, -17.869, -17.826,            // 7.1
        -17.795, -17.771, -17.753, -17.740, -17.730, -17.722, -17.717, -17.713},  // 7.1
       {-18.387, -18.191, -18.053, -17.952, -17.878, -17.823, -17.782,            // 7.2
        -17.752, -17.729, -17.711, -17.698, -17.689, -17.681, -17.676, -17.672},  // 7.2
       {-18.161, -17.990, -17.874, -17.793, -17.736, -17.696, -17.668,            // 7.3
        -17.648, -17.634, -17.625, -17.619, -17.616, -17.614, -17.614, -17.615},  // 7.3
       {-17.908, -17.774, -17.690, -17.637, -17.604, -17.583, -17.572,            // 7.4
        -17.567, -17.566, -17.568, -17.571, -17.576, -17.581, -17.587, -17.593},  // 7.4
       {-17.681, -17.589, -17.540, -17.515, -17.506, -17.505, -17.511,            // 7.5
        -17.520, -17.530, -17.542, -17.554, -17.566, -17.578, -17.589, -17.600},  // 7.5
       {-17.647, -17.606, -17.584, -17.575, -17.573, -17.576, -17.582,            // 7.6
        -17.589, -17.597, -17.605, -17.614, -17.623, -17.631, -17.639, -17.646},  // 7.6
       {-17.300, -17.291, -17.291, -17.297, -17.307, -17.319, -17.333,            // 7.7
        -17.347, -17.361, -17.375, -17.389, -17.402, -17.415, -17.427, -17.438},  // 7.7
       {-16.786, -16.802, -16.825, -16.853, -16.883, -16.914, -16.944,            // 7.8
        -16.974, -17.003, -17.030, -17.055, -17.079, -17.101, -17.122, -17.141},  // 7.8
       {-16.489, -16.533, -16.579, -16.625, -16.670, -16.713, -16.754,            // 7.9
        -16.793, -16.830, -16.864, -16.896, -16.925, -16.952, -16.977, -17.000},  // 7.9
       {-16.694, -16.724, -16.756, -16.789, -16.823, -16.856, -16.888,            // 8.0
        -16.919, -16.949, -16.976, -17.002, -17.026, -17.048, -17.069, -17.088},  // 8.0
       {-16.935, -16.951, -16.971, -16.993, -17.016, -17.040, -17.064,            // 8.1
        -17.088, -17.111, -17.132, -17.153, -17.172, -17.190, -17.206, -17.222},  // 8.1
       {-17.200, -17.208, -17.220, -17.235, -17.251, -17.269, -17.286,            // 8.2
        -17.304, -17.322, -17.338, -17.354, -17.369, -17.384, -17.397, -17.409},  // 8.2
       {-17.597, -17.591, -17.589, -17.590, -17.594, -17.600, -17.608,            // 8.3
        -17.617, -17.626, -17.635, -17.645, -17.654, -17.662, -17.671, -17.679},  // 8.3
       {-18.166, -18.134, -18.107, -18.085, -18.068, -18.056, -18.047,            // 8.4
        -18.041, -18.038, -18.036, -18.035, -18.035, -18.036, -18.038, -18.039},  // 8.4
       {-19.000, -18.917, -18.838, -18.770, -18.714, -18.669, -18.632,            // 8.5
        -18.603, -18.579, -18.560, -18.545, -18.532, -18.522, -18.514, -18.507},  // 8.5
       {-20.313, -19.982, -19.754, -19.592, -19.472, -19.380, -19.309,            // 8.6
        -19.253, -19.208, -19.172, -19.143, -19.119, -19.099, -19.083, -19.069},  // 8.6
       {-19.751, -19.611, -19.520, -19.461, -19.423, -19.398, -19.382,            // 8.7
        -19.372, -19.366, -19.364, -19.363, -19.364, -19.366, -19.368, -19.371},  // 8.7
       {-19.581, -19.431, -19.337, -19.277, -19.240, -19.218, -19.207,            // 8.8
        -19.202, -19.203, -19.207, -19.212, -19.220, -19.228, -19.236, -19.245},  // 8.8
       {-19.685, -19.506, -19.389, -19.311, -19.258, -19.222, -19.199,            // 8.9
        -19.184, -19.175, -19.170, -19.168, -19.169, -19.171, -19.174, -19.177},  // 8.9
       {-19.977, -19.756, -19.606, -19.501, -19.425, -19.370, -19.330,            // 9.0
        -19.300, -19.278, -19.262, -19.250, -19.241, -19.235, -19.230, -19.227},  // 9.0
       {-20.445, -20.158, -19.958, -19.815, -19.711, -19.633, -19.574,            // 9.1
        -19.528, -19.493, -19.465, -19.442, -19.425, -19.410, -19.398, -19.389},  // 9.1
       {-20.980, -20.625, -20.391, -20.229, -20.110, -20.020, -19.949,            // 9.2
        -19.892, -19.846, -19.807, -19.775, -19.748, -19.724, -19.704, -19.687},  // 9.2
       {-21.404, -21.023, -20.771, -20.594, -20.461, -20.358, -20.274,            // 9.3
        -20.205, -20.148, -20.099, -20.058, -20.022, -19.991, -19.965, -19.942},  // 9.3
       {-21.309, -20.970, -20.753, -20.603, -20.495, -20.412, -20.348,            // 9.4
        -20.295, -20.252, -20.215, -20.185, -20.158, -20.135, -20.115, -20.098},  // 9.4
       {-21.221, -20.906, -20.707, -20.574, -20.480, -20.412, -20.361,            // 9.5
        -20.322, -20.292, -20.268, -20.249, -20.233, -20.221, -20.210, -20.201},  // 9.5
       {-21.441, -21.097, -20.878, -20.728, -20.623, -20.546, -20.489,            // 9.6
        -20.446, -20.413, -20.387, -20.368, -20.352, -20.340, -20.330, -20.322},  // 9.6
       {-21.668, -21.305, -21.071, -20.911, -20.797, -20.713, -20.650,            // 9.7
        -20.602, -20.565, -20.536, -20.514, -20.496, -20.481, -20.470, -20.460},  // 9.7
       {-21.926, -21.556, -21.316, -21.150, -21.031, -20.942, -20.874,            // 9.8
        -20.822, -20.782, -20.750, -20.724, -20.704, -20.687, -20.674, -20.663},  // 9.8
       {-22.319, -21.937, -21.686, -21.510, -21.380, -21.282, -21.206,            // 9.9
        -21.147, -21.099, -21.061, -21.031, -21.006, -20.985, -20.968, -20.954},  // 9.9
       {-22.969, -22.561, -22.288, -22.092, -21.945, -21.832, -21.743,            //10.0
        -21.672, -21.616, -21.570, -21.533, -21.503, -21.477, -21.457, -21.439},  //10.0
       {-24.001, -23.527, -23.199, -22.957, -22.772, -22.629, -22.516,            //10.1
        -22.427, -22.355, -22.297, -22.250, -22.212, -22.180, -22.153, -22.131},  //10.1
       {-24.233, -23.774, -23.477, -23.273, -23.128, -23.022, -22.943,            //10.2
        -22.883, -22.837, -22.802, -22.774, -22.752, -22.735, -22.721, -22.710},  //10.2
       {-24.550, -23.913, -23.521, -23.266, -23.094, -22.976, -22.893,            //10.3
        -22.836, -22.796, -22.768, -22.750, -22.737, -22.730, -22.726, -22.725},  //10.3
       {-24.301, -23.665, -23.274, -23.019, -22.848, -22.730, -22.648,            //10.4
        -22.591, -22.552, -22.525, -22.507, -22.495, -22.489, -22.485, -22.485},  //10.4
       {-24.519, -23.883, -23.491, -23.237, -23.065, -22.948, -22.866,            //10.5
        -22.809, -22.770, -22.743, -22.724, -22.713, -22.706, -22.703, -22.702}}; //10.5

  double WAVENO, EVOLT, EN, TN, CROSSCHT[15], CHop;
  int N, IT;

  WAVENO = state->FREQ / CLIGHTcm;
  EVOLT = WAVENO / 8065.479e0;
  N = EVOLT * 10.;
  if (N < 20 || N >= 105)
    return 0.;
  if (state->T[J] >= 9000.)
    return 0.;

  EN = N * 0.1;
  for (IT = 0; IT < 15; IT++)
    CROSSCHT[IT] = CROSSCH[N - 1][IT] + (CROSSCH[N][IT] - CROSSCH[N - 1][IT]) * (EVOLT - EN) / 0.1;
  IT = (state->T[J] - 2000.) / 500.;
  IT = max(IT, 0);
  TN = (IT + 1) * 500. + 1500.;
  CHop = pow10(CROSSCHT[IT] + (CROSSCHT[IT + 1] - CROSSCHT[IT]) * (state->T[J] - TN) / 500.);
  return CHop * state->PARTITION_FUNCTIONS[J][state->IXCH];
}

double NHOP(int J, GlobalState *state)
{
  /*
  Cross-sections of Fe 1 photoionization time
  This routine is based on data provided by Phillip Stancil
  */
  static double WL0, WLSTEP;
  static int n_WL = 4701, n_Temp = 15, first = 1;
  static float T_TBL[15];
  static double GCROSS[4701][15][3];
  double WAVE, factor_wl, factor_temp, f1, f2, NHop;
  int i_wl, i_temp, i;

  if (first)
  {
    char path[512];
    FILE *NHop_data;
    int headlen, n_etrans, ii;
    char head[2048];
    float gauss_fwhm;

    strncpy(path, state->PATH, state->PATHLEN + 1);
    strcat(path, DATAFILE_NH);
    NHop_data = fopen(path, "rb");

    i = fread(&headlen, sizeof(int), 1, NHop_data);
    if (state->change_byte_order)
      headlen = *(int *)ByteSwap((char *)&headlen, 4);

    i = fread(head, 1, headlen, NHop_data);

    i = fread(&gauss_fwhm, sizeof(float), 1, NHop_data);
    if (state->change_byte_order)
      gauss_fwhm = *(float *)ByteSwap((char *)&gauss_fwhm, 4);

    i = fread(&n_etrans, sizeof(int), 1, NHop_data);
    if (state->change_byte_order)
      n_etrans = *(int *)ByteSwap((char *)&n_etrans, 4);

    i = fread(&n_Temp, sizeof(int), 1, NHop_data);
    if (state->change_byte_order)
      n_Temp = *(int *)ByteSwap((char *)&n_Temp, 4);

    i = fread(&n_WL, sizeof(int), 1, NHop_data);
    if (state->change_byte_order)
      n_WL = *(int *)ByteSwap((char *)&n_WL, 4);

    i = fread(&WL0, sizeof(double), 1, NHop_data);
    if (state->change_byte_order)
      WL0 = *(double *)ByteSwap((char *)&WL0, 8);

    i = fread(&WLSTEP, sizeof(double), 1, NHop_data);
    if (state->change_byte_order)
      WLSTEP = *(double *)ByteSwap((char *)&WLSTEP, 8);

    i = fread(T_TBL, sizeof(float), n_Temp, NHop_data);
    if (state->change_byte_order)
    {
      for (i_temp = 0; i_temp < n_Temp; i_temp++)
        T_TBL[i_temp] = *(float *)ByteSwap((char *)(T_TBL + i_temp), 4);
    }
    i = fread(GCROSS, sizeof(double), n_etrans * n_Temp * n_WL, NHop_data);
    if (state->change_byte_order)
    {
      ii = 0;
      for (i_wl = 0; i_wl < n_WL; i_wl++)
        for (i_temp = 0; i_temp < n_Temp; i_temp++)
          for (i = 0; i < 3; i++)
          {
            GCROSS[i_wl][i_temp][i] = *(double *)ByteSwap((char *)(GCROSS + ii), 4);
            ii++;
          }
    }
    fclose(NHop_data);
    first = 0;
  }

  WAVE = CLIGHT / state->FREQ;
  if (WAVE < WL0 || WAVE > WL0 + WLSTEP * (n_WL - 1))
    return 0.;
  if (state->T[J] < T_TBL[0] || state->T[J] > T_TBL[n_Temp - 1])
    return 0.;

  i_wl = (WAVE - WL0) / WLSTEP;
  factor_wl = (WAVE - WL0 - i_wl * WLSTEP) / WLSTEP;

  for (i_temp = 0; i_temp < n_Temp - 1; i_temp++)
    if (T_TBL[i_temp + 1] > state->T[J])
      break;
  factor_temp = (state->T[J] - T_TBL[i_temp]) / (T_TBL[i_temp + 1] - T_TBL[i_temp]);

  f1 = (GCROSS[i_wl][i_temp + 1][0] - GCROSS[i_wl][i_temp][0]) * factor_temp + GCROSS[i_wl][i_temp][0];
  f2 = (GCROSS[i_wl + 1][i_temp + 1][0] - GCROSS[i_wl + 1][i_temp][0]) * factor_temp + GCROSS[i_wl + 1][i_temp][0];
  NHop = (f2 - f1) * factor_wl + f1;

  f1 = (GCROSS[i_wl][i_temp + 1][1] - GCROSS[i_wl][i_temp][1]) * factor_temp + GCROSS[i_wl][i_temp][1];
  f2 = (GCROSS[i_wl + 1][i_temp + 1][1] - GCROSS[i_wl + 1][i_temp][1]) * factor_temp + GCROSS[i_wl + 1][i_temp][1];
  NHop += (f2 - f1) * factor_wl + f1;

  factor_temp = (1. / state->T[J] - 1. / T_TBL[i_temp]) / (1. / T_TBL[i_temp + 1] - 1. / T_TBL[i_temp]);
  f1 = (GCROSS[i_wl][i_temp + 1][2] - GCROSS[i_wl][i_temp][2]) * factor_temp + GCROSS[i_wl][i_temp][2];
  f2 = (GCROSS[i_wl + 1][i_temp + 1][2] - GCROSS[i_wl + 1][i_temp][2]) * factor_temp + GCROSS[i_wl + 1][i_temp][2];
  NHop += pow10((f2 - f1) * factor_wl + f1);

  return NHop * state->PARTITION_FUNCTIONS[J][state->IXNH];
}

double OHOP(int J, GlobalState *state)
{
  static double CROSSOH[130][15] =
      {{-30.855, -29.121, -27.976, -27.166, -26.566, -26.106, -25.742,              // 2.1
        -25.448, -25.207, -25.006, -24.836, -24.691, -24.566, -24.457, -24.363},    // 2.1
       {-30.494, -28.760, -27.615, -26.806, -26.206, -25.745, -25.381,              // 2.2
        -25.088, -24.846, -24.645, -24.475, -24.330, -24.205, -24.097, -24.002},    // 2.2
       {-30.157, -28.425, -27.280, -26.472, -25.872, -25.411, -25.048,              // 2.3
        -24.754, -24.513, -24.312, -24.142, -23.997, -23.872, -23.764, -23.669},    // 2.3
       {-29.848, -28.117, -26.974, -26.165, -25.566, -25.105, -24.742,              // 2.4
        -24.448, -24.207, -24.006, -23.836, -23.692, -23.567, -23.458, -23.364},    // 2.4
       {-29.567, -27.837, -26.693, -25.885, -25.286, -24.826, -24.462,              // 2.5
        -24.169, -23.928, -23.727, -23.557, -23.412, -23.287, -23.179, -23.084},    // 2.5
       {-29.307, -27.578, -26.436, -25.628, -25.029, -24.569, -24.205,              // 2.6
        -23.912, -23.671, -23.470, -23.300, -23.155, -23.031, -22.922, -22.828},    // 2.6
       {-29.068, -27.341, -26.199, -25.391, -24.792, -24.332, -23.969,              // 2.7
        -23.676, -23.435, -23.234, -23.064, -22.920, -22.795, -22.687, -22.592},    // 2.7
       {-28.820, -27.115, -25.978, -25.172, -24.574, -24.115, -23.752,              // 2.8
        -23.459, -23.218, -23.017, -22.848, -22.703, -22.579, -22.470, -22.376},    // 2.8
       {-28.540, -26.891, -25.768, -24.968, -24.372, -23.914, -23.552,              // 2.9
        -23.259, -23.019, -22.818, -22.649, -22.504, -22.380, -22.272, -22.177},    // 2.9
       {-28.275, -26.681, -25.574, -24.779, -24.186, -23.729, -23.368,              // 3.0
        -23.076, -22.836, -22.636, -22.467, -22.322, -22.198, -22.090, -21.996},    // 3.0
       {-27.993, -26.470, -25.388, -24.602, -24.014, -23.560, -23.200,              // 3.1
        -22.909, -22.669, -22.470, -22.301, -22.157, -22.033, -21.925, -21.831},    // 3.1
       {-27.698, -26.252, -25.204, -24.433, -23.851, -23.401, -23.043,              // 3.2
        -22.754, -22.515, -22.316, -22.148, -22.005, -21.881, -21.773, -21.679},    // 3.2
       {-27.398, -26.026, -25.019, -24.267, -23.696, -23.251, -22.896,              // 3.3
        -22.609, -22.372, -22.174, -22.007, -21.864, -21.741, -21.634, -21.540},    // 3.3
       {-27.100, -25.791, -24.828, -24.102, -23.543, -23.106, -22.756,              // 3.4
        -22.472, -22.238, -22.041, -21.875, -21.733, -21.611, -21.504, -21.411},    // 3.4
       {-26.807, -25.549, -24.631, -23.933, -23.391, -22.964, -22.621,              // 3.5
        -22.341, -22.109, -21.915, -21.751, -21.610, -21.488, -21.383, -21.290},    // 3.5
       {-26.531, -25.310, -24.431, -23.761, -23.238, -22.823, -22.488,              // 3.6
        -22.214, -21.986, -21.795, -21.633, -21.494, -21.374, -21.269, -21.178},    // 3.6
       {-26.239, -25.066, -24.225, -23.585, -23.082, -22.681, -22.356,              // 3.7
        -22.089, -21.866, -21.679, -21.520, -21.383, -21.265, -21.162, -21.072},    // 3.7
       {-25.945, -24.824, -24.017, -23.405, -22.923, -22.538, -22.223,              // 3.8
        -21.964, -21.748, -21.565, -21.410, -21.276, -21.160, -21.059, -20.970},    // 3.8
       {-25.663, -24.587, -23.810, -23.222, -22.761, -22.391, -22.088,              // 3.9
        -21.838, -21.629, -21.452, -21.300, -21.170, -21.057, -20.958, -20.872},    // 3.9
       {-25.372, -24.350, -23.603, -23.038, -22.596, -22.241, -21.950,              // 4.0
        -21.710, -21.508, -21.337, -21.190, -21.064, -20.954, -20.858, -20.774},    // 4.0
       {-25.076, -24.111, -23.396, -22.853, -22.429, -22.088, -21.809,              // 4.1
        -21.578, -21.384, -21.220, -21.078, -20.957, -20.851, -20.758, -20.676},    // 4.1
       {-24.779, -23.870, -23.189, -22.669, -22.261, -21.934, -21.667,              // 4.2
        -21.445, -21.259, -21.101, -20.965, -20.848, -20.746, -20.656, -20.578},    // 4.2
       {-24.486, -23.629, -22.983, -22.486, -22.095, -21.781, -21.524,              // 4.3
        -21.311, -21.132, -20.980, -20.850, -20.737, -20.639, -20.553, -20.478},    // 4.3
       {-24.183, -23.382, -22.774, -22.302, -21.928, -21.627, -21.381,              // 4.4
        -21.177, -21.005, -20.859, -20.734, -20.625, -20.531, -20.449, -20.376},    // 4.4
       {-23.867, -23.127, -22.561, -22.116, -21.761, -21.474, -21.238,              // 4.5
        -21.043, -20.878, -20.738, -20.617, -20.513, -20.423, -20.344, -20.274},    // 4.5
       {-23.538, -22.862, -22.340, -21.926, -21.592, -21.320, -21.096,              // 4.6
        -20.909, -20.751, -20.617, -20.502, -20.402, -20.315, -20.239, -20.172},    // 4.6
       {-23.234, -22.604, -22.120, -21.734, -21.422, -21.166, -20.953,              // 4.7
        -20.776, -20.625, -20.497, -20.387, -20.291, -20.208, -20.135, -20.071},    // 4.7
       {-22.934, -22.347, -21.898, -21.541, -21.250, -21.010, -20.811,              // 4.8
        -20.643, -20.500, -20.378, -20.273, -20.182, -20.102, -20.033, -19.971},    // 4.8
       {-22.637, -22.092, -21.676, -21.345, -21.075, -20.853, -20.666,              // 4.9
        -20.508, -20.374, -20.259, -20.159, -20.073, -19.997, -19.931, -19.872},    // 4.9
       {-22.337, -21.835, -21.452, -21.147, -20.899, -20.693, -20.520,              // 5.0
        -20.373, -20.247, -20.139, -20.046, -19.964, -19.892, -19.830, -19.774},    // 5.0
       {-22.049, -21.584, -21.230, -20.950, -20.721, -20.531, -20.372,              // 5.1
        -20.236, -20.119, -20.019, -19.931, -19.855, -19.788, -19.729, -19.676},    // 5.1
       {-21.768, -21.337, -21.011, -20.754, -20.544, -20.370, -20.223,              // 5.2
        -20.098, -19.991, -19.898, -19.817, -19.746, -19.683, -19.628, -19.579},    // 5.2
       {-21.494, -21.096, -20.796, -20.559, -20.367, -20.208, -20.074,              // 5.3
        -19.960, -19.861, -19.776, -19.701, -19.636, -19.578, -19.527, -19.482},    // 5.3
       {-21.233, -20.861, -20.585, -20.368, -20.193, -20.048, -19.926,              // 5.4
        -19.821, -19.732, -19.654, -19.586, -19.526, -19.473, -19.426, -19.384},    // 5.4
       {-20.983, -20.635, -20.380, -20.181, -20.021, -19.889, -19.778,              // 5.5
        -19.683, -19.602, -19.531, -19.469, -19.415, -19.367, -19.324, -19.286},    // 5.5
       {-20.743, -20.418, -20.182, -19.999, -19.853, -19.733, -19.633,              // 5.6
        -19.547, -19.474, -19.410, -19.354, -19.305, -19.261, -19.223, -19.189},    // 5.6
       {-20.515, -20.210, -19.991, -19.824, -19.690, -19.581, -19.490,              // 5.7
        -19.413, -19.347, -19.290, -19.240, -19.196, -19.157, -19.122, -19.092},    // 5.7
       {-20.297, -20.011, -19.808, -19.654, -19.532, -19.434, -19.352,              // 5.8
        -19.282, -19.223, -19.172, -19.127, -19.088, -19.054, -19.023, -18.996},    // 5.8
       {-20.090, -19.822, -19.633, -19.491, -19.381, -19.291, -19.218,              // 5.9
        -19.156, -19.103, -19.057, -19.018, -18.983, -18.952, -18.925, -18.901},    // 5.9
       {-19.893, -19.642, -19.467, -19.337, -19.236, -19.155, -19.089,              // 6.0
        -19.034, -18.987, -18.946, -18.912, -18.881, -18.854, -18.831, -18.810},    // 6.0
       {-19.705, -19.472, -19.309, -19.190, -19.098, -19.025, -18.966,              // 6.1
        -18.917, -18.876, -18.840, -18.810, -18.783, -18.760, -18.739, -18.721},    // 6.1
       {-19.527, -19.310, -19.161, -19.051, -18.968, -18.903, -18.851,              // 6.2
        -18.807, -18.771, -18.740, -18.713, -18.690, -18.670, -18.653, -18.637},    // 6.2
       {-19.357, -19.159, -19.022, -18.922, -18.847, -18.789, -18.743,              // 6.3
        -18.704, -18.673, -18.646, -18.623, -18.603, -18.586, -18.571, -18.558},    // 6.3
       {-19.195, -19.016, -18.892, -18.803, -18.736, -18.684, -18.643,              // 6.4
        -18.610, -18.583, -18.560, -18.540, -18.523, -18.509, -18.496, -18.485},    // 6.4
       {-19.042, -18.883, -18.772, -18.693, -18.634, -18.589, -18.553,              // 6.5
        -18.525, -18.501, -18.481, -18.465, -18.451, -18.438, -18.428, -18.419},    // 6.5
       {-18.894, -18.758, -18.662, -18.593, -18.542, -18.503, -18.473,              // 6.6
        -18.448, -18.428, -18.412, -18.398, -18.386, -18.376, -18.367, -18.359},    // 6.6
       {-18.752, -18.639, -18.559, -18.501, -18.458, -18.426, -18.400,              // 6.7
        -18.380, -18.363, -18.350, -18.338, -18.328, -18.320, -18.313, -18.306},    // 6.7
       {-18.611, -18.523, -18.460, -18.415, -18.381, -18.355, -18.334,              // 6.8
        -18.318, -18.304, -18.293, -18.284, -18.276, -18.269, -18.263, -18.258},    // 6.8
       {-18.471, -18.408, -18.362, -18.329, -18.304, -18.285, -18.269,              // 6.9
        -18.257, -18.247, -18.238, -18.231, -18.224, -18.219, -18.214, -18.210},    // 6.9
       {-18.330, -18.290, -18.261, -18.239, -18.223, -18.211, -18.201,              // 7.0
        -18.192, -18.185, -18.179, -18.174, -18.169, -18.165, -18.162, -18.159},    // 7.0
       {-18.190, -18.168, -18.154, -18.143, -18.135, -18.129, -18.124,              // 7.1
        -18.120, -18.116, -18.112, -18.109, -18.106, -18.104, -18.102, -18.100},    // 7.1
       {-18.055, -18.047, -18.043, -18.042, -18.040, -18.039, -18.039,              // 7.2
        -18.038, -18.037, -18.036, -18.035, -18.034, -18.033, -18.033, -18.032},    // 7.2
       {-17.929, -17.931, -17.935, -17.939, -17.943, -17.946, -17.948,              // 7.3
        -17.950, -17.952, -17.953, -17.955, -17.956, -17.957, -17.958, -17.959},    // 7.3
       {-17.818, -17.826, -17.834, -17.842, -17.849, -17.855, -17.860,              // 7.4
        -17.865, -17.869, -17.872, -17.875, -17.878, -17.881, -17.883, -17.886},    // 7.4
       {-17.724, -17.736, -17.747, -17.758, -17.767, -17.775, -17.782,              // 7.5
        -17.788, -17.793, -17.798, -17.803, -17.807, -17.811, -17.815, -17.819},    // 7.5
       {-17.651, -17.665, -17.678, -17.690, -17.701, -17.710, -17.718,              // 7.6
        -17.725, -17.732, -17.738, -17.744, -17.749, -17.755, -17.760, -17.765},    // 7.6
       {-17.601, -17.615, -17.629, -17.642, -17.653, -17.663, -17.672,              // 7.7
        -17.680, -17.688, -17.695, -17.701, -17.708, -17.714, -17.720, -17.726},    // 7.7
       {-17.572, -17.587, -17.602, -17.614, -17.626, -17.636, -17.645,              // 7.8
        -17.654, -17.662, -17.670, -17.677, -17.684, -17.691, -17.698, -17.704},    // 7.8
       {-17.565, -17.581, -17.595, -17.607, -17.619, -17.629, -17.638,              // 7.9
        -17.647, -17.656, -17.664, -17.671, -17.679, -17.686, -17.693, -17.700},    // 7.9
       {-17.580, -17.594, -17.608, -17.620, -17.630, -17.640, -17.650,              // 8.0
        -17.658, -17.667, -17.675, -17.682, -17.690, -17.697, -17.704, -17.711},    // 8.0
       {-17.613, -17.626, -17.639, -17.649, -17.659, -17.669, -17.677,              // 8.1
        -17.686, -17.694, -17.701, -17.709, -17.716, -17.723, -17.730, -17.737},    // 8.1
       {-17.663, -17.675, -17.685, -17.695, -17.703, -17.711, -17.719,              // 8.2
        -17.727, -17.734, -17.741, -17.748, -17.755, -17.761, -17.768, -17.774},    // 8.2
       {-17.728, -17.737, -17.745, -17.752, -17.759, -17.766, -17.772,              // 8.3
        -17.778, -17.785, -17.791, -17.797, -17.803, -17.808, -17.814, -17.820},    // 8.3
       {-17.803, -17.809, -17.814, -17.818, -17.823, -17.828, -17.832,              // 8.4
        -17.837, -17.842, -17.847, -17.852, -17.856, -17.861, -17.866, -17.871},    // 8.4
       {-17.884, -17.886, -17.888, -17.889, -17.891, -17.893, -17.896,              // 8.5
        -17.899, -17.902, -17.905, -17.908, -17.912, -17.915, -17.919, -17.922},    // 8.5
       {-17.966, -17.964, -17.961, -17.959, -17.958, -17.958, -17.958,              // 8.6
        -17.959, -17.960, -17.961, -17.963, -17.964, -17.966, -17.968, -17.970},    // 8.6
       {-18.040, -18.034, -18.028, -18.023, -18.019, -18.016, -18.013,              // 8.7
        -18.012, -18.010, -18.010, -18.009, -18.009, -18.009, -18.009, -18.010},    // 8.7
       {-18.096, -18.087, -18.078, -18.071, -18.065, -18.059, -18.055,              // 8.8
        -18.051, -18.047, -18.045, -18.042, -18.040, -18.039, -18.037, -18.036},    // 8.8
       {-18.125, -18.115, -18.105, -18.097, -18.089, -18.082, -18.076,              // 8.9
        -18.070, -18.065, -18.061, -18.057, -18.053, -18.051, -18.048, -18.046},    // 8.9
       {-18.120, -18.112, -18.103, -18.095, -18.087, -18.079, -18.072,              // 9.0
        -18.066, -18.060, -18.055, -18.050, -18.046, -18.042, -18.039, -18.036},    // 9.0
       {-18.083, -18.078, -18.071, -18.064, -18.057, -18.050, -18.044,              // 9.1
        -18.037, -18.032, -18.026, -18.022, -18.017, -18.014, -18.010, -18.007},    // 9.1
       {-18.025, -18.022, -18.017, -18.012, -18.006, -18.000, -17.994,              // 9.2
        -17.989, -17.984, -17.979, -17.975, -17.971, -17.968, -17.965, -17.963},    // 9.2
       {-17.957, -17.955, -17.952, -17.948, -17.943, -17.938, -17.934,              // 9.3
        -17.929, -17.925, -17.922, -17.918, -17.916, -17.913, -17.911, -17.910},    // 9.3
       {-17.890, -17.889, -17.886, -17.882, -17.879, -17.875, -17.871,              // 9.4
        -17.867, -17.864, -17.862, -17.860, -17.858, -17.857, -17.856, -17.855},    // 9.4
       {-17.831, -17.829, -17.826, -17.822, -17.819, -17.815, -17.812,              // 9.5
        -17.810, -17.807, -17.806, -17.804, -17.803, -17.803, -17.803, -17.803},    // 9.5
       {-17.786, -17.782, -17.777, -17.773, -17.769, -17.766, -17.763,              // 9.6
        -17.761, -17.759, -17.758, -17.757, -17.757, -17.757, -17.758, -17.759},    // 9.6
       {-17.753, -17.747, -17.741, -17.735, -17.731, -17.727, -17.724,              // 9.7
        -17.722, -17.721, -17.720, -17.720, -17.720, -17.721, -17.722, -17.724},    // 9.7
       {-17.733, -17.724, -17.716, -17.709, -17.703, -17.699, -17.696,              // 9.8
        -17.694, -17.693, -17.692, -17.692, -17.693, -17.694, -17.695, -17.697},    // 9.8
       {-17.723, -17.711, -17.700, -17.691, -17.685, -17.680, -17.676,              // 9.9
        -17.674, -17.673, -17.672, -17.673, -17.673, -17.675, -17.676, -17.678},    // 9.9
       {-17.718, -17.702, -17.689, -17.679, -17.672, -17.667, -17.663,              //10.0
        -17.660, -17.659, -17.659, -17.659, -17.660, -17.661, -17.663, -17.665},    //10.0
       {-17.713, -17.695, -17.681, -17.670, -17.662, -17.656, -17.653,              //10.1
        -17.650, -17.649, -17.649, -17.649, -17.650, -17.651, -17.653, -17.655},    //10.1
       {-17.705, -17.686, -17.671, -17.660, -17.652, -17.647, -17.643,              //10.2
        -17.641, -17.640, -17.640, -17.640, -17.641, -17.643, -17.645, -17.647},    //10.2
       {-17.690, -17.671, -17.657, -17.647, -17.640, -17.635, -17.632,              //10.3
        -17.630, -17.630, -17.630, -17.631, -17.632, -17.634, -17.636, -17.639},    //10.3
       {-17.667, -17.649, -17.637, -17.629, -17.623, -17.619, -17.618,              //10.4
        -17.617, -17.617, -17.618, -17.619, -17.621, -17.623, -17.626, -17.628},    //10.4
       {-17.635, -17.621, -17.611, -17.605, -17.601, -17.600, -17.599,              //10.5
        -17.599, -17.601, -17.602, -17.604, -17.607, -17.609, -17.612, -17.615},    //10.5
       {-17.596, -17.585, -17.579, -17.576, -17.575, -17.575, -17.576,              //10.6
        -17.578, -17.580, -17.582, -17.585, -17.588, -17.591, -17.595, -17.598},    //10.6
       {-17.550, -17.544, -17.542, -17.542, -17.544, -17.546, -17.548,              //10.7
        -17.552, -17.555, -17.558, -17.562, -17.566, -17.570, -17.573, -17.577},    //10.7
       {-17.501, -17.500, -17.501, -17.504, -17.508, -17.513, -17.517,              //10.8
        -17.521, -17.526, -17.530, -17.535, -17.539, -17.544, -17.548, -17.553},    //10.8
       {-17.449, -17.452, -17.457, -17.463, -17.470, -17.476, -17.482,              //10.9
        -17.488, -17.493, -17.499, -17.504, -17.509, -17.514, -17.519, -17.524},    //10.9
       {-17.396, -17.403, -17.412, -17.420, -17.429, -17.437, -17.444,              //11.0
        -17.451, -17.458, -17.464, -17.470, -17.476, -17.481, -17.487, -17.492},    //11.0
       {-17.344, -17.355, -17.366, -17.377, -17.387, -17.396, -17.405,              //11.1
        -17.413, -17.420, -17.427, -17.434, -17.440, -17.446, -17.452, -17.458},    //11.1
       {-17.295, -17.307, -17.321, -17.333, -17.345, -17.355, -17.365,              //11.2
        -17.373, -17.382, -17.389, -17.397, -17.404, -17.410, -17.417, -17.423},    //11.2
       {-17.249, -17.264, -17.278, -17.292, -17.304, -17.316, -17.326,              //11.3
        -17.335, -17.344, -17.352, -17.360, -17.368, -17.375, -17.382, -17.389},    //11.3
       {-17.209, -17.225, -17.241, -17.255, -17.268, -17.280, -17.291,              //11.4
        -17.301, -17.310, -17.319, -17.327, -17.335, -17.343, -17.350, -17.357},    //11.4
       {-17.177, -17.194, -17.210, -17.225, -17.239, -17.251, -17.262,              //11.5
        -17.272, -17.282, -17.291, -17.300, -17.308, -17.316, -17.324, -17.331},    //11.5
       {-17.154, -17.172, -17.189, -17.204, -17.218, -17.230, -17.242,              //11.6
        -17.252, -17.262, -17.272, -17.280, -17.289, -17.298, -17.306, -17.314},    //11.6
       {-17.144, -17.162, -17.179, -17.194, -17.208, -17.220, -17.232,              //11.7
        -17.242, -17.253, -17.262, -17.271, -17.280, -17.289, -17.297, -17.306},    //11.7
       {-17.146, -17.164, -17.181, -17.196, -17.210, -17.222, -17.234,              //11.8
        -17.245, -17.255, -17.265, -17.274, -17.283, -17.292, -17.301, -17.309},    //11.8
       {-17.163, -17.180, -17.197, -17.212, -17.225, -17.237, -17.249,              //11.9
        -17.260, -17.270, -17.280, -17.289, -17.298, -17.307, -17.316, -17.325},    //11.9
       {-17.193, -17.211, -17.227, -17.241, -17.254, -17.266, -17.277,              //12.0
        -17.288, -17.298, -17.308, -17.317, -17.327, -17.336, -17.345, -17.353},    //12.0
       {-17.239, -17.256, -17.271, -17.284, -17.297, -17.309, -17.320,              //12.1
        -17.330, -17.340, -17.350, -17.359, -17.369, -17.378, -17.387, -17.395},    //12.1
       {-17.299, -17.315, -17.329, -17.342, -17.354, -17.365, -17.376,              //12.2
        -17.386, -17.396, -17.405, -17.415, -17.424, -17.433, -17.442, -17.451},    //12.2
       {-17.373, -17.388, -17.402, -17.414, -17.425, -17.436, -17.446,              //12.3
        -17.456, -17.466, -17.475, -17.484, -17.493, -17.502, -17.511, -17.520},    //12.3
       {-17.462, -17.476, -17.489, -17.500, -17.511, -17.521, -17.531,              //12.4
        -17.541, -17.550, -17.559, -17.569, -17.578, -17.587, -17.595, -17.604},    //12.4
       {-17.567, -17.581, -17.592, -17.603, -17.613, -17.623, -17.632,              //12.5
        -17.641, -17.651, -17.660, -17.669, -17.678, -17.686, -17.695, -17.704},    //12.5
       {-17.689, -17.701, -17.712, -17.722, -17.732, -17.741, -17.750,              //12.6
        -17.759, -17.768, -17.777, -17.786, -17.795, -17.803, -17.812, -17.821},    //12.6
       {-17.829, -17.840, -17.851, -17.860, -17.869, -17.878, -17.887,              //12.7
        -17.896, -17.904, -17.913, -17.922, -17.930, -17.939, -17.948, -17.956},    //12.7
       {-17.988, -18.000, -18.010, -18.019, -18.028, -18.036, -18.045,              //12.8
        -18.053, -18.062, -18.070, -18.079, -18.087, -18.096, -18.104, -18.112},    //12.8
       {-18.171, -18.183, -18.192, -18.201, -18.210, -18.218, -18.227,              //12.9
        -18.235, -18.243, -18.252, -18.260, -18.268, -18.277, -18.285, -18.293},    //12.9
       {-18.381, -18.393, -18.403, -18.413, -18.422, -18.430, -18.438,              //13.0
        -18.447, -18.455, -18.463, -18.471, -18.479, -18.487, -18.495, -18.503},    //13.0
       {-18.625, -18.638, -18.650, -18.660, -18.669, -18.678, -18.687,              //13.1
        -18.695, -18.703, -18.711, -18.719, -18.726, -18.734, -18.742, -18.750},    //13.1
       {-18.912, -18.929, -18.943, -18.955, -18.966, -18.975, -18.984,              //13.2
        -18.993, -19.001, -19.008, -19.016, -19.023, -19.031, -19.038, -19.045},    //13.2
       {-19.260, -19.283, -19.303, -19.320, -19.333, -19.345, -19.355,              //13.3
        -19.364, -19.372, -19.380, -19.387, -19.394, -19.400, -19.407, -19.413},    //13.3
       {-19.704, -19.740, -19.771, -19.796, -19.816, -19.832, -19.845,              //13.4
        -19.855, -19.863, -19.870, -19.876, -19.882, -19.887, -19.892, -19.897},    //13.4
       {-20.339, -20.386, -20.424, -20.454, -20.476, -20.492, -20.502,              //13.5
        -20.509, -20.513, -20.516, -20.518, -20.520, -20.521, -20.523, -20.524},    //13.5
       {-21.052, -21.075, -21.093, -21.105, -21.114, -21.120, -21.123,              //13.6
        -21.125, -21.126, -21.127, -21.128, -21.130, -21.131, -21.133, -21.135},    //13.6
       {-21.174, -21.203, -21.230, -21.255, -21.278, -21.299, -21.320,              //13.7
        -21.339, -21.357, -21.375, -21.392, -21.408, -21.424, -21.439, -21.454},    //13.7
       {-21.285, -21.317, -21.346, -21.372, -21.395, -21.416, -21.435,              //13.8
        -21.452, -21.468, -21.483, -21.497, -21.511, -21.524, -21.536, -21.548},    //13.8
       {-21.396, -21.429, -21.459, -21.486, -21.511, -21.532, -21.551,              //13.9
        -21.569, -21.585, -21.600, -21.614, -21.627, -21.640, -21.652, -21.663},    //13.9
       {-21.516, -21.549, -21.580, -21.609, -21.635, -21.658, -21.678,              //14.0
        -21.696, -21.713, -21.728, -21.742, -21.755, -21.767, -21.779, -21.790},    //14.0
       {-21.651, -21.681, -21.711, -21.738, -21.763, -21.785, -21.804,              //14.1
        -21.821, -21.837, -21.851, -21.864, -21.876, -21.887, -21.898, -21.908},    //14.1
       {-21.810, -21.831, -21.853, -21.874, -21.893, -21.910, -21.925,              //14.2
        -21.938, -21.950, -21.961, -21.971, -21.980, -21.989, -21.998, -22.006},    //14.2
       {-22.009, -22.016, -22.026, -22.037, -22.048, -22.058, -22.066,              //14.3
        -22.074, -22.081, -22.088, -22.094, -22.099, -22.105, -22.111, -22.117},    //14.3
       {-22.353, -22.317, -22.296, -22.284, -22.276, -22.270, -22.266,              //14.4
        -22.262, -22.260, -22.258, -22.257, -22.257, -22.257, -22.258, -22.259},    //14.4
       {-22.705, -22.609, -22.552, -22.515, -22.488, -22.468, -22.451,              //14.5
        5 - 22.438, -22.427, -22.418, -22.410, -22.405, -22.400, -22.397, -22.395}, //14.5
       {-22.889, -22.791, -22.731, -22.690, -22.659, -22.634, -22.612,              //14.6
        -22.594, -22.579, -22.566, -22.555, -22.546, -22.539, -22.533, -22.528},    //14.6
       {-23.211, -23.109, -23.041, -22.989, -22.945, -22.906, -22.872,              //14.7
        -22.842, -22.816, -22.793, -22.774, -22.757, -22.743, -22.732, -22.722},    //14.7
       {-25.312, -24.669, -24.250, -23.959, -23.746, -23.587, -23.463,              //14.8
        -23.366, -23.288, -23.225, -23.173, -23.131, -23.095, -23.066, -23.041},    //14.8
       {-25.394, -24.752, -24.333, -24.041, -23.829, -23.669, -23.546,              //14.9
        -23.449, -23.371, -23.308, -23.256, -23.214, -23.178, -23.149, -23.124},    //14.9
       {-25.430, -24.787, -24.369, -24.077, -23.865, -23.705, -23.582,              //15.0
        -23.484, -23.407, -23.344, -23.292, -23.249, -23.214, -23.185, -23.160}};   //15.0

  double WAVENO, EVOLT, EN, TN, CROSSOHT[15], OHop;
  int N, IT;

  WAVENO = state->FREQ / CLIGHTcm;
  EVOLT = WAVENO / 8065.479e0;
  N = EVOLT * 10. - 20.;
  if (N <= 0 || N >= 130)
    return 0.;
  if (state->T[J] >= 9000.)
    return 0.;

  EN = N * 0.1 + 2.;
  for (IT = 0; IT < 15; IT++)
    CROSSOHT[IT] = CROSSOH[N - 1][IT] + (CROSSOH[N][IT] - CROSSOH[N - 1][IT]) * (EVOLT - EN) / 0.1;
  IT = (state->T[J] - 2000.) / 500.;
  IT = max(IT, 0);
  TN = (IT + 1) * 500. + 1500.;
  OHop = pow10(CROSSOHT[IT] + (CROSSOHT[IT + 1] - CROSSOHT[IT]) * (state->T[J] - TN) / 500.);
  return OHop * state->PARTITION_FUNCTIONS[J][state->IXOH];
}

void COOLOP(double *acool, GlobalState *state) /* Si1, Mg1, Al1, C1, Fe1 */
{
  int J;

  if (state->PATHLEN > 0)
  {
    for (J = 0; J < state->NRHOX; J++)
    {
      acool[J] = (C1OP_new(J, state) * state->FRACT[J][state->IXC1] + MG1OP_new(J, state) * state->FRACT[J][state->IXMG1] + AL1OP_new(J, state) * state->FRACT[J][state->IXAL1] + SI1OP_new(J, state) * state->FRACT[J][state->IXSI1] + FE1OP_new(J, state) * state->FRACT[J][state->IXFE1] + CHOP(J, state) * state->FRACT[J][state->IXCH] + NHOP(J, state) * state->FRACT[J][state->IXNH] + OHOP(J, state) * state->FRACT[J][state->IXOH]) * state->STIM[J] / state->RHO[J];
    }
  }
  else
  {
    for (J = 0; J < state->NRHOX; J++)
    {
      acool[J] = (C1OP_new(J, state) * state->FRACT[J][state->IXC1] + MG1OP_new(J, state) * state->FRACT[J][state->IXMG1] + AL1OP_new(J, state) * state->FRACT[J][state->IXAL1] + SI1OP_new(J, state) * state->FRACT[J][state->IXSI1] + FE1OP(J, state) * state->FRACT[J][state->IXFE1] + CHOP(J, state) * state->FRACT[J][state->IXCH] + OHOP(J, state) * state->FRACT[J][state->IXOH]) * state->STIM[J] / state->RHO[J];
    }
  }
  return;
}

double N1OP(int J, GlobalState *state) /* Cross-section */
{
  double C1130, C1020, X1130, X1020, X853;

  C1130 = 6. * exp(-3.575 / state->TKEV[J]);
  C1020 = 10. * exp(-2.384 / state->TKEV[J]);
  X1130 = 0.;
  X1020 = 0.;
  X853 = 0.;
  if (state->FREQ >= 3.517915e15)
    X853 = SEATON(state->FREQ, 3.517915e15, 1.142e-17, 2.0, 4.29);
  if (state->FREQ >= 2.941534e15)
    X1020 = SEATON(state->FREQ, 2.941534e15, 4.410e-18, 1.5, 3.85);
  if (state->FREQ >= 2.653317e15)
    X1130 = SEATON(state->FREQ, 2.653317e15, 4.200e-18, 1.5, 4.34);
  return X853 * 4. + X1020 * C1020 + X1130 * C1130;
}

double O1OP(int J, GlobalState *state) /*  CROSS-SECTION TIMES PARTITION FUNCTION */
{
  return (state->FREQ >= 3.28805e15) ? 9. * SEATON(state->FREQ, 3.28805e15, 2.94e-18, 1., 2.66) : 0;
}

double MG2OP(int J, GlobalState *state) /* CROSS-SECTION TIMES PARTITION FUNCTION */
{
  double C1169, X1169, X824, XXX;

  C1169 = 6. * exp(-4.43 / state->TKEV[J]);
  X1169 = 0.;
  X824 = 0.;

  if (state->FREQ >= 3.635492E15)
    X824 = SEATON(state->FREQ, 3.635492E15, 1.40E-19, 4., 6.7);
  if (state->FREQ >= 2.564306E15)
  {
    XXX = (2.564306E15 / state->FREQ);
    XXX = XXX * XXX * XXX;
    X1169 = 5.11E-19 * XXX;
  }
  return X824 * 2. + X1169 * C1169;
}

double SI2OP(int J, GlobalState *state) /* CROSS-SECTION TIMES THE PARTITION FUNCTION */
{
  static double PEACH[14][6] =
      /*    10000     12000     14000     16000     18000     20000       WAVE(A) */
      {{-43.8941, -43.8941, -43.8941, -43.8941, -43.8941, -43.8941},  /*    500 */
       {-42.2444, -42.2444, -42.2444, -42.2444, -42.2444, -42.2444},  /*    600 */
       {-40.6054, -40.6054, -40.6054, -40.6054, -40.6054, -40.6054},  /*    759 */
       {-54.2389, -52.2906, -50.8799, -49.8033, -48.9485, -48.2490},  /*    760 */
       {-50.4108, -48.4892, -47.1090, -46.0672, -45.2510, -44.5933},  /*   1905 */
       {-52.0936, -50.0741, -48.5999, -47.4676, -46.5649, -45.8246},  /*   1906 */
       {-51.9548, -49.9371, -48.4647, -47.3340, -46.4333, -45.6947},  /*   1975 */
       {-54.2407, -51.7319, -49.9178, -48.5395, -47.4529, -46.5709},  /*   1976 */
       {-52.7355, -50.2218, -48.4059, -47.0267, -45.9402, -45.0592},  /*   3245 */
       {-53.5387, -50.9189, -49.0200, -47.5750, -46.4341, -45.5082},  /*   3246 */
       {-53.2417, -50.6234, -48.7252, -47.2810, -46.1410, -45.2153},  /*   3576 */
       {-53.5097, -50.8535, -48.9263, -47.4586, -46.2994, -45.3581},  /*   3577 */
       {-54.0561, -51.2365, -49.1980, -47.6497, -46.4302, -45.4414},  /*   3900 */
       {-53.8469, -51.0256, -48.9860, -47.4368, -46.2162, -45.2266}}; /*  4200 */
  static double FREQSI[7] = {4.9965417e15, 3.9466738e15, 1.5736321e15,
                             1.5171539e15, 9.2378947e14, 8.3825004e14,
                             7.6869872e14};
  /*     2P,2D,2P,2D,2P */
  static double FLOG[9] = {36.32984, 36.14752, 35.91165, 34.99216, 34.95561,
                           34.45941, 34.36234, 34.27572, 34.20161};
  static double TLG[6] = {9.21034, 9.39266, 9.54681, 9.68034, 9.79813, 9.90349};
  double DT, D, D1, XWL1, XWL2;
  int NT, N;

  NT = min(5, (int)floor(state->T[J] / 2000.) - 4);
  if (NT < 1)
    NT = 1;
  DT = (state->TLOG[J] - TLG[NT - 1]) / (TLG[NT] - TLG[NT - 1]);
  for (N = 0; N < 7; N++)
    if (state->FREQ > FREQSI[N])
      break;
  D = (state->FREQLG - FLOG[N]) / (FLOG[N + 1] - FLOG[N]);
  /* 24-11-2009 Eric Stempels noted a bug when porting this subroutine from FORTRAN
   The checks below should be against 1 and 13 and not 2 and 14 as N is smaller
   by one compared to it FOTRAN counterpart */
  if (N > 1)
    N = 2 * N - 2;
  if (N == 13)
    N = 12;
  D1 = 1. - D;
  XWL1 = PEACH[N + 1][NT - 1] * D + PEACH[N][NT - 1] * D1;
  XWL2 = PEACH[N + 1][NT] * D + PEACH[N][NT] * D1;
  return exp(XWL1 * (1. - DT) + XWL2 * DT) * 6.;
}

double CA2OP(int J, GlobalState *state) /* CROSS-SECTION TIMES THE PARTITION FUNCTION */
{
  double C1218, C1420, X1218, X1420, X1044, XXX;

  C1218 = 10. * exp(-1.697 / state->TKEV[J]);
  C1420 = 6. * exp(-3.142 / state->TKEV[J]);
  X1044 = 0.;
  X1218 = 0.;
  X1420 = 0.;
  if (state->FREQ >= 2.870454e15)
  {
    XXX = (2.870454e15 / state->FREQ);
    XXX = XXX * XXX * XXX;
    X1044 = 1.08e-19 * XXX;
  }
  if (state->FREQ >= 2.460127e15)
    X1218 = 1.64e-17 * sqrt(2.460127e15 / state->FREQ);
  if (state->FREQ >= 2.110779e15)
    X1420 = SEATON(state->FREQ, 2.110779e15, 4.13e-18, 3., 0.69);
  return X1044 + X1218 * C1218 + X1420 * C1420;
}

void LUKEOP(double *aluke, GlobalState *state) /*  SI2,MG2,CA2,N1,O1 */
{
  int J;

  for (J = 0; J < state->NRHOX; J++)
    aluke[J] = (N1OP(J, state) * state->FRACT[J][state->IXN1] + O1OP(J, state) * state->FRACT[J][state->IXO1] +
                MG2OP(J, state) * state->FRACT[J][state->IXMG2] + SI2OP(J, state) * state->FRACT[J][state->IXSI2] +
                CA2OP(J, state) * state->FRACT[J][state->IXCA2]) *
               state->STIM[J] / state->RHO[J];
  return;
}

void HOTOP(double *ahot, GlobalState *state)
{
  static int NUM = 60;
  static double A[420] = {
      4.149945E15, 6.90E-18, 1.000, 6., 6., 13.71, 2., //     6.01
      4.574341E15, 2.50E-18, 1.000, 4., 2., 11.96, 2., //     6.01
      5.220770E15, 1.08E-17, 1.000, 4., 10., 9.28, 2., //     6.01
      5.222307E15, 5.35E-18, 3.769, 2., 1., 0.00, 16., //    10.00
      5.892577E15, 4.60E-18, 1.950, 6., 6., 0.00, 2.,  //     6.01
      6.177022E15, 3.50E-18, 1.000, 4., 12., 5.33, 2., //     6.01
      6.181062E15, 6.75E-18, 3.101, 5., 1., 4.05, 6.,  //     7.01
      6.701879E15, 6.65E-18, 2.789, 5., 5., 1.90, 6.,  //     7.01
      7.158382E15, 6.65E-18, 2.860, 6., 9., 0.00, 6.,  //     7.01

      7.284488E15, 3.43E-18, 4.174, 5., 6., 5.02, 11.,  //     8.01
      7.693612E15, 3.53E-18, 3.808, 5., 10., 3.33, 11., //     8.01
      7.885955E15, 2.32E-18, 3.110, 5., 6., 5.02, 11.,  //     8.01
      8.295079E15, 3.97E-18, 3.033, 5., 10., 3.33, 11., //     8.01
      8.497686E15, 7.32E-18, 3.837, 5., 4., 0.00, 11.,  //     8.01
      8.509966E15, 2.00E-18, 1.750, 7., 3., 12.69, 3.,  //     6.02
      8.572854E15, 1.68E-18, 3.751, 5., 6., 5.02, 11.,  //     8.01
      9.906370E15, 4.16E-18, 2.717, 3., 6., 0.00, 17.,  //    10.01
      1.000693E16, 2.40E-18, 1.750, 7., 9., 6.50, 3.,   //     6.02

      1.046078E16, 4.80E-18, 1.000, 4., 10., 12.53, 7., //     7.02
      1.067157E16, 2.71E-18, 2.148, 3., 6., 0.00, 17.,  //    10.01
      1.146734E16, 2.06E-18, 1.626, 6., 6., 0.00, 7.,   //     7.02
      1.156813E16, 5.20E-19, 2.126, 3., 6., 0.00, 17.,  //    10.01
      1.157840E16, 9.10E-19, 4.750, 4., 1., 0.00, 3.,   //     6.02
      1.177220E16, 5.30E-18, 1.000, 4., 12., 7.10, 7.,  //     7.02
      1.198813E16, 3.97E-18, 2.780, 6., 1., 5.35, 12.,  //     8.02
      1.325920E16, 3.79E-18, 2.777, 6., 5., 2.51, 12.,  //     8.02
      1.327649E16, 3.65E-18, 2.014, 6., 9., 0.00, 12.,  //     8.02

      1.361466E16, 7.00E-18, 1.000, 2., 5., 7.48, 12., //     8.02
      1.365932E16, 9.30E-19, 1.500, 7., 6., 8.00, 4.,  //     6.03
      1.481487E16, 1.10E-18, 1.750, 7., 3., 16.20, 8., //     7.03
      1.490032E16, 5.49E-18, 3.000, 5., 1., 6.91, 18., //    10.02
      1.533389E16, 1.80E-18, 2.277, 4., 9., 0.00, 18., //    10.02
      1.559452E16, 8.70E-19, 3.000, 6., 2., 0.00, 4.,  //     6.03
      1.579688E16, 4.17E-18, 2.074, 4., 5., 3.20, 18., //    10.02
      1.643205E16, 1.39E-18, 2.792, 5., 5., 3.20, 18., //    10.02
      1.656208E16, 2.50E-18, 2.346, 5., 9., 0.00, 18., //    10.02

      1.671401E16, 1.30E-18, 1.750, 7., 9., 8.35, 8.,    //     7.03
      1.719725E16, 1.48E-18, 2.225, 5., 9., 0.00, 18.,   //    10.02
      1.737839E16, 2.70E-18, 1.000, 4., 10., 15.74, 13., //     8.03
      1.871079E16, 1.27E-18, .831, 6., 6., 0.00, 13.,    //     8.03
      1.873298E16, 9.10E-19, 3.000, 4., 1., 0.00, 8.,    //     7.03
      1.903597E16, 2.90E-18, 1.000, 4., 12., 8.88, 13.,  //     8.03
      2.060738E16, 4.60E-18, 1.000, 3., 12., 22.84, 19., //    10.03
      2.125492E16, 5.90E-19, 1.000, 6., 6., 9.99, 9.,    //     7.04
      2.162610E16, 1.69E-18, 1.937, 5., 6., 7.71, 19.,   //    10.03

      2.226127E16, 1.69E-18, 1.841, 5., 10., 5.08, 19., //    10.03
      2.251163E16, 9.30E-19, 2.455, 6., 6., 7.71, 19.,  //    10.03
      2.278001E16, 7.90E-19, 1.000, 6., 9., 10.20, 14., //     8.04
      2.317678E16, 1.65E-18, 2.277, 6., 10., 5.08, 19., //    10.03
      2.348946E16, 3.11E-18, 1.963, 6., 4., 0.00, 19.,  //    10.03
      2.351911E16, 7.30E-19, 1.486, 5., 6., 7.71, 19.,  //    10.03
      2.366973E16, 5.00E-19, 1.000, 4., 2., 0.00, 9.,   //     7.04
      2.507544E16, 6.90E-19, 1.000, 6., 3., 19.69, 14., //     8.04
      2.754065E16, 7.60E-19, 1.000, 2., 1., 0.00, 14.,  //     8.04

      2.864850E16, 1.54E-18, 2.104, 6., 1., 7.92, 20.,  //    10.04
      2.965598E16, 1.53E-18, 2.021, 6., 5., 3.76, 20.,  //    10.04
      3.054151E16, 1.40E-18, 1.471, 6., 9., 0.00, 20.,  //    10.04
      3.085141E16, 2.80E-18, 1.000, 4., 5., 11.01, 20., //    10.04
      3.339687E16, 3.60E-19, 1.000, 6., 2., 0.00, 15.,  //     8.05
      3.818757E16, 4.90E-19, 1.145, 6., 6., 0.00, 21.}; //    10.05
  double FREE, XSECT;
  float XX, TEMP, XNATOM, XNELEC, POTI[8];
  double XNFC[MOSIZE * 6], XNFN[MOSIZE * 6], XNFO[MOSIZE * 6], XNFNE[MOSIZE * 6],
      XNFMG[MOSIZE * 6], XNFSI[MOSIZE * 6], XNFS[MOSIZE * 6], XNFFE[MOSIZE * 6],
      XNFP[MOSIZE * 21];
  int I, J, L, ID, MAXION, IONSIZ, ITAU;

  for (ITAU = 0; ITAU < state->NRHOX; ITAU++)
  {
    TEMP = state->T[ITAU];
    XNELEC = state->XNE[ITAU];
    XNATOM = state->XNA[ITAU];
    J = 2;
    MAXION = IONSIZ = 6;
    I = 6;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFC + 6 * ITAU, J); /* C  */
    I = 7;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFN + 6 * ITAU, J); /* N  */
    I = 8;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFO + 6 * ITAU, J); /* O  */
    I = 10;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFNE + 6 * ITAU, J); /* Ne */
    I = 12;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFMG + 6 * ITAU, J); /* Mg */
    I = 14;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFSI + 6 * ITAU, J); /* Si */
    I = 16;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFS + 6 * ITAU, J); /* S  */
    MAXION = IONSIZ = 5;
    I = 26;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFFE + 6 * ITAU, J);

    J = 1;
    MAXION = IONSIZ = 4;
    I = 6;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFP + 21 * ITAU, J); /* C  */
    MAXION = IONSIZ = 5;
    I = 7;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFP + 21 * ITAU + 4, J); /* N  */
    MAXION = IONSIZ = 6;
    I = 8;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFP + 21 * ITAU + 9, J); /* O  */
    I = 10;
    xsaha_(I, TEMP, XNELEC, XNATOM, MAXION, POTI, XNFP + 21 * ITAU + 15, J); /* Ne */
  }
  /* FREE-FREE */

  for (J = 0; J < state->NRHOX; J++)
  {
    int J2, J3, J4, J5, J6;
    J2 = J * 6 + 1;
    J3 = J2 + 1;
    J4 = J3 + 1;
    J5 = J4 + 1;
    J6 = J5 + 1;
    FREE = COULFF(J, 1, state) * 1. * (XNFC[J2] + XNFN[J2] + XNFO[J2] + XNFNE[J2] + XNFMG[J2] + XNFSI[J2] + XNFS[J2] + XNFFE[J2]) +
           COULFF(J, 2, state) * 4. * (XNFC[J3] + XNFN[J3] + XNFO[J3] + XNFNE[J3] + XNFMG[J3] + XNFSI[J3] + XNFS[J3] + XNFFE[J3]) +
           COULFF(J, 3, state) * 9. * (XNFC[J4] + XNFN[J4] + XNFO[J4] + XNFNE[J4] + XNFMG[J4] + XNFSI[J4] + XNFS[J4] + XNFFE[J4]) +
           COULFF(J, 4, state) * 16. * (XNFC[J5] + XNFN[J5] + XNFO[J5] + XNFNE[J5] + XNFMG[J5] + XNFSI[J5] + XNFS[J5] + XNFFE[J5]) +
           COULFF(J, 5, state) * 25. * (XNFC[J6] + XNFN[J6] + XNFO[J6] + XNFNE[J6] + XNFMG[J6] + XNFSI[J6] + XNFS[J6]);
    ahot[J] = FREE * 3.6919e8 / state->FREQ / state->FREQ / state->FREQ * state->XNE[J] / sqrt(state->T[J]);
  }
  L = -7;
  for (I = 1; I <= NUM; I++)
  {
    L += 7;
    if (state->FREQ < A[L])
      continue;
    XSECT = A[L + 1] * (A[L + 2] + (A[L] / state->FREQ) - A[L + 2] * (A[L] / state->FREQ)) *
            sqrt(pow(A[L] / state->FREQ, ((int)A[L + 3])));
    ID = ((int)A[L + 6]) - 1;
    for (J = 0; J < state->NRHOX; J++)
    {
      XX = XSECT * XNFP[J * 21 + ID] * A[L + 4];
      if (XX > ahot[J] / 100.)
        ahot[J] += XX / exp(A[L + 5] / state->TKEV[J]);
    }
  }
  for (J = 0; J < state->NRHOX; J++)
  {
    ahot[J] *= state->STIM[J] / state->RHO[J];
    /*    printf("%d %f\n",J,ahot[J]); */
  }
}

void ELECOP(double *sigel, GlobalState *state)
{
  int J;

  for (J = 0; J < state->NRHOX; J++)
    sigel[J] = 0.6653e-24 * state->XNE[J] / state->RHO[J];
}

void H2RAOP(double *sigh2, int iH2mol, GlobalState *state)
{
  double WAVE, WW, SIG, ARG;
  int J;

  WAVE = CLIGHT / min(state->FREQ, 2.922e15);
  WW = WAVE * WAVE;
  SIG = (8.14e-13 + 1.28e-6 / WW + 1.61 / (WW * WW)) / (WW * WW);
  for (J = 0; J < state->NRHOX; J++)
  {
    sigh2[J] = state->FRACT[J][iH2mol] * state->PARTITION_FUNCTIONS[J][iH2mol] / state->RHO[J] * SIG;
  }
}

extern "C" char const *SME_DLL GetOpacity(int n, void *arg[], GlobalState *state) /* Returns specific cont. opacity */
{
  short i, j, nrhox, key;
  double *a1;
  IDL_STRING *species, *a4;

  if (n < 3)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  if (!state->flagCONTIN)
  {
    strcpy(state->result, "Opacity has not been calculated");
    return state->result;
  }
  j = *(short *)arg[0]; /* state->IFOP number */
  i = *(short *)arg[1]; /* Length of IDL arrays */
  nrhox = min(state->NRHOX, i);
  a1 = (double *)arg[2];
  switch (j)
  {
  case -3:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->COPSTD[i];
    return &OK_response;
  case -2:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->COPRED[i];
    return &OK_response;
  case -1:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->COPBLU[i];
    return &OK_response;
  case 0:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->AHYD[i];
    return &OK_response;
  case 1:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->AH2P[i];
    return &OK_response;
  case 2:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->AHMIN[i];
    return &OK_response;
  case 3:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->SIGH[i];
    return &OK_response;
  case 4:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->AHE1[i];
    return &OK_response;
  case 5:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->AHE2[i];
    return &OK_response;
  case 6:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->AHEMIN[i];
    return &OK_response;
  case 7:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->SIGHE[i];
    return &OK_response;
  case 8:
    if (n > 3)
    {
      species = (IDL_STRING *)arg[3];
      key = 0;
      if (n == 5)
      {
        a4 = (IDL_STRING *)arg[4];
        if (!strncmp(a4->s, "new", a4->slen))
          key = 1;
        if (!strncmp(a4->s, "old", a4->slen))
          key = 2;
        if (!strncmp(a4->s, "fraction", a4->slen))
          key = 3;
      }
      if (!strcmp(species->s, "C1"))
      {
        switch (key)
        {
        case 0:
          for (i = 0; i < nrhox; i++)
            a1[i] = C1OP_new(i, state) * state->FRACT[i][state->IXC1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        case 1:
          for (i = 0; i < nrhox; i++)
            a1[i] = C1OP_new(i, state);
          return &OK_response;
        case 2:
          for (i = 0; i < nrhox; i++)
            a1[i] = C1OP(i, state);
          return &OK_response;
        case 3:
          for (i = 0; i < nrhox; i++)
            a1[i] = state->FRACT[i][state->IXC1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        }
      }
      else if (!strcmp(species->s, "Mg1"))
      {
        switch (key)
        {
        case 0:
          for (i = 0; i < nrhox; i++)
            a1[i] = MG1OP_new(i, state) * state->FRACT[i][state->IXMG1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        case 1:
          for (i = 0; i < nrhox; i++)
            a1[i] = MG1OP_new(i, state);
          return &OK_response;
        case 2:
          for (i = 0; i < nrhox; i++)
            a1[i] = MG1OP(i, state);
          return &OK_response;
        case 3:
          for (i = 0; i < nrhox; i++)
            a1[i] = state->FRACT[i][state->IXMG1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        }
      }
      else if (!strcmp(species->s, "Al1"))
      {
        switch (key)
        {
        case 0:
          for (i = 0; i < nrhox; i++)
            a1[i] = AL1OP_new(i, state) * state->FRACT[i][state->IXAL1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        case 1:
          for (i = 0; i < nrhox; i++)
            a1[i] = AL1OP_new(i, state);
          return &OK_response;
        case 2:
          for (i = 0; i < nrhox; i++)
            a1[i] = AL1OP(i, state);
          return &OK_response;
        case 3:
          for (i = 0; i < nrhox; i++)
            a1[i] = state->FRACT[i][state->IXAL1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        }
      }
      else if (!strcmp(species->s, "Si1"))
      {
        switch (key)
        {
        case 0:
          for (i = 0; i < nrhox; i++)
            a1[i] = SI1OP_new(i, state) * state->FRACT[i][state->IXSI1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        case 1:
          for (i = 0; i < nrhox; i++)
            a1[i] = SI1OP_new(i, state);
          return &OK_response;
        case 2:
          for (i = 0; i < nrhox; i++)
            a1[i] = SI1OP(i, state);
          return &OK_response;
        case 3:
          for (i = 0; i < nrhox; i++)
            a1[i] = state->FRACT[i][state->IXSI1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        }
      }
      else if (!strcmp(species->s, "Fe1"))
      {
        switch (key)
        {
        case 0:
          for (i = 0; i < nrhox; i++)
            a1[i] = FE1OP_new(i, state) * state->FRACT[i][state->IXFE1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        case 1:
          for (i = 0; i < nrhox; i++)
            a1[i] = FE1OP_new(i, state);
          return &OK_response;
        case 2:
          for (i = 0; i < nrhox; i++)
            a1[i] = FE1OP(i, state);
          return &OK_response;
        case 3:
          for (i = 0; i < nrhox; i++)
            a1[i] = state->FRACT[i][state->IXFE1] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        }
      }
      else if (!strcmp(species->s, "CH"))
      {
        switch (key)
        {
        case 0:
          for (i = 0; i < nrhox; i++)
            a1[i] = CHOP(i, state) * state->FRACT[i][state->IXCH] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        case 1:
          for (i = 0; i < nrhox; i++)
            a1[i] = CHOP(i, state);
          return &OK_response;
        case 2:
          for (i = 0; i < nrhox; i++)
            a1[i] = CHOP(i, state);
          return &OK_response;
        case 3:
          for (i = 0; i < nrhox; i++)
            a1[i] = state->FRACT[i][state->IXCH] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        }
      }
      else if (!strcmp(species->s, "NH"))
      {
        switch (key)
        {
        case 0:
          for (i = 0; i < nrhox; i++)
            a1[i] = NHOP(i, state) * state->FRACT[i][state->IXNH] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        case 1:
          for (i = 0; i < nrhox; i++)
            a1[i] = NHOP(i, state);
          return &OK_response;
        case 2:
          for (i = 0; i < nrhox; i++)
            a1[i] = NHOP(i, state);
          return &OK_response;
        case 3:
          for (i = 0; i < nrhox; i++)
            a1[i] = state->FRACT[i][state->IXNH] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        }
      }
      else if (!strcmp(species->s, "OH"))
      {
        switch (key)
        {
        case 0:
          for (i = 0; i < nrhox; i++)
            a1[i] = OHOP(i, state) * state->FRACT[i][state->IXOH] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        case 1:
          for (i = 0; i < nrhox; i++)
            a1[i] = OHOP(i, state);
          return &OK_response;
        case 2:
          for (i = 0; i < nrhox; i++)
            a1[i] = OHOP(i, state);
          return &OK_response;
        case 3:
          for (i = 0; i < nrhox; i++)
            a1[i] = state->FRACT[i][state->IXOH] * state->STIM[i] / state->RHO[i];
          return &OK_response;
        }
      }
      else
      {
        sprintf(state->result, "SME cannot compute continuous opacity for %s", species->s);
        return state->result;
      }
    }
    else
    {
      for (i = 0; i < nrhox; i++)
        a1[i] = state->ACOOL[i];
      return &OK_response;
    }
  case 9:
    if (n > 3)
    {
      species = (IDL_STRING *)arg[3];
      if (!strcmp(species->s, "N1"))
      {
        for (i = 0; i < nrhox; i++)
          a1[i] = N1OP(i, state) * state->FRACT[i][state->IXN1] * state->STIM[i] / state->RHO[i];
        return &OK_response;
      }
      else if (!strcmp(species->s, "O1"))
      {
        for (i = 0; i < nrhox; i++)
          a1[i] = O1OP(i, state) * state->FRACT[i][state->IXO1] * state->STIM[i] / state->RHO[i];
        return &OK_response;
      }
      else if (!strcmp(species->s, "Mg2"))
      {
        for (i = 0; i < nrhox; i++)
          a1[i] = MG2OP(i, state) * state->FRACT[i][state->IXMG2] * state->STIM[i] / state->RHO[i];
        return &OK_response;
      }
      else if (!strcmp(species->s, "Si2"))
      {
        for (i = 0; i < nrhox; i++)
          a1[i] = SI2OP(i, state) * state->FRACT[i][state->IXSI2] * state->STIM[i] / state->RHO[i];
        return &OK_response;
      }
      else if (!strcmp(species->s, "Ca2"))
      {
        for (i = 0; i < nrhox; i++)
          a1[i] = CA2OP(i, state) * state->FRACT[i][state->IXCA2] * state->STIM[i] / state->RHO[i];
        return &OK_response;
      }
      else
      {
        sprintf(state->result, "SME cannot compute continuous opacity for %s", species->s);
        return state->result;
      }
    }
    else
    {
      for (i = 0; i < nrhox; i++)
        a1[i] = state->ALUKE[i];
      return &OK_response;
    }
  case 10:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->AHOT[i];
    return &OK_response;
  case 11:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->SIGEL[i];
    return &OK_response;
  case 12:
    for (i = 0; i < nrhox; i++)
      a1[i] = state->SIGH2[i];
    return &OK_response;
  default:
    strcpy(state->result, "Wrong opacity switch number");
    return state->result;
  }
}

void AutoIonization(GlobalState *state)
{
  /*  CHECK FOR AUTOIONIZATION LINES */
  int OPEN, LINE;
  double EXUP;
  FILE *file12;

  OPEN = 0;
  for (LINE = 0; LINE < state->NLINES; LINE++)
  {
    state->MARK[LINE] = 0;
    state->AUTOION[LINE] = 0;
    EXUP = state->EXCIT[LINE] + 1. / (state->WLCENT[LINE] * 8065.544e-8);
    if (EXUP >= state->POTION[state->SPINDEX[LINE]])
    {
      if (!OPEN)
      {
        file12 = fopen("syntherr.log", "wt");
        if (file12 != NULL)
          OPEN = 1;
        if (OPEN)
          fprintf(file12, "Lines are numbered from 0\n");
      }
      state->AUTOION[LINE] = 1;
      if (state->GAMQST[LINE] > 0.0 && state->GAMVW[LINE] > 0.0)
      {
        if (OPEN)
          fprintf(file12, "Autoionizing line \'%s\' #%d will be computed\n",
                  strtrim(Terminator(state->SPLIST + 8 * state->SPINDEX[LINE], 8)), LINE);
      }
      else
      {
        if (OPEN)
          fprintf(file12, "Autoionizing line \'%s\' #%d will not be computed\n",
                  strtrim(Terminator(state->SPLIST + 8 * state->SPINDEX[LINE], 8)), LINE);
        state->MARK[LINE] = 2;
      }
    }
  }
  if (OPEN)
    fclose(file12);

  /* IF YOU EVER REMEMBER SOMETHING THAT CAN BE PRECALCULATED,
   JUST PUT IT IN HERE!!! */
}

extern "C" char const *SME_DLL Ionization(int n, void *arg[], GlobalState *state)
{
  /*
   Interface routine between the C++ part of SME the FORTRAN 77 code
   eosmag that solves the equation of molecular equilibrium. All it does
   is to compile the list of species from the line list, pass them to
   the eqcount subroutine in eosmag. eqcount counts the number of
   different species state->N_SPLIST including the basic set defined in eosmag.
   ESO_count_species then allocates the arrays state->SPLIST[state->N_SPLIST] and
   state->SPINDEX[state->NLINES]
  */

  int LINE;
  char *species_list;
  int i, NITER, nelem, eos_mode, pf_mode, j;
  int use_electron_density_from_EOS, use_particle_density_from_EOS,
      use_gas_density_from_EOS;
  short switches;
  char *c, tmpname[13];
  float xna, xne, TEMP, XNATOM, XNELEC, XNA_estim, XNE_estim, RHO_estim,
      Pgas, Pelec, max_Ne_err;
  int dump01, dump02, return_pfs, return1, return2, return3, i_max_Ne_err;

  if (!state->flagMODEL)
  {
    strcpy(state->result, "Model atmosphere not set");
    return state->result;
  }
  if (!state->flagABUND)
  {
    strcpy(state->result, "Abundances not set");
    return state->result;
  }
  if (!state->flagLINELIST)
  {
    strcpy(state->result, "No line list set yet");
    return state->result;
  }
  if (state->SPLIST != NULL)
    FREE(state->SPLIST);

  species_list = NULL;
  CALLOC(species_list, state->NLINES * 8, char);
  if (species_list == NULL)
  {
    strcpy(state->result, "No enough space in EOS_count_species");
    return state->result;
  }

  /* The only allowed argument in call to Ionization contains switches
   indicating that electron and/or particle density
   must be substituted with number densities computed by EOS*/
  if (n > 0)
  {
    switches = *(short *)arg[0];
    use_particle_density_from_EOS = (switches & 0x01);
    use_electron_density_from_EOS = (switches & 0x02);
    use_gas_density_from_EOS = (switches & 0x04);
    dump01 = (switches & 0x08);
    dump02 = (switches & 0x10);
    return_pfs = (switches & 0x20);
  }
  else
  {
    use_particle_density_from_EOS = 0;
    use_electron_density_from_EOS = 0;
    use_gas_density_from_EOS = 0;
    dump01 = 0;
    dump02 = 0;
    return_pfs = 0;
  }

  for (LINE = 0; LINE < state->NLINES; LINE++)
  {
    strncpy(tmpname, state->spname + 8 * LINE, 8);
    tmpname[8] = '\0';
    c = strchr(tmpname, ' ');
    if (c != NULL)
      *c = '\0'; /* Cut the ionization stage */
    strcpy(species_list + 8 * LINE, tmpname);
    i = strlen(tmpname);
    if (i < 8)
      for (; i < 8; i++)
        species_list[8 * LINE + i] = ' ';
  }

  /* First determine the size of the complete list returned by eqcount in as state->N_SPLIST */

  state->N_SPLIST = 0; /* That is to indicate that no default list has been set yet */

  nelem = MAX_ELEM - 1;
  switch (i = eqcount_(ELEMEN + 1, species_list, state->ION, state->NLINES, state->N_SPLIST, nelem, 3, 8))
  {
  case 0:
    break;
  case 1:
    FREE(species_list);
    strcpy(state->result, "EOS_count_species found illegal species name");
    return state->result;
  default:
    FREE(species_list);
    sprintf(state->result, "EOS_count_species - SPLSIZ must be larger than %d", i);
    return state->result;
  }

  /* Now allocate space for the complete list of species and the index */

  CALLOC(state->SPLIST, state->N_SPLIST * 8, char);
  if (state->SPLIST == NULL)
  {
    strcpy(state->result, "Not enough space in EOS_count_species");
    return state->result;
  }

  /* Construct a complete list of species */

  i = 0;
  switch (eqlist_(state->ABUND + 1, ELEMEN + 1, species_list, state->ION, state->SPINDEX, state->SPLIST,
                  state->NLINES, i, state->N_SPLIST, nelem, 3, 8, 8))
  {
  case 0:
    break;
  case 1:
    FREE(species_list);
    FREE(state->SPLIST);
    strcpy(state->result, "EOS_list_species found illegal species name");
    return state->result;
  case 2:
    FREE(species_list);
    FREE(state->SPLIST);
    strcpy(state->result, "EOS_list_species received too small state->N_SPLIST");
    return state->result;
  case 3:
    FREE(species_list);
    FREE(state->SPLIST);
    strcpy(state->result, "EOS_list_species could not match ionization state");
    return state->result;
  case 4:
    FREE(species_list);
    FREE(state->SPLIST);
    strcpy(state->result, "EOS_list_species found e- in the middle of the list");
    return state->result;
  case 5:
    FREE(species_list);
    FREE(state->SPLIST);
    strcpy(state->result, "EOS_list_species - Unreasonable abundances");
    return state->result;
  default:
    FREE(species_list);
    FREE(state->SPLIST);
    strcpy(state->result, "EOS_list_species - this error should never happen");
    return state->result;
  }
  FREE(species_list);
  state->N_SPLIST = i;

  /* Now call the solver for molecular equilibrium eqstat. Parameters are:
   state->T         - temperature (var)
   state->XNA       - atomic number density (var)
   state->XNE       - electron number density (var)
   state->ABUND     - abundances (array)
   ELEMEN    - array of element names (char, should be converted to FORTRAN?)
   AMASS     - atomic masses (array)
   state->SPINDEX   - index for each sp. line to the EOS list of species (array)
   state->SPLIST    - EOS list of species(array of char, created by eqlist, so should
               already be in FORTRAN 77 format)
   state->FRACT     - number densities / partition functions (array of state->N_SPLIST*state->NRHOX)
   state->POTION    - ionization potential for each species (array)
   state->MOLWEIGHT - molecular weight of each species (array)
   state->H1FRACT   - number density of neutral Hydrogen (array of state->NRHOX elements)
   state->HE1FRACT  - number density of neutral Helium (array of state->NRHOX elements)
   state->NLINES    - the number of sp. lines (var)
   state->N_SPLIST  - the total number of species (var)
   xne       - number density of electrons computed by EOS
   xna       - number density of particles computed by EOS
  */

  if (state->FRACT != NULL)
  {
    for (i = 0; i < state->NRHOX_allocated; i++)
      FREE(state->FRACT[i]);
    FREE(state->FRACT);
  }
  if (state->PARTITION_FUNCTIONS != NULL)
  {
    for (i = 0; i < state->NRHOX_allocated; i++)
      FREE(state->PARTITION_FUNCTIONS[i]);
    FREE(state->PARTITION_FUNCTIONS);
  }
  state->flagIONIZ = 0;

  if (state->POTION != NULL)
    FREE(state->POTION);
  if (state->MOLWEIGHT != NULL)
    FREE(state->MOLWEIGHT);

  CALLOC(state->FRACT, state->NRHOX, float *);
  for (i = 0; i < state->NRHOX; i++)
  {
    CALLOC(state->FRACT[i], state->N_SPLIST, float);
    if (state->FRACT[i] == NULL)
    {
      strcpy(state->result, "Ionization: Not enough memory");
      return state->result;
    }
  }
  CALLOC(state->PARTITION_FUNCTIONS, state->NRHOX, float *);
  for (i = 0; i < state->NRHOX; i++)
  {
    CALLOC(state->PARTITION_FUNCTIONS[i], state->N_SPLIST, float);
    if (state->PARTITION_FUNCTIONS[i] == NULL)
    {
      strcpy(state->result, "Ionization: Not enough memory");
      return state->result;
    }
  }
  state->NRHOX_allocated = state->NRHOX;

  CALLOC(state->POTION, state->N_SPLIST, float);
  if (state->POTION == NULL)
  {
    strcpy(state->result, "Ionization: Not enough memory");
    return state->result;
  }

  CALLOC(state->MOLWEIGHT, state->N_SPLIST, float);
  if (state->MOLWEIGHT == NULL)
  {
    strcpy(state->result, "Ionization: Not enough memory");
    return state->result;
  }

  /* Find out the location of continuous absorbers */

  for (i = 0; i < state->N_SPLIST; i++)
  {
    if (!strncmp(state->SPLIST + 8 * i, "H ", 2))
      state->IXH1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "H+ ", 3))
      state->IXH2 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "H- ", 3))
      state->IXHMIN = i;
    else if (!strncmp(state->SPLIST + 8 * i, "H2 ", 3))
      state->IXH2mol = i;
    else if (!strncmp(state->SPLIST + 8 * i, "H2+ ", 4))
      state->IXH2pl = i;
    else if (!strncmp(state->SPLIST + 8 * i, "He ", 3))
      state->IXHE1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "He+ ", 4))
      state->IXHE2 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "He++ ", 5))
      state->IXHE3 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "C ", 2))
      state->IXC1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "Al ", 3))
      state->IXAL1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "Si ", 3))
      state->IXSI1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "Si+ ", 4))
      state->IXSI2 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "Ca ", 3))
      state->IXCA1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "Ca+ ", 4))
      state->IXCA2 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "Mg ", 3))
      state->IXMG1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "Mg+ ", 4))
      state->IXMG2 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "N ", 2))
      state->IXN1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "Fe ", 3))
      state->IXFE1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "O ", 2))
      state->IXO1 = i;
    else if (!strncmp(state->SPLIST + 8 * i, "CH ", 3))
      state->IXCH = i;
    else if (!strncmp(state->SPLIST + 8 * i, "NH ", 3))
      state->IXNH = i;
    else if (!strncmp(state->SPLIST + 8 * i, "OH ", 3))
      state->IXOH = i;
    state->POTION[i] = -1.;
    state->MOLWEIGHT[i] = -1.;
  }

  eos_mode = (use_electron_density_from_EOS) ? 0 : 10;
  if (return_pfs)
  {
    for (i = 0; i < state->NRHOX; i++)
    {
      TEMP = state->T[i];
      Pelec = state->XNE[i] * state->TK[i];
      Pgas = Pelec + state->XNA[i] * state->TK[i];
      eqpf_(TEMP, Pgas, Pelec, state->ABUND + 1, ELEMEN + 1, AMASS + 1,
            nelem, state->SPLIST, state->N_SPLIST, state->PARTITION_FUNCTIONS[i],
            3, 8);
    }
    return &OK_response;
  }

  i_max_Ne_err = -1;
  max_Ne_err = 0.;
  for (i = 0; i < state->NRHOX; i++)
  {
    TEMP = state->T[i];
    Pelec = state->XNE[i] * state->TK[i];
    Pgas = Pelec + state->XNA[i] * state->TK[i];

    eqstat_(eos_mode, TEMP, Pgas, Pelec, state->ABUND + 1, ELEMEN + 1, AMASS + 1,
            nelem, state->SPINDEX, state->SPLIST, state->FRACT[i], state->PARTITION_FUNCTIONS[i], state->POTION,
            state->MOLWEIGHT, state->NLINES, state->N_SPLIST, XNE_estim, XNA_estim, RHO_estim, NITER, 3, 8);

    if (fabs(state->XNE[i] - XNE_estim) / state->XNE[i] > max_Ne_err)
    {
      i_max_Ne_err = i;
      max_Ne_err = fabs(state->XNE[i] - XNE_estim) / state->XNE[i];
    }
    state->H1FRACT[i] = state->FRACT[i][state->IXH1] * state->PARTITION_FUNCTIONS[i][state->IXH1];
    state->HE1FRACT[i] = state->FRACT[i][state->IXHE1] * state->PARTITION_FUNCTIONS[i][state->IXHE1];
    state->H2molFRACT[i] = state->FRACT[i][state->IXH2mol] * state->PARTITION_FUNCTIONS[i][state->IXH2mol];
    state->XNE_eos[i] = XNE_estim;
    state->XNA_eos[i] = XNA_estim;
    state->RHO_eos[i] = RHO_estim;

    if (dump02)
    {
      printf("%f %d %d %s %f %f\n", TEMP, i, 79, Terminator(state->SPLIST + 8 * 79, 8),
             state->PARTITION_FUNCTIONS[i][79], // Fe
             log10(state->FRACT[i][79] * state->PARTITION_FUNCTIONS[i][79] / state->RHO[i]));
      printf("%f %d %d %s %f %f\n", TEMP, i, 80, Terminator(state->SPLIST + 8 * 80, 8),
             state->PARTITION_FUNCTIONS[i][80], // Fe+
             log10(state->FRACT[i][80] * state->PARTITION_FUNCTIONS[i][80] / state->RHO[i]));
      printf("%f %d %d %s %f %f\n", TEMP, i, 145, Terminator(state->SPLIST + 8 * 145, 8),
             state->PARTITION_FUNCTIONS[i][145], // CN
             log10(state->FRACT[i][145] * state->PARTITION_FUNCTIONS[i][145] / state->RHO[i]));
    }

    if (dump01 && i == state->NRHOX - 1)
    {
      printf("Atmospheric layer #%d out of %d (%g %g %g)\n", i, state->NRHOX - 1, state->T[i], state->XNE[i], state->XNA[i]);
      for (j = 0; j < state->N_SPLIST; j++)
        printf("%d %s %f %10.4g %f\n", j, Terminator(state->SPLIST + 8 * j, 8),
               state->PARTITION_FUNCTIONS[i][j],
               state->FRACT[i][j],
               state->FRACT[i][j] / state->RHO[i]);
    }
    state->FRACT[i][state->N_SPLIST - 1] = XNE_estim;
    if (use_electron_density_from_EOS)
      state->XNE[i] = XNE_estim;
    if (use_particle_density_from_EOS)
      state->XNA[i] = XNA_estim;
    if (use_gas_density_from_EOS)
      state->RHO[i] = RHO_estim;
  }
  for (i = 0; i < state->NLINES; i++)
    state->SPINDEX[i]--; /* Index in FORTRAN is 1-based */

  state->flagIONIZ = 1;
  if (max_Ne_err > 0.5)
  {
    sprintf(state->result, "WARNING: EOS-computed electron density differs from the model by %d%% in layer %d",
            round(max_Ne_err * 100), i_max_Ne_err + 1);
    return state->result;
  }

  return &OK_response;
}

extern "C" char const *SME_DLL GetFraction(int n, void *arg[], GlobalState *state)
{
  short i, l, mode;
  IDL_STRING *a0;
  char sp[9];
  int j;
  double *a;

  if (!state->flagMODEL)
  {
    strcpy(state->result, "No model atmosphere has been set");
    return state->result;
  }

  mode = *(short *)arg[1]; /* Return mode=0 - number densities
                                             =1 - partition functions
                                          other - number densities/pf */
  if (!state->flagIONIZ && mode != 1)
  {
    strcpy(state->result, "Molecular-ionization equilibrium was not computed");
    return state->result;
  }

  if (n < 4)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  a0 = (IDL_STRING *)arg[0]; /* Pointer to the name of species */

  if (!strncmp("e-", a0->s, a0->slen))
    mode = 10;          /* Ignore PF when dealing
                                                  with electrons */
  l = *(short *)arg[2]; /* Array length */
  a = (double *)arg[3]; /* Array */

  for (i = 0; i < state->N_SPLIST; i++) /* Search for requested species */
  {
    if (!strncmp(state->SPLIST + 8 * i, a0->s, a0->slen))
    {
      switch (mode)
      {
      case 0:
        for (j = 0; j < min(state->NRHOX, l); j++)
          a[j] = state->FRACT[j][i] *
                 state->PARTITION_FUNCTIONS[j][i];
        return &OK_response;
      case 1:
        for (j = 0; j < min(state->NRHOX, l); j++)
          a[j] = state->PARTITION_FUNCTIONS[j][i];
        return &OK_response;
      default:
        for (j = 0; j < min(state->NRHOX, l); j++)
          a[j] = state->FRACT[j][i];
        return &OK_response;
      }
    }
  }
  sprintf(state->result, "Requested species %s not found", Terminator(a0->s, a0->slen));
  return state->result;
}

extern "C" char const *SME_DLL GetDensity(int n, void *arg[], GlobalState *state)
{
  short l;
  char sp[9];
  int j;
  double *a;

  if (!state->flagMODEL)
  {
    strcpy(state->result, "No model atmosphere has been set");
    return state->result;
  }

  if (!state->flagIONIZ)
  {
    strcpy(state->result, "Molecular-ionization equilibrium was not computed");
    return state->result;
  }

  if (n < 2)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  l = *(short *)arg[0]; /* Array length */
  a = (double *)arg[1]; /* Array */
  for (j = 0; j < min(state->NRHOX, l); j++)
    a[j] = state->RHO_eos[j];
  return &OK_response;
}

extern "C" char const *SME_DLL GetNatom(int n, void *arg[], GlobalState *state)
{
  short l;
  int j;
  double *a;

  if (!state->flagMODEL)
  {
    strcpy(state->result, "No model atmosphere has been set");
    return state->result;
  }

  if (!state->flagIONIZ)
  {
    strcpy(state->result, "Molecular-ionization equilibrium was not computed");
    return state->result;
  }

  if (n < 2)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  l = *(short *)arg[0]; /* Array length */
  a = (double *)arg[1]; /* Array */
  for (j = 0; j < min(state->NRHOX, l); j++)
    a[j] = state->XNA_eos[j];
  return &OK_response;
}

extern "C" char const *SME_DLL GetNelec(int n, void *arg[], GlobalState *state)
{
  short l;
  int j;
  double *a;

  if (!state->flagMODEL)
  {
    strcpy(state->result, "No model atmosphere has been set");
    return state->result;
  }

  if (!state->flagIONIZ)
  {
    strcpy(state->result, "Molecular-ionization equilibrium was not computed");
    return state->result;
  }

  if (n < 2)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  l = *(short *)arg[0]; /* Array length */
  a = (double *)arg[1]; /* Array */
  for (j = 0; j < min(state->NRHOX, l); j++)
    a[j] = state->XNE_eos[j];
  return &OK_response;
}

extern "C" char const *SME_DLL Transf(int n, void *arg[], GlobalState *state)
{
  /*  THIS SUBROUTINE EXPLICITLY SOLVES THE TRANSFER EQUATION
    FOR A SET OF NODES ON THE STAR DISK. THE RESULTS ARE:
    AN ARRAY TABLE(WAVELENGTH) WITH SPECIFIC INTENSITIES
    (LINE OPACITY INCLUDED) AND FC* WITH CONTINUUM INTENSITIES
    AT BOTH ENDS OF SPECTRAL INTERVAL. THE RESULTS ARE
    WRITTEN TO THE FILE #11, AS WELL AS THE INFORMATION ABOUT
    THE NUMBER OF WAVELENGTHS, THE NUMBER OF NODES ON THE DISK,
    MODEL TEMPERATURE AND GRAVITY, THE ABUNDANCE AND
    THE WAVELENGTH RANGE.

    Author: N.Piskunov

    LAST UPDATE: September 13, 1993.
    C++ Version: October 26, 1994
  */

  double *TABLE, *WL, *FCBLUE, *FCRED, *MU, EPS1, EPS2;
  int NWSIZE, NWL;
  int imu, im;
  double MU_sph[MOSIZE], rhox[MUSIZE * MOSIZE], rhox_sph[MUSIZE][2 * MOSIZE],
      P_impact, WW, delta_lambda;
  double opacity_tot[MOSIZE], opacity_cont[MOSIZE], source[MOSIZE],
      source_cont[MOSIZE];
  short NMU, iret, keep_lineop, long_continuum;
  int line;

  /* Check if everything is set and pre-calculated */

  if (!state->flagMODEL)
  {
    strcpy(state->result, "No model atmosphere has been set");
    return state->result;
  }
  if (!state->flagWLRANGE)
  {
    strcpy(state->result, "No wavelength range has been set");
    return state->result;
  }
  if (!state->flagABUND)
  {
    strcpy(state->result, "No list of abundances has been set");
    return state->result;
  }
  if (!state->flagLINELIST)
  {
    strcpy(state->result, "No line list has been set");
    return state->result;
  }
  if (!state->flagIONIZ)
  {
    strcpy(state->result, "Molecular-ionization equilibrium was not computed");
    return state->result;
  }
  if (!state->flagCONTIN)
  {
    strcpy(state->result, "No arrays have been allocated for continous opacity calculations");
    return state->result;
  }
  if (!state->lineOPACITIES)
  {
    strcpy(state->result, "No memory has been allocated for storing line opacities");
    return state->result;
  }

  /* Get the arguments */

  if (n < 9)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  if (n > 10) /* New SME software capable of using predefined wavelength grid */
  {
    NMU = *(short *)arg[0];          /* Number of limb points */
    MU = (double *)arg[1];           /* Array of limb points */
    FCBLUE = (double *)arg[2];       /* Continuum specific intensity on the blue end */
    FCRED = (double *)arg[3];        /* Continuum specific intensity on the red end */
    NWSIZE = *(int *)arg[4];         /* Length of the arrays for synthesis */
    NWL = *(int *)arg[5];            /* Length of predefined wavelength vector */
    WL = (double *)arg[6];           /* Array for wavelengths */
    TABLE = (double *)arg[7];        /* Array for synthetic spectrum */
    EPS1 = *(double *)arg[8];        /* Accuracy of the radiative transfer integration */
    EPS2 = *(double *)arg[9];        /* Accuracy of the interpolation on wl grid */
    keep_lineop = *(short *)arg[10]; /* For several spectral segments there is no 
                                      point recomputing line opacities. This flag
                                      tells when recalculations are needed */
  }
  else /* Old SME software */
  {
    NMU = *(short *)arg[0];    /* Number of limb points */
    MU = (double *)arg[1];     /* Array of limb points */
    FCBLUE = (double *)arg[2]; /* Continuum specific intensity on the blue end */
    FCRED = (double *)arg[3];  /* Continuum specific intensity on the red end */
    NWSIZE = *(long *)arg[4];  /* Length of the arrays for synthesis */
    WL = (double *)arg[5];     /* Array for wavelengths */
    TABLE = (double *)arg[6];  /* Array for synthetic spectrum */
    EPS1 = *(double *)arg[7];  /* Accuracy of the radiative transfer integration */
    EPS2 = *(double *)arg[8];  /* Accuracy of the interpolation on wl grid */
    state->change_byte_order = 0;
  }

  if (NMU > MUSIZE)
  {
    snprintf(state->result, 511, "Specified number of limb angles (%d) exceeds MUSIZE (%d)", NMU, MUSIZE);
    return state->result;
  }

  if (n > 11) /* Check of continuum is needed at every wavelength */
  {           /* If this flag is true FCBLUE must be an arrays of */
              /* the size NWSIZE. On exit FCRED keeps its meaning */
    long_continuum = *(short *)arg[11];
  }
  else
    long_continuum = 0;

  if (!keep_lineop)
  {
    /* Allocate temporary arrays */
    CALLOC(state->YABUND, state->NLINES, double);
    CALLOC(state->XMASS, state->NLINES, double);
    CALLOC(state->EXCUP, state->NLINES, double);
    CALLOC(state->ENU4, state->NLINES, double);
    CALLOC(state->ENL4, state->NLINES, double);
    if (state->ENL4 == NULL)
    {
      strcpy(state->result, "Not enough memory");
      return state->result;
    }

    /* Check autoionization lines */

    AutoIonization(state);

    /* Initialize flags prepare central line opacities and the Voigt function parameters */

    for (line = 0; line < state->NLINES; line++)
    {
      LINEOPAC(line, state);
      if (NWL == 0)
      {
        state->MARK[line] = (state->ALMAX[line] < EPS1) ? 2 : -1;
        state->Wlim_left[line] = max(state->WLCENT[line] - 1000., 0.); /* Initialize line contribution limits */
        state->Wlim_right[line] = min(state->WLCENT[line] + 1000., 2000000.);
      }
      state->ALMAX[line] = 0.;
    }
    FREE(state->YABUND);
    FREE(state->XMASS);
    FREE(state->EXCUP);
    FREE(state->ENU4);
    FREE(state->ENL4);

    // Line contribution limits
    for (line = 0; line < state->NLINES; line++) // Check the line contribution at various detunings
    {
      delta_lambda = 0.2;
      WW = state->WLCENT[line];
      if (state->MARK[line] == -1)
      {
        state->MARK[line] = 0;
        do
        {
          delta_lambda = delta_lambda * 1.5;
          OPMTRX(WW + delta_lambda, opacity_tot, opacity_cont,
                 source, source_cont, line, line, state); // Assess line contribution at a given offset
        } while (state->ALMAX[line] > EPS1);
        state->Wlim_left[line] = max(WW - delta_lambda, 0.);
        state->Wlim_right[line] = min(WW + delta_lambda, 2000000.);
      }
    }
  }

  if (state->MOTYPE == 3) /* If things get spherical initialize a 2D array of MUs and do the RT */
  {
    double sintheta, deltaR, meanR, meanZ, path;
    int nrhox, grazing[MUSIZE], NRHOXs[MUSIZE];
    /*
    The main idea here is that we simply scale up delta m (or delta tau) by the ratio of
    geometrical path along the ray and along the radius. Rays are characterized by the impact
    parameter P that is derived from Mu at the outer surface. Z distance along the ray is
    measured from the plane perpendicular to the line-of-sight and crossing the stellar center.
    The main relation is: Z^2 = R^2 - P^2.
            Z2 - Z1   (Z2^2 - Z1^2)   R2 + R1   R2 + R1
    dZ/dR = ------- = ------------- * ------- = -------.
            R2 - R1   (R2^2 - R1^2)   Z2 + Z1   Z2 + Z1
    The corresponding change in dm is then:
                      dZ            Rmean
    dm_sph = dm_rad * -- = dm_rad * -----
                      dR            Zmean
    */
    for (imu = 0; imu < NMU; imu++)
    {
      P_impact = (state->RADIUS + state->RAD_ATMO[0]) * sqrt(1. - MU[imu] * MU[imu]);
      grazing[imu] = (P_impact > state->RADIUS + state->RAD_ATMO[state->NRHOX - 1]) ? 1 : 0;
      if (grazing[imu]) /* Dealing with grazing rays that do not penetrate optically thick layers */
      {
        for (nrhox = 1; nrhox < state->NRHOX; nrhox++)
          if (P_impact >= state->RADIUS + state->RAD_ATMO[nrhox])
            break;
        deltaR = state->RAD_ATMO[nrhox - 1] - state->RAD_ATMO[nrhox]; // The layer where we do not cross both
        path = state->RAD_ATMO[nrhox - 1] + state->RADIUS;            // boundaries gets special treatment
        path = 2. * sqrt(path * path - P_impact * P_impact);          // Geometrical path through the inner ring
        rhox_sph[imu][0] = state->RHOX[0] / MU[imu];                  // Scale the top mass value by projected path
        for (im = 1; im < nrhox; im++)                                // Loop from the surface to the deepest layer
        {
          meanR = state->RAD_ATMO[im] + state->RAD_ATMO[im - 1] + 2 * state->RADIUS;
          meanZ = sqrt((state->RAD_ATMO[im] + state->RADIUS) * (state->RAD_ATMO[im] + state->RADIUS) - P_impact * P_impact) +
                  sqrt((state->RAD_ATMO[im - 1] + state->RADIUS) * (state->RAD_ATMO[im - 1] + state->RADIUS) - P_impact * P_impact);
          rhox_sph[imu][im] = rhox_sph[imu][im - 1] + (state->RHOX[im] - state->RHOX[im - 1]) * meanR / meanZ;
        }
        rhox_sph[imu][nrhox] = rhox_sph[imu][nrhox - 1] + // Column mass across the deepest layer
                               path * (state->RHOX[nrhox] - state->RHOX[nrhox - 1]) / (state->RAD_ATMO[nrhox - 1] - state->RAD_ATMO[nrhox]);
        for (im = nrhox + 1; im < 2 * nrhox; im++) // The rest of the grazing ray back to the surface
        {                                          // We have column mass chunks stored in rhox_sph already
          rhox_sph[imu][im] = rhox_sph[imu][im - 1] + (rhox_sph[imu][2 * nrhox - im] - rhox_sph[imu][2 * nrhox - im - 1]);
        }
        NRHOXs[imu] = 2 * nrhox;
      }
      else /* Normal rays are treated as in plane parallel case except for variable Mu */
      {
        rhox_sph[imu][0] = state->RHOX[0] / MU[imu]; // Scale the top mass value by projected path
        for (im = 1; im < state->NRHOX; im++)
        {
          meanR = state->RAD_ATMO[im] + state->RAD_ATMO[im - 1] + 2 * state->RADIUS;
          meanZ = sqrt((state->RAD_ATMO[im] + state->RADIUS) * (state->RAD_ATMO[im] + state->RADIUS) - P_impact * P_impact) +
                  sqrt((state->RAD_ATMO[im - 1] + state->RADIUS) * (state->RAD_ATMO[im - 1] + state->RADIUS) - P_impact * P_impact);
          rhox_sph[imu][im] = rhox_sph[imu][im - 1] + (state->RHOX[im] - state->RHOX[im - 1]) * meanR / meanZ;
        }
        NRHOXs[imu] = state->NRHOX;
      }
    }
    iret = RKINTS_sph(rhox_sph, NMU, NRHOXs, EPS1, EPS2, FCBLUE, FCRED, TABLE, NWSIZE, NWL,
                      WL, long_continuum, grazing, state);
  }
  else /* Plane-parallel case is handled by simpler routine RKINTS which
          is responsible for the adaptive wavelength grid */
  {
    for (imu = 0; imu < NMU; imu++) /* Prepare state->RHOX arrays for each Mu */
    {
      for (im = 0; im < state->NRHOX; im++)
        rhox[imu * state->NRHOX + im] = state->RHOX[im] / MU[imu];
    }
    iret = RKINTS(rhox, NMU, EPS1, EPS2, FCBLUE, FCRED, TABLE, NWSIZE, NWL,
                  WL, long_continuum, state);
  }

  *((int *)arg[5]) = NWL;

  return iret ? "Not enough array length to store all the points" : "";
}

extern "C" char const *SME_DLL GetLineRange(int n, void *arg[], GlobalState *state) /* Get importance range for every line */
{
  int nlines, line;
  double *b;

  if (!state->flagMODEL)
  {
    strcpy(state->result, "No model atmosphere has been set");
    return state->result;
  }
  if (!state->flagWLRANGE)
  {
    strcpy(state->result, "No wavelength range has been set");
    return state->result;
  }
  if (!state->flagABUND)
  {
    strcpy(state->result, "No list of abundances has been set");
    return state->result;
  }
  if (!state->flagLINELIST)
  {
    strcpy(state->result, "No line list has been set");
    return state->result;
  }
  if (!state->flagIONIZ)
  {
    strcpy(state->result, "Molecular-ionization equilibrium was not computed");
    return state->result;
  }
  if (!state->flagCONTIN)
  {
    strcpy(state->result, "No arrays have been allocated for continous opacity calculations");
    return state->result;
  }
  if (!state->lineOPACITIES)
  {
    strcpy(state->result, "No memory has been allocated for storing line opacities");
    return state->result;
  }

  if (n < 2) // Check if arguments are present
  {
    strcpy(state->result, "GetLineRange: Requires an double array pointer and its length");
    return state->result;
  }

  b = (double *)arg[0];
  nlines = *(int *)arg[1];

  for (line = 0; line < min(nlines, state->NLINES); line++)
  {
    if (state->MARK[line])
    {
      b[2 * line] = b[2 * line + 1] = state->WLCENT[line];
    }
    else
    {
      b[2 * line] = state->Wlim_left[line];
      b[2 * line + 1] = state->Wlim_right[line];
    }
  }

  return &OK_response;
}

extern "C" char const *SME_DLL CentralDepth(int n, void *arg[], GlobalState *state)
{
  /*
  THIS SUBROUTINE EXPLICITLY SOLVES THE TRANSFER EQUATION
  FOR A SET OF NODES ON THE STAR DISK IN THE CENTERS OF SPETRAL
  LINES. THE RESULTS ARE SPECIFIC INTENSITIES

  Author: N.Piskunov

  LAST UPDATE: September 13, 1993.
  C++ Version: January 15, 1999
  */

  double TBL[81], WEIGHTS[81], *MU, EPS1, FC, s0, s1, opacity[MOSIZE], wlstd;
  float *TABLE;
  int NMU, IMU, line, im, IM, NWSIZE;

  /* Check if everything is set and pre-calculated */

  if (!state->flagMODEL)
  {
    strcpy(state->result, "No model atmosphere has been set");
    return state->result;
  }
  if (!state->flagWLRANGE)
  {
    strcpy(state->result, "No wavelength range has been set");
    return state->result;
  }
  if (!state->flagABUND)
  {
    strcpy(state->result, "No list of abundances has been set");
    return state->result;
  }
  if (!state->flagLINELIST)
  {
    strcpy(state->result, "No line list has been set");
    return state->result;
  }
  if (!state->flagIONIZ)
  {
    strcpy(state->result, "Molecular-ionization equilibrium was not computed");
    return state->result;
  }
  if (!state->flagCONTIN)
  {
    strcpy(state->result, "No arrays have been allocated for continous opacity calculations");
    return state->result;
  }
  if (!state->lineOPACITIES)
  {
    strcpy(state->result, "No memory has been allocated for storing line opacities");
    return state->result;
  }

  /* Get the arguments */

  if (n < 5)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  NMU = *(int *)arg[0]; /* Number of limb points */
  if (NMU > 81)
  {
    strcpy(state->result, "SME library is limited to maximum 81 mu angles");
    return state->result;
  }
  MU = (double *)arg[1];    /* Array of limb points */
  NWSIZE = *(int *)arg[2];  /* Length of the arrays for synthesis */
  TABLE = (float *)arg[3];  /* Array for synthetic spectrum */
  EPS1 = *(double *)arg[4]; /* Accuracy of the radiative transfer integration */
  if (NWSIZE < state->NLINES)
  {
    strcpy(state->result, "Array size is smaller than the number of sp.lines");
    return state->result;
  }

  /* Check autoionization lines */

  AutoIonization(state);

  /* Initialize intensity vector */

  for (line = 0; line < state->NLINES; line++)
  {
    TABLE[line] = 0.;
  }

  /* Calculate weights for combining intensities into fluxes. The normalized
   weights are proportional to the projected area represented by each mu
   value. The annular area between consecutive mu values is divided equally
   between the two mu values. The first mu value in the list is assumed to
   be the largest, and the corresponding region extends all the way to disk
   center. The final mu value is assumed to be the smallest, and the region
   extends all the way to the limb. */

  s1 = 0.0;
  for (IMU = 0; IMU < NMU; IMU++)
  {
    s0 = s1;
    s1 = (IMU < NMU - 1) ? 1.0 - 0.5 * (MU[IMU] * MU[IMU] + MU[IMU + 1] * MU[IMU + 1]) : 1.0;
    WEIGHTS[IMU] = s1 - s0;
  }

  /* INTEGRATE TRANSFER EQUATION FOR SPECIFIC INTENSITIES */

  CONTOP(state->WLSTD, state->COPSTD, state);
  for (line = 0; line < state->NLINES; line++)
  {
    FC = 0.0;
    CONTOP(state->WLCENT[line], opacity, state); /* Compute continuous opacity at the line center */

    CENTERINTG(MU, NMU, line, opacity, TBL, state);
    for (IMU = 0; IMU < NMU; IMU++)
    {

      TABLE[line] = TABLE[line] + WEIGHTS[IMU] * TBL[IMU];
      FC = FC + WEIGHTS[IMU] * FCINTG(MU[IMU], state->WLCENT[line], opacity, state);
    }

    TABLE[line] = (TABLE[line] < FC) ? 1.0 - TABLE[line] / FC : 0.0;
  }

  return &OK_response;
}

#define EPS3 6.
#define DVEL_MIN 3.e4 // minimum wavelength points spacing in velocity scale [cm/s] \
                      // corresponding to R=1000000 with 2 point sampling

int RKINTS_sph(double rhox[][2 * MOSIZE], int NMU, int NRHOXs[], double EPS1, double EPS2,
               double *FCBLUE, double *FCRED, double *TABLE, int NWSIZE, int &NWL,
               double *WL, short long_continuum, int grazing[], GlobalState *state)
{
  /* 
  THIS SUBROUTINE CALLS SUBROUTINE FCINTG TO INTEGRATE THE EMMERGING
  SPECIFIC INTENSITIES FOR CONTINUUM AT THE EDGES OF SPECTRAL
  INTERVAL (returned as "FC*") AND SUBROUTINE TBINTG FOR THE LINE
  (returned as "TABLE").

  Author: N.Piskunov

  UPDATES: 13-Sep-1993 written.
            26-Oct-1994 C++ Version
            25-Sep-2010 Modified to allow for spherical geometry in 1D models
            12-Jan-2015 Modified the loop limits according to the new approximation
                        for grazing rays
  */
  double WW, FCL, FNORM;
  double opacity_tot[2 * MOSIZE], opacity_cont[2 * MOSIZE],
      source[2 * MOSIZE], source_cont[2 * MOSIZE];
  double DWL_MIN;
  int nrhox;
  int line, line_first, line_last, i, IMU, IM, IWL;

  /* If the wavelength grid is pre-set, just do the calculations */

  if (NWL > 0 && NWL <= NWSIZE)
  {
    for (IWL = 0; IWL < NWL; IWL++)
    {
      OPMTRX(WL[IWL], opacity_tot, opacity_cont,
             source, source_cont, 0, state->NLINES - 1, state);

      for (IMU = 0; IMU < NMU; IMU++)
      {
        nrhox = NRHOXs[IMU];
        if (grazing[IMU])
        {
          for (IM = 0; IM < nrhox / 2; IM++)
          {
            opacity_tot[nrhox - IM - 1] = opacity_tot[IM];
            opacity_cont[nrhox - IM - 1] = opacity_cont[IM];
            source[nrhox - IM - 1] = source[IM];
            source_cont[nrhox - IM - 1] = source_cont[IM];
          }
        }
        TBINTG_sph(nrhox, rhox[IMU], opacity_tot, source, TABLE + IWL * NMU + IMU, grazing[IMU], state);
        if (long_continuum)
        {
          TBINTG_sph(nrhox, rhox[IMU], opacity_cont, source_cont, FCBLUE + IWL * NMU + IMU, grazing[IMU], state);
          if (IMU == 0)
            FNORM = FCBLUE[IWL * NMU];
        }
        else if (fabs(WL[IWL] - state->WFIRST) < 1.e-4)
        {
          TBINTG_sph(nrhox, rhox[IMU], opacity_cont, source_cont, FCBLUE + IMU, grazing[IMU], state);
        }
        if (fabs(WL[IWL] - state->WLAST) < 1.e-4)
        {
          TBINTG_sph(nrhox, rhox[IMU], opacity_cont, source_cont, FCRED + IMU, grazing[IMU], state);
        }
      }
    }
    return 0;
  }

  /* Wavelength grid is not pre-set. Construct an adaptive grid starting from the blue */

  WL[0] = state->WFIRST;
  OPMTRX(state->WFIRST, opacity_tot, opacity_cont,
         source, source_cont, 0, state->NLINES - 1, state);
  for (IMU = 0; IMU < NMU; IMU++)
  {
    nrhox = NRHOXs[IMU];
    if (grazing[IMU])
    {
      for (IM = 0; IM < nrhox / 2; IM++)
      {
        opacity_tot[nrhox - IM - 1] = opacity_tot[IM];
        opacity_cont[nrhox - IM - 1] = opacity_cont[IM];
        source[nrhox - IM - 1] = source[IM];
        source_cont[nrhox - IM - 1] = source_cont[IM];
      }
    }
    TBINTG_sph(nrhox, rhox[IMU], opacity_tot, source, TABLE + IMU, grazing[IMU], state);
    TBINTG_sph(nrhox, rhox[IMU], opacity_cont, source_cont, FCBLUE + IMU, grazing[IMU], state);
    if (IMU == 0)
      FNORM = FCBLUE[IMU];
  }

  /* Now add one line point at each line center. Check line contribution. */

  IWL = 0;
  for (line = 0; line < state->NLINES; line++)
  {
    WW = state->WLCENT[line];
    DWL_MIN = WW * DVEL_MIN / CLIGHTcm;
    if (WW > state->WFIRST && WW < state->WLAST && WW - WL[IWL] > DWL_MIN && !state->MARK[line])
    {
      // Next pair of wavelength points associated with the next line
      IWL++;
      if (IWL > NWSIZE - 1)
        return 1;
      WL[IWL] = (WW + WL[IWL - 1]) * 0.5; // Intermediate wavelength step

      OPMTRX(WL[IWL], opacity_tot, opacity_cont,
             source, source_cont, 0, state->NLINES - 1, state);
      if (state->Wlim_right[line] > WL[IWL] && state->WLCENT[line] <= WL[IWL] &&
          state->ALMAX[line] < EPS1)
        state->Wlim_right[line] = WL[IWL];
      if (state->Wlim_left[line] < WL[IWL] && state->WLCENT[line] > WL[IWL] &&
          state->ALMAX[line] < EPS1)
        state->Wlim_left[line] = WL[IWL];

      for (IMU = 0; IMU < NMU; IMU++)
      {
        nrhox = NRHOXs[IMU];
        if (grazing[IMU])
        {
          for (IM = 0; IM < nrhox / 2; IM++)
          {
            opacity_tot[nrhox - IM - 1] = opacity_tot[IM];
            opacity_cont[nrhox - IM - 1] = opacity_cont[IM];
            source[nrhox - IM - 1] = source[IM];
            source_cont[nrhox - IM - 1] = source_cont[IM];
          }
        }
        TBINTG_sph(nrhox, rhox[IMU], opacity_tot, source, TABLE + IWL * NMU + IMU, grazing[IMU], state);
        if (long_continuum)
        {
          TBINTG_sph(nrhox, rhox[IMU], opacity_cont, source_cont, FCBLUE + IWL * NMU + IMU, grazing[IMU], state);
          if (IMU == 0)
            FNORM = FCBLUE[IWL * NMU];
        }
      }

      // 2nd point in the pair
      IWL++;
      if (IWL >= NWSIZE - 1)
        return 1;
      WL[IWL] = WW; // Put a point in the line center

      OPMTRX(WL[IWL], opacity_tot, opacity_cont,
             source, source_cont, 0, state->NLINES - 1, state);
      if (state->Wlim_right[line] > WL[IWL] && state->WLCENT[line] <= WL[IWL] &&
          state->ALMAX[line] < EPS1)
        state->Wlim_right[line] = WL[IWL];
      if (state->Wlim_left[line] < WL[IWL] && state->WLCENT[line] > WL[IWL] &&
          state->ALMAX[line] < EPS1)
        state->Wlim_left[line] = WL[IWL];

      for (IMU = 0; IMU < NMU; IMU++)
      {
        nrhox = NRHOXs[IMU];
        if (grazing[IMU])
        {
          for (IM = 0; IM < nrhox / 2; IM++)
          {
            opacity_tot[nrhox - IM - 1] = opacity_tot[IM];
            opacity_cont[nrhox - IM - 1] = opacity_cont[IM];
            source[nrhox - IM - 1] = source[IM];
            source_cont[nrhox - IM - 1] = source_cont[IM];
          }
        }
        TBINTG_sph(nrhox, rhox[IMU], opacity_tot, source, TABLE + IWL * NMU + IMU, grazing[IMU], state);
        if (long_continuum)
        {
          TBINTG_sph(nrhox, rhox[IMU], opacity_cont, source_cont, FCBLUE + IWL * NMU + IMU, grazing[IMU], state);
          if (IMU == 0)
            FNORM = FCBLUE[IWL * NMU];
        }
      }
    }
  }

  /* One more point at the red end of the spectral interval */

  DWL_MIN = state->WLAST * DVEL_MIN / CLIGHTcm;
  if (state->WLAST - WL[IWL] > DWL_MIN)
    IWL++;
  if (IWL > NWSIZE - 1)
    return 1;
  WL[IWL] = state->WLAST;
  OPMTRX(WL[IWL], opacity_tot, opacity_cont, source, source_cont, 0, state->NLINES - 1, state);
  for (IMU = 0; IMU < NMU; IMU++)
  {
    nrhox = NRHOXs[IMU];
    if (grazing[IMU])
    {
      for (IM = 0; IM < nrhox / 2; IM++)
      {
        opacity_tot[nrhox - IM - 1] = opacity_tot[IM];
        opacity_cont[nrhox - IM - 1] = opacity_cont[IM];
        source[nrhox - IM - 1] = source[IM];
        source_cont[nrhox - IM - 1] = source_cont[IM];
      }
    }
    TBINTG_sph(nrhox, rhox[IMU], opacity_tot, source, TABLE + IWL * NMU + IMU, grazing[IMU], state);
    TBINTG_sph(nrhox, rhox[IMU], opacity_cont, source_cont, FCRED + IMU, grazing[IMU], state);
    if (long_continuum)
      FCBLUE[IWL * NMU + IMU] = FCRED[IMU];
    FNORM = (FCBLUE[0] + FCRED[0]) * 0.5;
  }
  NWL = IWL + 1;

  /*  AND NOW ADJUST STEP SIZE OF ABS(TABLE(IWL)-TABLE(IWL-1)) IS TOO BIG */

  IWL = 1;
  line_first = 0;
  line_last = state->NLINES - 1;
  while (IWL < NWL)
  {
    if (NWL >= NWSIZE - 1)
      return 1;
    for (i = NWL; i > IWL; i--)
    {
      WL[i] = WL[i - 1];
      for (IMU = 0; IMU < NMU; IMU++)
        TABLE[i * NMU + IMU] = TABLE[(i - 1) * NMU + IMU];
      if (long_continuum)
      {
        for (IMU = 0; IMU < NMU; IMU++)
          FCBLUE[i * NMU + IMU] = FCBLUE[(i - 1) * NMU + IMU];
      }
    }
    WL[IWL] = (WL[IWL] + WL[IWL - 1]) * 0.5;
    NWL++;

    /* Get the value of the middle point */

    OPMTRX(WL[IWL], opacity_tot, opacity_cont,
           source, source_cont, line_first, line_last, state);
    for (IMU = 0; IMU < NMU; IMU++)
    {
      nrhox = NRHOXs[IMU];
      if (grazing[IMU])
      {
        for (IM = 0; IM < nrhox / 2; IM++)
        {
          opacity_tot[nrhox - IM - 1] = opacity_tot[IM];
          opacity_cont[nrhox - IM - 1] = opacity_cont[IM];
          source[nrhox - IM - 1] = source[IM];
          source_cont[nrhox - IM - 1] = source_cont[IM];
        }
      }
      TBINTG_sph(nrhox, rhox[IMU], opacity_tot, source, TABLE + IWL * NMU + IMU, grazing[IMU], state);
      if (long_continuum)
      {
        TBINTG_sph(nrhox, rhox[IMU], opacity_cont, source_cont, FCBLUE + IWL * NMU + IMU, grazing[IMU], state);
        if (IMU == 0)
          FNORM = FCBLUE[IWL * NMU];
      }
    }

    FCL = fabs(TABLE[IWL * NMU] - 0.5 * (TABLE[(IWL - 1) * NMU] + TABLE[(IWL + 1) * NMU])) +
          0.005 * fabs(TABLE[(IWL - 1) * NMU] - TABLE[(IWL + 1) * NMU]);
    FCL /= FNORM;

    /* Here is a new version that I hope is fiinally robust */

    DWL_MIN = WL[IWL - 1] * DVEL_MIN / CLIGHTcm;
    if (FCL < EPS2 || WL[IWL] - WL[IWL - 1] <= DWL_MIN) /* Check if linear approx. is OK */
    {
      /*
      Now we will move right of the WL(IWL) and will never come back, mark
      permanently all weak lines left of this wavelength. Unmark all
      temporary marked lines.
      */

      /* Here is a new version that I hope is finally robust */

      for (line = state->NLINES - 1; line >= line_last; line--)
      {
        if (state->Wlim_left[line] < WL[IWL + 2] && state->MARK[line] == 0)
        {
          line_last = line;
          break;
        }
      }
      for (line = line_first; line <= line_last; line++)
      {
        if (state->Wlim_right[line] > WL[IWL] && state->MARK[line] == 0)
        {
          line_first = line;
          break;
        }
      }
      IWL += 2; /* Advance to the next point */
    }
    else
    {
      /* At this point we are about to add more points to the left, so we can
      ignore all weak lines to the right of this wavelength. */

      for (line = 0; line <= line_first; line++)
      {
        if (state->Wlim_right[line] > WL[IWL - 1] && state->MARK[line] == 0)
        {
          line_first = line;
          break;
        }
      }
      for (line = line_last; line >= line_first; line--)
      {
        if (state->Wlim_left[line] < WL[IWL] && state->MARK[line] == 0)
        {
          line_last = line;
          break;
        }
      }
    }
  }
  return 0;
}

int RKINTS(double *rhox, int NMU, double EPS1, double EPS2,
           double *FCBLUE, double *FCRED, double *TABLE,
           int NWSIZE, int &NWL, double *WL,
           short long_continuum, GlobalState *state)
{
  /*
  THIS SUBROUTINE CALLS SUBROUTINE FCINTG TO INTEGRATE THE EMERGING
  SPECIFIC INTENSITIES FOR CONTINUUM AT THE EDGES OF SPECTRAL
  INTERVAL (returned as "FC*") AND SUBROUTINE TBINTG FOR THE LINE
  (returned as "TABLE").

  Author: N.Piskunov

  UPDATES:  13-Sep-1993 written.
            26-Oct-1994 C++ Version
            25-Sep-2010 Modified to allow for spherical geometry in 1D models
  */
  double WW, FCL, FNORM;
  double opacity_tot[MOSIZE], opacity_cont[MOSIZE], source[MOSIZE],
      source_cont[MOSIZE];
  double ddd, opacity_tot_n[MOSIZE], opacity_cont_n[MOSIZE];
  double DWL_MIN;
  int line, line_first, line_last, i, IMU, IM, IWL, NNWL;

  if (NWL > 0 && NWL <= NWSIZE) // If the wavelength grid is preset, just do it
  {                             // No adaptive grid in this case
    if (!long_continuum)
    {
      OPMTRX(state->WFIRST, opacity_tot, opacity_cont, source, source_cont, 0, state->NLINES - 1, state);
      TBINTG(NMU, rhox, opacity_cont, source_cont, FCBLUE, state);
    }

    line_first = 0;
    line_last = state->NLINES - 1;
    while (state->Wlim_right[line_first] < WL[0] && line_first < line_last)
      line_first++;
    while (state->Wlim_left[line_last] > WL[NWL - 1] && line_first < line_last)
      line_last--;

    NNWL = NWL;
    for (IWL = 0; IWL < NNWL; IWL++)
    {
      OPMTRX(WL[IWL], opacity_tot, opacity_cont, source, source_cont, line_first, line_last, state);
      TBINTG(NMU, rhox, opacity_tot, source, TABLE + IWL * NMU, state);
      if (long_continuum)
      {
        TBINTG(NMU, rhox, opacity_cont, source_cont, FCBLUE + IWL * NMU, state);
      }
    }
    OPMTRX(state->WLAST, opacity_tot, opacity_cont, source, source_cont, 0, state->NLINES - 1, state);
    TBINTG(NMU, rhox, opacity_cont, source_cont, FCRED, state);
    return 0;
  }

  /* CALCULATE CONTINUUM FLUX FOR BOTH ENDS OF THE INTERVAL
   FIRST WE CALCULATE FLUX AT THE BLUE END OF SPECTRAL INTERVAL */

  WL[0] = state->WFIRST;
  OPMTRX(state->WFIRST, opacity_tot, opacity_cont, source, source_cont, 0, state->NLINES - 1, state);

  TBINTG(NMU, rhox, opacity_tot, source, TABLE, state);
  TBINTG(NMU, rhox, opacity_cont, source_cont, FCBLUE, state);
  FNORM = FCBLUE[0];

  /*  Add one point at each line center and one in between */

  IWL = 0;
  for (line = 0; line < state->NLINES; line++)
  {
    WW = state->WLCENT[line];
    DWL_MIN = WW * DVEL_MIN / CLIGHTcm;
    if (WW > state->WFIRST && WW < state->WLAST && WW - WL[IWL] > DWL_MIN && !state->MARK[line])
    {
      IWL++;
      if (IWL > NWSIZE - 1)
        return 1;
      // Add one point between the previous point and the next line center
      WL[IWL] = (WW + WL[IWL - 1]) * 0.5; // Half-way between the next line center and the previous wavelength point
      OPMTRX(WL[IWL], opacity_tot, opacity_cont, source, source_cont, 0, state->NLINES - 1, state);
      if (state->Wlim_right[line] > WL[IWL] && state->WLCENT[line] <= WL[IWL] &&
          state->ALMAX[line] < EPS1)
        state->Wlim_right[line] = WL[IWL];
      if (state->Wlim_left[line] < WL[IWL] && state->WLCENT[line] > WL[IWL] &&
          state->ALMAX[line] < EPS1)
        state->Wlim_left[line] = WL[IWL];
      TBINTG(NMU, rhox, opacity_tot, source, TABLE + IWL * NMU, state);
      if (long_continuum)
      {
        TBINTG(NMU, rhox, opacity_cont, source_cont, FCBLUE + IWL * NMU, state);
        FNORM = FCBLUE[IWL * NMU];
      }

      // Add one point at the line center and test if line is at all important
      IWL++;
      if (IWL > NWSIZE - 1)
        return 1;
      WL[IWL] = WW; // Smack in the next line center
      OPMTRX(WL[IWL], opacity_tot, opacity_cont, source, source_cont, 0, state->NLINES - 1, state);
      if (state->Wlim_right[line] > WL[IWL] && state->WLCENT[line] <= WL[IWL] &&
          state->ALMAX[line] < EPS1)
        state->Wlim_right[line] = WL[IWL];
      if (state->Wlim_left[line] < WL[IWL] && state->WLCENT[line] > WL[IWL] &&
          state->ALMAX[line] < EPS1)
        state->Wlim_left[line] = WL[IWL];
      TBINTG(NMU, rhox, opacity_tot, source, TABLE + IWL * NMU, state);
      if (long_continuum)
      {
        state->debug_print = 0;
        TBINTG(NMU, rhox, opacity_cont, source_cont, FCBLUE + IWL * NMU, state);
        FNORM = FCBLUE[IWL * NMU];
      }

      if (1. - TABLE[IWL * NMU] / FNORM < EPS2)
        state->MARK[line] = 2;
    }
  }

  /* ... and finally add one more point at the red end of the spectral interval */

  DWL_MIN = state->WLAST * DVEL_MIN / CLIGHTcm;
  if (state->WLAST - WL[IWL] > DWL_MIN)
    IWL++;
  if (IWL > NWSIZE - 1)
    return 1;
  WL[IWL] = state->WLAST;
  OPMTRX(WL[IWL], opacity_tot, opacity_cont, source, source_cont, 0, state->NLINES - 1, state);
  TBINTG(NMU, rhox, opacity_tot, source, TABLE + IWL * NMU, state);
  state->debug_print = 1;
  TBINTG(NMU, rhox, opacity_cont, source_cont, FCRED, state);
  state->debug_print = 0;
  if (long_continuum)
  {
    for (IMU = 0; IMU < NMU; IMU++)
      FCBLUE[IWL * NMU + IMU] = FCRED[IMU];
  }
  else
  {
    FNORM = (FCBLUE[0] + FCRED[0]) * 0.5;
  }
  NWL = IWL + 1;

  /* Now we go on refining the wavelength grid based on comparing the actual value
   disk center intensity with linear interpolation between adjacent points */

  IWL = 1;
  line_first = 0;
  line_last = state->NLINES - 1;
  while (IWL < NWL)
  {
    if (NWL >= NWSIZE - 1)
      return 1;
    for (i = NWL; i > IWL; i--)
    {
      WL[i] = WL[i - 1];
      for (IMU = 0; IMU < NMU; IMU++)
        TABLE[i * NMU + IMU] = TABLE[(i - 1) * NMU + IMU];
      if (long_continuum)
      {
        for (IMU = 0; IMU < NMU; IMU++)
          FCBLUE[i * NMU + IMU] = FCBLUE[(i - 1) * NMU + IMU];
      }
    }
    WL[IWL] = (WL[IWL] + WL[IWL - 1]) * 0.5;
    NWL++;

    /* Get the value of the middle point */

    OPMTRX(WL[IWL], opacity_tot, opacity_cont, source, source_cont,
           line_first, line_last, state);

    TBINTG(NMU, rhox, opacity_tot, source, TABLE + IWL * NMU, state);
    if (long_continuum)
    {
      TBINTG(NMU, rhox, opacity_cont, source_cont, FCBLUE + IWL * NMU, state);
      FNORM = FCBLUE[IWL * NMU];
    }

    FCL = fabs(TABLE[IWL * NMU] - 0.5 * (TABLE[(IWL - 1) * NMU] + TABLE[(IWL + 1) * NMU])) +
          0.005 * fabs(TABLE[(IWL - 1) * NMU] - TABLE[(IWL + 1) * NMU]);
    FCL /= FNORM;

    DWL_MIN = WL[IWL] * DVEL_MIN / CLIGHTcm;
    if (FCL < EPS2 || WL[IWL] - WL[IWL - 1] <= DWL_MIN) /* Check if linear approx. is OK */
    {
      /*  Now we will move right of the WL(IWL) and will never comeback, mark
      permanently all weak lines left of this wavelength. Unmark all
      temporary marked lines. Here is a new version that I hope is fiinally robust */

      for (line = state->NLINES - 1; line >= line_last; line--)
      {
        if (state->Wlim_left[line] < WL[IWL + 2])
        {
          line_last = line;
          break;
        }
      }
      for (line = line_first; line <= line_last; line++)
      {
        if (state->Wlim_right[line] > WL[IWL])
        {
          line_first = line;
          break;
        }
      }

      IWL += 2; /* Advance to the next point */
    }
    else
    {
      /* At this point we are about to add more points to the left, so we can
      ignore all weak lines to the right of this wavelength. */

      for (line = 0; line <= line_first; line++)
      {
        if (state->Wlim_right[line] > WL[IWL - 1])
        {
          line_first = line;
          break;
        }
      }
      for (line = line_last; line >= line_first; line--)
      {
        if (state->Wlim_left[line] < WL[IWL])
        {
          line_last = line;
          break;
        }
      }
    }
  }
  return 0;
}

#undef EPS3
#undef DVEL_MIN

#define FLUX_SCALE 1.0686475e5

double FCINTG(double MU, double WAVE, double *COPWL, GlobalState *state)
{
  /*
  Quadratic DELO with Bezier spline RT solver
  AUTHOR: N.Piskunov
  LAST UPDATE: May 4, 2009
  */
  double OPC_A, OPC_B, OPC_C, SRC_A, SRC_B, SRC_C, INTENSITY;
  double CNTR_AB, CNTR_BC, SPRIME_A, SPRIME_B;
  double STEP_AB, STEP_BC, DER, DER1, DELTA, DELTA1;
  double ALPHA, BETA, GAMMA, EPS, B, LAMBDA, SPRIME_SAVE, DBNU;
  double CONWL5, HNUK;
  int IM;

  /* Useful things for the Planck function */

  CONWL5 = exp(50.7649141 - 5. * log(WAVE));
  HNUK = 1.43868e8 / WAVE;

  SRC_B = CONWL5 / (exp(HNUK / state->T[state->NRHOX - 1]) - 1.); // Source function
  SRC_C = CONWL5 / (exp(HNUK / state->T[state->NRHOX - 2]) - 1.);
  OPC_B = (state->MOTYPE == 0) ? COPWL[state->NRHOX - 1] / state->COPSTD[state->NRHOX - 1] : COPWL[state->NRHOX - 1]; // Opacities
  OPC_C = (state->MOTYPE == 0) ? COPWL[state->NRHOX - 2] / state->COPSTD[state->NRHOX - 2] : COPWL[state->NRHOX - 2];

  DBNU = 2.0 * (SRC_B - SRC_C) / ((state->RHOX[state->NRHOX - 1] - state->RHOX[state->NRHOX - 2]) * (OPC_B + OPC_C)) * MU;
  INTENSITY = 0.5 * (SRC_B + SRC_C) + DBNU; // Intensity at the bottom

  SPRIME_SAVE = 0.0; // Initialize S'

  for (IM = state->NRHOX - 2; IM > 0; IM--) // Work your way from the deepest
  {                                         // layer to the surface
    SRC_A = SRC_B;                          // Shift source functions and opacities
    OPC_A = OPC_B;
    SRC_B = SRC_C;
    OPC_B = OPC_C;
    SRC_C = CONWL5 / (exp(HNUK / state->T[IM - 1]) - 1.); // Downwind point
    OPC_C = (state->MOTYPE == 0) ? COPWL[IM - 1] / state->COPSTD[IM - 1] : COPWL[IM - 1];
    /*
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! New version based on monotoneous quadratic Bezier splines
    !
    ! If we define for points A and B along a ray:
    !    u = (tau - tau_a)/(tau_b - tau_a)
    ! then any function can be fit with a Bezier spline as
    !    f(u) = f(tau_a) * (1 - u)^2 + f(tau_b) * u^2 + 2*C*u*(1-u)
    ! where C is the local control parameter.
    !
    ! We solve RT using short characteristics method in order to get the intensity
    ! propagating through point IM in the direction IM+1->IM->IM+1:
    ! I_b = eps * I_a + b
    !  where: b       = alpha * S_a + beta * S_b + gamma * Cont_ab
    !         eps     = exp(-delta)
    !         delta   = tau_b - tau_a
    !         delta'  = tau_c - tau_b
    !         alpha   = (1 - 2/delta) + 2/delta^2 * (1- eps)
    !         beta    = 2/delta^2 * (1 - eps) - eps * (1 + 2/delta)
    !         gamma   = 2/delta * (1 + eps) - 4/delta^2 * (1 - eps)
    !         S_a     - source function in the upwind point A
    !         S_b     - source function in the central point B
    !         Cont_ab - local control parameter
    !
    !  Control parameter for interval [x_a, x_b] can be computed in two ways
    !    C' = f(x_a) + delta/2*S'_a
    !  and
    !    C" = f(x_b) - delta/2*S'_b
    !
    !  We take the mean for all intermediate steps: Cont_ab = (C' + C") / 2
    !  For the first step:                          Cont_ab = C"
    !  For the last step:                           Cont_ab = C'
    !
    !  If D(b-1/2)*D(b+1/2) > 0 then
    !    S'_b  = D(b-1/2)*D(b+1/2) / (lambda*D(b+1/2) + (1-lambda)*D(b-1/2))
    !  Else
    !    S'_b  = 0
    !
    !         D(b-1/2) = (S_b - S_a) / delta
    !         D(b+1/2) = (S_c - S_b) / delta'
    !         lambda   = [1 + delta'/(delta + delta')]/3
    !
    ! A few additional notations:
    !         U_0   = 1 - eps
    !         U_1   = 2/delta
    !         U_2   = 2/delta^2 = U_1/delta
    !         U_3   = U_0 * U_1
    !         U_4   = U_3 / delta
    !         alpha = (1 - U_1) + U_4           = (delta^2 - 2*delta + 2 - 2*eps)/delta^2
    !         beta  = U_4 - eps * (1 + U_1)     = [2 - (2 + 2*delta + delta^2)*eps]/delta^2
    !         gamma = U_1 * (1 + eps) - 2 * U_4 = [2*delta - 4 + (2*delta + 4)*eps]/delta^2
    !
    ! Special care must be take when delta is small.
    ! In this case (using x instead of delta to make formulas shorter)
    !
    !         eps = exp(-x) = 1 - x + x^2/2 - x^3/6 + x^4/24 - x^5/120
    !         U_1 = 2/x
    !         1 - eps = 1 - exp(-x) = x - x^2/2 + x^3/6 - x^4/24 + x^5/120
    !         U_4 = (1 - eps)*2/x^2 = 2/x - 1 - x/3 - x^2/12 + x^3/60
    ! and
    !         alpha = 1 -U_1 + U_4 = x/3 - x^2/12 + x^3/60
    !         beta  = U_4 - eps*(1 + U_1) = x/3 - x^2/4 + x^3/10
    !         gamma = U_1 * (1 + eps) - 2 * U_4 = x/3 -x^2/6 + x^3/20
    !
    ! Note that we kept the 3rd order in x throughout the whole expansion.
    !
    ! In order to compute delta and delta' we approximate the opacity between
    ! points [A,B] and [B,C] with Bezier spline as explained above and integrate
    ! the optical path analytically. Note that the control parameters are different
    ! for [A,B] and [B,C]:
    !    delta   = L_ab/3*(k_a + k_b + C_ab)
    !    delta'  = L_bc/3*(k_b + k_c + C_bc)
    !
    !    C_ab = k_b - d_ab/2*S'_b
    !    C_bc = k_b + d_bc/2*S'_b
    !
    ! Now to the the actual computing. delta and delta' first (assuming equispaced
    ! geometrical grid lambda is 1/2):
    */
    STEP_AB = (state->RHOX[IM + 1] - state->RHOX[IM]) / MU;
    STEP_BC = (state->RHOX[IM] - state->RHOX[IM - 1]) / MU;
    DER = (OPC_B - OPC_A) / STEP_AB;
    DER1 = (OPC_C - OPC_B) / STEP_BC;
    LAMBDA = (1.0 + STEP_BC / (STEP_AB + STEP_BC)) / 3.0;
    SPRIME_A = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
    CNTR_AB = OPC_B - STEP_AB / 2.0 * SPRIME_A;
    CNTR_BC = OPC_B + STEP_BC / 2.0 * SPRIME_A;
    DELTA = STEP_AB / 3.0 * (OPC_A + OPC_B + CNTR_AB);
    DELTA1 = STEP_BC / 3.0 * (OPC_B + OPC_C + CNTR_BC);
    /*
    Next we switch to optical depth and compute the contribution
    from the source function:
    */
    EPS = (DELTA < 100.0) ? exp(-DELTA) : 0.0; // Avoiding underflow
    /*
    Calculate parabolic coefficients for the source function
    Special provision is taken for the case of a very small
    DELTA resulting in precision loss when evaluating EPS and differences.
    Here we do Taylor expansion up to delta^3 for ALPHA, BETA and GAMMA.
    */
    if (DELTA < 1.e-3) // Use analytical expansion for small DELTA
    {
      ALPHA = DELTA / 3.0 - DELTA * DELTA / 12.0 + DELTA * DELTA * DELTA / 60.0;
      BETA = DELTA / 3.0 - DELTA * DELTA / 4.0 + DELTA * DELTA * DELTA / 10.0;
      GAMMA = DELTA / 3.0 - DELTA * DELTA / 6.0 + DELTA * DELTA * DELTA / 20.0;
    }
    else // or accurate calculations otherwise
    {
      ALPHA = (DELTA * DELTA - 2.0 * DELTA + 2.0 - 2.0 * EPS) / (DELTA * DELTA);
      BETA = (2.0 - (2.0 + 2.0 * DELTA + DELTA * DELTA) * EPS) / (DELTA * DELTA);
      GAMMA = (2.0 * DELTA - 4.0 + (2.0 * DELTA + 4.0) * EPS) / (DELTA * DELTA);
    }
    /*
    The last thing is the control parameter in optical path:
    */
    DER = (SRC_B - SRC_A) / DELTA;
    DER1 = (SRC_C - SRC_B) / DELTA1;
    LAMBDA = (1.0 + DELTA1 / (DELTA + DELTA1)) / 3.0;
    SPRIME_A = SPRIME_SAVE;
    SPRIME_B = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
    SPRIME_SAVE = SPRIME_B;
    if (IM == state->NRHOX - 2)
    {
      CNTR_AB = SRC_B - DELTA / 2.0 * SPRIME_B;
    }
    else
    {
      CNTR_AB = (SRC_A + DELTA * 0.5 * SPRIME_A + SRC_B - DELTA * 0.5 * SPRIME_B) * 0.5;
    }
    /*
    Finally, we are ready to compute the intensity in point B
    */
    B = ALPHA * SRC_B + BETA * SRC_A + GAMMA * CNTR_AB;
    INTENSITY = EPS * INTENSITY + B;
  }

  /* Continuum intensity at the surface */

  return INTENSITY * FLUX_SCALE;
}

void TBINTG_sph(int nrhox, double rhox[], double opacity[], double source[],
                double *RESULT, int grazing, GlobalState *state)
{
  /*
  RT solver
  AUTHOR: N.Piskunov
  UPDATES:  May  4, 2009 Re-written as quadratic DELO with Bezier splines
            Sep 26, 2010 Simplified the structure by moving the opacity and the
                         source function calculations to RKINTS which is the
                         caller of TBINTG. This version is for spherical models
  */
  double OPC_A, OPC_B, OPC_C, SRC_A, SRC_B, SRC_C, INTENSITY;
  double CNTR_AB, CNTR_BC, SPRIME_A, SPRIME_B;
  double STEP_AB, STEP_BC, DER, DER1, DELTA, DELTA1;
  double ALPHA, BETA, GAMMA, EPS, B, LAMBDA, SPRIME_SAVE, DBNU;
  int IM, IMU;

  /* Useful things for the Planck function */

  SRC_B = source[nrhox - 1]; // Source function
  SRC_C = source[nrhox - 2];
  OPC_B = opacity[nrhox - 1]; // Opacities
  OPC_C = opacity[nrhox - 2];
  DBNU = 2.0 * (SRC_B - SRC_C) / ((rhox[nrhox - 1] - rhox[nrhox - 2]) * (OPC_B + OPC_C));
  INTENSITY = (grazing) ? 0. : 0.5 * (SRC_B + SRC_C) + DBNU; // Line intensity at the boundary

  SPRIME_SAVE = 0.0; // Initialize S'

  for (IM = nrhox - 2; IM > 0; IM--) // Work your way from the deepest
  {                                  // layer to the surface
    SRC_A = SRC_B;                   // Shift source functions and opacities
    OPC_A = OPC_B;
    SRC_B = SRC_C;
    OPC_B = OPC_C;
    SRC_C = source[IM - 1]; // Downwind point
    OPC_C = opacity[IM - 1];
    /*
    Steps in monochromatic optical depth
    */
    STEP_AB = (rhox[IM + 1] - rhox[IM]);
    STEP_BC = (rhox[IM] - rhox[IM - 1]);
    DER = (OPC_B - OPC_A) / STEP_AB;
    DER1 = (OPC_C - OPC_B) / STEP_BC;
    LAMBDA = (1.0 + STEP_BC / (STEP_AB + STEP_BC)) / 3.0;
    SPRIME_A = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
    CNTR_AB = OPC_B - STEP_AB / 2.0 * SPRIME_A;
    CNTR_BC = OPC_B + STEP_BC / 2.0 * SPRIME_A;
    DELTA = STEP_AB / 3.0 * (OPC_A + OPC_B + CNTR_AB);
    DELTA1 = STEP_BC / 3.0 * (OPC_B + OPC_C + CNTR_BC);
    /*
    Next we switch to optical depth and compute the contribution
    from the source function:
    */
    EPS = (DELTA < 100.0) ? exp(-DELTA) : 0.0; // Avoiding underflow
    /*
    Calculate parabolic coefficients for the source function
    Special provision is taken for the case of a very small
    DELTA resulting in precision loss when evaluating EPS and differences.
    Here we do Taylor expansion up to delta^3 for ALPHA, BETA and GAMMA.
    */
    if (DELTA < 1.e-3) // Use analytical expansion for small DELTA
    {
      ALPHA = DELTA / 3.0 - DELTA * DELTA / 12.0 + DELTA * DELTA * DELTA / 60.0;
      BETA = DELTA / 3.0 - DELTA * DELTA / 4.0 + DELTA * DELTA * DELTA / 10.0;
      GAMMA = DELTA / 3.0 - DELTA * DELTA / 6.0 + DELTA * DELTA * DELTA / 20.0;
    }
    else // or accurate calculations otherwise
    {
      ALPHA = (DELTA * DELTA - 2.0 * DELTA + 2.0 - 2.0 * EPS) / (DELTA * DELTA);
      BETA = (2.0 - (2.0 + 2.0 * DELTA + DELTA * DELTA) * EPS) / (DELTA * DELTA);
      GAMMA = (2.0 * DELTA - 4.0 + (2.0 * DELTA + 4.0) * EPS) / (DELTA * DELTA);
    }
    /*
    The last thing is the control parameter in optical path:
    */
    DER = (SRC_B - SRC_A) / DELTA;
    DER1 = (SRC_C - SRC_B) / DELTA1;
    LAMBDA = (1.0 + DELTA1 / (DELTA + DELTA1)) / 3.0;
    SPRIME_A = SPRIME_SAVE;
    SPRIME_B = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
    SPRIME_SAVE = SPRIME_B;
    if (IM == nrhox - 2)
    {
      CNTR_AB = SRC_B - DELTA / 2.0 * SPRIME_B;
    }
    else
    {
      CNTR_AB = (SRC_A + DELTA * 0.5 * SPRIME_A + SRC_B - DELTA * 0.5 * SPRIME_B) * 0.5;
    }
    /*
    Finally, we are ready to compute the intensity in point B
    */
    B = ALPHA * SRC_B + BETA * SRC_A + GAMMA * CNTR_AB;
    INTENSITY = EPS * INTENSITY + B;
  }
  *RESULT = INTENSITY * FLUX_SCALE;
}

void TBINTG1(double rhox[], double opacity[], double source[], double *RESULT, GlobalState *state)
{
  /*
  RT solver
  AUTHOR: N.Piskunov
  UPDATES:  May  4, 2009 Re-written as quadratic DELO with Bezier splines
            Sep 26, 2010 Simplified the structure by moving the opacity and the
                         source function calculations to RKINTS which is the
                         caller of TBINTG
  */
  double OPC_A, OPC_B, OPC_C, SRC_A, SRC_B, SRC_C, INTENSITY;
  double CNTR_AB, CNTR_BC, SPRIME_A, SPRIME_B;
  double STEP_AB, STEP_BC, DER, DER1, DELTA, DELTA1;
  double ALPHA, BETA, GAMMA, EPS, B, LAMBDA, SPRIME_SAVE, DBNU;
  int IM;

  /* Useful things for the Planck function */

  SRC_B = source[state->NRHOX - 1]; // Source function
  SRC_C = source[state->NRHOX - 2];
  OPC_B = opacity[state->NRHOX - 1]; // Opacities
  OPC_C = opacity[state->NRHOX - 2];
  DBNU = 2.0 * (SRC_B - SRC_C) / ((rhox[state->NRHOX - 1] - rhox[state->NRHOX - 2]) * (OPC_B + OPC_C));
  INTENSITY = 0.5 * (SRC_B + SRC_C) + DBNU; // Line intensity at the bottom

  SPRIME_SAVE = 0.0; // Initialize S'

  for (IM = state->NRHOX - 2; IM > 0; IM--) // Work your way from the deepest
  {                                         // layer to the surface
    SRC_A = SRC_B;                          // Shift source functions and opacities
    OPC_A = OPC_B;
    SRC_B = SRC_C;
    OPC_B = OPC_C;
    SRC_C = source[IM - 1]; // Downwind point
    OPC_C = opacity[IM - 1];
    /*
    Steps in monochromatic optical depth
    */
    STEP_AB = (rhox[IM + 1] - rhox[IM]);
    STEP_BC = (rhox[IM] - rhox[IM - 1]);
    DER = (OPC_B - OPC_A) / STEP_AB;
    DER1 = (OPC_C - OPC_B) / STEP_BC;
    LAMBDA = (1.0 + STEP_BC / (STEP_AB + STEP_BC)) / 3.0;
    SPRIME_A = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
    CNTR_AB = OPC_B - STEP_AB / 2.0 * SPRIME_A;
    CNTR_BC = OPC_B + STEP_BC / 2.0 * SPRIME_A;
    DELTA = STEP_AB / 3.0 * (OPC_A + OPC_B + CNTR_AB);
    DELTA1 = STEP_BC / 3.0 * (OPC_B + OPC_C + CNTR_BC);
    /*
    Next we switch to optical depth and compute the contribution
    from the source function:
    */
    EPS = (DELTA < 100.0) ? exp(-DELTA) : 0.0; // Avoiding underflow
    /*
    Calculate parabolic coefficients for the source function
    Special provision is taken for the case of a very small
    DELTA resulting in precision loss when evaluating EPS and differences.
    Here we do Taylor expansion up to delta^3 for ALPHA, BETA and GAMMA.
    */
    if (DELTA < 1.e-3) // Use analytical expansion for small DELTA
    {
      ALPHA = DELTA / 3.0 - DELTA * DELTA / 12.0 + DELTA * DELTA * DELTA / 60.0;
      BETA = DELTA / 3.0 - DELTA * DELTA / 4.0 + DELTA * DELTA * DELTA / 10.0;
      GAMMA = DELTA / 3.0 - DELTA * DELTA / 6.0 + DELTA * DELTA * DELTA / 20.0;
    }
    else // or accurate calculations otherwise
    {
      ALPHA = (DELTA * DELTA - 2.0 * DELTA + 2.0 - 2.0 * EPS) / (DELTA * DELTA);
      BETA = (2.0 - (2.0 + 2.0 * DELTA + DELTA * DELTA) * EPS) / (DELTA * DELTA);
      GAMMA = (2.0 * DELTA - 4.0 + (2.0 * DELTA + 4.0) * EPS) / (DELTA * DELTA);
    }
    /*
    The last thing is the control parameter in optical path:
    */
    DER = (SRC_B - SRC_A) / DELTA;
    DER1 = (SRC_C - SRC_B) / DELTA1;
    LAMBDA = (1.0 + DELTA1 / (DELTA + DELTA1)) / 3.0;
    SPRIME_A = SPRIME_SAVE;
    SPRIME_B = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
    SPRIME_SAVE = SPRIME_B;
    if (IM == state->NRHOX - 2)
    {
      CNTR_AB = SRC_B - DELTA / 2.0 * SPRIME_B;
    }
    else
    {
      CNTR_AB = (SRC_A + DELTA * 0.5 * SPRIME_A + SRC_B - DELTA * 0.5 * SPRIME_B) * 0.5;
    }
    /*
    Finally, we are ready to compute the intensity in point B
    */
    B = ALPHA * SRC_B + BETA * SRC_A + GAMMA * CNTR_AB;
    INTENSITY = EPS * INTENSITY + B;
  }
  *RESULT = INTENSITY * FLUX_SCALE;
}

void TBINTG(int Nmu, double rhox[], double opacity[], double source[],
            double RESULT[], GlobalState *state)
{
  /*
  RT solver for plane parallel geometry
  AUTHOR: N.Piskunov
  UPDATES:  May  4, 2009 Re-written as quadratic DELO with Bezier splines
            Sep 26, 2010 Simplified the structure by moving the opacity and the
                         source function calculations to RKINTS which is the
                         caller of TBINTG
            Feb 14, 2011 Move the mu loop inside TBINTG to speed up things
  */
  double OPC_A, OPC_B, OPC_C, SRC_A, SRC_B, SRC_C;
  double CNTR_AB, CNTR_BC, SPRIME_A, SPRIME_B;
  double STEP_AB, STEP_BC, DER, DER1, DELTA, DELTA1;
  double ALPHA, BETA, GAMMA, EPS, B, LAMBDA, DBNU;
  double SPRIME_SAVE[MUSIZE], INTENSITY[MUSIZE];
  int IM, imu;

  /* Useful things for the Planck function */

  SRC_B = source[state->NRHOX - 1]; // Source function
  SRC_C = source[state->NRHOX - 2];
  OPC_B = opacity[state->NRHOX - 1]; // Opacities
  OPC_C = opacity[state->NRHOX - 2];
  for (imu = 0; imu < Nmu; imu++)
  {
    DBNU = 2.0 * (SRC_B - SRC_C) / ((rhox[imu * state->NRHOX + state->NRHOX - 1] - rhox[imu * state->NRHOX + state->NRHOX - 2]) * (OPC_B + OPC_C));
    INTENSITY[imu] = 0.5 * (SRC_B + SRC_C) + DBNU; // Line intensity at the bottom
    SPRIME_SAVE[imu] = 0.0;                        // Initialize S'
  }

  for (IM = state->NRHOX - 2; IM > 0; IM--) // Work your way from the deepest
  {                                         // layer to the surface
    SRC_A = SRC_B;                          // Shift source functions and opacities
    OPC_A = OPC_B;
    SRC_B = SRC_C;
    OPC_B = OPC_C;
    SRC_C = source[IM - 1]; // Downwind point
    OPC_C = opacity[IM - 1];
    /*
    Steps in monochromatic optical depth
    */
    for (imu = 0; imu < Nmu; imu++)
    {
      STEP_AB = (rhox[imu * state->NRHOX + IM + 1] - rhox[imu * state->NRHOX + IM]);
      STEP_BC = (rhox[imu * state->NRHOX + IM] - rhox[imu * state->NRHOX + IM - 1]);
      DER = (OPC_B - OPC_A) / STEP_AB;
      DER1 = (OPC_C - OPC_B) / STEP_BC;
      LAMBDA = (1.0 + STEP_BC / (STEP_AB + STEP_BC)) / 3.0;
      SPRIME_A = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
      CNTR_AB = OPC_B - STEP_AB / 2.0 * SPRIME_A;
      CNTR_BC = OPC_B + STEP_BC / 2.0 * SPRIME_A;
      DELTA = STEP_AB / 3.0 * (OPC_A + OPC_B + CNTR_AB);
      DELTA1 = STEP_BC / 3.0 * (OPC_B + OPC_C + CNTR_BC);
      /*
      Next we switch to optical depth and compute the contribution
      from the source function:
      */
      EPS = (DELTA < 100.0) ? exp(-DELTA) : 0.0; // Avoiding underflow
      /*
      Calculate parabolic coefficients for the source function
      Special provision is taken for the case of a very small
      DELTA resulting in precision loss when evaluating EPS and differences.
      Here we do Taylor expansion up to delta^3 for ALPHA, BETA and GAMMA.
      */
      if (DELTA < 1.e-3) // Use analytical expansion for small DELTA
      {
        ALPHA = DELTA / 3.0 - DELTA * DELTA / 12.0 + DELTA * DELTA * DELTA / 60.0;
        BETA = DELTA / 3.0 - DELTA * DELTA / 4.0 + DELTA * DELTA * DELTA / 10.0;
        GAMMA = DELTA / 3.0 - DELTA * DELTA / 6.0 + DELTA * DELTA * DELTA / 20.0;
      }
      else // or accurate calculations otherwise
      {
        ALPHA = (DELTA * DELTA - 2.0 * DELTA + 2.0 - 2.0 * EPS) / (DELTA * DELTA);
        BETA = (2.0 - (2.0 + 2.0 * DELTA + DELTA * DELTA) * EPS) / (DELTA * DELTA);
        GAMMA = (2.0 * DELTA - 4.0 + (2.0 * DELTA + 4.0) * EPS) / (DELTA * DELTA);
      }
      /*
      The last thing is the control parameter in optical path:
      */
      DER = (SRC_B - SRC_A) / DELTA;
      DER1 = (SRC_C - SRC_B) / DELTA1;
      LAMBDA = (1.0 + DELTA1 / (DELTA + DELTA1)) / 3.0;
      SPRIME_A = SPRIME_SAVE[imu];
      SPRIME_B = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
      SPRIME_SAVE[imu] = SPRIME_B;
      if (IM == state->NRHOX - 2)
      {
        CNTR_AB = SRC_B - DELTA / 2.0 * SPRIME_B;
      }
      else
      {
        CNTR_AB = (SRC_A + DELTA * 0.5 * SPRIME_A + SRC_B - DELTA * 0.5 * SPRIME_B) * 0.5;
      }
      /*
      Finally, we are ready to compute the intensity in point B
      */
      B = ALPHA * SRC_B + BETA * SRC_A + GAMMA * CNTR_AB;
      INTENSITY[imu] = EPS * INTENSITY[imu] + B;
    }
  }
  for (imu = 0; imu < Nmu; imu++)
    RESULT[imu] = INTENSITY[imu] * FLUX_SCALE;
}

void CENTERINTG(double *MUs, int NMU, int LINE, double *contop, double *RESULT, GlobalState *state)
{
  /*
  Quadratic DELO with Bezier spline RT solver
  AUTHOR: N.Piskunov
  LAST UPDATE: May 4, 2009
  */
  double OPC_A, OPC_B, OPC_C, SRC_A, SRC_B, SRC_C, INTENSITY;
  double CNTR_AB, CNTR_BC, SPRIME_A, SPRIME_B;
  double STEP_AB, STEP_BC, DER, DER1, DELTA, DELTA1;
  double ALPHA, BETA, GAMMA, EPS, B, LAMBDA, SPRIME_SAVE, DBNU;
  double CONWL5, HNUK, MU, XK[MOSIZE];
  int IM, IMU;

  /* Useful things for the Planck function */

  CONWL5 = exp(50.7649141 - 5. * log(state->WLCENT[LINE]));
  HNUK = 1.43868e8 / state->WLCENT[LINE];

  OPMTRX1(LINE, XK, state);

  if (state->MOTYPE)
    for (IM = 0; IM < state->NRHOX; IM++)
      XK[IM] = XK[IM] + contop[IM];
  else
    for (IM = 0; IM < state->NRHOX; IM++)
      XK[IM] = XK[IM] + contop[IM] / state->COPSTD[IM];

  for (IMU = 0; IMU < NMU; IMU++)
  {
    MU = MUs[IMU];
    SRC_B = CONWL5 / (exp(HNUK / state->T[state->NRHOX - 1]) - 1.); // Source function
    SRC_C = CONWL5 / (exp(HNUK / state->T[state->NRHOX - 2]) - 1.);
    OPC_B = XK[state->NRHOX - 1]; // Opacities
    OPC_C = XK[state->NRHOX - 2];
    DBNU = 2.0 * (SRC_B - SRC_C) / ((state->RHOX[state->NRHOX - 1] - state->RHOX[state->NRHOX - 2]) * (OPC_B + OPC_C)) * MU;
    INTENSITY = 0.5 * (SRC_B + SRC_C) + DBNU; // Intensity at the bottom

    SPRIME_SAVE = 0.0; // Initialize S'

    for (IM = state->NRHOX - 2; IM > 0; IM--) // Work your way from the deepest
    {                                         // layer to the surface
      SRC_A = SRC_B;                          // Shift source functions and opacities
      OPC_A = OPC_B;
      SRC_B = SRC_C;
      OPC_B = OPC_C;
      SRC_C = CONWL5 / (exp(HNUK / state->T[IM - 1]) - 1.); // Downwind point
      OPC_C = XK[IM - 1];
      /*
      Steps in monochromatic optical depth
      */
      STEP_AB = (state->RHOX[IM + 1] - state->RHOX[IM]) / MU;
      STEP_BC = (state->RHOX[IM] - state->RHOX[IM - 1]) / MU;
      DER = (OPC_B - OPC_A) / STEP_AB;
      DER1 = (OPC_C - OPC_B) / STEP_BC;
      LAMBDA = (1.0 + STEP_BC / (STEP_AB + STEP_BC)) / 3.0;
      SPRIME_A = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
      CNTR_AB = OPC_B - STEP_AB / 2.0 * SPRIME_A;
      CNTR_BC = OPC_B + STEP_BC / 2.0 * SPRIME_A;
      DELTA = STEP_AB / 3.0 * (OPC_A + OPC_B + CNTR_AB);
      DELTA1 = STEP_BC / 3.0 * (OPC_B + OPC_C + CNTR_BC);
      /*
      Next we switch to optical depth and compute the contribution
      from the source function:
      */
      EPS = (DELTA < 100.0) ? exp(-DELTA) : 0.0; // Avoiding underflow
      /*
      Calculate parabolic coefficients for the source function
      Special provision is taken for the case of a very small
      DELTA resulting in precision loss when evaluating EPS and differences.
      Here we do Taylor expansion up to delta^3 for ALPHA, BETA and GAMMA.
      */
      if (DELTA < 1.e-3) // Use analytical expansion for small DELTA
      {
        ALPHA = DELTA / 3.0 - DELTA * DELTA / 12.0 + DELTA * DELTA * DELTA / 60.0;
        BETA = DELTA / 3.0 - DELTA * DELTA / 4.0 + DELTA * DELTA * DELTA / 10.0;
        GAMMA = DELTA / 3.0 - DELTA * DELTA / 6.0 + DELTA * DELTA * DELTA / 20.0;
      }
      else // or accurate calculations otherwise
      {
        ALPHA = (DELTA * DELTA - 2.0 * DELTA + 2.0 - 2.0 * EPS) / (DELTA * DELTA);
        BETA = (2.0 - (2.0 + 2.0 * DELTA + DELTA * DELTA) * EPS) / (DELTA * DELTA);
        GAMMA = (2.0 * DELTA - 4.0 + (2.0 * DELTA + 4.0) * EPS) / (DELTA * DELTA);
      }
      /*
      The last thing is the control parameter in optical path:
      */
      DER = (SRC_B - SRC_A) / DELTA;
      DER1 = (SRC_C - SRC_B) / DELTA1;
      LAMBDA = (1.0 + DELTA1 / (DELTA + DELTA1)) / 3.0;
      SPRIME_A = SPRIME_SAVE;
      SPRIME_B = (DER * DER1 > 0.0) ? DER / (LAMBDA * DER1 + (1.0 - LAMBDA) * DER) * DER1 : 0.0;
      SPRIME_SAVE = SPRIME_B;
      if (IM == state->NRHOX - 2)
      {
        CNTR_AB = SRC_B - DELTA / 2.0 * SPRIME_B;
      }
      else
      {
        CNTR_AB = (SRC_A + DELTA * 0.5 * SPRIME_A + SRC_B - DELTA * 0.5 * SPRIME_B) * 0.5;
      }
      /*
      Finally, we are ready to compute the intensity in point B
      */
      B = ALPHA * SRC_B + BETA * SRC_A + GAMMA * CNTR_AB;
      INTENSITY = EPS * INTENSITY + B;
    }
    RESULT[IMU] = INTENSITY * FLUX_SCALE;
  }
}

#undef FLUX_SCALE

extern "C" char const *SME_DLL GetLineOpacity(int n, void *arg[], GlobalState *state) /* Returns specific line opacity */
{
  int MOTYPE_orig;
  short i, j, nrhox;
  double *a1, *a2, *a3, *a4, *a5, WAVE, *XK, *XC, *SRC, *SRC_CONT;

  if (n < 3)
  {
    strcpy(state->result, "Not enough arguments");
    return state->result;
  }
  WAVE = *(double *)arg[0]; /* Wavelength */
  i = *(short *)arg[1];     /* Length of IDL opacity array */
  nrhox = min(state->NRHOX, i);
  a1 = (double *)arg[2];       /* Line opacity */
  a2 = (double *)arg[3];       /* Continuum opacity including scatter */
  a3 = (double *)arg[4];       /* Scatter */
  a4 = (double *)arg[5];       /* Total source function */
  a5 = (double *)arg[6];       /* Continuum source function */
  MOTYPE_orig = state->MOTYPE; /* Save state->MOTYPE */
  state->MOTYPE = -1;          /* Set state->MOTYPE to return only line opacity */

  /* Allocate temporary arrays */

  CALLOC(XK, state->NRHOX, double);
  CALLOC(XC, state->NRHOX, double);
  CALLOC(SRC, state->NRHOX, double);
  CALLOC(SRC_CONT, state->NRHOX, double);

  AutoIonization(state);
  OPMTRX(WAVE, XK, XC, SRC, SRC_CONT, 0, state->NLINES - 1, state);

  for (i = 0; i < nrhox; i++)
  {
    a1[i] = XK[i];
    a2[i] = XC[i];
    a3[i] = state->SIGH[i] + state->SIGEL[i] + state->SIGH2[i] + state->SIGHE[i];
    a4[i] = SRC[i];
    a5[i] = SRC_CONT[i];
  }

  FREE(XK);
  FREE(XC);
  FREE(SRC);
  FREE(SRC_CONT);

  state->MOTYPE = MOTYPE_orig;
  return &OK_response;
}

#define Z 4.9946686e-21
#define C4PI CLIGHT * 4. * PI
#define PI4 4. * PI
#define K 1.380658e-23
#define M0 1.660540e-27
#define A0 5.29177249e-11

void LINEOPAC(int LINE, GlobalState *state)
{
  /*
   This function computes central line opacity without the
   profile and the line width. The exception is the Hydrogen
   lines that are treated inside OPMTRX. Line opacity is per gram
   of matter in cm^2/g.

   Author: N.Piskunov

                    pi*e^2
   Line opacity is: ------ * gf * N_absorb * state->STIM
                     m*c
  
   The Hydrogen line profiles are computed externally by Kurucz
   approximation (HLINOP) or by interpolation in Stehle's tables (HTABLE)
   and are area normalized!

   Therefore the normalization factor Z=PI*e^2/(m*c) with speed
   of light in cm/s. The net state->result is that Z is in cm^2/s !!!

   Other constants: K  - Boltzmann's constant J/K,
                    M0 - unit atomic mass kg (Carbon 12 scale),
                    A0 - Bohr radius m

   Author: N.Piskunov

   C++ Version: October 26, 1994
   UPDATES: May 26, 1999
                 Consistent interface to HLINOP (same as in SYNTH)
            Jan 20, 2010
                Temperature dependent van der Waals if ALPHA and SIGMA are
                available and reduced mass of perturbers by Paul Barklem
            Aug 26, 2010
                Added calculations of continuum opacity and the source
                function
  */

  double HNUXXX, DDWL, WAVE;
  double OPCONB, OPCONR, OPCON, DNDOPL, DLDOPL, A, UAV, V4, W4, VOIGT,
      XXRHO, XNELEC, XNATOM, XTK, XSTIM, VH, H1FRC, HE1FRC, H2molFRC,
      GVWPRT, TEMP3, TEMP6, ALINE, WLC, FR, EFRACT, SHFT, TEMPER,
      DOPL, GQST, GVW, CW, GAMTOT, Vmicro, VTURB2, ALINE1,
      SIGMA, ALPHA, GX, X, GAMMAF, VBAR, CONWL5, HNUK;
  double opac[MOSIZE];
  short ion, ITAU;
  int i_cont;

  WAVE = state->WLCENT[LINE];
  CONTOP(WAVE, opac, state);
  state->ALMAX[LINE] = 0.;
  for (ITAU = 0; ITAU < state->NRHOX; ITAU++)
  {
    TEMPER = state->T[ITAU];
    HNUXXX = CLIGHT * 6.6256e-27 / WAVE;
    XXRHO = state->RHO[ITAU];  /* Density                 */
    XNELEC = state->XNE[ITAU]; /* Electron number density */
    XNATOM = state->XNA[ITAU]; /* Atom number density     */
    Vmicro = state->VTURB[ITAU];
    OPCON = opac[ITAU];

    /* Fractions of H I and He I */

    H1FRC = state->H1FRACT[ITAU];
    HE1FRC = state->HE1FRACT[ITAU];
    H2molFRC = state->H2molFRACT[ITAU];

    /* Some other useful things */

    XTK = TEMPER * 1.38054e-16;
    XSTIM = 1. - exp(-HNUXXX / XTK);
    TEMP6 = pow(TEMPER / 10000., 1. / 6.) * XNELEC;
    TEMP3 = pow(TEMPER / 10000., 0.3) * (H1FRC + 0.413 * HE1FRC +
                                         (state->flagH2broad ? 0.876 * H2molFRC : 0.));

    /* state->VTURB is in km/s, 1.E13 converts C to km/s, so VTURB2 is dimensionless */
    VTURB2 = 1.e26 / CLIGHT / CLIGHT * Vmicro * Vmicro;

    /* Loop through spectral lines */

    if (state->AUTOION[LINE] && (state->GAMVW[LINE] <= 0.0 || state->GAMQST[LINE] <= 0.0))
    {
      state->AVOIGT[ITAU][LINE] = 1.;
      state->VVOIGT[ITAU][LINE] = 1.;
      state->LINEOP[ITAU][LINE] = 0.;
      state->MARK[LINE] = 2;
    }
    else
    {
      WLC = state->WLCENT[LINE];
      ion = state->ION[LINE]; /* ion==1 for neutrals */

      /* The fraction number of absorbing atoms */

      FR = state->FRACT[ITAU][state->SPINDEX[LINE]];
      EFRACT = FR * exp(-state->EXCIT[LINE] / (8.6171e-5 * TEMPER));

      /* Wavelength independent things for a given line */

      state->YABUND[LINE] = Z * state->GF[LINE];
      state->XMASS[LINE] = 1.66355e24 / CLIGHT / CLIGHT / state->MOLWEIGHT[state->SPINDEX[LINE]];
      state->EXCUP[LINE] = state->EXCIT[LINE] + 1. / (WLC * 8065.544e-8);
      if (!state->AUTOION[LINE] && (state->GAMVW[LINE] == 0. || state->GAMQST[LINE] == 0.))
      {
        state->ENU4[LINE] = (ion * 13.598 * ion / (state->POTION[state->SPINDEX[LINE]] - state->EXCUP[LINE]));
        state->ENU4[LINE] = state->ENU4[LINE] * state->ENU4[LINE];
        state->ENL4[LINE] = (ion * 13.598 * ion / (state->POTION[state->SPINDEX[LINE]] - state->EXCIT[LINE]));
        state->ENL4[LINE] = state->ENL4[LINE] * state->ENL4[LINE];
      }

      /* Radiative damping */

      state->GAMRAD[LINE] = (state->GAMRAD[LINE] > 0.0) ? state->GAMRAD[LINE] : 0.222e16 / (WLC * WLC);

      /* Identify Helium lines included in Dimitrijevic & Sahal-Brechot table;
      Stark damping for those will be computed in subroutine GAMHE */

      state->IDHEL[LINE] = -1;
      if (!strncmp(state->spname + 8 * LINE, "He ", 3) && !state->MARK[LINE])
      {
        switch ((int)floor(WLC))
        {
        case 3819:
          state->IDHEL[LINE] = 0;
          break;
        case 3867:
          state->IDHEL[LINE] = 1;
          break;
        case 3871:
          state->IDHEL[LINE] = 2;
          break;
        case 3888:
          state->IDHEL[LINE] = 3;
          break;
        case 3926:
          state->IDHEL[LINE] = 4;
          break;
        case 3964:
          state->IDHEL[LINE] = 5;
          break;
        case 4009:
          state->IDHEL[LINE] = 6;
          break;
        case 4120:
        case 4121:
          state->IDHEL[LINE] = 7;
          break;
        case 4143:
          state->IDHEL[LINE] = 8;
          break;
        case 4168:
        case 4169:
          state->IDHEL[LINE] = 9;
          break;
        case 4437:
          state->IDHEL[LINE] = 10;
          break;
        case 4471:
          state->IDHEL[LINE] = 11;
          break;
        case 4713:
          state->IDHEL[LINE] = 12;
          break;
        case 4921:
        case 4922:
          state->IDHEL[LINE] = 13;
          break;
        case 5015:
        case 5016:
          state->IDHEL[LINE] = 14;
          break;
        case 5047:
          state->IDHEL[LINE] = 15;
          break;
        case 5875:
          state->IDHEL[LINE] = 16;
          break;
        case 6678:
          state->IDHEL[LINE] = 17;
          break;
        case 4026:
          state->IDHEL[LINE] = 18;
          break;
        case 4387:
        case 4388:
          state->IDHEL[LINE] = 19;
          break;
        default:
          break;
        }
      }

      /*  Doppler broadening: DOPL is in fact delta_lambda/lambda
      DLDOPL is delta_lambda in Angstroems
      DNDOPL is delta_nu in Hz. */

      DOPL = sqrt(TEMPER * state->XMASS[LINE] + VTURB2);
      DLDOPL = WAVE * DOPL;
      state->VVOIGT[ITAU][LINE] = 1. / DLDOPL;
      DNDOPL = DOPL / WAVE;

      if (!strncmp(state->spname + 8 * LINE, "H ", 2)) // This is a hydrogen line
      {
        double HNORM;

        HNORM = SQRTPI * EFRACT * state->YABUND[LINE] * XSTIM / XXRHO;
        state->VVOIGT[ITAU][LINE] = DOPL;
        state->LINEOP[ITAU][LINE] = HNORM;
        state->ALMAX[LINE] = 1.e6;
      }
      else // Non-hydrogen line
      {

        /*  Qudratic Stark effect (if the constant is available, compute according
        to D.Gray, otherwise - follow C.Cowley). For Helium - Dimitrijevich
        tables are used. */

        if (state->IDHEL[LINE] < 0) /* If not Helium */
        {
          if (state->GAMQST[LINE] > 0.0 || state->AUTOION[LINE])
            GQST = state->GAMQST[LINE] * TEMP6;
          else
          {
            GQST = (ion - 1) ? 5.42e-7 * state->ENU4[LINE] * XNELEC / ((ion + 1) * (ion + 1)) : 2.26e-7 * state->ENU4[LINE] * XNELEC;
          }
        }
        else /* Compute Stark broadenning for Helium separately */
        {
          GAMHE(state->IDHEL[LINE], TEMPER, XNELEC, XNATOM, GQST, SHFT, state);
        }

        /*  Van der Waals damping parameter */
        if (state->ANSTEE[LINE])
        {
          /*
          This van der Waals part is written by Paul Barklem
          Compute the broadening by hydrogen from cross-section data which is in m^2
          Unpack the temperature dependent van der Waals parameters:
          integer part is SIGMA and decimal part is ALPHA.
          */
          SIGMA = ((int)state->GAMVW[LINE]) * A0 * A0;
          ALPHA = state->GAMVW[LINE] - (int)state->GAMVW[LINE];

          //  Compute the Gamma function of X, this function is valid over the range 1<X<2

          X = 2.e0 - ALPHA * 0.5e0;
          GX = X - 1.e0;
          GAMMAF = 1.e0 + (-0.5748646e0 + (0.9512363e0 + (-0.6998588e0 + (0.4245549e0 - 0.1010678e0 * GX) * GX) * GX) * GX) * GX;

          //  Compute the halfwidth per unit perturber density for vbar

          GVW = pow(4. / PI, ALPHA * 0.5) * GAMMAF * 1.e4 * SIGMA;
          VBAR = sqrt(8. * K * TEMPER / PI / M0 * (1. / 1.008 + 1. / state->MOLWEIGHT[state->SPINDEX[LINE]]));
          GVW = GVW * pow(VBAR / 1.e4, 1. - ALPHA);

          //  Fullwidth given H1FRC perturbers per cm^3 with approximate HeI and
          //  molecular H2 contributions. The factor of 2 in the end comes from
          //  converting to the full width.

          GVW = GVW * (H1FRC + 0.413 * HE1FRC + (state->flagH2broad ? 0.876 * H2molFRC : 0.)) * 1.e6 * 2.;
        }
        else if ((!state->ANSTEE[LINE] && state->GAMVW[LINE] > 0.0) || state->AUTOION[LINE])
        { // Input was log line width per unit density (rad/s cm^3)
          GVW = state->GAMVW[LINE] * TEMP3 * state->VW_scale;
        }
        else
        { // Input was zero and so we use Unsold theory
          CW = 1.61e-33 * (state->ENU4[LINE] - state->ENL4[LINE]) / (ion * ion);
          state->GAMVW[LINE] = 78654.213 * pow(CW, 0.4);
          GVW = state->GAMVW[LINE] * TEMP3 * state->VW_scale;
        }

        /*  Total broadening and VOIGT function parameters */

        GAMTOT = state->GAMRAD[LINE] + GQST + GVW;
        state->AVOIGT[ITAU][LINE] = GAMTOT / (DNDOPL * C4PI);
        A = state->AVOIGT[ITAU][LINE];

        /*  VOIGT function calculation: Humlicek, J. 1982, J.Q.S.R.state->T. 27, 437
        stripted for the case of line center (V==0) */

        UAV = A * A;
        if (A >= 15.)
          W4 = A * 0.5641896 / (0.5 + UAV);
        else if (A >= 5.5)
          W4 = A * (1.410474 + UAV * 0.5641896) / (0.75 + UAV * (3. + UAV));
        else if (A >= -0.176)
          W4 = (16.4955 + A * (20.20933 + A * (11.96482 + A * (3.778987 + A * 0.5642236)))) /
               (16.4955 + A * (38.82363 + A * (39.27121 + A * (21.69274 + A * (6.699398 + A)))));
        else
        {
          W4 = A * (36183.31 - UAV * (3321.9905 - UAV * (1540.787 - UAV * (219.0313 - UAV * (35.76683 - UAV * (1.320522 - UAV * .56419))))));
          V4 = (32066.6 - UAV * (24322.84 - UAV * (9022.228 - UAV * (2186.181 - UAV * (364.2191 - UAV * (61.57037 - UAV * (1.841439 - UAV)))))));
          W4 = exp(UAV) - W4 / V4;
        }
        VOIGT = W4;

        /*  Line absorption without the VOIGT function */

        state->LINEOP[ITAU][LINE] = EFRACT * state->YABUND[LINE] * XSTIM / (XXRHO * DNDOPL);
        if (state->LINEOP[ITAU][LINE] * VOIGT / OPCON > state->ALMAX[LINE])
          state->ALMAX[LINE] = state->LINEOP[ITAU][LINE] * VOIGT / OPCON;
      }
    }
  }
}

void OPMTRX(double WAVE, double *XK, double *XC, double *source_line,
            double *source_cont, int LINE_START, int LINE_FINISH, GlobalState *state)
{
  /*
   THIS FUNCTION CALCULATES THE OPACITY OR OPACITY RATIO (OPACWL/OPACSTD)
   PER GRAMM OF STELLAR MATER (CM**2/GM) PER ANGSTROEM AT DEPTH #IM
   OF THE STANDARD MODEL DEPTH SCALE. WAVELENGTH IS TAKEN EITHER FROM
   WAVE (ICODE=0) OR FROM EDGES OF SPECTRAL INTERVAL (ICODE=1,2).

   Author: N.Piskunov

                    pi*e^2
   Line opacity is: ------ * gf * N_absorb * state->STIM * f(wl-wl0)
                     m*c
 
   where the line profile f(wl) is assumed to be nomalized so that:
 
   \integ f(wl-wl0) d wl = 1
 
   This is true for Voigt, Hydrogen and (I hope) Fano profiles.
                                                      1
   E.g., in case of Voigt profile f(wl-wl0)= -------------------- * H(a,v)
                                             sqrt(pi)*del_nu_Dopp
   where del_Dopp = DNDOPL  is in Hz,

   where H(a,v) is the Voigt function with normalization:
   \integ H(a,v) d v = sqrt(pi)
 
   Two Hydrogen line profiles are computed externally by Kurucz
   approximation (HLINOP) or by interpolation in Stehle's tables (HTABLE)
   and are area normalized!

   Therefore the normalization factor Z=PI*e^2/(m*c) with speed
   of light in cm/s. The net state->result is that Z is in cm^2/s !!!

   Other constants: K  - Boltzmann's constant J/K,
                    M0 - unit atomic mass kg (Carbon 12 scale),
                    A0 - Bohr radius m

   Author: N.Piskunov

   C++ Version: October 26, 1994
   UPDATES: May 26, 1999
                 Consistent interface to HLINOP (same as in SYNTH)
            Jan 20, 2010
                Temperature dependent van der Waals if ALPHA and SIGMA are
                available and reduced mass of perturbers by Paul Barklem
            Aug 26, 2010
                Added calculations of continuum opacity and the source
                function
  */

  double HNUXXX, DDWL;
  double OPCONB, OPCONR, OPCON, DNDOPL, DLDOPL, A, V,
      XNELEC, XNATOM, H1FRC, HE1FRC,
      ALINE, WLC, GQST, SHFT, VOIGT, TEMPER,
      DOPL, ALINE1, CONWL5, HNUK, EHNUKT, XNLTE, SRC_cont, SRC_line;
  double opcon[MOSIZE];
  short ion, ITAU;
  int i_cont;
  int LINE;

  CONWL5 = exp(50.7649141 - 5. * log(WAVE));
  HNUK = 1.43868e8 / WAVE;
  for (LINE = LINE_START; LINE <= LINE_FINISH; LINE++)
    state->ALMAX[LINE] = 0.;

  CONTOP(WAVE, opcon, state);
  for (ITAU = 0; ITAU < state->NRHOX; ITAU++)
  {
    TEMPER = state->T[ITAU];
    OPCON = opcon[ITAU];
    XNELEC = state->XNE[ITAU]; /* Electron number density */
    XNATOM = state->XNA[ITAU]; /* Atom number density     */

    EHNUKT = exp(HNUK / TEMPER);
    if (state->initNLTE)
    {
      SRC_cont = CONWL5 / (EHNUKT - 1.); // LTE source function used for continuum
      source_cont[ITAU] = SRC_cont;
      source_line[ITAU] = 0.;
    }
    else
    {
      source_cont[ITAU] = CONWL5 / (EHNUKT - 1.);
      source_line[ITAU] = source_cont[ITAU];
    }

    /* Loop through spectral lines */

    ALINE = 0.;
    for (LINE = LINE_START; LINE <= LINE_FINISH; LINE++)
    {
      if (state->MARK[LINE] || WAVE <= state->Wlim_left[LINE] || WAVE >= state->Wlim_right[LINE])
        continue;
      if (state->AUTOION[LINE] && (state->GAMVW[LINE] <= 0.0 || state->GAMQST[LINE] <= 0.0))
        continue;
      WLC = state->WLCENT[LINE];

      if (state->initNLTE) // NLTE correction
      {
        XNLTE = state->BNLTE_low[LINE][ITAU] / (EHNUKT - 1.) *
                (EHNUKT - state->BNLTE_upp[LINE][ITAU] / state->BNLTE_low[LINE][ITAU]);
        SRC_line = CONWL5 / // NLTE source function for line
                   (state->BNLTE_low[LINE][ITAU] / state->BNLTE_upp[LINE][ITAU] * EHNUKT - 1.);
      }

      if (!strncmp(state->spname + 8 * LINE, "H ", 2)) // This is a hydrogen line
      {
        int NBLO, NBUP;
        double HNORM;
        float temper, xnelec, h1frc, he1frc, dopl, aline1, aline2;
        double wave, wlcent;

        NBLO = (int)(state->GAMQST[LINE] + 0.1);
        NBUP = (int)(state->GAMVW[LINE] + 0.1);

        temper = TEMPER;
        xnelec = state->XNE[ITAU];
        h1frc = state->H1FRACT[ITAU];
        he1frc = state->HE1FRACT[ITAU];
        wave = WAVE;
        wlcent = state->WLCENT[LINE];
        dopl = state->VVOIGT[ITAU][LINE];
        hlinprof_(wave, wlcent, temper, xnelec, NBLO, NBUP,
                  h1frc, he1frc, dopl, aline1, state->PATH, &state->PATHLEN, &state->change_byte_order);
        ALINE1 = aline1 * state->LINEOP[ITAU][LINE] * wave * wave;
        if (state->initNLTE)
        {
          ALINE1 *= XNLTE; // NLTE correction to the line opacity
          source_line[ITAU] += ALINE1 * SRC_line;
        }
        state->ALMAX[LINE] = ALINE1 / OPCON;
      }
      else // Non-hydrogen line
      {
        double TR, TI, UR, UI, SAV, XX, YY, X1, Y1, X2, Y2, UU, VV;

        if (state->IDHEL[LINE] > 0)
        {
          GAMHE(state->IDHEL[LINE], TEMPER, XNELEC, state->FRACT[ITAU][1], GQST, SHFT, state);
          WLC = WLC + SHFT;
        }

        A = state->AVOIGT[ITAU][LINE];
        V = (WAVE - WLC) * state->VVOIGT[ITAU][LINE];

        /*  VOIGT function calculation: Humlicek, J. 1982, J.Q.S.R.state->T. 27, 437 */

        TR = A;
        TI = -V;
        UR = A * A - V * V;
        UI = -2 * A * V;
        SAV = fabs(V) + A;
        if (SAV >= 15.)
        {
          UR = UR + 0.5;
          XX = max(A * A, V * V);
          TR = TR / XX;
          TI = TI / XX;
          UR = UR / XX;
          UI = UI / XX;
          VOIGT = 0.5641896 * (TR * UR + TI * UI) / (UR * UR + UI * UI);
        }
        else if (SAV >= 5.5)
        {
          X1 = UR * 0.5641896 + 1.410474;
          Y1 = UI * 0.5641896;
          XX = X1 * TR - Y1 * TI;
          YY = X1 * TI + Y1 * TR;
          X1 = UR + 3.;
          Y1 = UI;
          UU = X1 * UR - Y1 * UI + 0.75;
          VV = X1 * UI + Y1 * UR;
          VOIGT = (XX * UU + YY * VV) / (UU * UU + VV * VV);
        }
        else if (A >= 0.195 * fabs(V) - 0.176)
        {
          X1 = 3.778987 + TR * 0.5642236;
          Y1 = TI * 0.5642236;
          X2 = X1 * TR - Y1 * TI + 11.96482;
          Y2 = X1 * TI + Y1 * TR;
          X1 = X2 * TR - Y2 * TI + 20.20933;
          Y1 = X2 * TI + Y2 * TR;
          XX = X1 * TR - Y1 * TI + 16.4955;
          YY = X1 * TI + Y1 * TR;
          X1 = TR + 6.699398;
          Y1 = TI;
          X2 = X1 * TR - Y1 * TI + 21.69274;
          Y2 = X1 * TI + Y1 * TR;
          X1 = X2 * TR - Y2 * TI + 39.27121;
          Y1 = X2 * TI + Y2 * TR;
          X2 = X1 * TR - Y1 * TI + 38.82363;
          Y2 = X1 * TI + Y1 * TR;
          UU = X2 * TR - Y2 * TI + 16.4955;
          VV = X2 * TI + Y2 * TR;
          VOIGT = (XX * UU + YY * VV) / (UU * UU + VV * VV);
        }
        else
        {
          X1 = 1.320522 - UR * 0.56419;
          Y1 = -UI * 0.56419;
          X2 = 35.76683 - (X1 * UR - Y1 * UI);
          Y2 = -(X1 * UI + Y1 * UR);
          X1 = 219.0313 - (X2 * UR - Y2 * UI);
          Y1 = -(X2 * UI + Y2 * UR);
          X2 = 1540.787 - (X1 * UR - Y1 * UI);
          Y2 = -(X1 * UI + Y1 * UR);
          X1 = 3321.9905 - (X2 * UR - Y2 * UI);
          Y1 = -(X2 * UI + Y2 * UR);
          X2 = 36183.31 - (X1 * UR - Y1 * UI);
          Y2 = -(X1 * UI + Y1 * UR);
          XX = X2 * TR - Y2 * TI;
          YY = X2 * TI + Y2 * TR;
          X1 = 1.841439 - UR;
          Y1 = -UI;
          X2 = 61.57037 - (X1 * UR - Y1 * UI);
          Y2 = -(X1 * UI + Y1 * UR);
          X1 = 364.2191 - (X2 * UR - Y2 * UI);
          Y1 = -(X2 * UI + Y2 * UR);
          X2 = 2186.181 - (X1 * UR - Y1 * UI);
          Y2 = -(X1 * UI + Y1 * UR);
          X1 = 9022.228 - (X2 * UR - Y2 * UI);
          Y1 = -(X2 * UI + Y2 * UR);
          X2 = 24322.84 - (X1 * UR - Y1 * UI);
          Y2 = -(X1 * UI + Y1 * UR);
          UU = 32066.6 - (X2 * UR - Y2 * UI);
          VV = -(X2 * UI + Y2 * UR);
          VOIGT = exp(UR) * cos(UI) - (XX * UU + YY * VV) / (UU * UU + VV * VV);
        }

        /*  Line absorption with the VOIGT function */

        ALINE1 = VOIGT * state->LINEOP[ITAU][LINE];
        if (state->initNLTE)
        {
          ALINE1 *= XNLTE; // NLTE correction to the line opacity
          source_line[ITAU] += ALINE1 * SRC_line;
        }
        if (ALINE1 / OPCON > state->ALMAX[LINE])
          state->ALMAX[LINE] = ALINE1 / OPCON;
      }
      ALINE += ALINE1;
    }

    /* Compute total opacity */

    if (state->MOTYPE > 0) // state->RHOX model
    {
      XK[ITAU] = ALINE + OPCON;
      XC[ITAU] = OPCON;
    }
    else if (state->MOTYPE == 0) // TAU model
    {
      XK[ITAU] = (ALINE + OPCON) / state->COPSTD[ITAU];
      XC[ITAU] = OPCON / state->COPSTD[ITAU];
    }
    else if (state->MOTYPE == -1)
    {
      XK[ITAU] = ALINE;
      XC[ITAU] = OPCON;
    }
    if (state->initNLTE)
      source_line[ITAU] = (source_line[ITAU] + OPCON * SRC_cont) / (ALINE + OPCON);
  }
}

#undef Z
#undef PI4
#undef K
#undef M0
#undef A0

void OPMTRX1(int LINE, double *XK, GlobalState *state)
{
  /*
   THIS FUNCTION CALCULATES THE OPACITY OR OPACITY RATIO (OPACWL/OPACSTD)
   PER GRAMM OF STELLAR MATER (CM**2/GM) PER ANGSTROEM AT DEPTH #IM
   OF THE STANDARD MODEL DEPTH SCALE. WAVELENGTH IS THE CENTRAL
   WAVELENGTH OF LINE "LINE".
   
   For comments and constants description see OPMTRX above.

   Author: N.Piskunov
 
   C++ Version: January 15, 1999
   LAST UPDATE: See OPMTRX above
  */

#define Z 0.026540045e0
#define PI4 4. * PI
#define K 1.380658e-23
#define M0 1.660540e-27
#define A0 5.29177249e-11

  double OPCON, A, UAV, W4, V4,
      XNELEC, XNATOM, ALINE, VOIGT,
      TEMPER, DOPL;
  short ITAU;

  for (ITAU = 0; ITAU < state->NRHOX; ITAU++)
  {
    TEMPER = state->T[ITAU];
    XNELEC = state->XNE[ITAU]; /* Electron number density */
    XNATOM = state->XNA[ITAU]; /* Atom number density     */

    /* Loop through spectral lines */

    ALINE = 0.;
    {
      if (!strncmp(state->spname + 8 * LINE, "H ", 2)) // This is a hydrogen line
      {
        int NBLO, NBUP;
        float temper, xnelec, h1frc, he1frc, dopl, aline;
        double wave, wlcent;

        NBLO = (int)(state->GAMQST[LINE] + 0.1);
        NBUP = (int)(state->GAMVW[LINE] + 0.1);
        temper = TEMPER;
        xnelec = XNELEC;
        h1frc = state->H1FRACT[ITAU];
        he1frc = state->HE1FRACT[ITAU];
        dopl = state->VVOIGT[ITAU][LINE];
        wave = state->WLCENT[LINE];
        wlcent = state->WLCENT[LINE];

        hlinprof_(wave, wlcent, temper, xnelec, NBLO, NBUP,
                  h1frc, he1frc, dopl, aline, state->PATH, &state->PATHLEN, &state->change_byte_order);
        ALINE = aline * state->LINEOP[ITAU][LINE];
      }
      else // Non-hydrogen line
      {

        /*  VOIGT function calculation: Humlicek, J. 1982, J.Q.S.R.state->T. 27, 437
        stripted for the case of line center (V==0) */

        A = state->AVOIGT[ITAU][LINE] * state->WLCENT[LINE];
        UAV = A * A;
        if (A >= 15.)
          W4 = A * 0.5641896 / (0.5 + UAV);
        else if (A >= 5.5)
          W4 = A * (1.410474 + UAV * 0.5641896) / (0.75 + UAV * (3. + UAV));
        else if (A >= -0.176)
          W4 = (16.4955 + A * (20.20933 + A * (11.96482 + A * (3.778987 + A * 0.5642236)))) /
               (16.4955 + A * (38.82363 + A * (39.27121 + A * (21.69274 + A * (6.699398 + A)))));
        else
        {
          W4 = A * (36183.31 - UAV * (3321.9905 - UAV * (1540.787 - UAV * (219.0313 - UAV * (35.76683 - UAV * (1.320522 - UAV * .56419))))));
          V4 = (32066.6 - UAV * (24322.84 - UAV * (9022.228 - UAV * (2186.181 - UAV * (364.2191 - UAV * (61.57037 - UAV * (1.841439 - UAV)))))));
          W4 = exp(UAV) - W4 / V4;
        }
        VOIGT = W4;

        /*  Line absorption with the VOIGT function */

        ALINE = VOIGT * state->LINEOP[ITAU][LINE] * state->WLCENT[LINE];
      }
    }

    /* Compute total opacity */

    if (state->MOTYPE > 0)
      XK[ITAU] = ALINE;
    else if (state->MOTYPE == 0)
      XK[ITAU] = ALINE / state->COPSTD[ITAU];
    else if (state->MOTYPE == -1)
      XK[ITAU] = ALINE;
  }
}

#undef Z
#undef PI4
#undef K
#undef M0
#undef A0

void GAMHE(short IND, double temp, double ANE, double ANP,
           double &GAM, double &SHIFT, GlobalState *state)
{
  /*   NEUTRAL HELIUM STARK BROADENING PARAMETERS
     AFTER DIMITRIJEVIC AND SAHAL-BRECHOT, 1984, J.Q.S.R.state->T. 31, 301
     OR FREUDENSTEIN AND COOPER, 1978, AP.J. 224, 1079  (FOR C(IND)>0)
  */
  static double W[20][5] =
      /*   ELECTRONS state->T= 5000   10000   20000   40000     LAMBDA */
      {{5.990, 6.650, 6.610, 6.210, 3819.60},
       {2.950, 3.130, 3.230, 3.300, 3867.50},
       {109.000, 94.400, 79.500, 65.700, 3871.79},
       {0.142, 0.166, 0.182, 0.190, 3888.65},
       {70.700, 60.700, 50.900, 41.900, 3926.53},
       {1.540, 1.480, 1.400, 1.290, 3964.73},
       {41.600, 50.500, 57.400, 65.800, 4009.27},
       {1.320, 1.350, 1.380, 1.460, 4120.80},
       {7.830, 8.750, 8.690, 8.040, 4143.76},
       {5.830, 6.370, 6.820, 6.990, 4168.97},
       {2.280, 2.320, 2.360, 2.430, 4437.55},
       {2.470, 2.200, 1.910, 1.650, 4471.50},
       {0.588, 0.620, 0.641, 0.659, 4713.20},
       {2.600, 2.480, 2.240, 1.960, 4921.93},
       {0.627, 0.597, 0.568, 0.532, 5015.68},
       {1.050, 1.090, 1.110, 1.140, 5047.74},
       {0.277, 0.298, 0.296, 0.293, 5875.70},
       {0.714, 0.666, 0.602, 0.538, 6678.15},
       {3.490, 3.630, 3.470, 3.190, 4026.20},
       {4.970, 5.100, 4.810, 4.310, 4387.93}};
  static double V[20][4] =
      /*   PROTONS   state->T= 5000   10000   20000   40000 */
      {{1.520, 4.540, 9.140, 10.200},
       {0.607, 0.710, 0.802, 0.901},
       {0.000, 0.000, 0.000, 0.000},
       {0.0396, 0.0434, 0.0476, 0.0526},
       {0.000, 0.000, 0.000, 0.000},
       {0.507, 0.585, 0.665, 0.762},
       {0.930, 1.710, 13.600, 27.200},
       {0.288, 0.325, 0.365, 0.410},
       {1.330, 6.800, 12.900, 14.300},
       {1.100, 1.370, 1.560, 1.760},
       {0.516, 0.579, 0.650, 0.730},
       {1.520, 1.730, 1.830, 1.630},
       {0.128, 0.143, 0.161, 0.181},
       {2.040, 2.740, 2.950, 2.740},
       {0.187, 0.210, 0.237, 0.270},
       {0.231, 0.260, 0.291, 0.327},
       {0.0591, 0.0650, 0.0719, 0.0799},
       {0.231, 0.260, 0.295, 0.339},
       {2.180, 3.760, 4.790, 4.560},
       {1.860, 5.320, 7.070, 7.150}};
  static double SHIFTE[20][4] =
      /*  Shifts due to electrons */
      {{-0.698, -0.558, -0.354, -0.216},
       {1.800, 1.930, 1.810, 1.670},
       {8.510, 5.340, 2.560, 1.560},
       {0.075, 0.061, 0.049, 0.035},
       {7.130, 4.270, 1.960, 0.560},
       {-0.459, -0.345, -0.249, -0.179},
       {10.400, 20.700, 29.700, 38.000},
       {0.890, 0.931, 0.851, 0.677},
       {0.924, 0.856, 0.775, 0.656},
       {3.120, 3.430, 3.490, 3.500},
       {1.690, 1.600, 1.270, 0.906},
       {0.062, -0.064, -0.015, -0.006},
       {0.409, 0.456, 0.439, 0.349},
       {0.436, 0.368, 0.298, 0.221},
       {-0.236, -0.179, -0.132, -0.095},
       {0.730, 0.745, 0.668, 0.528},
       {-0.073, -0.040, -0.012, -0.005},
       {0.249, 0.222, 0.180, 0.144},
       {-0.425, -0.315, -0.209, -0.136},
       {0.665, 0.558, 0.450, 0.336}};
  static double SHIFTP[20][4] =
      /*  Shifts due to protons */
      {{0.000, 0.055, 1.790, 6.100},
       {0.243, 0.422, 0.579, 0.725},
       {0.000, 0.000, 0.000, 0.000},
       {0.028, 0.033, 0.039, 0.044},
       {0.000, 0.000, 0.000, 0.000},
       {-0.232, -0.367, -0.488, -0.602},
       {0.000, 0.000, 0.089, 4.630},
       {0.170, 0.234, 0.294, 0.351},
       {0.000, 0.028, 1.540, 6.750},
       {0.280, 0.676, 1.030, 1.340},
       {0.465, 0.532, 0.604, 0.684},
       {1.350, 1.560, 1.840, 2.110},
       {0.094, 0.117, 0.139, 0.161},
       {0.261, 1.140, 2.010, 2.650},
       {-0.131, -0.164, -0.197, -0.231},
       {0.158, 0.203, 0.246, 0.288},
       {-0.045, -0.052, -0.060, -0.069},
       {0.171, 0.211, 0.250, 0.292},
       {0.002, 0.544, 2.200, 3.680},
       {0.001, 0.359, 2.770, 5.140}};
  static double C[20] = {0., 0., 1.83e-4, 0., 1.13e-4, 0., 0., 0., 0., 0., 1.6e-4,
                         0., 0., 0., 0., 0., 0., 0., 0., 0.};
  static double TT1 = 3.699, TT2 = 4., TT3 = 4.301, TT4 = 4.602;
  double TLG, TJ, TJ0, TJ1, TJ2;
  short J;

  if (W[IND][0] != 0.0)
  {

    /* CUBIC INTERPOLATION OVER state->T=5000,10000,20000,40000 IN LOG SCALE */

    TLG = log10(temp);
    if (TLG <= TT3)
    {
      J = 3;
      TJ = (TT3 - TT2) * (TT3 - TT1) * (TT2 - TT1);
      TJ0 = (TLG - TT1) * (TLG - TT2) * (TT2 - TT1) / TJ;
      TJ1 = (TLG - TT1) * (TT3 - TLG) * (TT3 - TT1) / TJ;
      TJ2 = (TLG - TT2) * (TLG - TT3) * (TT3 - TT2) / TJ;
    }
    else
    {
      J = 4;
      TJ = (TT4 - TT3) * (TT4 - TT2) * (TT3 - TT2);
      TJ0 = (TLG - TT2) * (TLG - TT3) * (TT3 - TT2) / TJ;
      TJ1 = (TLG - TT2) * (TT4 - TLG) * (TT4 - TT2) / TJ;
      TJ2 = (TLG - TT3) * (TLG - TT4) * (TT4 - TT3) / TJ;
    }
    GAM = ((TJ0 * W[IND][J] + TJ1 * W[IND][J - 1] + TJ2 * W[IND][J - 2]) * ANE + (TJ0 * V[IND][J] + TJ1 * V[IND][J - 1] + TJ2 * V[IND][J - 2]) * ANP) * 1.884e3 / (W[IND][4] * W[IND][4]);
    if (GAM < 0.)
      GAM = 0.;
    SHIFT = (TJ0 * SHIFTE[IND][J] + TJ1 * SHIFTE[IND][J - 1] + TJ2 * SHIFTE[IND][J - 2]) * (ANE / 1.e16) +
            (TJ0 * SHIFTP[IND][J] + TJ1 * SHIFTP[IND][J - 1] + TJ2 * SHIFTP[IND][J - 2]) * (ANP / 1.e16);
  }
  else
  {
    GAM = C[IND] * pow(temp, 0.16667) * ANE;
    SHIFT = 0;
  }
}

double VACAIR(double W)
{
  //  W IS VACUUM WAVELENGTH IN Angstroms

  double WAVEN;

  WAVEN = 1.e8 / W;
  WAVEN *= WAVEN;
  return W / (1.00008342130 + 2406030.0 / (1.30e10 - WAVEN) + 15997.0 / (3.89e9 - WAVEN));
}
