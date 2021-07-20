#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "platform.h"
#include "sme_idl_interface.h"
#include "sme_synth_parallel.h"

#define pow10(x) exp(2.30258509299405e0 * (x))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define round(x) (x >= 0) ? (int)(x + 0.5) : (int)(x - 0.5)

// IDL only has a single global state
// Leave it at NULL until requested, at which point it will be initialized
GlobalState *idl_state = NULL;

GlobalState *IDL_GetState()
{
  if (idl_state == NULL)
  {
    idl_state = _NewState();
  }
  return idl_state;
}

extern "C" int SME_DLL GetNLINES(int n, void *args[])
{
  GlobalState *state = IDL_GetState();
  return state->NLINES;
}

extern "C" short SME_DLL GetNRHOX(int n, void *args[])
{
  GlobalState *state = IDL_GetState();
  return state->NRHOX;
}

extern "C" char *SME_DLL GetSPNAME(int n, void *args[])
{
  GlobalState *state = IDL_GetState();
  return state->spname;
}

extern "C" char const *SME_DLL SMELibraryVersion(int n, void *arg[]) /* Return SME library version */
{
  return _SMELibraryVersion();
}

extern "C" char const *SME_DLL GetDataFiles(int n, void *arg[]) /* Return SME library version */
{
  return _GetDataFiles();
}

extern "C" char const *SME_DLL GetLibraryPath(int n, void *arg[])
{
  GlobalState *state = IDL_GetState();
  return _GetLibraryPath(state);
}

/*
  Set SME library datafile location
  If smelib was installed using make install the default location should point to the data files already
*/
extern "C" char const *SME_DLL SetLibraryPath(int n, void *arg[])
{
  if (n != 1)
  {
    return "No path was specified";
  }
  GlobalState *state = IDL_GetState();
  IDL_STRING s = *(IDL_STRING *)arg[0];
  return _SetLibraryPath(s.s, s.slen, state);
}

extern "C" char const *SME_DLL InputWaveRange(int n, void *arg[]) /* Read in Wavelength range */
{
  if (n < 2)
  {
    return "Only one argument found";
  }

  GlobalState *state = IDL_GetState();
  double wfirst = *(double *)arg[0];
  double wlast = *(double *)arg[1];

  return _InputWaveRange(wfirst, wlast, state);
}

extern "C" char const *SME_DLL SetVWscale(int n, void *arg[]) /* Set van der Waals scaling factor */
{
  if (n < 1)
  {
    return "Not enough arguments";
  }
  GlobalState *state = IDL_GetState();
  double vw_scale = *(double *)arg[0];
  return _SetVWscale(vw_scale, state);
}

extern "C" char const *SME_DLL SetH2broad(int n, void *arg[]) /* Set flag for H2 molecule */
{
  GlobalState *state = IDL_GetState();
  return _SetH2broad(1, state);
}

extern "C" char const *SME_DLL ClearH2broad(int n, void *arg[]) /* Clear flag for H2 molecule */
{
  GlobalState *state = IDL_GetState();
  return _SetH2broad(0, state);
}

extern "C" char const *SME_DLL InputLineList(int n, void *arg[]) /* Read in line list */
{

  if (n < 2)
  {
    return "Not enough arguments";
  }

  GlobalState *state = IDL_GetState();
  int slen = 8;
  int nlines = *(int *)arg[0];
  IDL_STRING *a0 = (IDL_STRING *)arg[1];
  double *linelist = (double *)arg[2];
  char *species = (char *)calloc(slen * nlines, sizeof(char));

  for (int i = 0; i < nlines; i++)
  {
    memcpy(species + slen * i, a0[i].s, a0[i].slen);
    if (a0[i].slen < slen)
      for (int l = a0[i].slen; l < slen; l++)
        species[slen * i + l] = ' ';
  }
  const char *response = _InputLineList(nlines, slen, species, linelist, state);

  free(species);

  return response;
}

extern "C" char const *SME_DLL OutputLineList(int n, void *arg[]) /* Return line list */
{
  /*
   state->NLINES - NUMBERS OF SPECTRAL LINES;
   For each line:
   state->GAMRAD - RADIATION DAMPING (C1);
   state->GAMQST - QUADRATIC STARK DUMPING (C4);
   state->GAMVW  - VAN DER WAALS DUMPING (C6);
  */

  if (n < 2)
  {
    return "Not enough arguments";
  }

  GlobalState *state = IDL_GetState();
  int nlines = *(int *)arg[0];
  double *a1 = (double *)arg[1];

  return _OutputLineList(nlines, a1, state);
}

extern "C" char const *SME_DLL UpdateLineList(int n, void *arg[]) /* Change line list parameters */
{
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
    return "Not enough arguments";
  }

  GlobalState *state = IDL_GetState();
  short nupdate = *(short *)arg[0];
  IDL_STRING *a0 = (IDL_STRING *)arg[1]; /* Setup pointers for species        */
  double *a1 = (double *)arg[2];         /* Setup pointers to line parameters */
  short *index = (short *)arg[3];

  int slen = 8;
  char *species = (char *)calloc(slen * nupdate, sizeof(char));

  for (int i = 0; i < nupdate; i++)
  {
    memcpy(species + slen * i, a0[i].s, a0[i].slen);
    if (a0[i].slen < slen)
      for (int l = a0[i].slen; l < slen; l++)
        species[slen * i + l] = ' ';
  }

  const char *response = _UpdateLineList(nupdate, slen, index, species, a1, state);

  free(species);

  return response;
}

extern "C" char const *SME_DLL InputModel(int n, void *arg[]) /* Read in model atmosphere */
{
  if (n < 12)
  {
    return "Not enough arguments";
  }

  GlobalState *state = IDL_GetState();
  short nrhox = *(short *)arg[0];
  double teff = *(double *)arg[1];
  double grav = *(double *)arg[2];
  double wlstd = *(double *)arg[3];

  char motype[8];
  IDL_STRING *s = (IDL_STRING *)arg[4];
  short l = min(7, s->slen);
  strncpy(motype, s->s, l);
  motype[l] = 0;
  for (size_t i = 0; i < strlen(motype); i++)
  {
    motype[i] = toupper(motype[i]);
  }

  int arg_offset = 0;
  double radius = 0;
  double *height = NULL;
  // Adding provision for spherical models
  if (!strncmp(motype, "SPH", 3))
  {
    arg_offset = 1;
    radius = *(double *)arg[5];
    height = (double *)arg[12 + arg_offset];
  }

  short *ifop = (short *)arg[5 + arg_offset];
  double *depth = (double *)arg[6 + arg_offset];
  double *temp = (double *)arg[7 + arg_offset];
  double *xne = (double *)arg[8 + arg_offset];
  double *xna = (double *)arg[9 + arg_offset];
  double *rho = (double *)arg[10 + arg_offset];
  double *vt = (double *)arg[11 + arg_offset];

  return _InputModel(nrhox, teff, grav, wlstd, motype, 8, ifop, depth, temp, xne, xna, rho, vt, radius, height, state);
}

extern "C" char const *SME_DLL InputDepartureCoefficients(int n, void *arg[])
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

  // We assume that the caller will provide 2*state->NRHOX element array, so
  // be careful on the IDL side. The other argument is the line number.
  if (n < 2)
  {
    return "No arguments found";
  }

  GlobalState *state = IDL_GetState();
  double *b = (double *)arg[0];
  int line = *(int *)arg[1];

  return _InputDepartureCoefficients(b, line, state);
}

extern "C" char const *SME_DLL GetDepartureCoefficients(int n, void *arg[]) /* Get NLTE b's for specific line */
{
  if (n < 3) // Check if arguments are present
  {
    return "Requires an array pointer, its length and line number";
  }

  GlobalState *state = IDL_GetState();
  int line = *(int *)arg[2];
  double *b = (double *)arg[0];
  int nrhox = *(int *)arg[1];

  return _GetDepartureCoefficients(b, nrhox, line, state);
}

extern "C" char const *SME_DLL GetNLTEflags(int n, void *arg[]) /* Get NLTE flag for every line */
{
  if (n < 2) // Check if arguments are present
  {
    return "GetNLTELines: Requires an array pointer and its length";
  }

  GlobalState *state = IDL_GetState();
  short *b = (short *)arg[0];
  int nlines = *(int *)arg[1];

  return _GetNLTEflags(b, nlines, state);
}

extern "C" char const *SME_DLL ResetDepartureCoefficients(int n, void *arg[]) /* Reset LTE */
{
  GlobalState *state = IDL_GetState();
  return _ResetDepartureCoefficients(state);
}

extern "C" char const *SME_DLL InputAbund(int n, void *arg[]) /* Read in abundances */
{
  if (n < 1)
  {
    return "Not enough arguments";
  }
  GlobalState *state = IDL_GetState();
  double *a = (double *)arg[0];
  return _InputAbund(a, 99, state);
}

extern "C" char const *SME_DLL Opacity(int n, void *arg[]) /* Calculate opacities */
{
  short i, nrhox;
  double *a1 = NULL;
  double *a2 = NULL;
  double *a3 = NULL;
  short request_output = 0;
  GlobalState *state = IDL_GetState();

  if (n > 0)
  {
    if ((state->MOTYPE != 0 && n < 3) ||
        (state->MOTYPE == 0 && n < 4))
    {
      return "Opacity: Not enough arguments";
    }
  }

  if (n >= 3)
  {
    request_output = 1;
    i = *(short *)arg[0]; /* Length of IDL arrays */
    nrhox = min(state->NRHOX, i);
    a1 = (double *)arg[1];
    a2 = (double *)arg[2];
    if (state->MOTYPE == 0)
      a3 = (double *)arg[3];
  }

  // short request_output, short nout, double * out1, double * out2, double * out3
  const char *response = _Opacity(request_output, nrhox, a1, a2, a3, state);

  return response;
}

extern "C" char const *SME_DLL Ionization(int n, void *arg[])
{
  short switches;
  GlobalState * state = IDL_GetState();

  if (n > 0)
  {
    switches = *(short *)arg[0];
  }
  else
  {
    switches = 0;
  }

  return _Ionization(switches, state);
}

extern "C" char const *SME_DLL Transf(int n, void *arg[])
{
  short nmu, keep_lineop, long_continuum;
  double *mu, *fcblue, *fcred, *table, *wl;
  int nwsize, nwl;
  double eps1, eps2;
  GlobalState *state = IDL_GetState();

  if (n > 10) /* New SME software capable of using predefined wavelength grid */
  {
    nmu = *(short *)arg[0];          /* Number of limb points */
    mu = (double *)arg[1];           /* Array of limb points */
    fcblue = (double *)arg[2];       /* Continuum specific intensity on the blue end */
    fcred = (double *)arg[3];        /* Continuum specific intensity on the red end */
    nwsize = *(int *)arg[4];         /* Length of the arrays for synthesis */
    nwl = *(int *)arg[5];            /* Length of predefined wavelength vector */
    wl = (double *)arg[6];           /* Array for wavelengths */
    table = (double *)arg[7];        /* Array for synthetic spectrum */
    eps1 = *(double *)arg[8];        /* Accuracy of the radiative transfer integration */
    eps2 = *(double *)arg[9];        /* Accuracy of the interpolation on wl grid */
    keep_lineop = *(short *)arg[10]; /* For several spectral segments there is no 
                                      point recomputing line opacities. This flag
                                      tells when recalculations are needed */
  }
  else /* Old SME software */
  {
    nmu = *(short *)arg[0];    /* Number of limb points */
    mu = (double *)arg[1];     /* Array of limb points */
    fcblue = (double *)arg[2]; /* Continuum specific intensity on the blue end */
    fcred = (double *)arg[3];  /* Continuum specific intensity on the red end */
    nwsize = *(long *)arg[4];  /* Length of the arrays for synthesis */
    wl = (double *)arg[5];     /* Array for wavelengths */
    table = (double *)arg[6];  /* Array for synthetic spectrum */
    eps1 = *(double *)arg[7];  /* Accuracy of the radiative transfer integration */
    eps2 = *(double *)arg[8];  /* Accuracy of the interpolation on wl grid */
    keep_lineop = 0;
    nwl = 400000;
  }
  if (n > 11) /* Check of continuum is needed at every wavelength */
  {           /* If this flag is true FCBLUE must be an arrays of */
              /* the size NWSIZE. On exit FCRED keeps its meaning */
    long_continuum = *(short *)arg[11];
  }
  else
    long_continuum = 0;

  return _Transf(nmu, mu, fcblue, fcred, nwsize, nwl, wl, table, eps1, eps2, keep_lineop, long_continuum, state);
}