#include <stddef.h>
#include <string.h>
#include "sme_cython_interface.h"
#include "sme_synth_parallel.h"

// Cython only has a single global state
// Leave it at NULL until requested, at which point it will be initialized
GlobalState *cython_state = NULL;

GlobalState *Cython_GetState()
{
    if (cython_state == NULL)
    {
        cython_state = _NewState();
    }
    return cython_state;
}

/* Return SME library version */
const char *Cython_SMELibraryVersion()
{
    return _SMELibraryVersion();
}

/* Return the required data files */
const char *Cython_GetDataFiles()
{
    return _GetDataFiles();
}

/* Return the current data file directory */
const char *Cython_GetLibraryPath()
{
    GlobalState *state = Cython_GetState();
    return _GetLibraryPath(state);
}

/* Set the data file directory */
const char *Cython_SetLibraryPath(const char *path, int pathlen)
{
    GlobalState *state = Cython_GetState();
    return _SetLibraryPath(path, pathlen, state);
}

/* Read in Wavelength range */
const char *Cython_InputWaveRange(double wmin, double wmax)
{
    GlobalState *state = Cython_GetState();
    return _InputWaveRange(wmin, wmax, state);
}

/* Set van der Waals scaling factor */
const char *Cython_SetVWscale(double gamma6)
{
    GlobalState *state = Cython_GetState();
    return _SetVWscale(gamma6, state);
}

/* Set flag for H2 molecule */
const char *Cython_SetH2broad()
{
    GlobalState *state = Cython_GetState();
    return _SetH2broad(1, state);
}

/* Clear flag for H2 molecule */
const char *Cython_ClearH2broad()
{
    GlobalState *state = Cython_GetState();
    return _SetH2broad(0, state);
}

/* Read in line list */
const char *Cython_InputLineList(
    int nlines,
    int slen,
    const char *species,
    double *atomic)
{
    GlobalState *state = Cython_GetState();
    return _InputLineList(nlines, slen, species, atomic, state);
}

/* Return line list */
const char *Cython_OutputLineList(int nlines, double *atomic)
{
    GlobalState *state = Cython_GetState();
    return _OutputLineList(nlines, atomic, state);
}

/* Change line list parameters */
const char *Cython_UpdateLineList(
    short nlines,
    int slen,
    const char *species,
    double *atomic,
    short *index)
{
    GlobalState *state = Cython_GetState();
    return _UpdateLineList(nlines, slen, index, species, atomic, state);
}

/* Read in model atmosphere */
const char *Cython_InputModel(
    short ndepth,
    double teff,
    double grav,
    double wlstd,
    const char *motype,
    int mlen,
    short *opflag,
    double *depth,
    double *temp,
    double *xne,
    double *xna,
    double *rho,
    double *vt,
    double radius,
    double *height)
{
    GlobalState *state = Cython_GetState();
    return _InputModel(ndepth, teff, grav, wlstd, motype, mlen, opflag, depth, temp, xne, xna, rho, vt, radius, height, state);
}

const char *Cython_InputDepartureCoefficients(double *bmat, int lineindex)
{
    GlobalState *state = Cython_GetState();
    return _InputDepartureCoefficients(bmat, lineindex, state);
}

/* Get NLTE b's for specific line */
const char *Cython_GetDepartureCoefficients(double *bmat, int nrhox, int line)
{
    GlobalState *state = Cython_GetState();
    return _GetDepartureCoefficients(bmat, nrhox, line, state);
}

/* Get line list NLTE flags */
const char *Cython_GetNLTEflags(short *nlte_flags, int nlines)
{
    GlobalState *state = Cython_GetState();
    return _GetNLTEflags(nlte_flags, nlines, state);
}

/* Reset LTE */
const char *Cython_ResetDepartureCoefficients()
{
    GlobalState *state = Cython_GetState();
    return _ResetDepartureCoefficients(state);
}

/* Read in abundances */
const char *Cython_InputAbund(double *abund)
{
    GlobalState *state = Cython_GetState();
    return _InputAbund(abund, 99, state);
}

/* Calculate opacities */
const char *Cython_Opacity()
{
    GlobalState *state = Cython_GetState();
    return _Opacity(0, 0, NULL, NULL, NULL, state);
}

/* Returns specific cont. opacity */
const char *Cython_GetOpacity(short ifop, short length, double *result, const char *species, int slen, const char *key, int klen)
{
    GlobalState *state = Cython_GetState();
    return _GetOpacity(ifop, length, result, species, slen, key, klen, state);
}

/* Perfrom EOS calculations */
const char *Cython_Ionization(short ion)
{
    GlobalState *state = Cython_GetState();
    return _Ionization(ion, state);
}

/* Returns density in g/cm^3 */
const char *Cython_GetDensity(short length, double *result)
{
    GlobalState *state = Cython_GetState();
    return _GetDensity(length, result, state);
}

/* Returns atomic number density */
const char *Cython_GetNatom(short length, double *result)
{
    GlobalState *state = Cython_GetState();
    return _GetNatom(length, result, state);
}

/* Returns electron number density */
const char *Cython_GetNelec(short length, double *result)
{
    GlobalState *state = Cython_GetState();
    return _GetNelec(length, result, state);
}

/* Computes spectral synthesis */
const char *Cython_Transf(
    short nmu,
    double *mu,
    double *cint_seg,
    double *cintr_seg,
    int nwmax,
    int nw,
    double *wint_seg,
    double *sint_seg,
    double accrt,
    double accwi,
    short keep_lineop,
    short long_continuum)
{
    GlobalState *state = Cython_GetState();
    return _Transf(nmu, mu, cint_seg, cintr_seg, nwmax, nw, wint_seg, sint_seg, accrt, accwi, keep_lineop, long_continuum, state);
}

/* Computes line central depths */
const char *Cython_CentralDepth(int nmu, double *mu, int nwsize, float *table, double accrt)
{
    GlobalState *state = Cython_GetState();
    return _CentralDepth(nmu, mu, nwsize, table, accrt, state);
}

/* Returns specific line opacity */
const char *Cython_GetLineOpacity(double wave, short nrhox, double *lop, double *cop, double *scr, double *tsf, double *csf)
{
    GlobalState *state = Cython_GetState();
    return _GetLineOpacity(wave, nrhox, lop, cop, scr, tsf, csf, state);
}

/* Get validity range for every line */
const char *Cython_GetLineRange(double *linerange, int nlines)
{
    GlobalState *state = Cython_GetState();
    return _GetLineRange(linerange, nlines, state);
}
