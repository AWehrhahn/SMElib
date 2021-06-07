#include <stddef.h>
#include <string.h>
#include "sme_python_bridge.h"
#include "sme_synth_faster.h"

// define the external methods

/* Return SME library version */
const char * Python_SMELibraryVersion(){
    return SMELibraryVersion(0, NULL);
} 

/* Return the required data files */
const char * Python_GetDataFiles(){
    return GetDataFiles(0, NULL);
}

/* Return the current data file directory */
const char * Python_GetLibraryPath(){
    return GetLibraryPath(0, NULL);
}

/* Set the data file directory */
const char * Python_SetLibraryPath(IDL_STRING * path){
    void * args[] = {path};
    return SetLibraryPath(1, args);
}

/* Read in Wavelength range */
const char * Python_InputWaveRange(double wmin, double wmax){
    void * args[] = {&wmin, &wmax};
    return InputWaveRange(2, args); 
}

/* Set van der Waals scaling factor */
const char * Python_SetVWscale(double gamma6){
    void * args[] = {&gamma6};
    return SetVWscale(1, args);
}

/* Set flag for H2 molecule */
const char * Python_SetH2broad(){
    return SetH2broad(0, NULL);
}

/* Clear flag for H2 molecule */
const char * Python_ClearH2broad(){
    return ClearH2broad(0, NULL);
}

/* Read in line list */
const char * Python_InputLineList(
        int nlines, 
        IDL_STRING * species, 
        double * atomic){
    void * args[] = {&nlines, species, atomic};
    return InputLineList(3, args);
}

/* Return line list */
const char * Python_OutputLineList(int nlines, double * atomic){
    void * args[] = {&nlines, atomic};
    return OutputLineList(2, args);
}

/* Change line list parameters */
const char * Python_UpdateLineList(
        short nlines,
        IDL_STRING * species,
        double * atomic,
        short * index){
    void * args[] = {&nlines, species, atomic, index};
    return UpdateLineList(4, args);
}

/* Read in model atmosphere */
const char * Python_InputModel(
        short ndepth, 
        double teff, 
        double grav, 
        double wlstd, 
        IDL_STRING * motype, 
        short * opflag, 
        double * depth, 
        double * temp, 
        double * xne, 
        double * xna, 
        double * rho, 
        double * vt, 
        double radius, 
        double * height){
    if (strcmp(motype->s, "SPH")){
        void * args[] = {&ndepth, &teff, &grav, &wlstd, motype, &radius, opflag, depth, temp, xne, xna, rho, vt, height};
        return InputModel(14, args);
    } else {
        void * args[] = {&ndepth, &teff, &grav, &wlstd, motype, opflag, depth, temp, xne, xna, rho, vt};
        return InputModel(12, args);
    }
}

const char * Python_InputDepartureCoefficients(double * bmat, int lineindex){
    void * args[] = {bmat, &lineindex};
    return InputDepartureCoefficients(2, args);
}

/* Get NLTE b's for specific line */
const char * Python_GetDepartureCoefficients(double * bmat, int nrhox, int line){
    void * args[] = {bmat, &nrhox, &line};
    return GetDepartureCoefficients(3, args);
}

/* Get line list NLTE flags */
const char * Python_GetNLTEflags(short * nlte_flags, int nlines){
    void * args[] = {nlte_flags, &nlines};
    return GetNLTEflags(2, args);
}

/* Reset LTE */
const char * Python_ResetDepartureCoefficients(){
    return ResetDepartureCoefficients(0, NULL);
}

/* Read in abundances */
const char * Python_InputAbund(double * abund){
    void * args[] = {abund};
    return InputAbund(1, args);
}

/* Calculate opacities */
const char * Python_Opacity(){
    return Opacity(0, NULL);
}

/* Returns specific cont. opacity */
const char * Python_GetOpacity(short ifop, short length, double * result, IDL_STRING * species, IDL_STRING * key){
    if (ifop == 8){
        void * args[] = {&ifop, &length, result, species, key};
        return GetOpacity(5, args);
    } else if (ifop == 9)
    {
        void * args[] = {&ifop, &length, result, species};
        return GetOpacity(4, args);
    } else {
        void * args[] = {&ifop, &length, result};
        return GetOpacity(3, args);
    }
}

/* Perfrom EOS calculations */
const char * Python_Ionization(short ion){
    void * args[] = {&ion};
    return Ionization(1, args);
}

/* Returns density in g/cm^3 */
const char * Python_GetDensity(short length, double * result){
    void * args[] = {&length, result};
    return GetDensity(2, args);
}

/* Returns atomic number density */
const char * Python_GetNatom(short length, double * result){
    void * args[] = {&length, result};
    return GetNatom(2, args);
}

/* Returns electron number density */
const char * Python_GetNelec(short length, double * result){
    void * args[] = {&length, result};
    return GetNelec(2, args);
}

/* Computes spectral synthesis */
const char * Python_Transf(
    short nmu, 
    double * mu, 
    double * cint_seg, 
    double * cintr_seg, 
    int nwmax, 
    int nw, 
    double * wint_seg, 
    double * sint_seg, 
    double accrt, 
    double accwi, 
    short keep_lineop, 
    short long_continuum){
    void * args[] = {&nmu, mu, cint_seg, cintr_seg, &nwmax, &nw, wint_seg, sint_seg, &accrt, &accwi, &keep_lineop, & long_continuum};
    return Transf(12, args);
}

/* Computes line central depths */
const char * Python_CentralDepth(int nmu, double * mu, int nwsize, float * table, double accrt){
    void * args[] = {&nmu, mu, &nwsize, table, &accrt};
    return CentralDepth(5, args);
}

/* Returns specific line opacity */
const char * Python_GetLineOpacity(double wave, short nmu, double * lop, double * cop, double * scr, double * tsf, double * csf){
    void * args[] = {&wave, &nmu, lop, cop, scr, tsf, csf};
    return GetLineOpacity(7, args);
}

/* Get validity range for every line */
const char * Python_GetLineRange(double * linerange, int nlines){
    void * args[] = {linerange, &nlines};
    return GetLineRange(2, args);
}
