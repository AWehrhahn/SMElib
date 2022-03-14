C=========================================================================
C Kp_Q_spln: Returns equilibrium constant for a given molecule and temperature.
C
C Inputs:
C   SPNAME  [char] species name according to the table below.
C   TEMP    [real] temperature (in K) at which Kp is needed
C
C History:
C  28-jun-2007: First version written by N. Piskunov including 57 species.
C               Molecular equilibium tabulated by P. Barklem, resampled
C               for optimal spline interpolation and converted to Fortran
C               DATA statements by J. Valenti
C
C  15-dec-2007: Second version includes 58 molecular species.
C               Tabulated values are now alog10(Kp)+D0*5040/T vs alog10(T),
C               where Kp is an equilibrium constant in N/m^2, D0 is the
C               dissociation energy (eV) at 0 K, and T is temperature (K).
C               In this version, we start using a separate alog10(T) grid
C               for each species, rather than a common THETA=5040/T grid
C               for all species. We copied D0 from MOLCON in eos.f, except
C               for CH-, OH-, SiH-, SiN, and MgS, which we (JV) deduced from
C               Barklem data.
C
C  26-apr-2019: Subroutine data modified and the subroutine text generated
C               by IDL program qk_spl_nodes_f77.pro with errthr=0.000100
C
C Outputs:
C   K_spln [real*8] equilibrium constant (in dynes/cm^2) at temperature T,
C   Q_spln [real*8] partition functions at temperature T,
C          both interpolated from Paul Barklem's tables.
C
C To obtain molecular equilibrium constants, KP:
C
C   D2 = SPL_INIT(TK_<species>,K_<species>)
C   KP(T) = SPL_INTERP(TK_<species>,K_<species>,D2,TLOG)
C         - D0*5040/T
C
C To obtain partition functions,Q:
C
C   D2 = SPL_INIT(TQ_<species>,Q_<species>)
C   Q(T) = SPL_INTERP(TQ_<species>,Q_<species>,D2,TLOG)
C
C Note that KP_Q_SPLN returns log10(Q) and log10(Kp)+D0*5040/T
C
C Reference:
C   Paul Barklem, Remo Collet, 2016, A&A 588, 96.
C
      SUBROUTINE KP_Q_SPLN(SPNAME,TEMP,Q_spln,K_spln,BARKLEM)
C
      IMPLICIT NONE
      CHARACTER SPNAME*(*)
      REAL TEMP
      LOGICAL BARKLEM
      REAL*8 Q_spln,K_spln
C
C  Local variables
C
      LOGICAL FIRST
      INTEGER MSPEC,NTQ,NTK,KLO,KHI,I,II,ISPEC
      PARAMETER(MSPEC=291, NTQ=27, NTK=90)
      INTEGER MTQ(MSPEC),MTK(MSPEC)
      REAL*8 TLOG,A,U(90),SPL_INTERP
C
      CHARACTER SPLIST(MSPEC)*8
      REAL*8 TQ(NTQ,MSPEC),Q(NTQ,MSPEC),Q2(NTQ,MSPEC)
      REAL*8 TK(NTK,MSPEC),K(NTK,MSPEC),K2(NTK,MSPEC)
      REAL*8         TQ_H2   (NTQ),TQ_Li2  (NTQ),TQ_B2   (NTQ),
     * TQ_C2   (NTQ),TQ_N2   (NTQ),TQ_O2   (NTQ),TQ_F2   (NTQ),
     * TQ_Na2  (NTQ),TQ_Mg2  (NTQ),TQ_Al2  (NTQ),TQ_Si2  (NTQ),
     * TQ_P2   (NTQ),TQ_S2   (NTQ),TQ_Cl2  (NTQ),TQ_K2   (NTQ),
     * TQ_Cu2  (NTQ),TQ_As2  (NTQ),TQ_Se2  (NTQ),TQ_Sb2  (NTQ),
     * TQ_Te2  (NTQ),TQ_I2   (NTQ),TQ_Cs2  (NTQ),TQ_H2p  (NTQ),
     * TQ_He2p (NTQ),TQ_C2p  (NTQ),TQ_N2p  (NTQ),TQ_O2p  (NTQ),
     * TQ_Ne2p (NTQ),TQ_P2p  (NTQ),TQ_S2p  (NTQ),TQ_H2m  (NTQ),
     * TQ_C2m  (NTQ),TQ_LiH  (NTQ),TQ_BeH  (NTQ),TQ_BH   (NTQ),
     * TQ_CH   (NTQ),TQ_NH   (NTQ),TQ_OH   (NTQ),TQ_HF   (NTQ),
     * TQ_NaH  (NTQ),TQ_MgH  (NTQ),TQ_AlH  (NTQ),TQ_SiH  (NTQ),
     * TQ_PH   (NTQ),TQ_HS   (NTQ),TQ_HCl  (NTQ),TQ_KH   (NTQ),
     * TQ_CaH  (NTQ),TQ_TiH  (NTQ),TQ_CrH  (NTQ),TQ_MnH  (NTQ),
     * TQ_FeH  (NTQ),TQ_CoH  (NTQ),TQ_NiH  (NTQ),TQ_CuH  (NTQ),
     * TQ_ZnH  (NTQ),TQ_GaH  (NTQ),TQ_GeH  (NTQ),TQ_AsH  (NTQ),
     * TQ_SeH  (NTQ),TQ_HBr  (NTQ),TQ_RbH  (NTQ),TQ_SrH  (NTQ),
     * TQ_AgH  (NTQ),TQ_CdH  (NTQ),TQ_InH  (NTQ),TQ_SnH  (NTQ),
     * TQ_SbH  (NTQ),TQ_TeH  (NTQ),TQ_HI   (NTQ),TQ_CsH  (NTQ),
     * TQ_BaH  (NTQ),TQ_YbH  (NTQ),TQ_PtH  (NTQ),TQ_AuH  (NTQ),
     * TQ_HgH  (NTQ),TQ_TlH  (NTQ),TQ_PbH  (NTQ),TQ_BiH  (NTQ),
     * TQ_HeHp (NTQ),TQ_BeHp (NTQ),TQ_CHp  (NTQ),TQ_NHp  (NTQ),
     * TQ_OHp  (NTQ),TQ_HFp  (NTQ),TQ_NeHp (NTQ),TQ_MgHp (NTQ),
     * TQ_AlHp (NTQ),TQ_SiHp (NTQ),TQ_PHp  (NTQ),TQ_SHp  (NTQ),
     * TQ_HClp (NTQ),TQ_ZnHp (NTQ),TQ_HBrp (NTQ),TQ_CdHp (NTQ),
     * TQ_HgHp (NTQ),TQ_CHm  (NTQ),TQ_OHm  (NTQ),TQ_SiHm (NTQ),
     * TQ_HSm  (NTQ),TQ_CN   (NTQ),TQ_CO   (NTQ),TQ_CF   (NTQ),
     * TQ_SiC  (NTQ),TQ_CP   (NTQ),TQ_CS   (NTQ),TQ_CCl  (NTQ),
     * TQ_CSe  (NTQ),TQ_CBr  (NTQ),TQ_RhC  (NTQ),TQ_IrC  (NTQ),
     * TQ_PtC  (NTQ),TQ_CNp  (NTQ),TQ_COp  (NTQ),TQ_CNm  (NTQ),
     * TQ_CSm  (NTQ),TQ_BN   (NTQ),TQ_NO   (NTQ),TQ_NF   (NTQ),
     * TQ_AlN  (NTQ),TQ_SiN  (NTQ),TQ_PN   (NTQ),TQ_NS   (NTQ),
     * TQ_NCl  (NTQ),TQ_TiN  (NTQ),TQ_AsN  (NTQ),TQ_SeN  (NTQ),
     * TQ_ZrN  (NTQ),TQ_NOp  (NTQ),TQ_NSp  (NTQ),TQ_LiO  (NTQ),
     * TQ_BeO  (NTQ),TQ_BO   (NTQ),TQ_FO   (NTQ),TQ_NaO  (NTQ),
     * TQ_MgO  (NTQ),TQ_AlO  (NTQ),TQ_SiO  (NTQ),TQ_PO   (NTQ),
     * TQ_SO   (NTQ),TQ_ClO  (NTQ),TQ_KO   (NTQ),TQ_CaO  (NTQ),
     * TQ_ScO  (NTQ),TQ_TiO  (NTQ),TQ_VO   (NTQ),TQ_CrO  (NTQ),
     * TQ_MnO  (NTQ),TQ_FeO  (NTQ),TQ_NiO  (NTQ),TQ_CuO  (NTQ),
     * TQ_GaO  (NTQ),TQ_GeO  (NTQ),TQ_AsO  (NTQ),TQ_SeO  (NTQ),
     * TQ_BrO  (NTQ),TQ_RbO  (NTQ),TQ_SrO  (NTQ),TQ_YO   (NTQ),
     * TQ_ZrO  (NTQ),TQ_NbO  (NTQ),TQ_InO  (NTQ),TQ_SnO  (NTQ),
     * TQ_SbO  (NTQ),TQ_TeO  (NTQ),TQ_IO   (NTQ),TQ_BaO  (NTQ),
     * TQ_LaO  (NTQ),TQ_TbO  (NTQ),TQ_LuO  (NTQ),TQ_HfO  (NTQ),
     * TQ_TaO  (NTQ),TQ_WO   (NTQ),TQ_PtO  (NTQ),TQ_PbO  (NTQ),
     * TQ_BiO  (NTQ),TQ_ThO  (NTQ),TQ_BOp  (NTQ),TQ_SiOp (NTQ),
     * TQ_POp  (NTQ),TQ_SOp  (NTQ),TQ_AsOp (NTQ),TQ_TaOp (NTQ),
     * TQ_FeOm (NTQ),TQ_LiF  (NTQ),TQ_BeF  (NTQ),TQ_BF   (NTQ),
     * TQ_NaF  (NTQ),TQ_MgF  (NTQ),TQ_AlF  (NTQ),TQ_SiF  (NTQ),
     * TQ_PF   (NTQ),TQ_SF   (NTQ),TQ_KF   (NTQ),TQ_CaF  (NTQ),
     * TQ_ScF  (NTQ),TQ_MnF  (NTQ),TQ_NiF  (NTQ),TQ_CuF  (NTQ),
     * TQ_ZnF  (NTQ),TQ_GaF  (NTQ),TQ_GeF  (NTQ),TQ_AsF  (NTQ),
     * TQ_SeF  (NTQ),TQ_BrF  (NTQ),TQ_RbF  (NTQ),TQ_SrF  (NTQ),
     * TQ_YF   (NTQ),TQ_AgF  (NTQ),TQ_CdF  (NTQ),TQ_InF  (NTQ),
     * TQ_SnF  (NTQ),TQ_SbF  (NTQ),TQ_IF   (NTQ),TQ_CsF  (NTQ),
     * TQ_BaF  (NTQ),TQ_LaF  (NTQ),TQ_HoF  (NTQ),TQ_YbF  (NTQ),
     * TQ_LuF  (NTQ),TQ_HgF  (NTQ),TQ_TlF  (NTQ),TQ_PbF  (NTQ),
     * TQ_LiNa (NTQ),TQ_AsP  (NTQ),TQ_SbP  (NTQ),TQ_BeS  (NTQ),
     * TQ_BS   (NTQ),TQ_MgS  (NTQ),TQ_AlS  (NTQ),TQ_SiS  (NTQ),
     * TQ_PS   (NTQ),TQ_CaS  (NTQ),TQ_ScS  (NTQ),TQ_TiS  (NTQ),
     * TQ_CrS  (NTQ),TQ_CuS  (NTQ),TQ_GeS  (NTQ),TQ_AsS  (NTQ),
     * TQ_SeS  (NTQ),TQ_SrS  (NTQ),TQ_YS   (NTQ),TQ_SnS  (NTQ),
     * TQ_TeS  (NTQ),TQ_BaS  (NTQ),TQ_LaS  (NTQ),TQ_PbS  (NTQ),
     * TQ_BiS  (NTQ),TQ_LiCl (NTQ),TQ_BeCl (NTQ),TQ_BCl  (NTQ),
     * TQ_NaCl (NTQ),TQ_MgCl (NTQ),TQ_AlCl (NTQ),TQ_SiCl (NTQ),
     * TQ_PCl  (NTQ),TQ_KCl  (NTQ),TQ_CaCl (NTQ),TQ_ScCl (NTQ),
     * TQ_MnCl (NTQ),TQ_FeCl (NTQ),TQ_CuCl (NTQ),TQ_ZnCl (NTQ),
     * TQ_GaCl (NTQ),TQ_GeCl (NTQ),TQ_AsCl (NTQ),TQ_SeCl (NTQ),
     * TQ_BrCl (NTQ),TQ_RbCl (NTQ),TQ_SrCl (NTQ),TQ_YCl  (NTQ),
     * TQ_AgCl (NTQ),TQ_CdCl (NTQ),TQ_InCl (NTQ),TQ_SnCl (NTQ),
     * TQ_SbCl (NTQ),TQ_ICl  (NTQ),TQ_CsCl (NTQ),TQ_BaCl (NTQ),
     * TQ_YbCl (NTQ),TQ_AuCl (NTQ),TQ_HgCl (NTQ),TQ_TlCl (NTQ),
     * TQ_PbCl (NTQ),TQ_AlSe (NTQ),TQ_SiSe (NTQ),TQ_GeSe (NTQ),
     * TQ_KBr  (NTQ),TQ_SiTe (NTQ),TQ_GeTe (NTQ),TQ_KI   (NTQ)
      REAL*8          Q_H2   (NTQ), Q_Li2  (NTQ), Q_B2   (NTQ),
     *  Q_C2   (NTQ), Q_N2   (NTQ), Q_O2   (NTQ), Q_F2   (NTQ),
     *  Q_Na2  (NTQ), Q_Mg2  (NTQ), Q_Al2  (NTQ), Q_Si2  (NTQ),
     *  Q_P2   (NTQ), Q_S2   (NTQ), Q_Cl2  (NTQ), Q_K2   (NTQ),
     *  Q_Cu2  (NTQ), Q_As2  (NTQ), Q_Se2  (NTQ), Q_Sb2  (NTQ),
     *  Q_Te2  (NTQ), Q_I2   (NTQ), Q_Cs2  (NTQ), Q_H2p  (NTQ),
     *  Q_He2p (NTQ), Q_C2p  (NTQ), Q_N2p  (NTQ), Q_O2p  (NTQ),
     *  Q_Ne2p (NTQ), Q_P2p  (NTQ), Q_S2p  (NTQ), Q_H2m  (NTQ),
     *  Q_C2m  (NTQ), Q_LiH  (NTQ), Q_BeH  (NTQ), Q_BH   (NTQ),
     *  Q_CH   (NTQ), Q_NH   (NTQ), Q_OH   (NTQ), Q_HF   (NTQ),
     *  Q_NaH  (NTQ), Q_MgH  (NTQ), Q_AlH  (NTQ), Q_SiH  (NTQ),
     *  Q_PH   (NTQ), Q_HS   (NTQ), Q_HCl  (NTQ), Q_KH   (NTQ),
     *  Q_CaH  (NTQ), Q_TiH  (NTQ), Q_CrH  (NTQ), Q_MnH  (NTQ),
     *  Q_FeH  (NTQ), Q_CoH  (NTQ), Q_NiH  (NTQ), Q_CuH  (NTQ),
     *  Q_ZnH  (NTQ), Q_GaH  (NTQ), Q_GeH  (NTQ), Q_AsH  (NTQ),
     *  Q_SeH  (NTQ), Q_HBr  (NTQ), Q_RbH  (NTQ), Q_SrH  (NTQ),
     *  Q_AgH  (NTQ), Q_CdH  (NTQ), Q_InH  (NTQ), Q_SnH  (NTQ),
     *  Q_SbH  (NTQ), Q_TeH  (NTQ), Q_HI   (NTQ), Q_CsH  (NTQ),
     *  Q_BaH  (NTQ), Q_YbH  (NTQ), Q_PtH  (NTQ), Q_AuH  (NTQ),
     *  Q_HgH  (NTQ), Q_TlH  (NTQ), Q_PbH  (NTQ), Q_BiH  (NTQ),
     *  Q_HeHp (NTQ), Q_BeHp (NTQ), Q_CHp  (NTQ), Q_NHp  (NTQ),
     *  Q_OHp  (NTQ), Q_HFp  (NTQ), Q_NeHp (NTQ), Q_MgHp (NTQ),
     *  Q_AlHp (NTQ), Q_SiHp (NTQ), Q_PHp  (NTQ), Q_SHp  (NTQ),
     *  Q_HClp (NTQ), Q_ZnHp (NTQ), Q_HBrp (NTQ), Q_CdHp (NTQ),
     *  Q_HgHp (NTQ), Q_CHm  (NTQ), Q_OHm  (NTQ), Q_SiHm (NTQ),
     *  Q_HSm  (NTQ), Q_CN   (NTQ), Q_CO   (NTQ), Q_CF   (NTQ),
     *  Q_SiC  (NTQ), Q_CP   (NTQ), Q_CS   (NTQ), Q_CCl  (NTQ),
     *  Q_CSe  (NTQ), Q_CBr  (NTQ), Q_RhC  (NTQ), Q_IrC  (NTQ),
     *  Q_PtC  (NTQ), Q_CNp  (NTQ), Q_COp  (NTQ), Q_CNm  (NTQ),
     *  Q_CSm  (NTQ), Q_BN   (NTQ), Q_NO   (NTQ), Q_NF   (NTQ),
     *  Q_AlN  (NTQ), Q_SiN  (NTQ), Q_PN   (NTQ), Q_NS   (NTQ),
     *  Q_NCl  (NTQ), Q_TiN  (NTQ), Q_AsN  (NTQ), Q_SeN  (NTQ),
     *  Q_ZrN  (NTQ), Q_NOp  (NTQ), Q_NSp  (NTQ), Q_LiO  (NTQ),
     *  Q_BeO  (NTQ), Q_BO   (NTQ), Q_FO   (NTQ), Q_NaO  (NTQ),
     *  Q_MgO  (NTQ), Q_AlO  (NTQ), Q_SiO  (NTQ), Q_PO   (NTQ),
     *  Q_SO   (NTQ), Q_ClO  (NTQ), Q_KO   (NTQ), Q_CaO  (NTQ),
     *  Q_ScO  (NTQ), Q_TiO  (NTQ), Q_VO   (NTQ), Q_CrO  (NTQ),
     *  Q_MnO  (NTQ), Q_FeO  (NTQ), Q_NiO  (NTQ), Q_CuO  (NTQ),
     *  Q_GaO  (NTQ), Q_GeO  (NTQ), Q_AsO  (NTQ), Q_SeO  (NTQ),
     *  Q_BrO  (NTQ), Q_RbO  (NTQ), Q_SrO  (NTQ), Q_YO   (NTQ),
     *  Q_ZrO  (NTQ), Q_NbO  (NTQ), Q_InO  (NTQ), Q_SnO  (NTQ),
     *  Q_SbO  (NTQ), Q_TeO  (NTQ), Q_IO   (NTQ), Q_BaO  (NTQ),
     *  Q_LaO  (NTQ), Q_TbO  (NTQ), Q_LuO  (NTQ), Q_HfO  (NTQ),
     *  Q_TaO  (NTQ), Q_WO   (NTQ), Q_PtO  (NTQ), Q_PbO  (NTQ),
     *  Q_BiO  (NTQ), Q_ThO  (NTQ), Q_BOp  (NTQ), Q_SiOp (NTQ),
     *  Q_POp  (NTQ), Q_SOp  (NTQ), Q_AsOp (NTQ), Q_TaOp (NTQ),
     *  Q_FeOm (NTQ), Q_LiF  (NTQ), Q_BeF  (NTQ), Q_BF   (NTQ),
     *  Q_NaF  (NTQ), Q_MgF  (NTQ), Q_AlF  (NTQ), Q_SiF  (NTQ),
     *  Q_PF   (NTQ), Q_SF   (NTQ), Q_KF   (NTQ), Q_CaF  (NTQ),
     *  Q_ScF  (NTQ), Q_MnF  (NTQ), Q_NiF  (NTQ), Q_CuF  (NTQ),
     *  Q_ZnF  (NTQ), Q_GaF  (NTQ), Q_GeF  (NTQ), Q_AsF  (NTQ),
     *  Q_SeF  (NTQ), Q_BrF  (NTQ), Q_RbF  (NTQ), Q_SrF  (NTQ),
     *  Q_YF   (NTQ), Q_AgF  (NTQ), Q_CdF  (NTQ), Q_InF  (NTQ),
     *  Q_SnF  (NTQ), Q_SbF  (NTQ), Q_IF   (NTQ), Q_CsF  (NTQ),
     *  Q_BaF  (NTQ), Q_LaF  (NTQ), Q_HoF  (NTQ), Q_YbF  (NTQ),
     *  Q_LuF  (NTQ), Q_HgF  (NTQ), Q_TlF  (NTQ), Q_PbF  (NTQ),
     *  Q_LiNa (NTQ), Q_AsP  (NTQ), Q_SbP  (NTQ), Q_BeS  (NTQ),
     *  Q_BS   (NTQ), Q_MgS  (NTQ), Q_AlS  (NTQ), Q_SiS  (NTQ),
     *  Q_PS   (NTQ), Q_CaS  (NTQ), Q_ScS  (NTQ), Q_TiS  (NTQ),
     *  Q_CrS  (NTQ), Q_CuS  (NTQ), Q_GeS  (NTQ), Q_AsS  (NTQ),
     *  Q_SeS  (NTQ), Q_SrS  (NTQ), Q_YS   (NTQ), Q_SnS  (NTQ),
     *  Q_TeS  (NTQ), Q_BaS  (NTQ), Q_LaS  (NTQ), Q_PbS  (NTQ),
     *  Q_BiS  (NTQ), Q_LiCl (NTQ), Q_BeCl (NTQ), Q_BCl  (NTQ),
     *  Q_NaCl (NTQ), Q_MgCl (NTQ), Q_AlCl (NTQ), Q_SiCl (NTQ),
     *  Q_PCl  (NTQ), Q_KCl  (NTQ), Q_CaCl (NTQ), Q_ScCl (NTQ),
     *  Q_MnCl (NTQ), Q_FeCl (NTQ), Q_CuCl (NTQ), Q_ZnCl (NTQ),
     *  Q_GaCl (NTQ), Q_GeCl (NTQ), Q_AsCl (NTQ), Q_SeCl (NTQ),
     *  Q_BrCl (NTQ), Q_RbCl (NTQ), Q_SrCl (NTQ), Q_YCl  (NTQ),
     *  Q_AgCl (NTQ), Q_CdCl (NTQ), Q_InCl (NTQ), Q_SnCl (NTQ),
     *  Q_SbCl (NTQ), Q_ICl  (NTQ), Q_CsCl (NTQ), Q_BaCl (NTQ),
     *  Q_YbCl (NTQ), Q_AuCl (NTQ), Q_HgCl (NTQ), Q_TlCl (NTQ),
     *  Q_PbCl (NTQ), Q_AlSe (NTQ), Q_SiSe (NTQ), Q_GeSe (NTQ),
     *  Q_KBr  (NTQ), Q_SiTe (NTQ), Q_GeTe (NTQ), Q_KI   (NTQ)
      REAL*8         TK_H2   (NTK),TK_Li2  (NTK),TK_B2   (NTK),
     * TK_C2   (NTK),TK_N2   (NTK),TK_O2   (NTK),TK_F2   (NTK),
     * TK_Na2  (NTK),TK_Mg2  (NTK),TK_Al2  (NTK),TK_Si2  (NTK),
     * TK_P2   (NTK),TK_S2   (NTK),TK_Cl2  (NTK),TK_K2   (NTK),
     * TK_Cu2  (NTK),TK_As2  (NTK),TK_Se2  (NTK),TK_Sb2  (NTK),
     * TK_Te2  (NTK),TK_I2   (NTK),TK_Cs2  (NTK),TK_H2p  (NTK),
     * TK_He2p (NTK),TK_C2p  (NTK),TK_N2p  (NTK),TK_O2p  (NTK),
     * TK_Ne2p (NTK),TK_P2p  (NTK),TK_S2p  (NTK),TK_H2m  (NTK),
     * TK_C2m  (NTK),TK_LiH  (NTK),TK_BeH  (NTK),TK_BH   (NTK),
     * TK_CH   (NTK),TK_NH   (NTK),TK_OH   (NTK),TK_HF   (NTK),
     * TK_NaH  (NTK),TK_MgH  (NTK),TK_AlH  (NTK),TK_SiH  (NTK),
     * TK_PH   (NTK),TK_HS   (NTK),TK_HCl  (NTK),TK_KH   (NTK),
     * TK_CaH  (NTK),TK_TiH  (NTK),TK_CrH  (NTK),TK_MnH  (NTK),
     * TK_FeH  (NTK),TK_CoH  (NTK),TK_NiH  (NTK),TK_CuH  (NTK),
     * TK_ZnH  (NTK),TK_GaH  (NTK),TK_GeH  (NTK),TK_AsH  (NTK),
     * TK_SeH  (NTK),TK_HBr  (NTK),TK_RbH  (NTK),TK_SrH  (NTK),
     * TK_AgH  (NTK),TK_CdH  (NTK),TK_InH  (NTK),TK_SnH  (NTK),
     * TK_SbH  (NTK),TK_TeH  (NTK),TK_HI   (NTK),TK_CsH  (NTK),
     * TK_BaH  (NTK),TK_YbH  (NTK),TK_PtH  (NTK),TK_AuH  (NTK),
     * TK_HgH  (NTK),TK_TlH  (NTK),TK_PbH  (NTK),TK_BiH  (NTK),
     * TK_HeHp (NTK),TK_BeHp (NTK),TK_CHp  (NTK),TK_NHp  (NTK),
     * TK_OHp  (NTK),TK_HFp  (NTK),TK_NeHp (NTK),TK_MgHp (NTK),
     * TK_AlHp (NTK),TK_SiHp (NTK),TK_PHp  (NTK),TK_SHp  (NTK),
     * TK_HClp (NTK),TK_ZnHp (NTK),TK_HBrp (NTK),TK_CdHp (NTK),
     * TK_HgHp (NTK),TK_CHm  (NTK),TK_OHm  (NTK),TK_SiHm (NTK),
     * TK_HSm  (NTK),TK_CN   (NTK),TK_CO   (NTK),TK_CF   (NTK),
     * TK_SiC  (NTK),TK_CP   (NTK),TK_CS   (NTK),TK_CCl  (NTK),
     * TK_CSe  (NTK),TK_CBr  (NTK),TK_RhC  (NTK),TK_IrC  (NTK),
     * TK_PtC  (NTK),TK_CNp  (NTK),TK_COp  (NTK),TK_CNm  (NTK),
     * TK_CSm  (NTK),TK_BN   (NTK),TK_NO   (NTK),TK_NF   (NTK),
     * TK_AlN  (NTK),TK_SiN  (NTK),TK_PN   (NTK),TK_NS   (NTK),
     * TK_NCl  (NTK),TK_TiN  (NTK),TK_AsN  (NTK),TK_SeN  (NTK),
     * TK_ZrN  (NTK),TK_NOp  (NTK),TK_NSp  (NTK),TK_LiO  (NTK),
     * TK_BeO  (NTK),TK_BO   (NTK),TK_FO   (NTK),TK_NaO  (NTK),
     * TK_MgO  (NTK),TK_AlO  (NTK),TK_SiO  (NTK),TK_PO   (NTK),
     * TK_SO   (NTK),TK_ClO  (NTK),TK_KO   (NTK),TK_CaO  (NTK),
     * TK_ScO  (NTK),TK_TiO  (NTK),TK_VO   (NTK),TK_CrO  (NTK),
     * TK_MnO  (NTK),TK_FeO  (NTK),TK_NiO  (NTK),TK_CuO  (NTK),
     * TK_GaO  (NTK),TK_GeO  (NTK),TK_AsO  (NTK),TK_SeO  (NTK),
     * TK_BrO  (NTK),TK_RbO  (NTK),TK_SrO  (NTK),TK_YO   (NTK),
     * TK_ZrO  (NTK),TK_NbO  (NTK),TK_InO  (NTK),TK_SnO  (NTK),
     * TK_SbO  (NTK),TK_TeO  (NTK),TK_IO   (NTK),TK_BaO  (NTK),
     * TK_LaO  (NTK),TK_TbO  (NTK),TK_LuO  (NTK),TK_HfO  (NTK),
     * TK_TaO  (NTK),TK_WO   (NTK),TK_PtO  (NTK),TK_PbO  (NTK),
     * TK_BiO  (NTK),TK_ThO  (NTK),TK_BOp  (NTK),TK_SiOp (NTK),
     * TK_POp  (NTK),TK_SOp  (NTK),TK_AsOp (NTK),TK_TaOp (NTK),
     * TK_FeOm (NTK),TK_LiF  (NTK),TK_BeF  (NTK),TK_BF   (NTK),
     * TK_NaF  (NTK),TK_MgF  (NTK),TK_AlF  (NTK),TK_SiF  (NTK),
     * TK_PF   (NTK),TK_SF   (NTK),TK_KF   (NTK),TK_CaF  (NTK),
     * TK_ScF  (NTK),TK_MnF  (NTK),TK_NiF  (NTK),TK_CuF  (NTK),
     * TK_ZnF  (NTK),TK_GaF  (NTK),TK_GeF  (NTK),TK_AsF  (NTK),
     * TK_SeF  (NTK),TK_BrF  (NTK),TK_RbF  (NTK),TK_SrF  (NTK),
     * TK_YF   (NTK),TK_AgF  (NTK),TK_CdF  (NTK),TK_InF  (NTK),
     * TK_SnF  (NTK),TK_SbF  (NTK),TK_IF   (NTK),TK_CsF  (NTK),
     * TK_BaF  (NTK),TK_LaF  (NTK),TK_HoF  (NTK),TK_YbF  (NTK),
     * TK_LuF  (NTK),TK_HgF  (NTK),TK_TlF  (NTK),TK_PbF  (NTK),
     * TK_LiNa (NTK),TK_AsP  (NTK),TK_SbP  (NTK),TK_BeS  (NTK),
     * TK_BS   (NTK),TK_MgS  (NTK),TK_AlS  (NTK),TK_SiS  (NTK),
     * TK_PS   (NTK),TK_CaS  (NTK),TK_ScS  (NTK),TK_TiS  (NTK),
     * TK_CrS  (NTK),TK_CuS  (NTK),TK_GeS  (NTK),TK_AsS  (NTK),
     * TK_SeS  (NTK),TK_SrS  (NTK),TK_YS   (NTK),TK_SnS  (NTK),
     * TK_TeS  (NTK),TK_BaS  (NTK),TK_LaS  (NTK),TK_PbS  (NTK),
     * TK_BiS  (NTK),TK_LiCl (NTK),TK_BeCl (NTK),TK_BCl  (NTK),
     * TK_NaCl (NTK),TK_MgCl (NTK),TK_AlCl (NTK),TK_SiCl (NTK),
     * TK_PCl  (NTK),TK_KCl  (NTK),TK_CaCl (NTK),TK_ScCl (NTK),
     * TK_MnCl (NTK),TK_FeCl (NTK),TK_CuCl (NTK),TK_ZnCl (NTK),
     * TK_GaCl (NTK),TK_GeCl (NTK),TK_AsCl (NTK),TK_SeCl (NTK),
     * TK_BrCl (NTK),TK_RbCl (NTK),TK_SrCl (NTK),TK_YCl  (NTK),
     * TK_AgCl (NTK),TK_CdCl (NTK),TK_InCl (NTK),TK_SnCl (NTK),
     * TK_SbCl (NTK),TK_ICl  (NTK),TK_CsCl (NTK),TK_BaCl (NTK),
     * TK_YbCl (NTK),TK_AuCl (NTK),TK_HgCl (NTK),TK_TlCl (NTK),
     * TK_PbCl (NTK),TK_AlSe (NTK),TK_SiSe (NTK),TK_GeSe (NTK),
     * TK_KBr  (NTK),TK_SiTe (NTK),TK_GeTe (NTK),TK_KI   (NTK)
      REAL*8          K_H2   (NTK), K_Li2  (NTK), K_B2   (NTK),
     *  K_C2   (NTK), K_N2   (NTK), K_O2   (NTK), K_F2   (NTK),
     *  K_Na2  (NTK), K_Mg2  (NTK), K_Al2  (NTK), K_Si2  (NTK),
     *  K_P2   (NTK), K_S2   (NTK), K_Cl2  (NTK), K_K2   (NTK),
     *  K_Cu2  (NTK), K_As2  (NTK), K_Se2  (NTK), K_Sb2  (NTK),
     *  K_Te2  (NTK), K_I2   (NTK), K_Cs2  (NTK), K_H2p  (NTK),
     *  K_He2p (NTK), K_C2p  (NTK), K_N2p  (NTK), K_O2p  (NTK),
     *  K_Ne2p (NTK), K_P2p  (NTK), K_S2p  (NTK), K_H2m  (NTK),
     *  K_C2m  (NTK), K_LiH  (NTK), K_BeH  (NTK), K_BH   (NTK),
     *  K_CH   (NTK), K_NH   (NTK), K_OH   (NTK), K_HF   (NTK),
     *  K_NaH  (NTK), K_MgH  (NTK), K_AlH  (NTK), K_SiH  (NTK),
     *  K_PH   (NTK), K_HS   (NTK), K_HCl  (NTK), K_KH   (NTK),
     *  K_CaH  (NTK), K_TiH  (NTK), K_CrH  (NTK), K_MnH  (NTK),
     *  K_FeH  (NTK), K_CoH  (NTK), K_NiH  (NTK), K_CuH  (NTK),
     *  K_ZnH  (NTK), K_GaH  (NTK), K_GeH  (NTK), K_AsH  (NTK),
     *  K_SeH  (NTK), K_HBr  (NTK), K_RbH  (NTK), K_SrH  (NTK),
     *  K_AgH  (NTK), K_CdH  (NTK), K_InH  (NTK), K_SnH  (NTK),
     *  K_SbH  (NTK), K_TeH  (NTK), K_HI   (NTK), K_CsH  (NTK),
     *  K_BaH  (NTK), K_YbH  (NTK), K_PtH  (NTK), K_AuH  (NTK),
     *  K_HgH  (NTK), K_TlH  (NTK), K_PbH  (NTK), K_BiH  (NTK),
     *  K_HeHp (NTK), K_BeHp (NTK), K_CHp  (NTK), K_NHp  (NTK),
     *  K_OHp  (NTK), K_HFp  (NTK), K_NeHp (NTK), K_MgHp (NTK),
     *  K_AlHp (NTK), K_SiHp (NTK), K_PHp  (NTK), K_SHp  (NTK),
     *  K_HClp (NTK), K_ZnHp (NTK), K_HBrp (NTK), K_CdHp (NTK),
     *  K_HgHp (NTK), K_CHm  (NTK), K_OHm  (NTK), K_SiHm (NTK),
     *  K_HSm  (NTK), K_CN   (NTK), K_CO   (NTK), K_CF   (NTK),
     *  K_SiC  (NTK), K_CP   (NTK), K_CS   (NTK), K_CCl  (NTK),
     *  K_CSe  (NTK), K_CBr  (NTK), K_RhC  (NTK), K_IrC  (NTK),
     *  K_PtC  (NTK), K_CNp  (NTK), K_COp  (NTK), K_CNm  (NTK),
     *  K_CSm  (NTK), K_BN   (NTK), K_NO   (NTK), K_NF   (NTK),
     *  K_AlN  (NTK), K_SiN  (NTK), K_PN   (NTK), K_NS   (NTK),
     *  K_NCl  (NTK), K_TiN  (NTK), K_AsN  (NTK), K_SeN  (NTK),
     *  K_ZrN  (NTK), K_NOp  (NTK), K_NSp  (NTK), K_LiO  (NTK),
     *  K_BeO  (NTK), K_BO   (NTK), K_FO   (NTK), K_NaO  (NTK),
     *  K_MgO  (NTK), K_AlO  (NTK), K_SiO  (NTK), K_PO   (NTK),
     *  K_SO   (NTK), K_ClO  (NTK), K_KO   (NTK), K_CaO  (NTK),
     *  K_ScO  (NTK), K_TiO  (NTK), K_VO   (NTK), K_CrO  (NTK),
     *  K_MnO  (NTK), K_FeO  (NTK), K_NiO  (NTK), K_CuO  (NTK),
     *  K_GaO  (NTK), K_GeO  (NTK), K_AsO  (NTK), K_SeO  (NTK),
     *  K_BrO  (NTK), K_RbO  (NTK), K_SrO  (NTK), K_YO   (NTK),
     *  K_ZrO  (NTK), K_NbO  (NTK), K_InO  (NTK), K_SnO  (NTK),
     *  K_SbO  (NTK), K_TeO  (NTK), K_IO   (NTK), K_BaO  (NTK),
     *  K_LaO  (NTK), K_TbO  (NTK), K_LuO  (NTK), K_HfO  (NTK),
     *  K_TaO  (NTK), K_WO   (NTK), K_PtO  (NTK), K_PbO  (NTK),
     *  K_BiO  (NTK), K_ThO  (NTK), K_BOp  (NTK), K_SiOp (NTK),
     *  K_POp  (NTK), K_SOp  (NTK), K_AsOp (NTK), K_TaOp (NTK),
     *  K_FeOm (NTK), K_LiF  (NTK), K_BeF  (NTK), K_BF   (NTK),
     *  K_NaF  (NTK), K_MgF  (NTK), K_AlF  (NTK), K_SiF  (NTK),
     *  K_PF   (NTK), K_SF   (NTK), K_KF   (NTK), K_CaF  (NTK),
     *  K_ScF  (NTK), K_MnF  (NTK), K_NiF  (NTK), K_CuF  (NTK),
     *  K_ZnF  (NTK), K_GaF  (NTK), K_GeF  (NTK), K_AsF  (NTK),
     *  K_SeF  (NTK), K_BrF  (NTK), K_RbF  (NTK), K_SrF  (NTK),
     *  K_YF   (NTK), K_AgF  (NTK), K_CdF  (NTK), K_InF  (NTK),
     *  K_SnF  (NTK), K_SbF  (NTK), K_IF   (NTK), K_CsF  (NTK),
     *  K_BaF  (NTK), K_LaF  (NTK), K_HoF  (NTK), K_YbF  (NTK),
     *  K_LuF  (NTK), K_HgF  (NTK), K_TlF  (NTK), K_PbF  (NTK),
     *  K_LiNa (NTK), K_AsP  (NTK), K_SbP  (NTK), K_BeS  (NTK),
     *  K_BS   (NTK), K_MgS  (NTK), K_AlS  (NTK), K_SiS  (NTK),
     *  K_PS   (NTK), K_CaS  (NTK), K_ScS  (NTK), K_TiS  (NTK),
     *  K_CrS  (NTK), K_CuS  (NTK), K_GeS  (NTK), K_AsS  (NTK),
     *  K_SeS  (NTK), K_SrS  (NTK), K_YS   (NTK), K_SnS  (NTK),
     *  K_TeS  (NTK), K_BaS  (NTK), K_LaS  (NTK), K_PbS  (NTK),
     *  K_BiS  (NTK), K_LiCl (NTK), K_BeCl (NTK), K_BCl  (NTK),
     *  K_NaCl (NTK), K_MgCl (NTK), K_AlCl (NTK), K_SiCl (NTK),
     *  K_PCl  (NTK), K_KCl  (NTK), K_CaCl (NTK), K_ScCl (NTK),
     *  K_MnCl (NTK), K_FeCl (NTK), K_CuCl (NTK), K_ZnCl (NTK),
     *  K_GaCl (NTK), K_GeCl (NTK), K_AsCl (NTK), K_SeCl (NTK),
     *  K_BrCl (NTK), K_RbCl (NTK), K_SrCl (NTK), K_YCl  (NTK),
     *  K_AgCl (NTK), K_CdCl (NTK), K_InCl (NTK), K_SnCl (NTK),
     *  K_SbCl (NTK), K_ICl  (NTK), K_CsCl (NTK), K_BaCl (NTK),
     *  K_YbCl (NTK), K_AuCl (NTK), K_HgCl (NTK), K_TlCl (NTK),
     *  K_PbCl (NTK), K_AlSe (NTK), K_SiSe (NTK), K_GeSe (NTK),
     *  K_KBr  (NTK), K_SiTe (NTK), K_GeTe (NTK), K_KI   (NTK)
      EQUIVALENCE (TQ(1,  1),TQ_H2   ),(TQ(1,  2),TQ_Li2  )
      EQUIVALENCE (TQ(1,  3),TQ_B2   ),(TQ(1,  4),TQ_C2   )
      EQUIVALENCE (TQ(1,  5),TQ_N2   ),(TQ(1,  6),TQ_O2   )
      EQUIVALENCE (TQ(1,  7),TQ_F2   ),(TQ(1,  8),TQ_Na2  )
      EQUIVALENCE (TQ(1,  9),TQ_Mg2  ),(TQ(1, 10),TQ_Al2  )
      EQUIVALENCE (TQ(1, 11),TQ_Si2  ),(TQ(1, 12),TQ_P2   )
      EQUIVALENCE (TQ(1, 13),TQ_S2   ),(TQ(1, 14),TQ_Cl2  )
      EQUIVALENCE (TQ(1, 15),TQ_K2   ),(TQ(1, 16),TQ_Cu2  )
      EQUIVALENCE (TQ(1, 17),TQ_As2  ),(TQ(1, 18),TQ_Se2  )
      EQUIVALENCE (TQ(1, 19),TQ_Sb2  ),(TQ(1, 20),TQ_Te2  )
      EQUIVALENCE (TQ(1, 21),TQ_I2   ),(TQ(1, 22),TQ_Cs2  )
      EQUIVALENCE (TQ(1, 23),TQ_H2p  ),(TQ(1, 24),TQ_He2p )
      EQUIVALENCE (TQ(1, 25),TQ_C2p  ),(TQ(1, 26),TQ_N2p  )
      EQUIVALENCE (TQ(1, 27),TQ_O2p  ),(TQ(1, 28),TQ_Ne2p )
      EQUIVALENCE (TQ(1, 29),TQ_P2p  ),(TQ(1, 30),TQ_S2p  )
      EQUIVALENCE (TQ(1, 31),TQ_H2m  ),(TQ(1, 32),TQ_C2m  )
      EQUIVALENCE (TQ(1, 33),TQ_LiH  ),(TQ(1, 34),TQ_BeH  )
      EQUIVALENCE (TQ(1, 35),TQ_BH   ),(TQ(1, 36),TQ_CH   )
      EQUIVALENCE (TQ(1, 37),TQ_NH   ),(TQ(1, 38),TQ_OH   )
      EQUIVALENCE (TQ(1, 39),TQ_HF   ),(TQ(1, 40),TQ_NaH  )
      EQUIVALENCE (TQ(1, 41),TQ_MgH  ),(TQ(1, 42),TQ_AlH  )
      EQUIVALENCE (TQ(1, 43),TQ_SiH  ),(TQ(1, 44),TQ_PH   )
      EQUIVALENCE (TQ(1, 45),TQ_HS   ),(TQ(1, 46),TQ_HCl  )
      EQUIVALENCE (TQ(1, 47),TQ_KH   ),(TQ(1, 48),TQ_CaH  )
      EQUIVALENCE (TQ(1, 49),TQ_TiH  ),(TQ(1, 50),TQ_CrH  )
      EQUIVALENCE (TQ(1, 51),TQ_MnH  ),(TQ(1, 52),TQ_FeH  )
      EQUIVALENCE (TQ(1, 53),TQ_CoH  ),(TQ(1, 54),TQ_NiH  )
      EQUIVALENCE (TQ(1, 55),TQ_CuH  ),(TQ(1, 56),TQ_ZnH  )
      EQUIVALENCE (TQ(1, 57),TQ_GaH  ),(TQ(1, 58),TQ_GeH  )
      EQUIVALENCE (TQ(1, 59),TQ_AsH  ),(TQ(1, 60),TQ_SeH  )
      EQUIVALENCE (TQ(1, 61),TQ_HBr  ),(TQ(1, 62),TQ_RbH  )
      EQUIVALENCE (TQ(1, 63),TQ_SrH  ),(TQ(1, 64),TQ_AgH  )
      EQUIVALENCE (TQ(1, 65),TQ_CdH  ),(TQ(1, 66),TQ_InH  )
      EQUIVALENCE (TQ(1, 67),TQ_SnH  ),(TQ(1, 68),TQ_SbH  )
      EQUIVALENCE (TQ(1, 69),TQ_TeH  ),(TQ(1, 70),TQ_HI   )
      EQUIVALENCE (TQ(1, 71),TQ_CsH  ),(TQ(1, 72),TQ_BaH  )
      EQUIVALENCE (TQ(1, 73),TQ_YbH  ),(TQ(1, 74),TQ_PtH  )
      EQUIVALENCE (TQ(1, 75),TQ_AuH  ),(TQ(1, 76),TQ_HgH  )
      EQUIVALENCE (TQ(1, 77),TQ_TlH  ),(TQ(1, 78),TQ_PbH  )
      EQUIVALENCE (TQ(1, 79),TQ_BiH  ),(TQ(1, 80),TQ_HeHp )
      EQUIVALENCE (TQ(1, 81),TQ_BeHp ),(TQ(1, 82),TQ_CHp  )
      EQUIVALENCE (TQ(1, 83),TQ_NHp  ),(TQ(1, 84),TQ_OHp  )
      EQUIVALENCE (TQ(1, 85),TQ_HFp  ),(TQ(1, 86),TQ_NeHp )
      EQUIVALENCE (TQ(1, 87),TQ_MgHp ),(TQ(1, 88),TQ_AlHp )
      EQUIVALENCE (TQ(1, 89),TQ_SiHp ),(TQ(1, 90),TQ_PHp  )
      EQUIVALENCE (TQ(1, 91),TQ_SHp  ),(TQ(1, 92),TQ_HClp )
      EQUIVALENCE (TQ(1, 93),TQ_ZnHp ),(TQ(1, 94),TQ_HBrp )
      EQUIVALENCE (TQ(1, 95),TQ_CdHp ),(TQ(1, 96),TQ_HgHp )
      EQUIVALENCE (TQ(1, 97),TQ_CHm  ),(TQ(1, 98),TQ_OHm  )
      EQUIVALENCE (TQ(1, 99),TQ_SiHm ),(TQ(1,100),TQ_HSm  )
      EQUIVALENCE (TQ(1,101),TQ_CN   ),(TQ(1,102),TQ_CO   )
      EQUIVALENCE (TQ(1,103),TQ_CF   ),(TQ(1,104),TQ_SiC  )
      EQUIVALENCE (TQ(1,105),TQ_CP   ),(TQ(1,106),TQ_CS   )
      EQUIVALENCE (TQ(1,107),TQ_CCl  ),(TQ(1,108),TQ_CSe  )
      EQUIVALENCE (TQ(1,109),TQ_CBr  ),(TQ(1,110),TQ_RhC  )
      EQUIVALENCE (TQ(1,111),TQ_IrC  ),(TQ(1,112),TQ_PtC  )
      EQUIVALENCE (TQ(1,113),TQ_CNp  ),(TQ(1,114),TQ_COp  )
      EQUIVALENCE (TQ(1,115),TQ_CNm  ),(TQ(1,116),TQ_CSm  )
      EQUIVALENCE (TQ(1,117),TQ_BN   ),(TQ(1,118),TQ_NO   )
      EQUIVALENCE (TQ(1,119),TQ_NF   ),(TQ(1,120),TQ_AlN  )
      EQUIVALENCE (TQ(1,121),TQ_SiN  ),(TQ(1,122),TQ_PN   )
      EQUIVALENCE (TQ(1,123),TQ_NS   ),(TQ(1,124),TQ_NCl  )
      EQUIVALENCE (TQ(1,125),TQ_TiN  ),(TQ(1,126),TQ_AsN  )
      EQUIVALENCE (TQ(1,127),TQ_SeN  ),(TQ(1,128),TQ_ZrN  )
      EQUIVALENCE (TQ(1,129),TQ_NOp  ),(TQ(1,130),TQ_NSp  )
      EQUIVALENCE (TQ(1,131),TQ_LiO  ),(TQ(1,132),TQ_BeO  )
      EQUIVALENCE (TQ(1,133),TQ_BO   ),(TQ(1,134),TQ_FO   )
      EQUIVALENCE (TQ(1,135),TQ_NaO  ),(TQ(1,136),TQ_MgO  )
      EQUIVALENCE (TQ(1,137),TQ_AlO  ),(TQ(1,138),TQ_SiO  )
      EQUIVALENCE (TQ(1,139),TQ_PO   ),(TQ(1,140),TQ_SO   )
      EQUIVALENCE (TQ(1,141),TQ_ClO  ),(TQ(1,142),TQ_KO   )
      EQUIVALENCE (TQ(1,143),TQ_CaO  ),(TQ(1,144),TQ_ScO  )
      EQUIVALENCE (TQ(1,145),TQ_TiO  ),(TQ(1,146),TQ_VO   )
      EQUIVALENCE (TQ(1,147),TQ_CrO  ),(TQ(1,148),TQ_MnO  )
      EQUIVALENCE (TQ(1,149),TQ_FeO  ),(TQ(1,150),TQ_NiO  )
      EQUIVALENCE (TQ(1,151),TQ_CuO  ),(TQ(1,152),TQ_GaO  )
      EQUIVALENCE (TQ(1,153),TQ_GeO  ),(TQ(1,154),TQ_AsO  )
      EQUIVALENCE (TQ(1,155),TQ_SeO  ),(TQ(1,156),TQ_BrO  )
      EQUIVALENCE (TQ(1,157),TQ_RbO  ),(TQ(1,158),TQ_SrO  )
      EQUIVALENCE (TQ(1,159),TQ_YO   ),(TQ(1,160),TQ_ZrO  )
      EQUIVALENCE (TQ(1,161),TQ_NbO  ),(TQ(1,162),TQ_InO  )
      EQUIVALENCE (TQ(1,163),TQ_SnO  ),(TQ(1,164),TQ_SbO  )
      EQUIVALENCE (TQ(1,165),TQ_TeO  ),(TQ(1,166),TQ_IO   )
      EQUIVALENCE (TQ(1,167),TQ_BaO  ),(TQ(1,168),TQ_LaO  )
      EQUIVALENCE (TQ(1,169),TQ_TbO  ),(TQ(1,170),TQ_LuO  )
      EQUIVALENCE (TQ(1,171),TQ_HfO  ),(TQ(1,172),TQ_TaO  )
      EQUIVALENCE (TQ(1,173),TQ_WO   ),(TQ(1,174),TQ_PtO  )
      EQUIVALENCE (TQ(1,175),TQ_PbO  ),(TQ(1,176),TQ_BiO  )
      EQUIVALENCE (TQ(1,177),TQ_ThO  ),(TQ(1,178),TQ_BOp  )
      EQUIVALENCE (TQ(1,179),TQ_SiOp ),(TQ(1,180),TQ_POp  )
      EQUIVALENCE (TQ(1,181),TQ_SOp  ),(TQ(1,182),TQ_AsOp )
      EQUIVALENCE (TQ(1,183),TQ_TaOp ),(TQ(1,184),TQ_FeOm )
      EQUIVALENCE (TQ(1,185),TQ_LiF  ),(TQ(1,186),TQ_BeF  )
      EQUIVALENCE (TQ(1,187),TQ_BF   ),(TQ(1,188),TQ_NaF  )
      EQUIVALENCE (TQ(1,189),TQ_MgF  ),(TQ(1,190),TQ_AlF  )
      EQUIVALENCE (TQ(1,191),TQ_SiF  ),(TQ(1,192),TQ_PF   )
      EQUIVALENCE (TQ(1,193),TQ_SF   ),(TQ(1,194),TQ_KF   )
      EQUIVALENCE (TQ(1,195),TQ_CaF  ),(TQ(1,196),TQ_ScF  )
      EQUIVALENCE (TQ(1,197),TQ_MnF  ),(TQ(1,198),TQ_NiF  )
      EQUIVALENCE (TQ(1,199),TQ_CuF  ),(TQ(1,200),TQ_ZnF  )
      EQUIVALENCE (TQ(1,201),TQ_GaF  ),(TQ(1,202),TQ_GeF  )
      EQUIVALENCE (TQ(1,203),TQ_AsF  ),(TQ(1,204),TQ_SeF  )
      EQUIVALENCE (TQ(1,205),TQ_BrF  ),(TQ(1,206),TQ_RbF  )
      EQUIVALENCE (TQ(1,207),TQ_SrF  ),(TQ(1,208),TQ_YF   )
      EQUIVALENCE (TQ(1,209),TQ_AgF  ),(TQ(1,210),TQ_CdF  )
      EQUIVALENCE (TQ(1,211),TQ_InF  ),(TQ(1,212),TQ_SnF  )
      EQUIVALENCE (TQ(1,213),TQ_SbF  ),(TQ(1,214),TQ_IF   )
      EQUIVALENCE (TQ(1,215),TQ_CsF  ),(TQ(1,216),TQ_BaF  )
      EQUIVALENCE (TQ(1,217),TQ_LaF  ),(TQ(1,218),TQ_HoF  )
      EQUIVALENCE (TQ(1,219),TQ_YbF  ),(TQ(1,220),TQ_LuF  )
      EQUIVALENCE (TQ(1,221),TQ_HgF  ),(TQ(1,222),TQ_TlF  )
      EQUIVALENCE (TQ(1,223),TQ_PbF  ),(TQ(1,224),TQ_LiNa )
      EQUIVALENCE (TQ(1,225),TQ_AsP  ),(TQ(1,226),TQ_SbP  )
      EQUIVALENCE (TQ(1,227),TQ_BeS  ),(TQ(1,228),TQ_BS   )
      EQUIVALENCE (TQ(1,229),TQ_MgS  ),(TQ(1,230),TQ_AlS  )
      EQUIVALENCE (TQ(1,231),TQ_SiS  ),(TQ(1,232),TQ_PS   )
      EQUIVALENCE (TQ(1,233),TQ_CaS  ),(TQ(1,234),TQ_ScS  )
      EQUIVALENCE (TQ(1,235),TQ_TiS  ),(TQ(1,236),TQ_CrS  )
      EQUIVALENCE (TQ(1,237),TQ_CuS  ),(TQ(1,238),TQ_GeS  )
      EQUIVALENCE (TQ(1,239),TQ_AsS  ),(TQ(1,240),TQ_SeS  )
      EQUIVALENCE (TQ(1,241),TQ_SrS  ),(TQ(1,242),TQ_YS   )
      EQUIVALENCE (TQ(1,243),TQ_SnS  ),(TQ(1,244),TQ_TeS  )
      EQUIVALENCE (TQ(1,245),TQ_BaS  ),(TQ(1,246),TQ_LaS  )
      EQUIVALENCE (TQ(1,247),TQ_PbS  ),(TQ(1,248),TQ_BiS  )
      EQUIVALENCE (TQ(1,249),TQ_LiCl ),(TQ(1,250),TQ_BeCl )
      EQUIVALENCE (TQ(1,251),TQ_BCl  ),(TQ(1,252),TQ_NaCl )
      EQUIVALENCE (TQ(1,253),TQ_MgCl ),(TQ(1,254),TQ_AlCl )
      EQUIVALENCE (TQ(1,255),TQ_SiCl ),(TQ(1,256),TQ_PCl  )
      EQUIVALENCE (TQ(1,257),TQ_KCl  ),(TQ(1,258),TQ_CaCl )
      EQUIVALENCE (TQ(1,259),TQ_ScCl ),(TQ(1,260),TQ_MnCl )
      EQUIVALENCE (TQ(1,261),TQ_FeCl ),(TQ(1,262),TQ_CuCl )
      EQUIVALENCE (TQ(1,263),TQ_ZnCl ),(TQ(1,264),TQ_GaCl )
      EQUIVALENCE (TQ(1,265),TQ_GeCl ),(TQ(1,266),TQ_AsCl )
      EQUIVALENCE (TQ(1,267),TQ_SeCl ),(TQ(1,268),TQ_BrCl )
      EQUIVALENCE (TQ(1,269),TQ_RbCl ),(TQ(1,270),TQ_SrCl )
      EQUIVALENCE (TQ(1,271),TQ_YCl  ),(TQ(1,272),TQ_AgCl )
      EQUIVALENCE (TQ(1,273),TQ_CdCl ),(TQ(1,274),TQ_InCl )
      EQUIVALENCE (TQ(1,275),TQ_SnCl ),(TQ(1,276),TQ_SbCl )
      EQUIVALENCE (TQ(1,277),TQ_ICl  ),(TQ(1,278),TQ_CsCl )
      EQUIVALENCE (TQ(1,279),TQ_BaCl ),(TQ(1,280),TQ_YbCl )
      EQUIVALENCE (TQ(1,281),TQ_AuCl ),(TQ(1,282),TQ_HgCl )
      EQUIVALENCE (TQ(1,283),TQ_TlCl ),(TQ(1,284),TQ_PbCl )
      EQUIVALENCE (TQ(1,285),TQ_AlSe ),(TQ(1,286),TQ_SiSe )
      EQUIVALENCE (TQ(1,287),TQ_GeSe ),(TQ(1,288),TQ_KBr  )
      EQUIVALENCE (TQ(1,289),TQ_SiTe ),(TQ(1,290),TQ_GeTe )
      EQUIVALENCE (TQ(1,291),TQ_KI   )
      EQUIVALENCE ( Q(1,  1), Q_H2   ),( Q(1,  2), Q_Li2  )
      EQUIVALENCE ( Q(1,  3), Q_B2   ),( Q(1,  4), Q_C2   )
      EQUIVALENCE ( Q(1,  5), Q_N2   ),( Q(1,  6), Q_O2   )
      EQUIVALENCE ( Q(1,  7), Q_F2   ),( Q(1,  8), Q_Na2  )
      EQUIVALENCE ( Q(1,  9), Q_Mg2  ),( Q(1, 10), Q_Al2  )
      EQUIVALENCE ( Q(1, 11), Q_Si2  ),( Q(1, 12), Q_P2   )
      EQUIVALENCE ( Q(1, 13), Q_S2   ),( Q(1, 14), Q_Cl2  )
      EQUIVALENCE ( Q(1, 15), Q_K2   ),( Q(1, 16), Q_Cu2  )
      EQUIVALENCE ( Q(1, 17), Q_As2  ),( Q(1, 18), Q_Se2  )
      EQUIVALENCE ( Q(1, 19), Q_Sb2  ),( Q(1, 20), Q_Te2  )
      EQUIVALENCE ( Q(1, 21), Q_I2   ),( Q(1, 22), Q_Cs2  )
      EQUIVALENCE ( Q(1, 23), Q_H2p  ),( Q(1, 24), Q_He2p )
      EQUIVALENCE ( Q(1, 25), Q_C2p  ),( Q(1, 26), Q_N2p  )
      EQUIVALENCE ( Q(1, 27), Q_O2p  ),( Q(1, 28), Q_Ne2p )
      EQUIVALENCE ( Q(1, 29), Q_P2p  ),( Q(1, 30), Q_S2p  )
      EQUIVALENCE ( Q(1, 31), Q_H2m  ),( Q(1, 32), Q_C2m  )
      EQUIVALENCE ( Q(1, 33), Q_LiH  ),( Q(1, 34), Q_BeH  )
      EQUIVALENCE ( Q(1, 35), Q_BH   ),( Q(1, 36), Q_CH   )
      EQUIVALENCE ( Q(1, 37), Q_NH   ),( Q(1, 38), Q_OH   )
      EQUIVALENCE ( Q(1, 39), Q_HF   ),( Q(1, 40), Q_NaH  )
      EQUIVALENCE ( Q(1, 41), Q_MgH  ),( Q(1, 42), Q_AlH  )
      EQUIVALENCE ( Q(1, 43), Q_SiH  ),( Q(1, 44), Q_PH   )
      EQUIVALENCE ( Q(1, 45), Q_HS   ),( Q(1, 46), Q_HCl  )
      EQUIVALENCE ( Q(1, 47), Q_KH   ),( Q(1, 48), Q_CaH  )
      EQUIVALENCE ( Q(1, 49), Q_TiH  ),( Q(1, 50), Q_CrH  )
      EQUIVALENCE ( Q(1, 51), Q_MnH  ),( Q(1, 52), Q_FeH  )
      EQUIVALENCE ( Q(1, 53), Q_CoH  ),( Q(1, 54), Q_NiH  )
      EQUIVALENCE ( Q(1, 55), Q_CuH  ),( Q(1, 56), Q_ZnH  )
      EQUIVALENCE ( Q(1, 57), Q_GaH  ),( Q(1, 58), Q_GeH  )
      EQUIVALENCE ( Q(1, 59), Q_AsH  ),( Q(1, 60), Q_SeH  )
      EQUIVALENCE ( Q(1, 61), Q_HBr  ),( Q(1, 62), Q_RbH  )
      EQUIVALENCE ( Q(1, 63), Q_SrH  ),( Q(1, 64), Q_AgH  )
      EQUIVALENCE ( Q(1, 65), Q_CdH  ),( Q(1, 66), Q_InH  )
      EQUIVALENCE ( Q(1, 67), Q_SnH  ),( Q(1, 68), Q_SbH  )
      EQUIVALENCE ( Q(1, 69), Q_TeH  ),( Q(1, 70), Q_HI   )
      EQUIVALENCE ( Q(1, 71), Q_CsH  ),( Q(1, 72), Q_BaH  )
      EQUIVALENCE ( Q(1, 73), Q_YbH  ),( Q(1, 74), Q_PtH  )
      EQUIVALENCE ( Q(1, 75), Q_AuH  ),( Q(1, 76), Q_HgH  )
      EQUIVALENCE ( Q(1, 77), Q_TlH  ),( Q(1, 78), Q_PbH  )
      EQUIVALENCE ( Q(1, 79), Q_BiH  ),( Q(1, 80), Q_HeHp )
      EQUIVALENCE ( Q(1, 81), Q_BeHp ),( Q(1, 82), Q_CHp  )
      EQUIVALENCE ( Q(1, 83), Q_NHp  ),( Q(1, 84), Q_OHp  )
      EQUIVALENCE ( Q(1, 85), Q_HFp  ),( Q(1, 86), Q_NeHp )
      EQUIVALENCE ( Q(1, 87), Q_MgHp ),( Q(1, 88), Q_AlHp )
      EQUIVALENCE ( Q(1, 89), Q_SiHp ),( Q(1, 90), Q_PHp  )
      EQUIVALENCE ( Q(1, 91), Q_SHp  ),( Q(1, 92), Q_HClp )
      EQUIVALENCE ( Q(1, 93), Q_ZnHp ),( Q(1, 94), Q_HBrp )
      EQUIVALENCE ( Q(1, 95), Q_CdHp ),( Q(1, 96), Q_HgHp )
      EQUIVALENCE ( Q(1, 97), Q_CHm  ),( Q(1, 98), Q_OHm  )
      EQUIVALENCE ( Q(1, 99), Q_SiHm ),( Q(1,100), Q_HSm  )
      EQUIVALENCE ( Q(1,101), Q_CN   ),( Q(1,102), Q_CO   )
      EQUIVALENCE ( Q(1,103), Q_CF   ),( Q(1,104), Q_SiC  )
      EQUIVALENCE ( Q(1,105), Q_CP   ),( Q(1,106), Q_CS   )
      EQUIVALENCE ( Q(1,107), Q_CCl  ),( Q(1,108), Q_CSe  )
      EQUIVALENCE ( Q(1,109), Q_CBr  ),( Q(1,110), Q_RhC  )
      EQUIVALENCE ( Q(1,111), Q_IrC  ),( Q(1,112), Q_PtC  )
      EQUIVALENCE ( Q(1,113), Q_CNp  ),( Q(1,114), Q_COp  )
      EQUIVALENCE ( Q(1,115), Q_CNm  ),( Q(1,116), Q_CSm  )
      EQUIVALENCE ( Q(1,117), Q_BN   ),( Q(1,118), Q_NO   )
      EQUIVALENCE ( Q(1,119), Q_NF   ),( Q(1,120), Q_AlN  )
      EQUIVALENCE ( Q(1,121), Q_SiN  ),( Q(1,122), Q_PN   )
      EQUIVALENCE ( Q(1,123), Q_NS   ),( Q(1,124), Q_NCl  )
      EQUIVALENCE ( Q(1,125), Q_TiN  ),( Q(1,126), Q_AsN  )
      EQUIVALENCE ( Q(1,127), Q_SeN  ),( Q(1,128), Q_ZrN  )
      EQUIVALENCE ( Q(1,129), Q_NOp  ),( Q(1,130), Q_NSp  )
      EQUIVALENCE ( Q(1,131), Q_LiO  ),( Q(1,132), Q_BeO  )
      EQUIVALENCE ( Q(1,133), Q_BO   ),( Q(1,134), Q_FO   )
      EQUIVALENCE ( Q(1,135), Q_NaO  ),( Q(1,136), Q_MgO  )
      EQUIVALENCE ( Q(1,137), Q_AlO  ),( Q(1,138), Q_SiO  )
      EQUIVALENCE ( Q(1,139), Q_PO   ),( Q(1,140), Q_SO   )
      EQUIVALENCE ( Q(1,141), Q_ClO  ),( Q(1,142), Q_KO   )
      EQUIVALENCE ( Q(1,143), Q_CaO  ),( Q(1,144), Q_ScO  )
      EQUIVALENCE ( Q(1,145), Q_TiO  ),( Q(1,146), Q_VO   )
      EQUIVALENCE ( Q(1,147), Q_CrO  ),( Q(1,148), Q_MnO  )
      EQUIVALENCE ( Q(1,149), Q_FeO  ),( Q(1,150), Q_NiO  )
      EQUIVALENCE ( Q(1,151), Q_CuO  ),( Q(1,152), Q_GaO  )
      EQUIVALENCE ( Q(1,153), Q_GeO  ),( Q(1,154), Q_AsO  )
      EQUIVALENCE ( Q(1,155), Q_SeO  ),( Q(1,156), Q_BrO  )
      EQUIVALENCE ( Q(1,157), Q_RbO  ),( Q(1,158), Q_SrO  )
      EQUIVALENCE ( Q(1,159), Q_YO   ),( Q(1,160), Q_ZrO  )
      EQUIVALENCE ( Q(1,161), Q_NbO  ),( Q(1,162), Q_InO  )
      EQUIVALENCE ( Q(1,163), Q_SnO  ),( Q(1,164), Q_SbO  )
      EQUIVALENCE ( Q(1,165), Q_TeO  ),( Q(1,166), Q_IO   )
      EQUIVALENCE ( Q(1,167), Q_BaO  ),( Q(1,168), Q_LaO  )
      EQUIVALENCE ( Q(1,169), Q_TbO  ),( Q(1,170), Q_LuO  )
      EQUIVALENCE ( Q(1,171), Q_HfO  ),( Q(1,172), Q_TaO  )
      EQUIVALENCE ( Q(1,173), Q_WO   ),( Q(1,174), Q_PtO  )
      EQUIVALENCE ( Q(1,175), Q_PbO  ),( Q(1,176), Q_BiO  )
      EQUIVALENCE ( Q(1,177), Q_ThO  ),( Q(1,178), Q_BOp  )
      EQUIVALENCE ( Q(1,179), Q_SiOp ),( Q(1,180), Q_POp  )
      EQUIVALENCE ( Q(1,181), Q_SOp  ),( Q(1,182), Q_AsOp )
      EQUIVALENCE ( Q(1,183), Q_TaOp ),( Q(1,184), Q_FeOm )
      EQUIVALENCE ( Q(1,185), Q_LiF  ),( Q(1,186), Q_BeF  )
      EQUIVALENCE ( Q(1,187), Q_BF   ),( Q(1,188), Q_NaF  )
      EQUIVALENCE ( Q(1,189), Q_MgF  ),( Q(1,190), Q_AlF  )
      EQUIVALENCE ( Q(1,191), Q_SiF  ),( Q(1,192), Q_PF   )
      EQUIVALENCE ( Q(1,193), Q_SF   ),( Q(1,194), Q_KF   )
      EQUIVALENCE ( Q(1,195), Q_CaF  ),( Q(1,196), Q_ScF  )
      EQUIVALENCE ( Q(1,197), Q_MnF  ),( Q(1,198), Q_NiF  )
      EQUIVALENCE ( Q(1,199), Q_CuF  ),( Q(1,200), Q_ZnF  )
      EQUIVALENCE ( Q(1,201), Q_GaF  ),( Q(1,202), Q_GeF  )
      EQUIVALENCE ( Q(1,203), Q_AsF  ),( Q(1,204), Q_SeF  )
      EQUIVALENCE ( Q(1,205), Q_BrF  ),( Q(1,206), Q_RbF  )
      EQUIVALENCE ( Q(1,207), Q_SrF  ),( Q(1,208), Q_YF   )
      EQUIVALENCE ( Q(1,209), Q_AgF  ),( Q(1,210), Q_CdF  )
      EQUIVALENCE ( Q(1,211), Q_InF  ),( Q(1,212), Q_SnF  )
      EQUIVALENCE ( Q(1,213), Q_SbF  ),( Q(1,214), Q_IF   )
      EQUIVALENCE ( Q(1,215), Q_CsF  ),( Q(1,216), Q_BaF  )
      EQUIVALENCE ( Q(1,217), Q_LaF  ),( Q(1,218), Q_HoF  )
      EQUIVALENCE ( Q(1,219), Q_YbF  ),( Q(1,220), Q_LuF  )
      EQUIVALENCE ( Q(1,221), Q_HgF  ),( Q(1,222), Q_TlF  )
      EQUIVALENCE ( Q(1,223), Q_PbF  ),( Q(1,224), Q_LiNa )
      EQUIVALENCE ( Q(1,225), Q_AsP  ),( Q(1,226), Q_SbP  )
      EQUIVALENCE ( Q(1,227), Q_BeS  ),( Q(1,228), Q_BS   )
      EQUIVALENCE ( Q(1,229), Q_MgS  ),( Q(1,230), Q_AlS  )
      EQUIVALENCE ( Q(1,231), Q_SiS  ),( Q(1,232), Q_PS   )
      EQUIVALENCE ( Q(1,233), Q_CaS  ),( Q(1,234), Q_ScS  )
      EQUIVALENCE ( Q(1,235), Q_TiS  ),( Q(1,236), Q_CrS  )
      EQUIVALENCE ( Q(1,237), Q_CuS  ),( Q(1,238), Q_GeS  )
      EQUIVALENCE ( Q(1,239), Q_AsS  ),( Q(1,240), Q_SeS  )
      EQUIVALENCE ( Q(1,241), Q_SrS  ),( Q(1,242), Q_YS   )
      EQUIVALENCE ( Q(1,243), Q_SnS  ),( Q(1,244), Q_TeS  )
      EQUIVALENCE ( Q(1,245), Q_BaS  ),( Q(1,246), Q_LaS  )
      EQUIVALENCE ( Q(1,247), Q_PbS  ),( Q(1,248), Q_BiS  )
      EQUIVALENCE ( Q(1,249), Q_LiCl ),( Q(1,250), Q_BeCl )
      EQUIVALENCE ( Q(1,251), Q_BCl  ),( Q(1,252), Q_NaCl )
      EQUIVALENCE ( Q(1,253), Q_MgCl ),( Q(1,254), Q_AlCl )
      EQUIVALENCE ( Q(1,255), Q_SiCl ),( Q(1,256), Q_PCl  )
      EQUIVALENCE ( Q(1,257), Q_KCl  ),( Q(1,258), Q_CaCl )
      EQUIVALENCE ( Q(1,259), Q_ScCl ),( Q(1,260), Q_MnCl )
      EQUIVALENCE ( Q(1,261), Q_FeCl ),( Q(1,262), Q_CuCl )
      EQUIVALENCE ( Q(1,263), Q_ZnCl ),( Q(1,264), Q_GaCl )
      EQUIVALENCE ( Q(1,265), Q_GeCl ),( Q(1,266), Q_AsCl )
      EQUIVALENCE ( Q(1,267), Q_SeCl ),( Q(1,268), Q_BrCl )
      EQUIVALENCE ( Q(1,269), Q_RbCl ),( Q(1,270), Q_SrCl )
      EQUIVALENCE ( Q(1,271), Q_YCl  ),( Q(1,272), Q_AgCl )
      EQUIVALENCE ( Q(1,273), Q_CdCl ),( Q(1,274), Q_InCl )
      EQUIVALENCE ( Q(1,275), Q_SnCl ),( Q(1,276), Q_SbCl )
      EQUIVALENCE ( Q(1,277), Q_ICl  ),( Q(1,278), Q_CsCl )
      EQUIVALENCE ( Q(1,279), Q_BaCl ),( Q(1,280), Q_YbCl )
      EQUIVALENCE ( Q(1,281), Q_AuCl ),( Q(1,282), Q_HgCl )
      EQUIVALENCE ( Q(1,283), Q_TlCl ),( Q(1,284), Q_PbCl )
      EQUIVALENCE ( Q(1,285), Q_AlSe ),( Q(1,286), Q_SiSe )
      EQUIVALENCE ( Q(1,287), Q_GeSe ),( Q(1,288), Q_KBr  )
      EQUIVALENCE ( Q(1,289), Q_SiTe ),( Q(1,290), Q_GeTe )
      EQUIVALENCE ( Q(1,291), Q_KI   )
      EQUIVALENCE (TK(1,  1),TK_H2   ),(TK(1,  2),TK_Li2  )
      EQUIVALENCE (TK(1,  3),TK_B2   ),(TK(1,  4),TK_C2   )
      EQUIVALENCE (TK(1,  5),TK_N2   ),(TK(1,  6),TK_O2   )
      EQUIVALENCE (TK(1,  7),TK_F2   ),(TK(1,  8),TK_Na2  )
      EQUIVALENCE (TK(1,  9),TK_Mg2  ),(TK(1, 10),TK_Al2  )
      EQUIVALENCE (TK(1, 11),TK_Si2  ),(TK(1, 12),TK_P2   )
      EQUIVALENCE (TK(1, 13),TK_S2   ),(TK(1, 14),TK_Cl2  )
      EQUIVALENCE (TK(1, 15),TK_K2   ),(TK(1, 16),TK_Cu2  )
      EQUIVALENCE (TK(1, 17),TK_As2  ),(TK(1, 18),TK_Se2  )
      EQUIVALENCE (TK(1, 19),TK_Sb2  ),(TK(1, 20),TK_Te2  )
      EQUIVALENCE (TK(1, 21),TK_I2   ),(TK(1, 22),TK_Cs2  )
      EQUIVALENCE (TK(1, 23),TK_H2p  ),(TK(1, 24),TK_He2p )
      EQUIVALENCE (TK(1, 25),TK_C2p  ),(TK(1, 26),TK_N2p  )
      EQUIVALENCE (TK(1, 27),TK_O2p  ),(TK(1, 28),TK_Ne2p )
      EQUIVALENCE (TK(1, 29),TK_P2p  ),(TK(1, 30),TK_S2p  )
      EQUIVALENCE (TK(1, 31),TK_H2m  ),(TK(1, 32),TK_C2m  )
      EQUIVALENCE (TK(1, 33),TK_LiH  ),(TK(1, 34),TK_BeH  )
      EQUIVALENCE (TK(1, 35),TK_BH   ),(TK(1, 36),TK_CH   )
      EQUIVALENCE (TK(1, 37),TK_NH   ),(TK(1, 38),TK_OH   )
      EQUIVALENCE (TK(1, 39),TK_HF   ),(TK(1, 40),TK_NaH  )
      EQUIVALENCE (TK(1, 41),TK_MgH  ),(TK(1, 42),TK_AlH  )
      EQUIVALENCE (TK(1, 43),TK_SiH  ),(TK(1, 44),TK_PH   )
      EQUIVALENCE (TK(1, 45),TK_HS   ),(TK(1, 46),TK_HCl  )
      EQUIVALENCE (TK(1, 47),TK_KH   ),(TK(1, 48),TK_CaH  )
      EQUIVALENCE (TK(1, 49),TK_TiH  ),(TK(1, 50),TK_CrH  )
      EQUIVALENCE (TK(1, 51),TK_MnH  ),(TK(1, 52),TK_FeH  )
      EQUIVALENCE (TK(1, 53),TK_CoH  ),(TK(1, 54),TK_NiH  )
      EQUIVALENCE (TK(1, 55),TK_CuH  ),(TK(1, 56),TK_ZnH  )
      EQUIVALENCE (TK(1, 57),TK_GaH  ),(TK(1, 58),TK_GeH  )
      EQUIVALENCE (TK(1, 59),TK_AsH  ),(TK(1, 60),TK_SeH  )
      EQUIVALENCE (TK(1, 61),TK_HBr  ),(TK(1, 62),TK_RbH  )
      EQUIVALENCE (TK(1, 63),TK_SrH  ),(TK(1, 64),TK_AgH  )
      EQUIVALENCE (TK(1, 65),TK_CdH  ),(TK(1, 66),TK_InH  )
      EQUIVALENCE (TK(1, 67),TK_SnH  ),(TK(1, 68),TK_SbH  )
      EQUIVALENCE (TK(1, 69),TK_TeH  ),(TK(1, 70),TK_HI   )
      EQUIVALENCE (TK(1, 71),TK_CsH  ),(TK(1, 72),TK_BaH  )
      EQUIVALENCE (TK(1, 73),TK_YbH  ),(TK(1, 74),TK_PtH  )
      EQUIVALENCE (TK(1, 75),TK_AuH  ),(TK(1, 76),TK_HgH  )
      EQUIVALENCE (TK(1, 77),TK_TlH  ),(TK(1, 78),TK_PbH  )
      EQUIVALENCE (TK(1, 79),TK_BiH  ),(TK(1, 80),TK_HeHp )
      EQUIVALENCE (TK(1, 81),TK_BeHp ),(TK(1, 82),TK_CHp  )
      EQUIVALENCE (TK(1, 83),TK_NHp  ),(TK(1, 84),TK_OHp  )
      EQUIVALENCE (TK(1, 85),TK_HFp  ),(TK(1, 86),TK_NeHp )
      EQUIVALENCE (TK(1, 87),TK_MgHp ),(TK(1, 88),TK_AlHp )
      EQUIVALENCE (TK(1, 89),TK_SiHp ),(TK(1, 90),TK_PHp  )
      EQUIVALENCE (TK(1, 91),TK_SHp  ),(TK(1, 92),TK_HClp )
      EQUIVALENCE (TK(1, 93),TK_ZnHp ),(TK(1, 94),TK_HBrp )
      EQUIVALENCE (TK(1, 95),TK_CdHp ),(TK(1, 96),TK_HgHp )
      EQUIVALENCE (TK(1, 97),TK_CHm  ),(TK(1, 98),TK_OHm  )
      EQUIVALENCE (TK(1, 99),TK_SiHm ),(TK(1,100),TK_HSm  )
      EQUIVALENCE (TK(1,101),TK_CN   ),(TK(1,102),TK_CO   )
      EQUIVALENCE (TK(1,103),TK_CF   ),(TK(1,104),TK_SiC  )
      EQUIVALENCE (TK(1,105),TK_CP   ),(TK(1,106),TK_CS   )
      EQUIVALENCE (TK(1,107),TK_CCl  ),(TK(1,108),TK_CSe  )
      EQUIVALENCE (TK(1,109),TK_CBr  ),(TK(1,110),TK_RhC  )
      EQUIVALENCE (TK(1,111),TK_IrC  ),(TK(1,112),TK_PtC  )
      EQUIVALENCE (TK(1,113),TK_CNp  ),(TK(1,114),TK_COp  )
      EQUIVALENCE (TK(1,115),TK_CNm  ),(TK(1,116),TK_CSm  )
      EQUIVALENCE (TK(1,117),TK_BN   ),(TK(1,118),TK_NO   )
      EQUIVALENCE (TK(1,119),TK_NF   ),(TK(1,120),TK_AlN  )
      EQUIVALENCE (TK(1,121),TK_SiN  ),(TK(1,122),TK_PN   )
      EQUIVALENCE (TK(1,123),TK_NS   ),(TK(1,124),TK_NCl  )
      EQUIVALENCE (TK(1,125),TK_TiN  ),(TK(1,126),TK_AsN  )
      EQUIVALENCE (TK(1,127),TK_SeN  ),(TK(1,128),TK_ZrN  )
      EQUIVALENCE (TK(1,129),TK_NOp  ),(TK(1,130),TK_NSp  )
      EQUIVALENCE (TK(1,131),TK_LiO  ),(TK(1,132),TK_BeO  )
      EQUIVALENCE (TK(1,133),TK_BO   ),(TK(1,134),TK_FO   )
      EQUIVALENCE (TK(1,135),TK_NaO  ),(TK(1,136),TK_MgO  )
      EQUIVALENCE (TK(1,137),TK_AlO  ),(TK(1,138),TK_SiO  )
      EQUIVALENCE (TK(1,139),TK_PO   ),(TK(1,140),TK_SO   )
      EQUIVALENCE (TK(1,141),TK_ClO  ),(TK(1,142),TK_KO   )
      EQUIVALENCE (TK(1,143),TK_CaO  ),(TK(1,144),TK_ScO  )
      EQUIVALENCE (TK(1,145),TK_TiO  ),(TK(1,146),TK_VO   )
      EQUIVALENCE (TK(1,147),TK_CrO  ),(TK(1,148),TK_MnO  )
      EQUIVALENCE (TK(1,149),TK_FeO  ),(TK(1,150),TK_NiO  )
      EQUIVALENCE (TK(1,151),TK_CuO  ),(TK(1,152),TK_GaO  )
      EQUIVALENCE (TK(1,153),TK_GeO  ),(TK(1,154),TK_AsO  )
      EQUIVALENCE (TK(1,155),TK_SeO  ),(TK(1,156),TK_BrO  )
      EQUIVALENCE (TK(1,157),TK_RbO  ),(TK(1,158),TK_SrO  )
      EQUIVALENCE (TK(1,159),TK_YO   ),(TK(1,160),TK_ZrO  )
      EQUIVALENCE (TK(1,161),TK_NbO  ),(TK(1,162),TK_InO  )
      EQUIVALENCE (TK(1,163),TK_SnO  ),(TK(1,164),TK_SbO  )
      EQUIVALENCE (TK(1,165),TK_TeO  ),(TK(1,166),TK_IO   )
      EQUIVALENCE (TK(1,167),TK_BaO  ),(TK(1,168),TK_LaO  )
      EQUIVALENCE (TK(1,169),TK_TbO  ),(TK(1,170),TK_LuO  )
      EQUIVALENCE (TK(1,171),TK_HfO  ),(TK(1,172),TK_TaO  )
      EQUIVALENCE (TK(1,173),TK_WO   ),(TK(1,174),TK_PtO  )
      EQUIVALENCE (TK(1,175),TK_PbO  ),(TK(1,176),TK_BiO  )
      EQUIVALENCE (TK(1,177),TK_ThO  ),(TK(1,178),TK_BOp  )
      EQUIVALENCE (TK(1,179),TK_SiOp ),(TK(1,180),TK_POp  )
      EQUIVALENCE (TK(1,181),TK_SOp  ),(TK(1,182),TK_AsOp )
      EQUIVALENCE (TK(1,183),TK_TaOp ),(TK(1,184),TK_FeOm )
      EQUIVALENCE (TK(1,185),TK_LiF  ),(TK(1,186),TK_BeF  )
      EQUIVALENCE (TK(1,187),TK_BF   ),(TK(1,188),TK_NaF  )
      EQUIVALENCE (TK(1,189),TK_MgF  ),(TK(1,190),TK_AlF  )
      EQUIVALENCE (TK(1,191),TK_SiF  ),(TK(1,192),TK_PF   )
      EQUIVALENCE (TK(1,193),TK_SF   ),(TK(1,194),TK_KF   )
      EQUIVALENCE (TK(1,195),TK_CaF  ),(TK(1,196),TK_ScF  )
      EQUIVALENCE (TK(1,197),TK_MnF  ),(TK(1,198),TK_NiF  )
      EQUIVALENCE (TK(1,199),TK_CuF  ),(TK(1,200),TK_ZnF  )
      EQUIVALENCE (TK(1,201),TK_GaF  ),(TK(1,202),TK_GeF  )
      EQUIVALENCE (TK(1,203),TK_AsF  ),(TK(1,204),TK_SeF  )
      EQUIVALENCE (TK(1,205),TK_BrF  ),(TK(1,206),TK_RbF  )
      EQUIVALENCE (TK(1,207),TK_SrF  ),(TK(1,208),TK_YF   )
      EQUIVALENCE (TK(1,209),TK_AgF  ),(TK(1,210),TK_CdF  )
      EQUIVALENCE (TK(1,211),TK_InF  ),(TK(1,212),TK_SnF  )
      EQUIVALENCE (TK(1,213),TK_SbF  ),(TK(1,214),TK_IF   )
      EQUIVALENCE (TK(1,215),TK_CsF  ),(TK(1,216),TK_BaF  )
      EQUIVALENCE (TK(1,217),TK_LaF  ),(TK(1,218),TK_HoF  )
      EQUIVALENCE (TK(1,219),TK_YbF  ),(TK(1,220),TK_LuF  )
      EQUIVALENCE (TK(1,221),TK_HgF  ),(TK(1,222),TK_TlF  )
      EQUIVALENCE (TK(1,223),TK_PbF  ),(TK(1,224),TK_LiNa )
      EQUIVALENCE (TK(1,225),TK_AsP  ),(TK(1,226),TK_SbP  )
      EQUIVALENCE (TK(1,227),TK_BeS  ),(TK(1,228),TK_BS   )
      EQUIVALENCE (TK(1,229),TK_MgS  ),(TK(1,230),TK_AlS  )
      EQUIVALENCE (TK(1,231),TK_SiS  ),(TK(1,232),TK_PS   )
      EQUIVALENCE (TK(1,233),TK_CaS  ),(TK(1,234),TK_ScS  )
      EQUIVALENCE (TK(1,235),TK_TiS  ),(TK(1,236),TK_CrS  )
      EQUIVALENCE (TK(1,237),TK_CuS  ),(TK(1,238),TK_GeS  )
      EQUIVALENCE (TK(1,239),TK_AsS  ),(TK(1,240),TK_SeS  )
      EQUIVALENCE (TK(1,241),TK_SrS  ),(TK(1,242),TK_YS   )
      EQUIVALENCE (TK(1,243),TK_SnS  ),(TK(1,244),TK_TeS  )
      EQUIVALENCE (TK(1,245),TK_BaS  ),(TK(1,246),TK_LaS  )
      EQUIVALENCE (TK(1,247),TK_PbS  ),(TK(1,248),TK_BiS  )
      EQUIVALENCE (TK(1,249),TK_LiCl ),(TK(1,250),TK_BeCl )
      EQUIVALENCE (TK(1,251),TK_BCl  ),(TK(1,252),TK_NaCl )
      EQUIVALENCE (TK(1,253),TK_MgCl ),(TK(1,254),TK_AlCl )
      EQUIVALENCE (TK(1,255),TK_SiCl ),(TK(1,256),TK_PCl  )
      EQUIVALENCE (TK(1,257),TK_KCl  ),(TK(1,258),TK_CaCl )
      EQUIVALENCE (TK(1,259),TK_ScCl ),(TK(1,260),TK_MnCl )
      EQUIVALENCE (TK(1,261),TK_FeCl ),(TK(1,262),TK_CuCl )
      EQUIVALENCE (TK(1,263),TK_ZnCl ),(TK(1,264),TK_GaCl )
      EQUIVALENCE (TK(1,265),TK_GeCl ),(TK(1,266),TK_AsCl )
      EQUIVALENCE (TK(1,267),TK_SeCl ),(TK(1,268),TK_BrCl )
      EQUIVALENCE (TK(1,269),TK_RbCl ),(TK(1,270),TK_SrCl )
      EQUIVALENCE (TK(1,271),TK_YCl  ),(TK(1,272),TK_AgCl )
      EQUIVALENCE (TK(1,273),TK_CdCl ),(TK(1,274),TK_InCl )
      EQUIVALENCE (TK(1,275),TK_SnCl ),(TK(1,276),TK_SbCl )
      EQUIVALENCE (TK(1,277),TK_ICl  ),(TK(1,278),TK_CsCl )
      EQUIVALENCE (TK(1,279),TK_BaCl ),(TK(1,280),TK_YbCl )
      EQUIVALENCE (TK(1,281),TK_AuCl ),(TK(1,282),TK_HgCl )
      EQUIVALENCE (TK(1,283),TK_TlCl ),(TK(1,284),TK_PbCl )
      EQUIVALENCE (TK(1,285),TK_AlSe ),(TK(1,286),TK_SiSe )
      EQUIVALENCE (TK(1,287),TK_GeSe ),(TK(1,288),TK_KBr  )
      EQUIVALENCE (TK(1,289),TK_SiTe ),(TK(1,290),TK_GeTe )
      EQUIVALENCE (TK(1,291),TK_KI   )
      EQUIVALENCE ( K(1,  1), K_H2   ),( K(1,  2), K_Li2  )
      EQUIVALENCE ( K(1,  3), K_B2   ),( K(1,  4), K_C2   )
      EQUIVALENCE ( K(1,  5), K_N2   ),( K(1,  6), K_O2   )
      EQUIVALENCE ( K(1,  7), K_F2   ),( K(1,  8), K_Na2  )
      EQUIVALENCE ( K(1,  9), K_Mg2  ),( K(1, 10), K_Al2  )
      EQUIVALENCE ( K(1, 11), K_Si2  ),( K(1, 12), K_P2   )
      EQUIVALENCE ( K(1, 13), K_S2   ),( K(1, 14), K_Cl2  )
      EQUIVALENCE ( K(1, 15), K_K2   ),( K(1, 16), K_Cu2  )
      EQUIVALENCE ( K(1, 17), K_As2  ),( K(1, 18), K_Se2  )
      EQUIVALENCE ( K(1, 19), K_Sb2  ),( K(1, 20), K_Te2  )
      EQUIVALENCE ( K(1, 21), K_I2   ),( K(1, 22), K_Cs2  )
      EQUIVALENCE ( K(1, 23), K_H2p  ),( K(1, 24), K_He2p )
      EQUIVALENCE ( K(1, 25), K_C2p  ),( K(1, 26), K_N2p  )
      EQUIVALENCE ( K(1, 27), K_O2p  ),( K(1, 28), K_Ne2p )
      EQUIVALENCE ( K(1, 29), K_P2p  ),( K(1, 30), K_S2p  )
      EQUIVALENCE ( K(1, 31), K_H2m  ),( K(1, 32), K_C2m  )
      EQUIVALENCE ( K(1, 33), K_LiH  ),( K(1, 34), K_BeH  )
      EQUIVALENCE ( K(1, 35), K_BH   ),( K(1, 36), K_CH   )
      EQUIVALENCE ( K(1, 37), K_NH   ),( K(1, 38), K_OH   )
      EQUIVALENCE ( K(1, 39), K_HF   ),( K(1, 40), K_NaH  )
      EQUIVALENCE ( K(1, 41), K_MgH  ),( K(1, 42), K_AlH  )
      EQUIVALENCE ( K(1, 43), K_SiH  ),( K(1, 44), K_PH   )
      EQUIVALENCE ( K(1, 45), K_HS   ),( K(1, 46), K_HCl  )
      EQUIVALENCE ( K(1, 47), K_KH   ),( K(1, 48), K_CaH  )
      EQUIVALENCE ( K(1, 49), K_TiH  ),( K(1, 50), K_CrH  )
      EQUIVALENCE ( K(1, 51), K_MnH  ),( K(1, 52), K_FeH  )
      EQUIVALENCE ( K(1, 53), K_CoH  ),( K(1, 54), K_NiH  )
      EQUIVALENCE ( K(1, 55), K_CuH  ),( K(1, 56), K_ZnH  )
      EQUIVALENCE ( K(1, 57), K_GaH  ),( K(1, 58), K_GeH  )
      EQUIVALENCE ( K(1, 59), K_AsH  ),( K(1, 60), K_SeH  )
      EQUIVALENCE ( K(1, 61), K_HBr  ),( K(1, 62), K_RbH  )
      EQUIVALENCE ( K(1, 63), K_SrH  ),( K(1, 64), K_AgH  )
      EQUIVALENCE ( K(1, 65), K_CdH  ),( K(1, 66), K_InH  )
      EQUIVALENCE ( K(1, 67), K_SnH  ),( K(1, 68), K_SbH  )
      EQUIVALENCE ( K(1, 69), K_TeH  ),( K(1, 70), K_HI   )
      EQUIVALENCE ( K(1, 71), K_CsH  ),( K(1, 72), K_BaH  )
      EQUIVALENCE ( K(1, 73), K_YbH  ),( K(1, 74), K_PtH  )
      EQUIVALENCE ( K(1, 75), K_AuH  ),( K(1, 76), K_HgH  )
      EQUIVALENCE ( K(1, 77), K_TlH  ),( K(1, 78), K_PbH  )
      EQUIVALENCE ( K(1, 79), K_BiH  ),( K(1, 80), K_HeHp )
      EQUIVALENCE ( K(1, 81), K_BeHp ),( K(1, 82), K_CHp  )
      EQUIVALENCE ( K(1, 83), K_NHp  ),( K(1, 84), K_OHp  )
      EQUIVALENCE ( K(1, 85), K_HFp  ),( K(1, 86), K_NeHp )
      EQUIVALENCE ( K(1, 87), K_MgHp ),( K(1, 88), K_AlHp )
      EQUIVALENCE ( K(1, 89), K_SiHp ),( K(1, 90), K_PHp  )
      EQUIVALENCE ( K(1, 91), K_SHp  ),( K(1, 92), K_HClp )
      EQUIVALENCE ( K(1, 93), K_ZnHp ),( K(1, 94), K_HBrp )
      EQUIVALENCE ( K(1, 95), K_CdHp ),( K(1, 96), K_HgHp )
      EQUIVALENCE ( K(1, 97), K_CHm  ),( K(1, 98), K_OHm  )
      EQUIVALENCE ( K(1, 99), K_SiHm ),( K(1,100), K_HSm  )
      EQUIVALENCE ( K(1,101), K_CN   ),( K(1,102), K_CO   )
      EQUIVALENCE ( K(1,103), K_CF   ),( K(1,104), K_SiC  )
      EQUIVALENCE ( K(1,105), K_CP   ),( K(1,106), K_CS   )
      EQUIVALENCE ( K(1,107), K_CCl  ),( K(1,108), K_CSe  )
      EQUIVALENCE ( K(1,109), K_CBr  ),( K(1,110), K_RhC  )
      EQUIVALENCE ( K(1,111), K_IrC  ),( K(1,112), K_PtC  )
      EQUIVALENCE ( K(1,113), K_CNp  ),( K(1,114), K_COp  )
      EQUIVALENCE ( K(1,115), K_CNm  ),( K(1,116), K_CSm  )
      EQUIVALENCE ( K(1,117), K_BN   ),( K(1,118), K_NO   )
      EQUIVALENCE ( K(1,119), K_NF   ),( K(1,120), K_AlN  )
      EQUIVALENCE ( K(1,121), K_SiN  ),( K(1,122), K_PN   )
      EQUIVALENCE ( K(1,123), K_NS   ),( K(1,124), K_NCl  )
      EQUIVALENCE ( K(1,125), K_TiN  ),( K(1,126), K_AsN  )
      EQUIVALENCE ( K(1,127), K_SeN  ),( K(1,128), K_ZrN  )
      EQUIVALENCE ( K(1,129), K_NOp  ),( K(1,130), K_NSp  )
      EQUIVALENCE ( K(1,131), K_LiO  ),( K(1,132), K_BeO  )
      EQUIVALENCE ( K(1,133), K_BO   ),( K(1,134), K_FO   )
      EQUIVALENCE ( K(1,135), K_NaO  ),( K(1,136), K_MgO  )
      EQUIVALENCE ( K(1,137), K_AlO  ),( K(1,138), K_SiO  )
      EQUIVALENCE ( K(1,139), K_PO   ),( K(1,140), K_SO   )
      EQUIVALENCE ( K(1,141), K_ClO  ),( K(1,142), K_KO   )
      EQUIVALENCE ( K(1,143), K_CaO  ),( K(1,144), K_ScO  )
      EQUIVALENCE ( K(1,145), K_TiO  ),( K(1,146), K_VO   )
      EQUIVALENCE ( K(1,147), K_CrO  ),( K(1,148), K_MnO  )
      EQUIVALENCE ( K(1,149), K_FeO  ),( K(1,150), K_NiO  )
      EQUIVALENCE ( K(1,151), K_CuO  ),( K(1,152), K_GaO  )
      EQUIVALENCE ( K(1,153), K_GeO  ),( K(1,154), K_AsO  )
      EQUIVALENCE ( K(1,155), K_SeO  ),( K(1,156), K_BrO  )
      EQUIVALENCE ( K(1,157), K_RbO  ),( K(1,158), K_SrO  )
      EQUIVALENCE ( K(1,159), K_YO   ),( K(1,160), K_ZrO  )
      EQUIVALENCE ( K(1,161), K_NbO  ),( K(1,162), K_InO  )
      EQUIVALENCE ( K(1,163), K_SnO  ),( K(1,164), K_SbO  )
      EQUIVALENCE ( K(1,165), K_TeO  ),( K(1,166), K_IO   )
      EQUIVALENCE ( K(1,167), K_BaO  ),( K(1,168), K_LaO  )
      EQUIVALENCE ( K(1,169), K_TbO  ),( K(1,170), K_LuO  )
      EQUIVALENCE ( K(1,171), K_HfO  ),( K(1,172), K_TaO  )
      EQUIVALENCE ( K(1,173), K_WO   ),( K(1,174), K_PtO  )
      EQUIVALENCE ( K(1,175), K_PbO  ),( K(1,176), K_BiO  )
      EQUIVALENCE ( K(1,177), K_ThO  ),( K(1,178), K_BOp  )
      EQUIVALENCE ( K(1,179), K_SiOp ),( K(1,180), K_POp  )
      EQUIVALENCE ( K(1,181), K_SOp  ),( K(1,182), K_AsOp )
      EQUIVALENCE ( K(1,183), K_TaOp ),( K(1,184), K_FeOm )
      EQUIVALENCE ( K(1,185), K_LiF  ),( K(1,186), K_BeF  )
      EQUIVALENCE ( K(1,187), K_BF   ),( K(1,188), K_NaF  )
      EQUIVALENCE ( K(1,189), K_MgF  ),( K(1,190), K_AlF  )
      EQUIVALENCE ( K(1,191), K_SiF  ),( K(1,192), K_PF   )
      EQUIVALENCE ( K(1,193), K_SF   ),( K(1,194), K_KF   )
      EQUIVALENCE ( K(1,195), K_CaF  ),( K(1,196), K_ScF  )
      EQUIVALENCE ( K(1,197), K_MnF  ),( K(1,198), K_NiF  )
      EQUIVALENCE ( K(1,199), K_CuF  ),( K(1,200), K_ZnF  )
      EQUIVALENCE ( K(1,201), K_GaF  ),( K(1,202), K_GeF  )
      EQUIVALENCE ( K(1,203), K_AsF  ),( K(1,204), K_SeF  )
      EQUIVALENCE ( K(1,205), K_BrF  ),( K(1,206), K_RbF  )
      EQUIVALENCE ( K(1,207), K_SrF  ),( K(1,208), K_YF   )
      EQUIVALENCE ( K(1,209), K_AgF  ),( K(1,210), K_CdF  )
      EQUIVALENCE ( K(1,211), K_InF  ),( K(1,212), K_SnF  )
      EQUIVALENCE ( K(1,213), K_SbF  ),( K(1,214), K_IF   )
      EQUIVALENCE ( K(1,215), K_CsF  ),( K(1,216), K_BaF  )
      EQUIVALENCE ( K(1,217), K_LaF  ),( K(1,218), K_HoF  )
      EQUIVALENCE ( K(1,219), K_YbF  ),( K(1,220), K_LuF  )
      EQUIVALENCE ( K(1,221), K_HgF  ),( K(1,222), K_TlF  )
      EQUIVALENCE ( K(1,223), K_PbF  ),( K(1,224), K_LiNa )
      EQUIVALENCE ( K(1,225), K_AsP  ),( K(1,226), K_SbP  )
      EQUIVALENCE ( K(1,227), K_BeS  ),( K(1,228), K_BS   )
      EQUIVALENCE ( K(1,229), K_MgS  ),( K(1,230), K_AlS  )
      EQUIVALENCE ( K(1,231), K_SiS  ),( K(1,232), K_PS   )
      EQUIVALENCE ( K(1,233), K_CaS  ),( K(1,234), K_ScS  )
      EQUIVALENCE ( K(1,235), K_TiS  ),( K(1,236), K_CrS  )
      EQUIVALENCE ( K(1,237), K_CuS  ),( K(1,238), K_GeS  )
      EQUIVALENCE ( K(1,239), K_AsS  ),( K(1,240), K_SeS  )
      EQUIVALENCE ( K(1,241), K_SrS  ),( K(1,242), K_YS   )
      EQUIVALENCE ( K(1,243), K_SnS  ),( K(1,244), K_TeS  )
      EQUIVALENCE ( K(1,245), K_BaS  ),( K(1,246), K_LaS  )
      EQUIVALENCE ( K(1,247), K_PbS  ),( K(1,248), K_BiS  )
      EQUIVALENCE ( K(1,249), K_LiCl ),( K(1,250), K_BeCl )
      EQUIVALENCE ( K(1,251), K_BCl  ),( K(1,252), K_NaCl )
      EQUIVALENCE ( K(1,253), K_MgCl ),( K(1,254), K_AlCl )
      EQUIVALENCE ( K(1,255), K_SiCl ),( K(1,256), K_PCl  )
      EQUIVALENCE ( K(1,257), K_KCl  ),( K(1,258), K_CaCl )
      EQUIVALENCE ( K(1,259), K_ScCl ),( K(1,260), K_MnCl )
      EQUIVALENCE ( K(1,261), K_FeCl ),( K(1,262), K_CuCl )
      EQUIVALENCE ( K(1,263), K_ZnCl ),( K(1,264), K_GaCl )
      EQUIVALENCE ( K(1,265), K_GeCl ),( K(1,266), K_AsCl )
      EQUIVALENCE ( K(1,267), K_SeCl ),( K(1,268), K_BrCl )
      EQUIVALENCE ( K(1,269), K_RbCl ),( K(1,270), K_SrCl )
      EQUIVALENCE ( K(1,271), K_YCl  ),( K(1,272), K_AgCl )
      EQUIVALENCE ( K(1,273), K_CdCl ),( K(1,274), K_InCl )
      EQUIVALENCE ( K(1,275), K_SnCl ),( K(1,276), K_SbCl )
      EQUIVALENCE ( K(1,277), K_ICl  ),( K(1,278), K_CsCl )
      EQUIVALENCE ( K(1,279), K_BaCl ),( K(1,280), K_YbCl )
      EQUIVALENCE ( K(1,281), K_AuCl ),( K(1,282), K_HgCl )
      EQUIVALENCE ( K(1,283), K_TlCl ),( K(1,284), K_PbCl )
      EQUIVALENCE ( K(1,285), K_AlSe ),( K(1,286), K_SiSe )
      EQUIVALENCE ( K(1,287), K_GeSe ),( K(1,288), K_KBr  )
      EQUIVALENCE ( K(1,289), K_SiTe ),( K(1,290), K_GeTe )
      EQUIVALENCE ( K(1,291), K_KI   )
C
      SAVE
C
      DATA SPLIST/
     * 'H2    ','Li2   ','B2    ','C2    ','N2    ','O2    ',
     * 'F2    ','Na2   ','Mg2   ','Al2   ','Si2   ','P2    ',
     * 'S2    ','Cl2   ','K2    ','Cu2   ','As2   ','Se2   ',
     * 'Sb2   ','Te2   ','I2    ','Cs2   ','H2+   ','He2+  ',
     * 'C2+   ','N2+   ','O2+   ','Ne2+  ','P2+   ','S2+   ',
     * 'H2-   ','C2-   ','LiH   ','BeH   ','BH    ','CH    ',
     * 'NH    ','OH    ','HF    ','NaH   ','MgH   ','AlH   ',
     * 'SiH   ','PH    ','HS    ','HCl   ','KH    ','CaH   ',
     * 'TiH   ','CrH   ','MnH   ','FeH   ','CoH   ','NiH   ',
     * 'CuH   ','ZnH   ','GaH   ','GeH   ','AsH   ','SeH   ',
     * 'HBr   ','RbH   ','SrH   ','AgH   ','CdH   ','InH   ',
     * 'SnH   ','SbH   ','TeH   ','HI    ','CsH   ','BaH   ',
     * 'YbH   ','PtH   ','AuH   ','HgH   ','TlH   ','PbH   ',
     * 'BiH   ','HeH+  ','BeH+  ','CH+   ','NH+   ','OH+   ',
     * 'HF+   ','NeH+  ','MgH+  ','AlH+  ','SiH+  ','PH+   ',
     * 'SH+   ','HCl+  ','ZnH+  ','HBr+  ','CdH+  ','HgH+  ',
     * 'CH-   ','OH-   ','SiH-  ','HS-   ','CN    ','CO    ',
     * 'CF    ','SiC   ','CP    ','CS    ','CCl   ','CSe   ',
     * 'CBr   ','RhC   ','IrC   ','PtC   ','CN+   ','CO+   ',
     * 'CN-   ','CS-   ','BN    ','NO    ','NF    ','AlN   ',
     * 'SiN   ','PN    ','NS    ','NCl   ','TiN   ','AsN   ',
     * 'SeN   ','ZrN   ','NO+   ','NS+   ','LiO   ','BeO   ',
     * 'BO    ','FO    ','NaO   ','MgO   ','AlO   ','SiO   ',
     * 'PO    ','SO    ','ClO   ','KO    ','CaO   ','ScO   ',
     * 'TiO   ','VO    ','CrO   ','MnO   ','FeO   ','NiO   ',
     * 'CuO   ','GaO   ','GeO   ','AsO   ','SeO   ','BrO   ',
     * 'RbO   ','SrO   ','YO    ','ZrO   ','NbO   ','InO   ',
     * 'SnO   ','SbO   ','TeO   ','IO    ','BaO   ','LaO   ',
     * 'TbO   ','LuO   ','HfO   ','TaO   ','WO    ','PtO   ',
     * 'PbO   ','BiO   ','ThO   ','BO+   ','SiO+  ','PO+   ',
     * 'SO+   ','AsO+  ','TaO+  ','FeO-  ','LiF   ','BeF   ',
     * 'BF    ','NaF   ','MgF   ','AlF   ','SiF   ','PF    ',
     * 'SF    ','KF    ','CaF   ','ScF   ','MnF   ','NiF   ',
     * 'CuF   ','ZnF   ','GaF   ','GeF   ','AsF   ','SeF   ',
     * 'BrF   ','RbF   ','SrF   ','YF    ','AgF   ','CdF   ',
     * 'InF   ','SnF   ','SbF   ','IF    ','CsF   ','BaF   ',
     * 'LaF   ','HoF   ','YbF   ','LuF   ','HgF   ','TlF   ',
     * 'PbF   ','LiNa  ','AsP   ','SbP   ','BeS   ','BS    ',
     * 'MgS   ','AlS   ','SiS   ','PS    ','CaS   ','ScS   ',
     * 'TiS   ','CrS   ','CuS   ','GeS   ','AsS   ','SeS   ',
     * 'SrS   ','YS    ','SnS   ','TeS   ','BaS   ','LaS   ',
     * 'PbS   ','BiS   ','LiCl  ','BeCl  ','BCl   ','NaCl  ',
     * 'MgCl  ','AlCl  ','SiCl  ','PCl   ','KCl   ','CaCl  ',
     * 'ScCl  ','MnCl  ','FeCl  ','CuCl  ','ZnCl  ','GaCl  ',
     * 'GeCl  ','AsCl  ','SeCl  ','BrCl  ','RbCl  ','SrCl  ',
     * 'YCl   ','AgCl  ','CdCl  ','InCl  ','SnCl  ','SbCl  ',
     * 'ICl   ','CsCl  ','BaCl  ','YbCl  ','AuCl  ','HgCl  ',
     * 'TlCl  ','PbCl  ','AlSe  ','SiSe  ','GeSe  ','KBr   ',
     * 'SiTe  ','GeTe  ','KI    '/
C
C Molecular partition functions
C
      DATA TQ_H2/                                                       071215
     1  0.699999789529, 1.154399984723, 1.254400107040, 1.360900114337, H2   
     2  1.444399939753, 1.589900165736, 1.659399847868, 1.732099970393, H2   
     3  1.858100092901, 1.992199822742, 2.126600019701, 2.267200392272, H2   
     4  2.407999801624, 2.560700362716, 2.767900424219, 2.961800082309, H2   
     5  3.136900278615, 3.368899870082, 3.527099940601, 3.681300206711, H2   
     6  3.871899935943, 3.944299699681, 4.000000000000,      4*0.0D+00/ H2   
      DATA  Q_H2/                                                       071215
     1 -6.02059991D-01,-6.02034590D-01,-6.01763800D-01,-5.99733997D-01, H2   
     2 -5.93597533D-01,-5.55742686D-01,-5.17412974D-01,-4.61269255D-01, H2   
     3 -3.34511304D-01,-1.84441164D-01,-4.08644555D-02, 9.59462410D-02, H2   
     4  2.25207019D-01, 3.65408102D-01, 5.61513561D-01, 7.51273629D-01, H2   
     5  9.31227296D-01, 1.20108375D+00, 1.42119260D+00, 1.67404299D+00, H2   
     6  2.03451619D+00, 2.17855848D+00, 2.28974774D+00,      4*0.0D+00/ H2   
      DATA TQ_Li2/                                                      071215
     1  0.699999789529, 0.756300150457, 0.843900054273, 1.077999972882, Li2  
     2  1.360300129731, 1.659499845207, 1.803000107677, 1.939700032897, Li2  
     3  2.274600127090, 2.439699606371, 2.593900243314, 2.957299978718, Li2  
     4  3.138200308408, 3.295600094378, 3.421300118962, 3.559300390627, Li2  
     5  3.679900294873, 3.861699703316, 4.000000000000,      8*0.0D+00/ Li2  
      DATA  Q_Li2/                                                      071215
     1  4.43603986D-01, 4.96492474D-01, 5.79600491D-01, 8.05386198D-01, Li2  
     2  1.08224422D+00, 1.37872346D+00, 1.52176868D+00, 1.65945749D+00, Li2  
     3  2.02582494D+00, 2.23854702D+00, 2.46181616D+00, 3.07491463D+00, Li2  
     4  3.42138459D+00, 3.75180949D+00, 4.04337502D+00, 4.38604918D+00, Li2  
     5  4.68884851D+00, 5.13616014D+00, 5.47403023D+00,      8*0.0D+00/ Li2  
      DATA TQ_B2/                                                       071215
     1  0.699999789529, 0.727699981235, 0.770299989902, 0.881900068180, B2   
     2  1.210000014136, 1.543500048264, 1.912199884965, 2.134500194441, B2   
     3  2.336500081247, 2.482999603020, 2.639200321380, 2.832399996440, B2   
     4  3.083600014869, 3.330899937219, 3.475199879001, 3.599500369651, B2   
     5  3.824399791343, 3.929800318889, 3.972500091209, 4.000000000000, B2   
     6       7*0.0D+00/                                                 B2   
      DATA  Q_B2/                                                       071215
     1  6.88622485D-01, 7.12867698D-01, 7.50707679D-01, 8.52245200D-01, B2   
     2  1.16264599D+00, 1.48780626D+00, 1.85238814D+00, 2.07346438D+00, B2   
     3  2.27523860D+00, 2.42427341D+00, 2.59169384D+00, 2.82262601D+00, B2   
     4  3.17582347D+00, 3.58046076D+00, 3.83877238D+00, 4.07335751D+00, B2   
     5  4.53049798D+00, 4.76152545D+00, 4.85778181D+00, 4.92038505D+00, B2   
     6       7*0.0D+00/                                                 B2   
      DATA TQ_C2/                                                       071215
     1  0.699999789529, 0.718899824635, 0.749799982595, 0.824300057081, C2   
     2  0.928799873439, 1.031800173386, 1.291399978054, 1.672600007894, C2   
     3  1.854900019495, 2.021400361632, 2.141500216060, 2.262900298345, C2   
     4  2.367099817682, 2.433500039087, 2.498299947815, 2.581899969961, C2   
     5  2.675300201907, 2.793200050520, 2.910399840958, 3.112499710784, C2   
     6  3.279399786984, 3.416099990304, 3.540299964532, 3.837300085731, C2   
     7  3.936099883681, 3.974899918760, 4.000000000000/                 C2   
      DATA  Q_C2/                                                       071215
     1  8.67908445D-02, 9.78644669D-02, 1.17474831D-01, 1.71399820D-01, C2   
     2  2.58168689D-01, 3.50405485D-01, 5.93847753D-01, 9.63742050D-01, C2   
     3  1.14332994D+00, 1.30884692D+00, 1.43248606D+00, 1.57144975D+00, C2   
     4  1.71387281D+00, 1.81914170D+00, 1.93219633D+00, 2.08929475D+00, C2   
     5  2.27259569D+00, 2.50542488D+00, 2.73286731D+00, 3.11693686D+00, C2   
     6  3.43555388D+00, 3.70531252D+00, 3.96200836D+00, 4.63046056D+00, C2   
     7  4.87060205D+00, 4.96776543D+00, 5.03158435D+00/                 C2   
      DATA TQ_N2/                                                       071215
     1  0.699999789529, 0.726499952497, 0.767800040614, 0.868299917667, N2   
     2  1.016100067701, 1.163899929845, 1.422500014137, 1.681199865672, N2   
     3  2.124299971461, 2.392400165093, 2.620699890222, 2.813400172329, N2   
     4  3.030699761808, 3.158199790218, 3.289199955475, 3.558000363497, N2   
     5  3.719100231678, 3.863399742159, 3.947699776449, 3.979999552305, N2   
     6  4.000000000000,      6*0.0D+00/                                 N2   
      DATA  Q_N2/                                                       071215
     1  4.00389476D-02, 5.86348388D-02, 8.91017256D-02, 1.69600625D-01, N2   
     2  2.98849020D-01, 4.34744803D-01, 6.80462886D-01, 9.32110588D-01, N2   
     3  1.36974884D+00, 1.63651438D+00, 1.86442008D+00, 2.05937318D+00, N2   
     4  2.29400914D+00, 2.44705567D+00, 2.62011872D+00, 3.02703117D+00, N2   
     5  3.30062538D+00, 3.56173752D+00, 3.72162019D+00, 3.78486931D+00, N2   
     6  3.82480726D+00,      6*0.0D+00/                                 N2   
      DATA TQ_O2/                                                       071215
     1  0.699999789529, 0.730000036317, 0.776600143547, 0.891699906491, O2   
     2  1.164199936638, 1.369899883439, 1.599099916238, 1.800500166492, O2   
     3  2.001100025259, 2.255800137213, 2.466600175078, 2.645999905632, O2   
     4  2.864299750209, 3.041200005385, 3.288699943982, 3.533199802608, O2   
     5  3.702199857799, 3.832299976826, 3.933300085488, 3.973799997799, O2   
     6  4.000000000000,      6*0.0D+00/                                 O2   
      DATA  Q_O2/                                                       071215
     1  9.70305183D-01, 9.75289822D-01, 9.85082684D-01, 1.02107012D+00, O2   
     2  1.16838265D+00, 1.31985107D+00, 1.51237459D+00, 1.69422217D+00, O2   
     3  1.88259687D+00, 2.12813366D+00, 2.33482853D+00, 2.51460753D+00, O2   
     4  2.74978400D+00, 2.96689870D+00, 3.32291916D+00, 3.73976750D+00, O2   
     5  4.07141694D+00, 4.35364685D+00, 4.59199282D+00, 4.69343229D+00, O2   
     6  4.76107718D+00,      6*0.0D+00/                                 O2   
      DATA TQ_F2/                                                       071215
     1  0.699999789529, 0.778700194762, 0.895799998258, 1.222899862578, F2   
     2  1.544300064674, 1.883800029456, 2.099400375556, 2.302400251608, F2   
     3  2.591700196167, 2.705099934369, 2.831999987435, 3.214300123350, F2   
     4  3.356199921013, 3.500099992405, 3.589700152706, 3.675900199664, F2   
     5  3.861099689607, 3.945999738065, 3.979199609788, 4.000000000000, F2   
     6       7*0.0D+00/                                                 F2   
      DATA  Q_F2/                                                       071215
     1  3.31920651D-01, 4.04554960D-01, 5.14348743D-01, 8.28998832D-01, F2   
     2  1.14466755D+00, 1.48142017D+00, 1.69622117D+00, 1.89957614D+00, F2   
     3  2.20493631D+00, 2.33838239D+00, 2.50177103D+00, 3.09149017D+00, F2   
     4  3.34637577D+00, 3.62923548D+00, 3.82139788D+00, 4.01808718D+00, F2   
     5  4.46249947D+00, 4.66476510D+00, 4.74206501D+00, 4.78982234D+00, F2   
     6       7*0.0D+00/                                                 F2   
      DATA TQ_Na2/                                                      071215
     1  0.699999789529, 0.801100153082, 0.962899979950, 1.312500009149, Na2  
     2  1.460400016380, 1.605300023181, 1.791600005572, 1.973499868954, Na2  
     3  2.204899919770, 2.425900205283, 2.630200108135, 2.831699980682, Na2  
     4  3.008900198602, 3.211700062255, 3.285199863526, 3.357699807757, Na2  
     5  3.458800006393, 3.567499866216, 3.669400046982, 3.821999733321, Na2  
     6  3.914299954705, 4.000000000000,      5*0.0D+00/                 Na2  
      DATA  Q_Na2/                                                      071215
     1  1.05909841D+00, 1.15888275D+00, 1.31914114D+00, 1.66697762D+00, Na2  
     2  1.81471052D+00, 1.96086727D+00, 2.15697364D+00, 2.36900510D+00, Na2  
     3  2.68233774D+00, 3.02844353D+00, 3.38265229D+00, 3.75834771D+00, Na2  
     4  4.10984692D+00, 4.54805207D+00, 4.72041453D+00, 4.89759969D+00, Na2  
     5  5.15275217D+00, 5.42941704D+00, 5.68576570D+00, 6.07000268D+00, Na2  
     6  6.30887997D+00, 6.53496303D+00,      5*0.0D+00/                 Na2  
      DATA TQ_Mg2/                                                      071215
     1  0.699999789529, 0.915799991597, 1.022600038735, 1.127800068818, Mg2  
     2  1.396400008494, 1.555600045179, 1.740199784555, 1.933199893611, Mg2  
     3  2.096900312686, 2.341500197952, 2.543700040326, 2.737199713844, Mg2  
     4  2.923100142787, 3.244799883827, 3.365099769324, 3.496099899052, Mg2  
     5  3.701999853245, 3.847500313265, 3.941799643234, 3.977299746311, Mg2  
     6  4.000000000000,      6*0.0D+00/                                 Mg2  
      DATA  Q_Mg2/                                                      071215
     1  1.28625096D+00, 1.50096197D+00, 1.60801565D+00, 1.71512018D+00, Mg2  
     2  2.01138533D+00, 2.21592929D+00, 2.48967561D+00, 2.82283451D+00, Mg2  
     3  3.14439630D+00, 3.66433015D+00, 4.08178841D+00, 4.43526289D+00, Mg2  
     4  4.72070594D+00, 5.11041913D+00, 5.23626577D+00, 5.36872268D+00, Mg2  
     5  5.57481659D+00, 5.72372122D+00, 5.82612993D+00, 5.86692464D+00, Mg2  
     6  5.89383649D+00,      6*0.0D+00/                                 Mg2  
      DATA TQ_Al2/                                                      071215
     1  0.699999789529, 0.871499909176, 1.189500033715, 1.647200031946, Al2  
     2  1.796000096025, 1.935699947183, 2.145299950581, 2.285399853441, Al2  
     3  2.450199807466, 2.600500337538, 2.840100169699, 3.032799811414, Al2  
     4  3.128600085589, 3.216200167996, 3.411099858738, 3.505200111303, Al2  
     5  3.593400234862, 3.703199880566, 3.805700312655, 3.925600225576, Al2  
     6  3.970700220546, 4.000000000000,      5*0.0D+00/                 Al2  
      DATA  Q_Al2/                                                      071215
     1  1.41532052D+00, 1.58404536D+00, 1.89909806D+00, 2.35513418D+00, Al2  
     2  2.50388216D+00, 2.64471055D+00, 2.86559266D+00, 3.02778273D+00, Al2  
     3  3.24084818D+00, 3.45813467D+00, 3.84621638D+00, 4.18936153D+00, Al2  
     4  4.36895572D+00, 4.53890859D+00, 4.94791994D+00, 5.16819328D+00, Al2  
     5  5.38858885D+00, 5.67423455D+00, 5.94155444D+00, 6.24330749D+00, Al2  
     6  6.35207677D+00, 6.42114202D+00,      5*0.0D+00/                 Al2  
      DATA TQ_Si2/                                                      071215
     1  0.699999789529, 0.762900169857, 0.860700118031, 1.123099955810, Si2  
     2  1.423999974422, 1.764900106541, 1.923099997971, 2.071699725854, Si2  
     3  2.378500090781, 2.502800056197, 2.637300276362, 2.818999772428, Si2  
     4  3.027599914967, 3.267400398120, 3.490099759020, 3.601800250432, Si2  
     5  3.702199857799, 3.872399947254, 3.949299812574, 3.980599565776, Si2  
     6  4.000000000000,      6*0.0D+00/                                 Si2  
      DATA  Q_Si2/                                                      071215
     1  1.41029498D+00, 1.46384328D+00, 1.54954562D+00, 1.79023495D+00, Si2  
     2  2.07808059D+00, 2.41193359D+00, 2.56843849D+00, 2.71674248D+00, Si2  
     3  3.04261834D+00, 3.19220764D+00, 3.37033521D+00, 3.63854330D+00, Si2  
     4  3.98127780D+00, 4.41174536D+00, 4.83940607D+00, 5.06395680D+00, Si2  
     5  5.27391534D+00, 5.66334697D+00, 5.86008412D+00, 5.94396163D+00, Si2  
     6  5.99688744D+00,      6*0.0D+00/                                 Si2  
      DATA TQ_P2/                                                       071215
     1  0.699999789529, 0.866299970394, 1.144199969376, 1.844000051860, P2   
     2  2.047100131936, 2.241699825934, 2.554600287830, 2.676000217757, P2   
     3  2.814000129482, 3.041400009367, 3.360499647354, 3.540599971707, P2   
     4  3.640800271676, 3.730999578812, 3.888600323313, 3.955999959393, P2   
     5  3.983299626393, 4.000000000000,      9*0.0D+00/                 P2   
      DATA  Q_P2/                                                       071215
     1  7.72347201D-01, 9.34625905D-01, 1.20846506D+00, 1.90470727D+00, P2   
     2  2.10754099D+00, 2.30272630D+00, 2.63487945D+00, 2.78021367D+00, P2   
     3  2.96189146D+00, 3.30064397D+00, 3.84547943D+00, 4.18097240D+00, P2   
     4  4.37826500D+00, 4.56898434D+00, 4.96914711D+00, 5.18225094D+00, P2   
     5  5.27682750D+00, 5.33693046D+00,      9*0.0D+00/                 P2   
      DATA TQ_S2/                                                       071215
     1  0.699999789529, 0.725799935732, 0.765800093366, 0.868199920303, S2   
     2  1.139699884138, 1.475999920713, 1.853899996556, 2.052100246684, S2   
     3  2.230199552170, 2.481199556889, 2.584700035562, 2.708900028743, S2   
     4  2.851400280608, 2.985499664822, 3.327899869236, 3.512200274585, S2   
     5  3.682400130099, 3.785399869602, 3.883600204632, 3.953399902621, S2   
     6  3.982399606187, 4.000000000000,      5*0.0D+00/                 S2   
      DATA  Q_S2/                                                       071215
     1  1.33461514D+00, 1.35546896D+00, 1.38835863D+00, 1.47528039D+00, S2   
     2  1.71973382D+00, 2.03926864D+00, 2.40888170D+00, 2.60497773D+00, S2   
     3  2.78280998D+00, 3.04642089D+00, 3.16578265D+00, 3.32130030D+00, S2   
     4  3.51798955D+00, 3.72059946D+00, 4.30237947D+00, 4.64529520D+00, S2   
     5  4.97730901D+00, 5.18638325D+00, 5.39463833D+00, 5.55099974D+00, S2   
     6  5.61867573D+00, 5.66063229D+00,      5*0.0D+00/                 S2   
      DATA TQ_Cl2/                                                      071215
     1  0.699999789529, 0.843100073575, 1.080899942284, 1.746799921018, Cl2  
     2  1.927499887696, 2.095000264905, 2.368099841909, 2.526100003753, Cl2  
     3  2.691899626942, 2.852100228671, 3.047600132816, 3.269000436539, Cl2  
     4  3.356399905912, 3.448399760290, 3.519300440391, 3.599400367441, Cl2  
     5  3.667100000032, 3.733499630460, 3.826999854201, 3.894800002835, Cl2  
     6  3.959500035818, 3.984499653334, 4.000000000000,      4*0.0D+00/ Cl2  
      DATA  Q_Cl2/                                                      071215
     1  8.64822074D-01, 1.00507502D+00, 1.23980868D+00, 1.90248836D+00, Cl2  
     2  2.08296879D+00, 2.25103828D+00, 2.53802947D+00, 2.72435106D+00, Cl2  
     3  2.94435181D+00, 3.18198109D+00, 3.50195691D+00, 3.89688034D+00, Cl2  
     4  4.06114909D+00, 4.24123790D+00, 4.39033766D+00, 4.58175080D+00, Cl2  
     5  4.77402219D+00, 4.99417114D+00, 5.34349201D+00, 5.60679999D+00, Cl2  
     6  5.85298232D+00, 5.94524733D+00, 6.00146652D+00,      4*0.0D+00/ Cl2  
      DATA TQ_K2/                                                       071215
     1  0.699999789529, 1.212699956342, 1.391699896194, 1.591000140557, K2   
     2  1.766000076771, 2.007300167625, 2.252900069694, 2.501200013622, K2   
     3  2.733799640323, 2.954499914209, 3.131300150279, 3.219200238490, K2   
     4  3.308800382734, 3.454599904831, 3.521100376774, 3.593000226024, K2   
     5  3.695499714589, 3.837500090088, 3.928800296672, 4.000000000000, K2   
     6       7*0.0D+00/                                                 K2   
      DATA  Q_K2/                                                       071215
     1  1.49013601D+00, 2.00152464D+00, 2.18237477D+00, 2.39482765D+00, K2   
     2  2.60320055D+00, 2.93699423D+00, 3.33024362D+00, 3.77067695D+00, K2   
     3  4.21300325D+00, 4.65832155D+00, 5.03995045D+00, 5.24190112D+00, K2   
     4  5.45740526D+00, 5.82013310D+00, 5.98411414D+00, 6.15681971D+00, K2   
     5  6.39397214D+00, 6.71363439D+00, 6.92139260D+00, 7.08631501D+00, K2   
     6       7*0.0D+00/                                                 K2   
      DATA TQ_Cu2/                                                      071215
     1  0.699999789529, 0.823100086123, 1.022800043674, 1.455299924070, Cu2  
     2  1.628500078000, 1.783900123297, 2.035499860196, 2.180000306579, Cu2  
     3  2.346800317164, 2.506300149329, 2.789599967919, 3.083900020695, Cu2  
     4  3.334500027018, 3.454299897577, 3.569899693730, 3.685999879367, Cu2  
     5  3.802400240317, 3.919300083174, 4.000000000000,      8*0.0D+00/ Cu2  
      DATA  Q_Cu2/                                                      071215
     1  1.21033534D+00, 1.33232737D+00, 1.53078891D+00, 1.96200346D+00, Cu2  
     2  2.13505440D+00, 2.29117834D+00, 2.55558976D+00, 2.72418714D+00, Cu2  
     3  2.94203871D+00, 3.17529147D+00, 3.64262538D+00, 4.18287763D+00, Cu2  
     4  4.67784856D+00, 4.92740218D+00, 5.17984631D+00, 5.44962735D+00, Cu2  
     5  5.74092726D+00, 6.05432474D+00, 6.27887477D+00,      8*0.0D+00/ Cu2  
      DATA TQ_As2/                                                      071215
     1  0.699999789529, 0.829599928813, 1.046500098516, 1.662099881427, As2  
     2  1.805000060625, 1.961300021608, 2.102700198199, 2.332099976947, As2  
     3  2.481999577392, 2.642800137474, 2.934100004463, 3.243299849688, As2  
     4  3.348500355402, 3.450199798433, 3.553000259149, 3.655199743098, As2  
     5  3.752600048802, 3.848600338011, 3.940399611624, 3.976999767867, As2  
     6  4.000000000000,      6*0.0D+00/                                 As2  
      DATA  Q_As2/                                                      071215
     1  1.23821305D+00, 1.36672170D+00, 1.58238852D+00, 2.19654730D+00, As2  
     2  2.33934614D+00, 2.49604278D+00, 2.64026691D+00, 2.89185994D+00, As2  
     3  3.07776060D+00, 3.30071079D+00, 3.76253962D+00, 4.31414233D+00, As2  
     4  4.51220191D+00, 4.70835390D+00, 4.91322810D+00, 5.12891442D+00, As2  
     5  5.35453179D+00, 5.60553389D+00, 5.87569455D+00, 5.99059661D+00, As2  
     6  6.06441691D+00,      6*0.0D+00/                                 As2  
      DATA TQ_Se2/                                                      071215
     1  0.699999789529, 0.751500026057, 0.831599955811, 1.045300071891, Se2  
     2  1.612500077608, 1.789399987495, 1.969399810715, 2.111299704322, Se2  
     3  2.244299881142, 2.347600335158, 2.452099849878, 2.608399760790, Se2  
     4  2.791500012023, 2.933900018499, 3.105599975342, 3.296200107291, Se2  
     5  3.488599725638, 3.623299973012, 3.745399884860, 3.897299818612, Se2  
     6  3.959700040185, 3.984499653334, 4.000000000000,      4*0.0D+00/ Se2  
      DATA  Q_Se2/                                                      071215
     1  1.79106970D+00, 1.83965001D+00, 1.91584440D+00, 2.12205479D+00, Se2  
     2  2.68071073D+00, 2.85663733D+00, 3.03725002D+00, 3.18509339D+00, Se2  
     3  3.33605932D+00, 3.46742888D+00, 3.61629631D+00, 3.86964154D+00, Se2  
     4  4.20382649D+00, 4.48114823D+00, 4.82607207D+00, 5.21592417D+00, Se2  
     5  5.61467924D+00, 5.89815762D+00, 6.16133528D+00, 6.50710047D+00, Se2  
     6  6.65838306D+00, 6.72016996D+00, 6.75924860D+00,      4*0.0D+00/ Se2  
      DATA TQ_Sb2/                                                      071215
     1  0.699999789529, 1.087700109529, 1.541700011340, 1.684799949927, Sb2  
     2  1.820500178517, 2.158999835072, 2.317599882949, 2.472400077612, Sb2  
     3  2.819099765287, 2.986999698536, 3.165099955392, 3.327699864899, Sb2  
     4  3.474299943517, 3.598100338716, 3.705699937486, 3.799800183514, Sb2  
     5  3.880000119182, 3.953399902621, 3.982399606187, 4.000000000000, Sb2  
     6       7*0.0D+00/                                                 Sb2  
      DATA  Q_Sb2/                                                      071215
     1  1.54025288D+00, 1.92673076D+00, 2.38023116D+00, 2.52341747D+00, Sb2  
     2  2.66020904D+00, 3.02835948D+00, 3.22976701D+00, 3.44925878D+00, Sb2  
     3  4.01338713D+00, 4.31386516D+00, 4.64618900D+00, 4.95927771D+00, Sb2  
     4  5.24886673D+00, 5.50081617D+00, 5.72999725D+00, 5.94367363D+00, Sb2  
     5  6.13907080D+00, 6.32962212D+00, 6.40780295D+00, 6.45596508D+00, Sb2  
     6       7*0.0D+00/                                                 Sb2  
      DATA TQ_Te2/                                                      071215
     1  0.699999789529, 0.850699923150, 1.117399948762, 1.553100106080, Te2  
     2  1.690700053616, 1.820600175722, 2.127600040674, 2.318899789611, Te2  
     3  2.516900391962, 2.692499640332, 2.860399652536, 3.005800129426, Te2  
     4  3.131800161738, 3.381600144993, 3.542400014761, 3.712700093628, Te2  
     5  3.860999687322, 4.000000000000,      9*0.0D+00/                 Te2  
      DATA  Q_Te2/                                                      071215
     1  2.13163567D+00, 2.27895906D+00, 2.54193013D+00, 2.97489788D+00, Te2  
     2  3.11237970D+00, 3.24375600D+00, 3.58067132D+00, 3.82779386D+00, Te2  
     3  4.11987956D+00, 4.40747011D+00, 4.70753720D+00, 4.98991126D+00, Te2  
     4  5.25168841D+00, 5.80427401D+00, 6.17245775D+00, 6.57194942D+00, Te2  
     5  6.93617331D+00, 7.29248538D+00,      9*0.0D+00/                 Te2  
      DATA TQ_I2/                                                       071215
     1  0.699999789529, 1.074300070438, 1.419400066533, 1.584600032430, I2   
     2  1.731899975202, 1.939900037183, 2.103800119795, 2.238399751034, I2   
     3  2.369999887941, 2.659699853547, 2.886800271471, 3.102600192186, I2   
     4  3.201699836306, 3.298600158939, 3.471700129898, 3.552300244540, I2   
     5  3.628500088674, 3.723300022422, 3.811900269745, 3.875000006072, I2   
     6  3.927700272233, 3.971600155878, 3.988199736403, 4.000000000000, I2   
     7       3*0.0D+00/                                                 I2   
      DATA  Q_I2/                                                       071215
     1  1.67065503D+00, 2.04407127D+00, 2.38884131D+00, 2.55411377D+00, I2   
     2  2.70271139D+00, 2.92228442D+00, 3.11397716D+00, 3.28920025D+00, I2   
     3  3.47722280D+00, 3.94337981D+00, 4.34779221D+00, 4.75489499D+00, I2   
     4  4.94870773D+00, 5.14543344D+00, 5.54903550D+00, 5.78197911D+00, I2   
     5  6.03263621D+00, 6.37130332D+00, 6.69345453D+00, 6.91671638D+00, I2   
     6  7.09634971D+00, 7.24050509D+00, 7.29367113D+00, 7.33101035D+00, I2   
     7       3*0.0D+00/                                                 I2   
      DATA TQ_Cs2/                                                      071215
     1  0.699999789529, 0.850299913974, 1.011300182376, 1.225099909987, Cs2  
     2  1.429199836745, 1.586500080219, 1.735399891037, 2.053500278925, Cs2  
     3  2.318199839870, 2.597400318321, 2.703799902083, 2.813200186611, Cs2  
     4  3.007800174055, 3.104400062080, 3.202199847511, 3.368499859476, Cs2  
     5  3.447499741556, 3.530799747990, 3.673600144919, 3.831099950688, Cs2  
     6  3.938199732326, 3.975599868462, 4.000000000000,      4*0.0D+00/ Cs2  
      DATA  Q_Cs2/                                                      071215
     1  2.43893259D+00, 2.58920023D+00, 2.75128498D+00, 2.97609420D+00, Cs2  
     2  3.21704572D+00, 3.42838027D+00, 3.64979533D+00, 4.18126631D+00, Cs2  
     3  4.66740364D+00, 5.20970898D+00, 5.42351835D+00, 5.64850643D+00, Cs2  
     4  6.07212145D+00, 6.30182717D+00, 6.55013175D+00, 6.98839201D+00, Cs2  
     5  7.18951408D+00, 7.38924142D+00, 7.69798520D+00, 7.99745914D+00, Cs2  
     6  8.18505278D+00, 8.24845813D+00, 8.28932612D+00,      4*0.0D+00/ Cs2  
      DATA TQ_H2p/                                                      071215
     1  0.699999789529, 0.830599932887, 0.983499924737, 1.056100039702, H2p  
     2  1.127300056796, 1.218599830052, 1.305400055885, 1.380100096517, H2p  
     3  1.452799869735, 1.566500091055, 1.685699970990, 1.836900075409, H2p  
     4  1.968699828940, 2.096700307657, 2.216800187124, 2.541499990993, H2p  
     5  2.686599830628, 2.821799744468, 3.102800177730, 3.260500232439, H2p  
     6  3.463200104644, 3.569299736852, 3.678300256790, 3.897399811243, H2p  
     7  3.958100005248, 4.000000000000,      1*0.0D+00/                 H2p  
      DATA  Q_H2p/                                                      071215
     1 -3.01029825D-01,-3.01016082D-01,-3.00454420D-01,-2.98795782D-01, H2p  
     2 -2.94159387D-01,-2.78597577D-01,-2.46267675D-01,-2.00228977D-01, H2p  
     3 -1.39285243D-01,-2.15550072D-02, 1.11998437D-01, 2.72841711D-01, H2p  
     4  4.00066828D-01, 5.17031676D-01, 6.26231165D-01, 9.31227835D-01, H2p  
     5  1.07192364D+00, 1.20704010D+00, 1.52049852D+00, 1.72942914D+00, H2p  
     6  2.04187371D+00, 2.22601511D+00, 2.42926814D+00, 2.86712144D+00, H2p  
     7  2.99025654D+00, 3.07455288D+00,      1*0.0D+00/                 H2p  
      DATA TQ_He2p/                                                     071215
     1  0.699999789529, 0.807000018063, 0.929899846916, 1.024000073308, He2p 
     2  1.111000114888, 1.239399767535, 1.350699931060, 1.474899894956, He2p 
     3  1.629900108953, 1.824300072313, 2.049000175430, 2.268600422852, He2p 
     4  2.461700064562, 2.625700004832, 2.794800086754, 3.022300289709, He2p 
     5  3.291800012602, 3.386600271607, 3.488099714549, 3.581499960589, He2p 
     6  3.678800268691, 3.770400436996, 3.882000166654, 3.952799889519, He2p 
     7  3.982199601697, 4.000000000000,      1*0.0D+00/                 He2p 
      DATA  Q_He2p/                                                     071215
     1  3.01040679D-01, 3.01184077D-01, 3.02650335D-01, 3.07578958D-01, He2p 
     2  3.19529471D-01, 3.60450523D-01, 4.23741434D-01, 5.19002468D-01, He2p 
     3  6.56083610D-01, 8.37778230D-01, 1.05372827D+00, 1.26851615D+00, He2p 
     4  1.45952947D+00, 1.62437248D+00, 1.80251676D+00, 2.07336580D+00, He2p 
     5  2.46479835D+00, 2.62311199D+00, 2.80585558D+00, 2.98799640D+00, He2p 
     6  3.19357522D+00, 3.39988914D+00, 3.65795972D+00, 3.81865150D+00, He2p 
     7  3.88363152D+00, 3.92234990D+00,      1*0.0D+00/                 He2p 
      DATA TQ_C2p/                                                      071215
     1  0.699999789529, 0.730100033532, 0.776600143547, 0.895900000496, C2p  
     2  1.202899862312, 1.416700004458, 1.630400101202, 1.872199938448, C2p  
     3  2.125099988240, 2.308200386694, 2.479299576929, 2.814300108059, C2p  
     4  2.960200045521, 3.113999742871, 3.439399627900, 3.670600073512, C2p  
     5  3.831199952866, 3.936899826022, 3.975199897204, 4.000000000000, C2p  
     6       7*0.0D+00/                                                 C2p  
      DATA  Q_C2p/                                                      071215
     1  8.81468045D-01, 8.95510236D-01, 9.18955269D-01, 9.88006014D-01, C2p  
     2  1.21118396D+00, 1.39232356D+00, 1.58584708D+00, 1.81404081D+00, C2p  
     3  2.05893057D+00, 2.23866341D+00, 2.40843258D+00, 2.76372759D+00, C2p  
     4  2.94195355D+00, 3.15157414D+00, 3.66479403D+00, 4.07427860D+00, C2p  
     5  4.37469081D+00, 4.57872452D+00, 4.65397831D+00, 4.70312425D+00, C2p  
     6       7*0.0D+00/                                                 C2p  
      DATA TQ_N2p/                                                      071215
     1  0.699999789529, 0.726299947707, 0.767500048526, 0.867499938758, N2p  
     2  1.016400060534, 1.167000000044, 1.416199992962, 1.674499964166, N2p  
     3  2.117199824324, 2.374800002486, 2.597000309749, 2.773100242367, N2p  
     4  2.956799967199, 3.101500271695, 3.310100402031, 3.471900115561, N2p  
     5  3.612299712796, 3.716900184223, 3.851000299003, 3.943199674844, N2p  
     6  3.977799710384, 4.000000000000,      5*0.0D+00/                 N2p  
      DATA  Q_N2p/                                                      071215
     1  3.51449118D-01, 3.70329774D-01, 4.01326069D-01, 4.82533132D-01, N2p  
     2  6.13653921D-01, 7.52631047D-01, 9.89813342D-01, 1.24121109D+00, N2p  
     3  1.67856477D+00, 1.93490244D+00, 2.15674373D+00, 2.33479309D+00, N2p  
     4  2.53041287D+00, 2.69994225D+00, 2.98005508D+00, 3.23650170D+00, N2p  
     5  3.49832203D+00, 3.72081934D+00, 4.03760515D+00, 4.27271245D+00, N2p  
     6  4.36426628D+00, 4.42397744D+00,      5*0.0D+00/                 N2p  
      DATA TQ_O2p/                                                      071215
     1  0.699999789529, 0.733599936054, 0.785800067582, 0.918900060997, O2p  
     2  1.084500030825, 1.264700106388, 1.447899857958, 1.610900117616, O2p  
     3  1.912699897189, 2.046700122779, 2.184100012585, 2.439099648246, O2p  
     4  2.562100263074, 2.680700257665, 2.885100229904, 3.060600395421, O2p  
     5  3.240899795066, 3.404999728673, 3.663399924504, 3.754900103311, O2p  
     6  3.841900187284, 3.939699624215, 3.976299818164, 3.989399763344, O2p  
     7  4.000000000000,      2*0.0D+00/                                 O2p  
      DATA  Q_O2p/                                                      071215
     1  5.77966736D-01, 5.93513456D-01, 6.19872750D-01, 6.98075277D-01, O2p  
     2  8.13582025D-01, 9.56272505D-01, 1.11400789D+00, 1.26208211D+00, O2p  
     3  1.55991874D+00, 1.70888729D+00, 1.87327744D+00, 2.19366718D+00, O2p  
     4  2.34759894D+00, 2.49375616D+00, 2.74596448D+00, 2.97508543D+00, O2p  
     5  3.23364022D+00, 3.49240077D+00, 3.94234319D+00, 4.11364860D+00, O2p  
     6  4.28447510D+00, 4.49344196D+00, 4.57937971D+00, 4.61150809D+00, O2p  
     7  4.63809585D+00,      2*0.0D+00/                                 O2p  
      DATA TQ_Ne2p/                                                     071215
     1  0.699999789529, 0.723199873466, 0.759500233390, 0.852299959853, Ne2p 
     2  1.094700045881, 1.422900003546, 1.830999936543, 1.976999942815, Ne2p 
     3  2.121099904345, 2.343400240688, 2.555800315543, 2.831099967175, Ne2p 
     4  3.084000022637, 3.361599676521, 3.591800199508, 3.786099886106, Ne2p 
     5  4.000000000000,     10*0.0D+00/                                 Ne2p 
      DATA  Q_Ne2p/                                                     071215
     1  9.59787458D-01, 9.74591931D-01, 9.98623674D-01, 1.06435317D+00, Ne2p 
     2  1.25816744D+00, 1.55190866D+00, 1.94142875D+00, 2.08420633D+00, Ne2p 
     3  2.22739968D+00, 2.46142836D+00, 2.71700201D+00, 3.11040627D+00, Ne2p 
     4  3.52632288D+00, 4.02394592D+00, 4.45689475D+00, 4.83108264D+00, Ne2p 
     5  5.24859384D+00,     10*0.0D+00/                                 Ne2p 
      DATA TQ_P2p/                                                      071215
     1  0.699999789529, 0.749799982595, 0.827199986897, 1.029200201723, P2p  
     2  1.566100081636, 1.798700151530, 1.964799930481, 2.110899696187, P2p  
     3  2.539999957356, 2.736199692220, 2.935799885162, 3.111199682975, P2p  
     4  3.273800188258, 3.475699843158, 3.570499698488, 3.663199920421, P2p  
     5  3.764300328137, 3.866699817560, 3.946499749354, 3.979499588232, P2p  
     6  4.000000000000,      6*0.0D+00/                                 P2p  
      DATA  Q_P2p/                                                      071215
     1  1.44865192D+00, 1.49358521D+00, 1.56443035D+00, 1.75388474D+00, P2p  
     2  2.27571621D+00, 2.50689802D+00, 2.67821873D+00, 2.83977360D+00, P2p  
     3  3.40006312D+00, 3.70090338D+00, 4.03597281D+00, 4.35344588D+00, P2p  
     4  4.66432343D+00, 5.06736324D+00, 5.26216298D+00, 5.45705298D+00, P2p  
     5  5.67807394D+00, 5.91832831D+00, 6.12299745D+00, 6.21283460D+00, P2p  
     6  6.27015425D+00,      6*0.0D+00/                                 P2p  
      DATA TQ_S2p/                                                      071215
     1  0.699999789529, 0.766000088091, 0.867899928212, 1.147800048735, S2p  
     2  1.412399905596, 1.702699877846, 1.872899955765, 2.044200065549, S2p  
     3  2.193599661638, 2.338700133397, 2.447599750785, 2.557900364040, S2p  
     4  2.752600051197, 2.952899877347, 3.256100135312, 3.409799825748, S2p  
     5  3.547100127177, 3.667900016363, 3.780599756436, 3.912799916165, S2p  
     6  3.966000181200, 3.986599700481, 4.000000000000,      4*0.0D+00/ S2p  
      DATA  Q_S2p/                                                      071215
     1  1.14654707D+00, 1.20619342D+00, 1.29997094D+00, 1.56520780D+00, S2p  
     2  1.82238109D+00, 2.10835652D+00, 2.27713040D+00, 2.44835724D+00, S2p  
     3  2.60215346D+00, 2.76236416D+00, 2.89437279D+00, 3.04090802D+00, S2p  
     4  3.33145149D+00, 3.66546893D+00, 4.21470189D+00, 4.50536926D+00, S2p  
     5  4.76984550D+00, 5.00742384D+00, 5.23730214D+00, 5.52607990D+00, S2p  
     6  5.65027414D+00, 5.69970902D+00, 5.73227311D+00,      4*0.0D+00/ S2p  
      DATA TQ_H2m/                                                      071215
     1  0.699999789529, 1.032900143995, 1.328999957686, 1.473199855148, H2m  
     2  1.687100003756, 1.829799918596, 1.958100007507, 2.105300012880, H2m  
     3  2.246199921486, 2.377600069304, 2.506900165294, 2.668400046539, H2m  
     4  2.841900208750, 2.991999812772, 3.140500312349, 3.504000083327, H2m  
     5  3.609199714890, 3.717500197165, 3.882900188017, 3.951299856766, H2m  
     6  4.000000000000,      6*0.0D+00/                                 H2m  
      DATA  Q_H2m/                                                      071215
     1  1.76091259D-01, 1.76091319D-01, 1.76238392D-01, 1.77495707D-01, H2m  
     2  1.89032526D-01, 2.10700199D-01, 2.44777025D-01, 3.08035439D-01, H2m  
     3  3.97734309D-01, 5.02645054D-01, 6.17168038D-01, 7.67508056D-01, H2m  
     4  9.34210952D-01, 1.08254530D+00, 1.23610779D+00, 1.68778725D+00, H2m  
     5  1.85192638D+00, 2.04024517D+00, 2.35732361D+00, 2.49341851D+00, H2m  
     6  2.59050567D+00,      6*0.0D+00/                                 H2m  
      DATA TQ_C2m/                                                      071215
     1  0.699999789529, 0.720299804015, 0.753200070115, 0.832299971858, C2m  
     2  0.945899894469, 1.057400010618, 1.287800007289, 1.551800137749, C2m  
     3  2.070899705798, 2.324699812404, 2.554000273974, 2.695199700588, C2m  
     4  2.861599682590, 3.005300118268, 3.192599624783, 3.426700227116, C2m  
     5  3.560200390862, 3.690999621469, 3.782199794158, 3.866399810706, C2m  
     6  3.948999805801, 3.980599565776, 4.000000000000,      4*0.0D+00/ C2m  
      DATA  Q_C2m/                                                      071215
     1  3.98202493D-01, 4.10837922D-01, 4.32928875D-01, 4.92847298D-01, C2m  
     2  5.89909349D-01, 6.91287443D-01, 9.08407869D-01, 1.16386903D+00, C2m  
     3  1.67596886D+00, 1.92853431D+00, 2.15770129D+00, 2.30110623D+00, C2m  
     4  2.47867951D+00, 2.64661815D+00, 2.89385837D+00, 3.25077745D+00, C2m  
     5  3.47702069D+00, 3.71943037D+00, 3.90868231D+00, 4.10442020D+00, C2m  
     6  4.31697336D+00, 4.40281865D+00, 4.45648840D+00,      4*0.0D+00/ C2m  
      DATA TQ_LiH/                                                      071215
     1  0.699999789529, 0.730700016821, 0.778800197201, 0.895199984829, LiH  
     2  1.042400007546, 1.175399942663, 1.393299934424, 1.616899967588, LiH  
     3  1.876000032456, 2.146299880719, 2.341600200201, 2.515800368977, LiH  
     4  2.873399965124, 3.039299964956, 3.205799928186, 3.316499946864, LiH  
     5  3.422700147002, 3.551800234105, 3.636800262222, 3.719000229521, LiH  
     6  3.794200066108, 3.865599792427, 3.946599751612, 3.979499588232, LiH  
     7  4.000000000000,      2*0.0D+00/                                 LiH  
      DATA  Q_LiH/                                                      071215
     1  1.81871275D-02, 2.41476148D-02, 3.61007022D-02, 7.94297443D-02, LiH  
     2  1.61281721D-01, 2.54027634D-01, 4.29322890D-01, 6.27316160D-01, LiH  
     3  8.69754064D-01, 1.13110179D+00, 1.32323303D+00, 1.49726545D+00, LiH  
     4  1.88956420D+00, 2.10794320D+00, 2.35854290D+00, 2.54373390D+00, LiH  
     5  2.73678457D+00, 2.99410073D+00, 3.17636630D+00, 3.36010485D+00, LiH  
     6  3.53425890D+00, 3.70779339D+00, 3.92110094D+00, 4.01431925D+00, LiH  
     7  4.07446717D+00,      2*0.0D+00/                                 LiH  
      DATA TQ_BeH/                                                      071215
     1  0.699999789529, 0.734099922129, 0.789399968965, 0.930499854246, BeH  
     2  1.054800068786, 1.204599898664, 1.372699939755, 1.577499980917, BeH  
     3  1.783300138111, 1.992199822742, 2.280599756152, 2.519800452559, BeH  
     4  2.715000168729, 2.935399913233, 3.147099819858, 3.383100182977, BeH  
     5  3.490399766021, 3.598900356393, 3.710400044016, 3.820899706727, BeH  
     6  3.931100244051, 3.973100048097, 4.000000000000,      4*0.0D+00/ BeH  
      DATA  Q_BeH/                                                      071215
     1  3.79224761D-03, 5.87787992D-03, 1.11358995D-02, 4.02620528D-02, BeH  
     2  8.99158531D-02, 1.77352335D-01, 2.99718363D-01, 4.69414428D-01, BeH  
     3  6.53705390D-01, 8.49327176D-01, 1.12774319D+00, 1.36315701D+00, BeH  
     4  1.55893036D+00, 1.79469729D+00, 2.05603169D+00, 2.40514144D+00, BeH  
     5  2.58629797D+00, 2.78650944D+00, 3.01605732D+00, 3.27421672D+00, BeH  
     6  3.56309965D+00, 3.68006704D+00, 3.75649766D+00,      4*0.0D+00/ BeH  
      DATA TQ_BH/                                                       071215
     1  0.699999789529, 0.733199947194, 0.787100031971, 0.962799982554, BH   
     2  1.092700097037, 1.240299759776, 1.396700015662, 1.589300150645, BH   
     3  1.782800150457, 1.993599854557, 2.261200261212, 2.535399845476, BH   
     4  2.722200122283, 2.906499752472, 3.144999976559, 3.389400342512, BH   
     5  3.551700232019, 3.679700290113, 3.794600074494, 3.924600203359, BH   
     6  3.971000198990, 4.000000000000,      5*0.0D+00/                 BH   
      DATA  Q_BH/                                                       071215
     1  1.47489903D-03, 2.42692400D-03, 5.03778997D-03, 3.10120434D-02, BH   
     2  7.69785165D-02, 1.58106565D-01, 2.67889857D-01, 4.23463459D-01, BH   
     3  5.93664260D-01, 7.88807452D-01, 1.04527184D+00, 1.31404643D+00, BH   
     4  1.50025346D+00, 1.69155084D+00, 1.97071900D+00, 2.31805693D+00, BH   
     5  2.58974716D+00, 2.83196861D+00, 3.07496019D+00, 3.37720343D+00, BH   
     6  3.48994630D+00, 3.56117032D+00,      5*0.0D+00/                 BH   
      DATA TQ_CH/                                                       071215
     1  0.699999789529, 0.869299891303, 1.044300049703, 1.194699929527, CH   
     2  1.336400081666, 1.443699956112, 1.552900110953, 1.722399844607, CH   
     3  1.938500007183, 2.191699623576, 2.468700222442, 2.732699616537, CH   
     4  2.957299978718, 3.234199633717, 3.339100141763, 3.444899687433, CH   
     5  3.603900098454, 3.779999742290, 3.904599722954, 3.963600127415, CH   
     6  4.000000000000,      6*0.0D+00/                                 CH   
      DATA  Q_CH/                                                       071215
     1  1.07918131D+00, 1.07919302D+00, 1.07963725D+00, 1.08310999D+00, CH   
     2  1.09577727D+00, 1.11653631D+00, 1.15015128D+00, 1.22745964D+00, CH   
     3  1.36317588D+00, 1.55929023D+00, 1.80192879D+00, 2.04926865D+00, CH   
     4  2.27242936D+00, 2.59252337D+00, 2.73787615D+00, 2.90241811D+00, CH   
     5  3.18645666D+00, 3.55419639D+00, 3.84860706D+00, 3.99618659D+00, CH   
     6  4.08895880D+00,      6*0.0D+00/                                 CH   
      DATA TQ_NH/                                                       071215
     1  0.699999789529, 0.878200071439, 1.073100102078, 1.212799954201, NH   
     2  1.357900090828, 1.476199925396, 1.601799936098, 1.908199876296, NH   
     3  2.330699943760, 2.575099815375, 2.781699767765, 2.966500190371, NH   
     4  3.220800199242, 3.407699783278, 3.553400267497, 3.670900080652, NH   
     5  3.817299879894, 3.925500223354, 3.970900206175, 4.000000000000, NH   
     6       7*0.0D+00/                                                 NH   
      DATA  Q_NH/                                                       071215
     1  4.77231076D-01, 4.79695614D-01, 5.00976735D-01, 5.45018350D-01, NH   
     2  6.20703396D-01, 6.99680909D-01, 7.94879750D-01, 1.05678266D+00, NH   
     3  1.45326339D+00, 1.69141790D+00, 1.89574731D+00, 2.08314656D+00, NH   
     4  2.36707012D+00, 2.61372633D+00, 2.83618708D+00, 3.03921973D+00, NH   
     5  3.33084039D+00, 3.58008507D+00, 3.69308991D+00, 3.76773849D+00, NH   
     6       7*0.0D+00/                                                 NH   
      DATA TQ_OH/                                                       071215
     1  0.699999789529, 0.911899904287, 1.117399948762, 1.258300183474, OH   
     2  1.433099885940, 1.578999942405, 1.698499852916, 1.815900090198, OH   
     3  2.007700176810, 2.217900208498, 2.533599801697, 2.708800026260, OH   
     4  2.868899865413, 3.044400069101, 3.226199807431, 3.355499973866, OH   
     5  3.486799685718, 3.646799847820, 3.810500370817, 3.921900143372, OH   
     6  3.969700264120, 4.000000000000,      5*0.0D+00/                 OH   
      DATA  Q_OH/                                                       071215
     1  1.07918125D+00, 1.07918278D+00, 1.07939238D+00, 1.08119111D+00, OH   
     2  1.09316162D+00, 1.12145541D+00, 1.16139978D+00, 1.21541391D+00, OH   
     3  1.33065446D+00, 1.48598185D+00, 1.75538444D+00, 1.91641506D+00, OH   
     4  2.06851027D+00, 2.24239615D+00, 2.43925378D+00, 2.59670845D+00, OH   
     5  2.77559074D+00, 3.02223797D+00, 3.31092982D+00, 3.53268388D+00, OH   
     6  3.63467413D+00, 3.70130685D+00,      5*0.0D+00/                 OH   
      DATA TQ_HF/                                                       071215
     1  0.699999789529, 0.804700070698, 0.924499977120, 1.148300059757, HF   
     2  1.226599942311, 1.311399986641, 1.493999866319, 1.658699866494, HF   
     3  1.843500063924, 2.030899748561, 2.432600101900, 2.843900252141, HF   
     4  3.048200144763, 3.296000102987, 3.488999734509, 3.660699869389, HF   
     5  3.809500395954, 3.922300152259, 3.969600261879, 4.000000000000, HF   
     6       7*0.0D+00/                                                 HF   
      DATA  Q_HF/                                                       071215
     1  9.77535710D-06, 1.22437573D-04, 1.14337479D-03, 1.90562872D-02, HF   
     2  3.73617873D-02, 6.74616332D-02, 1.66561805D-01, 2.84307692D-01, HF   
     3  4.35107809D-01, 6.00521543D-01, 9.78327827D-01, 1.38136232D+00, HF   
     4  1.58732398D+00, 1.85979382D+00, 2.10881540D+00, 2.36641497D+00, HF   
     5  2.62012713D+00, 2.83406300D+00, 2.92985459D+00, 2.99327214D+00, HF   
     6       7*0.0D+00/                                                 HF   
      DATA TQ_NaH/                                                      071215
     1  0.699999789529, 0.725799935732, 0.765800093366, 0.866399967758, NaH  
     2  0.994399882231, 1.137999925775, 1.339300148162, 1.552100130441, NaH  
     3  1.789499985026, 2.042100017477, 2.237599731632, 2.409099828461, NaH  
     4  2.544500058265, 2.692899649259, 2.859599672197, 3.018300413371, NaH  
     5  3.194899678685, 3.328899890918, 3.519700449733, 3.613899751614, NaH  
     6  3.714700136768, 3.839800140184, 3.936999818815, 3.975399882833, NaH  
     7  4.000000000000,      2*0.0D+00/                                 NaH  
      DATA  Q_NaH/                                                      071215
     1  7.48799595D-02, 8.68509617D-02, 1.07246546D-01, 1.67176102D-01, NaH  
     2  2.57115597D-01, 3.70456723D-01, 5.43740034D-01, 7.38458595D-01, NaH  
     3  9.64057791D-01, 1.20982340D+00, 1.40257889D+00, 1.57364836D+00, NaH  
     4  1.71233635D+00, 1.87367534D+00, 2.07534738D+00, 2.29385435D+00, NaH  
     5  2.57023479D+00, 2.80398853D+00, 3.17810198D+00, 3.38337151D+00, NaH  
     6  3.61661171D+00, 3.92130706D+00, 4.16886848D+00, 4.26933463D+00, NaH  
     7  4.33442459D+00,      2*0.0D+00/                                 NaH  
      DATA TQ_MgH/                                                      071215
     1  0.699999789529, 0.727099966866, 0.769000008962, 0.872999945503, MgH  
     2  1.001000021343, 1.143199947332, 1.363000060461, 1.593500071323, MgH  
     3  1.846999979476, 2.113999759238, 2.318999782431, 2.502000034910, MgH  
     4  2.928700264041, 3.299400176155, 3.402999688226, 3.560500369301, MgH  
     5  3.662399904091, 3.792000019985, 3.914399957275, 4.000000000000, MgH  
     6       7*0.0D+00/                                                 MgH  
      DATA  Q_MgH/                                                      071215
     1  3.46993917D-01, 3.56563342D-01, 3.73567880D-01, 4.26602171D-01, MgH  
     2  5.08974267D-01, 6.15601736D-01, 7.99976531D-01, 1.00878850D+00, MgH  
     3  1.24885592D+00, 1.50838132D+00, 1.71043501D+00, 1.89282410D+00, MgH  
     4  2.36197962D+00, 2.90558688D+00, 3.08787196D+00, 3.39887738D+00, MgH  
     5  3.62678020D+00, 3.94565106D+00, 4.26808807D+00, 4.50022891D+00, MgH  
     6       7*0.0D+00/                                                 MgH  
      DATA TQ_AlH/                                                      071215
     1  0.699999789529, 0.729800031527, 0.776400138669, 0.889399883237, AlH  
     2  1.031800173386, 1.184799935855, 1.414099944681, 1.648400060249, AlH  
     3  1.900900059275, 2.166899995381, 2.362499706236, 2.534599826019, AlH  
     4  2.697899760844, 2.905299725487, 3.092800206632, 3.344700269929, AlH  
     5  3.451199822614, 3.547700141528, 3.697899764253, 3.787799926186, AlH  
     6  3.896699862825, 3.957699996514, 4.000000000000,      4*0.0D+00/ AlH  
      DATA  Q_AlH/                                                      071215
     1  3.36974449D-02, 4.24291080D-02, 5.88004720D-02, 1.11896303D-01, AlH  
     2  2.01029113D-01, 3.15272476D-01, 5.08074433D-01, 7.20949896D-01, AlH  
     3  9.60494293D-01, 1.21926704D+00, 1.41209412D+00, 1.58341693D+00, AlH  
     4  1.75005217D+00, 1.97936488D+00, 2.21824781D+00, 2.59949786D+00, AlH  
     5  2.78236857D+00, 2.95993656D+00, 3.26188398D+00, 3.45757904D+00, AlH  
     6  3.70347124D+00, 3.84131789D+00, 3.93550049D+00,      4*0.0D+00/ AlH  
      DATA TQ_SiH/                                                      071215
     1  0.699999789529, 0.879800110189, 1.078199967608, 1.233999910028, SiH  
     2  1.394599965485, 1.690600056189, 2.031799770402, 2.429800279553, SiH  
     3  2.585000042590, 2.723999992099, 2.926800222901, 3.090700156032, SiH  
     4  3.246599924794, 3.501700029706, 3.623599979685, 3.795600095460, SiH  
     5  3.912299903318, 3.966700196888, 4.000000000000,      8*0.0D+00/ SiH  
      DATA  Q_SiH/                                                      071215
     1  1.07933024D+00, 1.08182215D+00, 1.09956223D+00, 1.13762400D+00, SiH  
     2  1.20416772D+00, 1.38819508D+00, 1.66399650D+00, 2.02893899D+00, SiH  
     3  2.17809739D+00, 2.31522984D+00, 2.52898995D+00, 2.72426552D+00, SiH  
     4  2.93581599D+00, 3.34148317D+00, 3.56264277D+00, 3.91271234D+00, SiH  
     5  4.18102960D+00, 4.31355101D+00, 4.39624609D+00,      8*0.0D+00/ SiH  
      DATA TQ_PH/                                                       071215
     1  0.699999789529, 0.734399913773, 0.789099977183, 0.918600054281, PH   
     2  1.211699977747, 1.461099999160, 1.732899951155, 2.000000000000, PH   
     3  2.278999814935, 2.485799674780, 2.657699804533, 2.819899708158, PH   
     4  2.992699829154, 3.188899643623, 3.610399666698, 3.731799595339, PH   
     5  3.833600005141, 3.932000179184, 3.973600012170, 4.000000000000, PH   
     6       7*0.0D+00/                                                 PH   
      DATA  Q_PH/                                                       071215
     1  4.87411317D-01, 4.91909888D-01, 5.01917007D-01, 5.42584096D-01, PH   
     2  7.16982691D-01, 9.16960193D-01, 1.15996685D+00, 1.41216628D+00, PH   
     3  1.68314527D+00, 1.88688711D+00, 2.05779709D+00, 2.22240823D+00, PH   
     4  2.40884083D+00, 2.64867202D+00, 3.32685520D+00, 3.57635805D+00, PH   
     5  3.80820771D+00, 4.05167736D+00, 4.15922395D+00, 4.22845300D+00, PH   
     6       7*0.0D+00/                                                 PH   
      DATA TQ_HS/                                                       071215
     1  0.699999789529, 0.808299988314, 0.922300030166, 1.218599830052, HS   
     2  1.328399971347, 1.433199888210, 1.822200131004, 2.103800119795, HS   
     3  2.402199660122, 2.568499807571, 2.728299681105, 2.910599845847, HS   
     4  3.148399722852, 3.296200107291, 3.450699810524, 3.606799888579, HS   
     5  3.787599921470, 3.848700340261, 3.908599812783, 3.964600149825, HS   
     6  4.000000000000,      6*0.0D+00/                                 HS   
      DATA  Q_HS/                                                       071215
     1  7.78167493D-01, 7.78323696D-01, 7.79333212D-01, 8.05815320D-01, HS   
     2  8.34530520D-01, 8.74172999D-01, 1.10986912D+00, 1.33799074D+00, HS   
     3  1.60699334D+00, 1.76402382D+00, 1.91837111D+00, 2.10087116D+00, HS   
     4  2.36531977D+00, 2.55573216D+00, 2.78113386D+00, 3.03792565D+00, HS   
     5  3.37780741D+00, 3.50579649D+00, 3.63859667D+00, 3.76889168D+00, HS   
     6  3.85376876D+00,      6*0.0D+00/                                 HS   
      DATA TQ_HCl/                                                      071215
     1  0.699999789529, 0.736699849717, 0.796900108280, 0.946399882193, HCl  
     2  1.086600082475, 1.241699797451, 1.411399882606, 1.587700110402, HCl  
     3  1.774200073970, 1.965999899238, 2.263900320189, 2.547100116569, HCl  
     4  2.736299694382, 2.911299862961, 3.071899742099, 3.230599546283, HCl  
     5  3.376000022385, 3.531599766196, 3.677300232987, 3.831499959401, HCl  
     6  3.931500215221, 3.973300033726, 4.000000000000,      4*0.0D+00/ HCl  
      DATA  Q_HCl/                                                      071215
     1  3.23905174D-03, 5.25209836D-03, 1.06482061D-02, 4.15790456D-02, HCl  
     2  1.00105103D-01, 1.94536925D-01, 3.21057721D-01, 4.68117386D-01, HCl  
     3  6.34580255D-01, 8.13223801D-01, 1.09962143D+00, 1.37762937D+00, HCl  
     4  1.56543677D+00, 1.74246693D+00, 1.91412705D+00, 2.10156815D+00, HCl  
     5  2.29503623D+00, 2.52801686D+00, 2.77118108D+00, 3.05686630D+00, HCl  
     6  3.26044656D+00, 3.35022524D+00, 3.40899268D+00,      4*0.0D+00/ HCl  
      DATA TQ_KH/                                                       071215
     1  0.699999789529, 0.725599930943, 0.765300106554, 0.867199946667, KH   
     2  1.147700046531, 1.507900177925, 1.962899979950, 2.150899643511, KH   
     3  2.330099929537, 2.466800179589, 2.601300279133, 2.807100354048, KH   
     4  3.144999976559, 3.268900434138, 3.403899706427, 3.500500001730, KH   
     5  3.598900356393, 3.810800349159, 3.926100236685, 3.971100191804, KH   
     6  4.000000000000,      6*0.0D+00/                                 KH   
      DATA  Q_KH/                                                       071215
     1  4.61848231D-01, 4.78537510D-01, 5.05575953D-01, 5.80385899D-01, KH   
     2  8.13635067D-01, 1.14529198D+00, 1.58652136D+00, 1.77226124D+00, KH   
     3  1.95107353D+00, 2.09093703D+00, 2.23613422D+00, 2.48504837D+00, KH   
     4  2.98885498D+00, 3.20469817D+00, 3.46064541D+00, 3.65978952D+00, KH   
     5  3.87755012D+00, 4.37135572D+00, 4.62746352D+00, 4.72162633D+00, KH   
     6  4.78016966D+00,      6*0.0D+00/                                 KH   
      DATA TQ_CaH/                                                      071215
     1  0.699999789529, 0.725099918968, 0.763900143481, 0.862900060031, CaH  
     2  0.987399843714, 1.131200092324, 1.334800044978, 1.550500169418, CaH  
     3  1.796900114527, 2.054100292742, 2.262700293977, 2.454399901218, CaH  
     4  2.710800074089, 2.914899950976, 3.095100262051, 3.347700337407, CaH  
     5  3.496899917723, 3.646599861949, 3.744799871514, 3.840000144540, CaH  
     6  3.938799689081, 3.976499803794, 4.000000000000,      4*0.0D+00/ CaH  
      DATA  Q_CaH/                                                      071215
     1  4.04218507D-01, 4.17841276D-01, 4.40405080D-01, 5.05031783D-01, CaH  
     2  5.97249232D-01, 7.14308545D-01, 8.92845126D-01, 1.09233397D+00, CaH  
     3  1.32789842D+00, 1.57878387D+00, 1.78443186D+00, 1.97515303D+00, CaH  
     4  2.24269897D+00, 2.48533722D+00, 2.73310884D+00, 3.13613360D+00, CaH  
     5  3.40462993D+00, 3.70568387D+00, 3.93031277D+00, 4.17467812D+00, CaH  
     6  4.45432250D+00, 4.56608964D+00, 4.63662759D+00,      4*0.0D+00/ CaH  
      DATA TQ_TiH/                                                      071215
     1  0.699999789529, 0.894199962446, 1.091300132846, 1.317200105320, TiH  
     2  1.552600118261, 1.740099782487, 1.918600041436, 2.049400184587, TiH  
     3  2.246299923610, 2.393300098938, 2.530399723868, 2.734799661947, TiH  
     4  2.888300308148, 3.054700302595, 3.256600146225, 3.450099796015, TiH  
     5  3.569599715291, 3.688199726142, 3.826999854201, 3.930600280088, TiH  
     6  3.973000055282, 4.000000000000,      5*0.0D+00/                 TiH  
      DATA  Q_TiH/                                                      071215
     1  1.14613063D+00, 1.14634538D+00, 1.14993624D+00, 1.17528256D+00, TiH  
     2  1.25610823D+00, 1.37794808D+00, 1.54688381D+00, 1.70008016D+00, TiH  
     3  1.96176174D+00, 2.16824533D+00, 2.36428666D+00, 2.67147762D+00, TiH  
     4  2.92995283D+00, 3.24923567D+00, 3.69274144D+00, 4.17124766D+00, TiH  
     5  4.49304678D+00, 4.83017496D+00, 5.23401798D+00, 5.52899384D+00, TiH  
     6  5.64585639D+00, 5.71878724D+00,      5*0.0D+00/                 TiH  
      DATA TQ_CrH/                                                      071215
     1  0.699999789529, 0.727199969261, 0.769399998412, 0.873599960034, CrH  
     2  1.004400093911, 1.137899928224, 1.371299909224, 1.603799985860, CrH  
     3  1.866599969860, 2.143800055376, 2.363199723195, 2.547400123296, CrH  
     4  2.684499982624, 2.846000297702, 2.982099588403, 3.334300022030, CrH  
     5  3.466800182530, 3.607999801735, 3.680200283324, 3.753300065392, CrH  
     6  3.905399740920, 3.962300098281, 4.000000000000,      4*0.0D+00/ CrH  
      DATA  Q_CrH/                                                      071215
     1  8.15147006D-01, 8.23549195D-01, 8.38828212D-01, 8.87869009D-01, CrH  
     2  9.68462759D-01, 1.06594090D+00, 1.25906519D+00, 1.46863330D+00, CrH  
     3  1.71700714D+00, 1.98633167D+00, 2.20261881D+00, 2.38624206D+00, CrH  
     4  2.52635610D+00, 2.70180991D+00, 2.86502305D+00, 3.37892004D+00, CrH  
     5  3.61273488D+00, 3.90135470D+00, 4.07187121D+00, 4.26240978D+00, CrH  
     6  4.70246272D+00, 4.87351235D+00, 4.98639027D+00,      4*0.0D+00/ CrH  
      DATA TQ_MnH/                                                      071215
     1  0.699999789529, 0.727099966866, 0.769100006324, 0.873399955191, MnH  
     2  1.002700057627, 1.146700024486, 1.369199901398, 1.602299948539, MnH  
     3  1.855200026377, 2.122099925318, 2.326799857876, 2.513500320917, MnH  
     4  2.651099642787, 2.920700090821, 3.126800040721, 3.333399999580, MnH  
     5  3.458900008811, 3.572799753436, 3.651699661431, 3.731299585010, MnH  
     6  3.901199646599, 3.960600060182, 4.000000000000,      4*0.0D+00/ MnH  
      DATA  Q_MnH/                                                      071215
     1  8.94478193D-01, 9.04466149D-01, 9.22150496D-01, 9.76763653D-01, MnH  
     2  1.06133520D+00, 1.17044180D+00, 1.35833564D+00, 1.57038375D+00, MnH  
     3  1.81044828D+00, 2.07017910D+00, 2.27204462D+00, 2.45810315D+00, MnH  
     4  2.59857674D+00, 2.90141909D+00, 3.17751031D+00, 3.50208385D+00, MnH  
     5  3.72410087D+00, 3.94511951D+00, 4.11179614D+00, 4.29236786D+00, MnH  
     6  4.70979582D+00, 4.85975558D+00, 4.95863160D+00,      4*0.0D+00/ MnH  
      DATA TQ_FeH/                                                      071215
     1  0.699999789529, 0.906399944735, 1.218399834333, 1.329999934916, FeH  
     2  1.440500030897, 1.695599927535, 1.853099978204, 2.002300052814, FeH  
     3  2.378800097940, 2.507600183920, 2.642700144719, 2.950599824358, FeH  
     4  3.044000061136, 3.141500237729, 3.387400291866, 3.501900034369, FeH  
     5  3.618799870497, 3.763400306379, 3.884900235489, 3.955899957210, FeH  
     6  3.983199624148, 4.000000000000,      5*0.0D+00/                 FeH  
      DATA  Q_FeH/                                                      071215
     1  1.00000859D+00, 1.00058674D+00, 1.02040662D+00, 1.04384211D+00, FeH  
     2  1.07934353D+00, 1.20854380D+00, 1.32110396D+00, 1.45341492D+00, FeH  
     3  1.90036828D+00, 2.08601513D+00, 2.29780111D+00, 2.85846018D+00, FeH  
     4  3.05479569D+00, 3.27411555D+00, 3.88293754D+00, 4.18519990D+00, FeH  
     5  4.50173524D+00, 4.90213476D+00, 5.24150616D+00, 5.43752066D+00, FeH  
     6  5.51186318D+00, 5.55724606D+00,      5*0.0D+00/                 FeH  
      DATA TQ_CoH/                                                      071215
     1  0.699999789529, 0.730300027961, 0.777800172813, 0.892299919921, CoH  
     2  1.034500101246, 1.169900065715, 1.424899950594, 1.690700053616, CoH  
     3  1.936599966469, 2.203699887952, 2.362399703813, 2.518400423306, CoH  
     4  2.658899833941, 2.830699958170, 2.973000051951, 3.125400005824, CoH  
     5  3.519100435721, 3.816499937650, 3.923100170033, 4.000000000000, CoH  
     6       7*0.0D+00/                                                 CoH  
      DATA  Q_CoH/                                                      071215
     1  2.09745282D-02, 2.74608352D-02, 4.02633635D-02, 8.51696748D-02, CoH  
     2  1.65942068D-01, 2.61241280D-01, 4.69451261D-01, 7.09147735D-01, CoH  
     3  9.41946826D-01, 1.20143290D+00, 1.35764140D+00, 1.51269344D+00, CoH  
     4  1.65551677D+00, 1.84107104D+00, 2.01085622D+00, 2.21371777D+00, CoH  
     5  2.83829476D+00, 3.39086721D+00, 3.60492653D+00, 3.76360165D+00, CoH  
     6       7*0.0D+00/                                                 CoH  
      DATA TQ_NiH/                                                      071215
     1  0.699999789529, 0.846699986714, 1.006800145135, 1.249700012735, NiH  
     2  1.346600000080, 1.456899958845, 1.688000024819, 1.868999909269, NiH  
     3  2.093400224668, 2.321699747443, 2.539099935467, 2.717200218303, NiH  
     4  2.901499640034, 3.225799836454, 3.344700269929, 3.472700058213, NiH  
     5  3.622599957442, 3.739099746150, 3.836700072663, 3.918700067758, NiH  
     6  3.968700241710, 4.000000000000,      5*0.0D+00/                 NiH  
      DATA  Q_NiH/                                                      071215
     1  1.30103130D+00, 1.30108512D+00, 1.30200141D+00, 1.31624001D+00, NiH  
     2  1.33267774D+00, 1.36234593D+00, 1.46422226D+00, 1.57671715D+00, NiH  
     3  1.74526264D+00, 1.93989541D+00, 2.14694469D+00, 2.34196888D+00, NiH  
     4  2.57673310D+00, 3.07505125D+00, 3.28191948D+00, 3.51759856D+00, NiH  
     5  3.81270578D+00, 4.06062715D+00, 4.28195688D+00, 4.47447947D+00, NiH  
     6  4.59273397D+00, 4.66644676D+00,      5*0.0D+00/                 NiH  
      DATA TQ_CuH/                                                      071215
     1  0.699999789529, 0.730000036317, 0.777000153302, 0.891299897539, CuH  
     2  1.157599902898, 1.361500098944, 1.587800112917, 1.790899991182, CuH  
     3  1.991799813651, 2.265600357322, 2.514700345992, 2.628800075890, CuH  
     4  2.737899728980, 2.967500213363, 3.189299614614, 3.385900253881, CuH  
     5  3.517500398356, 3.635200229237, 3.843700227778, 3.936899826022, CuH  
     6  3.975399882833, 4.000000000000,      5*0.0D+00/                 CuH  
      DATA  Q_CuH/                                                      071215
     1  1.44227218D-02, 1.93513814D-02, 2.94040848D-02, 6.73877489D-02, CuH  
     2  2.23763224D-01, 3.82761202D-01, 5.79488945D-01, 7.66724751D-01, CuH  
     3  9.57952707D-01, 1.22439762D+00, 1.47050220D+00, 1.58466305D+00, CuH  
     4  1.69607102D+00, 1.94899170D+00, 2.23682739D+00, 2.53818572D+00, CuH  
     5  2.76619737D+00, 2.99155958D+00, 3.47088211D+00, 3.73352170D+00, CuH  
     6  3.85082661D+00, 3.92793744D+00,      5*0.0D+00/                 CuH  
      DATA TQ_ZnH/                                                      071215
     1  0.699999789529, 0.729900033922, 0.776700145986, 0.889999868442, ZnH  
     2  1.031600178729, 1.181799873391, 1.403700004197, 1.631100083770, ZnH  
     3  1.882100075068, 2.147799775924, 2.338000116803, 2.507500181260, ZnH  
     4  2.685399917483, 2.896499866497, 3.211700062255, 3.431400194121, ZnH  
     5  3.509800218545, 3.586000066019, 3.772600277792, 3.856499911235, ZnH  
     6  3.952099874234, 4.000000000000,      5*0.0D+00/                 ZnH  
      DATA  Q_ZnH/                                                      071215
     1  3.30304364D-01, 3.38304238D-01, 3.53531202D-01, 4.04049434D-01, ZnH  
     2  4.90163126D-01, 6.00461645D-01, 7.84983132D-01, 9.90038062D-01, ZnH  
     3  1.22707307D+00, 1.48494726D+00, 1.67222486D+00, 1.84092232D+00, ZnH  
     4  2.02327510D+00, 2.26088211D+00, 2.69922185D+00, 3.08422334D+00, ZnH  
     5  3.24079338D+00, 3.40247883D+00, 3.82182372D+00, 4.01131019D+00, ZnH  
     6  4.22468358D+00, 4.33152613D+00,      5*0.0D+00/                 ZnH  
      DATA TQ_GaH/                                                      071215
     1  0.699999789529, 0.726999964471, 0.769000008962, 0.872699938238, GaH  
     2  1.003400072568, 1.137199945369, 1.370399889597, 1.601999941074, GaH  
     3  1.862900063272, 2.138100277184, 2.352600195838, 2.539999957356, GaH  
     4  2.951399842789, 3.130300127362, 3.356199921013, 3.497199924725, GaH  
     5  3.629900119814, 3.714100123826, 3.801100211820, 3.922800163368, GaH  
     6  3.968700241710, 4.000000000000,      5*0.0D+00/                 GaH  
      DATA  Q_GaH/                                                      071215
     1  3.87788659D-02, 4.73656574D-02, 6.29453228D-02, 1.12549058D-01, GaH  
     2  1.93785569D-01, 2.91965377D-01, 4.85439972D-01, 6.94455135D-01, GaH  
     3  9.41150495D-01, 1.20856120D+00, 1.42000152D+00, 1.60684132D+00, GaH  
     4  2.05766078D+00, 2.29736564D+00, 2.65107218D+00, 2.90266934D+00, GaH  
     5  3.16944855D+00, 3.36165793D+00, 3.58315940D+00, 3.92478232D+00, GaH  
     6  4.05831429D+00, 4.14940430D+00,      5*0.0D+00/                 GaH  
      DATA TQ_GeH/                                                      071215
     1  0.699999789529, 0.866399967758, 1.036100058497, 1.191400009990, GeH  
     2  1.361200106641, 1.689800066946, 1.874499995348, 2.054800308863, GeH  
     3  2.277499921351, 2.467200188611, 2.624399975033, 2.776799964774, GeH  
     4  2.906899761467, 3.034399849209, 3.197799746650, 3.328599884413, GeH  
     5  3.426700227116, 3.523800180496, 3.640500292869, 3.758900198110, GeH  
     6  3.845900277270, 3.918200054911, 3.968600239469, 4.000000000000, GeH  
     7       3*0.0D+00/                                                 GeH  
      DATA  Q_GeH/                                                      071215
     1  1.07953964D+00, 1.08320616D+00, 1.10045049D+00, 1.13921745D+00, GeH  
     2  1.21122752D+00, 1.42271142D+00, 1.57001638D+00, 1.72608279D+00, GeH  
     3  1.92991052D+00, 2.10994451D+00, 2.26306303D+00, 2.41787761D+00, GeH  
     4  2.56001558D+00, 2.71309581D+00, 2.93423639D+00, 3.13321806D+00, GeH  
     5  3.29569737D+00, 3.46858248D+00, 3.69550609D+00, 3.95258059D+00, GeH  
     6  4.15909067D+00, 4.33845603D+00, 4.46508534D+00, 4.54379732D+00, GeH  
     7       3*0.0D+00/                                                 GeH  
      DATA TQ_AsH/                                                      071215
     1  0.699999789529, 0.731499994541, 0.781100196333, 0.899800087786, AsH  
     2  1.047300116266, 1.189500033715, 1.463199947499, 1.746499914815, AsH  
     3  1.983399926775, 2.282599796689, 2.609799658581, 2.764400339676, AsH  
     4  2.906499752472, 3.042200025296, 3.186999781414, 3.362599703036, AsH  
     5  3.526000020566, 3.682700109205, 3.838400109691, 3.934000035036, AsH  
     6  3.974299961872, 4.000000000000,      5*0.0D+00/                 AsH  
      DATA  Q_AsH/                                                      071215
     1  4.97527434D-01, 5.04169933D-01, 5.17487394D-01, 5.64471651D-01, AsH  
     2  6.49476237D-01, 7.50966731D-01, 9.77401929D-01, 1.23556915D+00, AsH  
     3  1.46120469D+00, 1.75297211D+00, 2.07709136D+00, 2.23327069D+00, AsH  
     4  2.38312222D+00, 2.53748077D+00, 2.71955343D+00, 2.96782248D+00, AsH  
     5  3.22510739D+00, 3.49306865D+00, 3.77828594D+00, 3.96292346D+00, AsH  
     6  4.04292255D+00, 4.09456564D+00,      5*0.0D+00/                 AsH  
      DATA TQ_SeH/                                                      071215
     1  0.699999789529, 0.875099996362, 1.087200097232, 1.191300012428, SeH  
     2  1.300200166544, 1.478699983936, 1.680599851630, 2.023900177064, SeH  
     3  2.377900076463, 2.516200377335, 2.650599630534, 2.950799828966, SeH  
     4  3.092900209041, 3.234399638575, 3.407999789345, 3.573899779715, SeH  
     5  3.854200073392, 3.941199629687, 4.000000000000,      8*0.0D+00/ SeH  
      DATA  Q_SeH/                                                      071215
     1  1.07927669D+00, 1.08102760D+00, 1.09746883D+00, 1.11863030D+00, SeH  
     2  1.15311481D+00, 1.23671016D+00, 1.36499744D+00, 1.63755689D+00, SeH  
     3  1.95799915D+00, 2.08912183D+00, 2.21959477D+00, 2.54302975D+00, SeH  
     4  2.72857275D+00, 2.93928428D+00, 3.22766364D+00, 3.52477988D+00, SeH  
     5  4.05394672D+00, 4.22242747D+00, 4.33711440D+00,      8*0.0D+00/ SeH  
      DATA TQ_HBr/                                                      071215
     1  0.699999789529, 0.731999980615, 0.782400160721, 0.904199995448, HBr  
     2  1.180099837995, 1.377300040070, 1.599699899621, 2.002100048221, HBr  
     3  2.565300035323, 2.764200334845, 2.960300047820, 3.103000163274, HBr  
     4  3.248999979416, 3.552700252888, 3.686899816684, 3.834000013854, HBr  
     5  3.932700128733, 3.973799997799, 4.000000000000,      8*0.0D+00/ HBr  
      DATA  Q_HBr/                                                      071215
     1  1.06701154D-02, 1.49255534D-02, 2.40633060D-02, 6.09830662D-02, HBr  
     2  2.19378780D-01, 3.72133138D-01, 5.64406983D-01, 9.40179501D-01, HBr  
     3  1.49156656D+00, 1.69014844D+00, 1.89339823D+00, 2.05404787D+00, HBr  
     4  2.23676019D+00, 2.69027499D+00, 2.92389839D+00, 3.20659031D+00, HBr  
     5  3.41495936D+00, 3.50640082D+00, 3.56605302D+00,      8*0.0D+00/ HBr  
      DATA TQ_RbH/                                                      071215
     1  0.699999789529, 0.726199945312, 0.766800066990, 0.871499909176, RbH  
     2  1.151400061434, 1.484799903037, 1.939200022183, 2.134900203635, RbH  
     3  2.317599882949, 2.578399889505, 2.790299984848, 3.117499817741, RbH  
     4  3.236399687149, 3.363499726899, 3.469600243109, 3.574199786882, RbH  
     5  3.797000124811, 3.920400110046, 3.968900246192, 4.000000000000, RbH  
     6       7*0.0D+00/                                                 RbH  
      DATA  Q_RbH/                                                      071215
     1  1.96028290D-01, 2.14319939D-01, 2.43682938D-01, 3.24219972D-01, RbH  
     2  5.63048360D-01, 8.72698076D-01, 1.31428125D+00, 1.50783068D+00, RbH  
     3  1.69036484D+00, 1.96485800D+00, 2.22151078D+00, 2.70900048D+00, RbH  
     4  2.91478860D+00, 3.15258411D+00, 3.36747332D+00, 3.59511411D+00, RbH  
     5  4.10463745D+00, 4.37149881D+00, 4.46966152D+00, 4.53040794D+00, RbH  
     6       7*0.0D+00/                                                 RbH  
      DATA TQ_SrH/                                                      071215
     1  0.699999789529, 0.725899938127, 0.765900090728, 0.869199893939, SrH  
     2  0.997299943219, 1.143699958354, 1.327599989563, 1.527100155270, SrH  
     3  1.769799973931, 2.019000441453, 2.236799712231, 2.434399976273, SrH  
     4  2.683300069479, 2.894300030047, 3.065500030056, 3.307100345321, SrH  
     5  3.455199919340, 3.594100250330, 3.693999683549, 3.787799926186, SrH  
     6  3.933100099903, 4.000000000000,      5*0.0D+00/                 SrH  
      DATA  Q_SrH/                                                      071215
     1  4.41523294D-01, 4.57569803D-01, 4.83658782D-01, 5.57173574D-01, SrH  
     2  6.57690965D-01, 7.81512834D-01, 9.46242806D-01, 1.13251569D+00, SrH  
     3  1.36566054D+00, 1.60950180D+00, 1.82488478D+00, 2.02237651D+00, SrH  
     4  2.28421507D+00, 2.53842862D+00, 2.77745341D+00, 3.16818840D+00, SrH  
     5  3.43921281D+00, 3.72134940D+00, 3.94732417D+00, 4.17893724D+00, SrH  
     6  4.56556559D+00, 4.74879780D+00,      5*0.0D+00/                 SrH  
      DATA TQ_AgH/                                                      071215
     1  0.699999789529, 0.730300027961, 0.777500165496, 0.891999913206, AgH  
     2  1.036800039794, 1.191500007551, 1.419400066533, 1.654199986231, AgH  
     3  1.906299923921, 2.172000106739, 2.368199844332, 2.539599947628, AgH  
     4  2.710700071835, 2.930200278154, 3.131900164029, 3.349500377895, AgH  
     5  3.459000011230, 3.575699822718, 3.666799993908, 3.762800291873, AgH  
     6  3.875300012858, 4.000000000000,      5*0.0D+00/                 AgH  
      DATA  Q_AgH/                                                      071215
     1  3.28155311D-02, 4.15553800D-02, 5.79573608D-02, 1.11520134D-01, AgH  
     2  2.02217398D-01, 3.17977897D-01, 5.09843552D-01, 7.23282288D-01, AgH  
     3  9.62509981D-01, 1.22100937D+00, 1.41441888D+00, 1.58494551D+00, AgH  
     4  1.75923881D+00, 2.00233256D+00, 2.26269657D+00, 2.59442022D+00, AgH  
     5  2.78266743D+00, 3.00103042D+00, 3.18638684D+00, 3.39674333D+00, AgH  
     6  3.65885271D+00, 3.95893120D+00,      5*0.0D+00/                 AgH  
      DATA TQ_CdH/                                                      071215
     1  0.699999789529, 0.726499952497, 0.767500048526, 0.869899875484, CdH  
     2  0.998199962146, 1.141999920879, 1.360800116903, 1.590300159943, CdH  
     3  1.840000148372, 2.103900112667, 2.300000195710, 2.478899605954, CdH  
     4  2.852700184153, 2.978599643352, 3.106199931974, 3.412799903470, CdH  
     5  3.514700332967, 3.625000010825, 3.847200306516, 3.938899681874, CdH  
     6  3.976099832535, 4.000000000000,      5*0.0D+00/                 CdH  
      DATA  Q_CdH/                                                      071215
     1  3.58611461D-01, 3.69276465D-01, 3.87837413D-01, 4.44076003D-01, CdH  
     2  5.30155032D-01, 6.40646229D-01, 8.26717303D-01, 1.03610245D+00, CdH  
     3  1.27341692D+00, 1.53036622D+00, 1.72380362D+00, 1.90222055D+00, CdH  
     4  2.30622115D+00, 2.46499285D+00, 2.64179052D+00, 3.12979144D+00, CdH  
     5  3.30981431D+00, 3.51462463D+00, 3.96863254D+00, 4.17854430D+00, CdH  
     6  4.26802438D+00, 4.32675566D+00,      5*0.0D+00/                 CdH  
      DATA TQ_InH/                                                      071215
     1  0.699999789529, 0.724899914179, 0.763500154031, 0.860700118031, InH  
     2  0.985099891497, 1.118599917613, 1.351099939936, 1.559899940427, InH  
     3  1.821600147773, 2.096700307657, 2.317899861409, 2.519000435843, InH  
     4  2.781899772832, 2.887300283697, 2.998599967235, 3.364499753415, InH  
     5  3.503400069339, 3.645599932592, 3.710700050487, 3.777799901493, InH  
     6  3.913499934150, 3.965500169995, 4.000000000000,      4*0.0D+00/ InH  
      DATA  Q_InH/                                                      071215
     1  7.14157862D-02, 8.26678243D-02, 1.01849617D-01, 1.58444523D-01, InH  
     2  2.44329030D-01, 3.48114738D-01, 5.47194613D-01, 7.38301936D-01, InH  
     3  9.87246925D-01, 1.25537800D+00, 1.47378116D+00, 1.67462371D+00, InH  
     4  1.95225000D+00, 2.07563529D+00, 2.21730948D+00, 2.77900103D+00, InH  
     5  3.03559913D+00, 3.34469336D+00, 3.51366674D+00, 3.71050015D+00, InH  
     6  4.16288078D+00, 4.34372925D+00, 4.46276818D+00,      4*0.0D+00/ InH  
      DATA TQ_SnH/                                                      071215
     1  0.699999789529, 0.784700097716, 0.938200004258, 1.073700086258, SnH  
     2  1.210100011995, 1.595100027013, 1.965099922671, 2.362599708658, SnH  
     3  2.527699884937, 2.679700301537, 2.832499998691, 2.992099815112, SnH  
     4  3.167300008958, 3.343100233940, 3.594800265797, 3.793000040950, SnH  
     5  3.916500011231, 3.967500214817, 4.000000000000,      8*0.0D+00/ SnH  
      DATA  Q_SnH/                                                      071215
     1  1.08080012D+00, 1.08393897D+00, 1.10012849D+00, 1.13223655D+00, SnH  
     2  1.18457754D+00, 1.42366294D+00, 1.73161523D+00, 2.10140667D+00, SnH  
     3  2.26170558D+00, 2.41484026D+00, 2.58022432D+00, 2.77274269D+00, SnH  
     4  3.01249642D+00, 3.28272832D+00, 3.71393874D+00, 4.08509942D+00, SnH  
     5  4.33056064D+00, 4.43493831D+00, 4.50219239D+00,      8*0.0D+00/ SnH  
      DATA TQ_SbH/                                                      071215
     1  0.699999789529, 0.724799911784, 0.763500154031, 0.859600127310, SbH  
     2  0.975999916507, 1.103699995613, 1.281500180861, 1.469599790058, SbH  
     3  1.688800043542, 1.978799980800, 2.104200091284, 2.227199743875, SbH  
     4  2.433200060025, 2.628800075890, 2.813200186611, 2.987499709774, SbH  
     5  3.313300174448, 3.496499908388, 3.662599908173, 3.869299876967, SbH  
     6  3.947899780964, 3.980099554550, 4.000000000000,      4*0.0D+00/ SbH  
      DATA  Q_SbH/                                                      071215
     1  5.24428426D-01, 5.33294765D-01, 5.49007101D-01, 5.97323324D-01, SbH  
     2  6.70513069D-01, 7.63913089D-01, 9.08965967D-01, 1.07474977D+00, SbH  
     3  1.27786431D+00, 1.55618225D+00, 1.67883238D+00, 1.80108115D+00, SbH  
     4  2.01666901D+00, 2.24592177D+00, 2.49098201D+00, 2.74904781D+00, SbH  
     5  3.29223192D+00, 3.62564222D+00, 3.94233938D+00, 4.36115687D+00, SbH  
     6  4.53209416D+00, 4.60449480D+00, 4.64994374D+00,      4*0.0D+00/ SbH  
      DATA TQ_TeH/                                                      071215
     1  0.699999789529, 0.795600078935, 0.976299922578, 1.130300114367, TeH  
     2  1.287700010044, 1.650900074038, 1.847199974650, 2.034299831074, TeH  
     3  2.250000002175, 2.431500178672, 2.587900110534, 2.751600026576, TeH  
     4  3.053600274043, 3.183000071501, 3.312600224232, 3.443699662454, TeH  
     5  3.685399921156, 3.862899730735, 4.000000000000,      8*0.0D+00/ TeH  
      DATA  Q_TeH/                                                      071215
     1  1.08040114D+00, 1.08347532D+00, 1.10337411D+00, 1.14440005D+00, TeH  
     2  1.21262967D+00, 1.45111691D+00, 1.61114867D+00, 1.77572255D+00, TeH  
     3  1.97490819D+00, 2.14775344D+00, 2.30063606D+00, 2.46901909D+00, TeH  
     4  2.83214113D+00, 3.02047478D+00, 3.23342887D+00, 3.47243923D+00, TeH  
     5  3.95618195D+00, 4.32728851D+00, 4.61533097D+00,      8*0.0D+00/ TeH  
      DATA TQ_HI/                                                       071215
     1  0.699999789529, 0.731000008466, 0.779500214272, 0.896800020640, HI   
     2  1.044800060797, 1.204799902941, 1.451099832787, 1.702499873141, HI   
     3  1.965399914860, 2.241099813193, 2.457299965951, 2.637700285839, HI   
     4  2.810400386561, 2.981099565926, 3.115999785654, 3.253900087296, HI   
     5  3.450399803269, 3.623399975236, 3.846400288519, 3.940299609366, HI   
     6  3.976799782237, 4.000000000000,      5*0.0D+00/                 HI   
      DATA  Q_HI/                                                       071215
     1  3.29344413D-02, 4.19132941D-02, 5.89007819D-02, 1.14404497D-01, HI   
     2  2.08111424D-01, 3.28928823D-01, 5.38280690D-01, 7.68820926D-01, HI   
     3  1.01983443D+00, 1.28903238D+00, 1.50260906D+00, 1.68217457D+00, HI   
     4  1.85740446D+00, 2.04132178D+00, 2.20144803D+00, 2.38353248D+00, HI   
     5  2.67884717D+00, 2.97467468D+00, 3.41388791D+00, 3.62286293D+00, HI   
     6  3.70789038D+00, 3.76295909D+00,      5*0.0D+00/                 HI   
      DATA TQ_CsH/                                                      071215
     1  0.699999789529, 0.727499976445, 0.770099985025, 0.879600105345, CsH  
     2  1.182999898376, 1.542500027751, 1.956199959560, 2.109799692137, CsH  
     3  2.259300218702, 2.411899896494, 2.652099667294, 2.805000309819, CsH  
     4  2.962600100702, 3.186399824927, 3.387500294398, 3.470200237426, CsH  
     5  3.551800234105, 3.629700115365, 3.696499735283, 3.757500164930, CsH  
     6  3.817099894333, 3.925200216689, 3.969800266362, 3.987799727422, CsH  
     7  4.000000000000,      2*0.0D+00/                                 CsH  
      DATA  Q_CsH/                                                      071215
     1  2.28832674D-01, 2.48986131D-01, 2.81171730D-01, 3.68452547D-01, CsH  
     2  6.33690859D-01, 9.72293657D-01, 1.37610679D+00, 1.52788960D+00, CsH  
     3  1.67656562D+00, 1.83137248D+00, 2.09585328D+00, 2.28852105D+00, CsH  
     4  2.51183896D+00, 2.87152847D+00, 3.23393307D+00, 3.39353862D+00, CsH  
     5  3.55822216D+00, 3.72744211D+00, 3.89435324D+00, 4.08117874D+00, CsH  
     6  4.30731958D+00, 4.80686557D+00, 5.02538307D+00, 5.11291249D+00, CsH  
     7  5.17174294D+00,      2*0.0D+00/                                 CsH  
      DATA TQ_BaH/                                                      071215
     1  0.699999789529, 0.724299899810, 0.762100190958, 0.859500125016, BaH  
     2  1.118499920209, 1.334300033514, 1.539000002062, 2.042800033501, BaH  
     3  2.230699564296, 2.411899896494, 2.745799897810, 2.878700077613, BaH  
     4  3.026100021026, 3.154799712782, 3.288199932488, 3.420200096931, BaH  
     5  3.496999920057, 3.574599796438, 3.646799847820, 3.714700136768, BaH  
     6  3.863899753583, 3.946999760643, 3.979699573861, 4.000000000000, BaH  
     7       3*0.0D+00/                                                 BaH  
      DATA  Q_BaH/                                                      071215
     1  1.62599277D-01, 1.78489722D-01, 2.04251207D-01, 2.75584825D-01, BaH  
     2  4.89384488D-01, 6.83973738D-01, 8.76450016D-01, 1.36677738D+00, BaH  
     3  1.55284543D+00, 1.73393259D+00, 2.09158368D+00, 2.25570938D+00, BaH  
     4  2.45867869D+00, 2.65483692D+00, 2.87806249D+00, 3.12789864D+00, BaH  
     5  3.29630764D+00, 3.49285104D+00, 3.70444152D+00, 3.92821236D+00, BaH  
     6  4.47135365D+00, 4.77688113D+00, 4.89373483D+00, 4.96494860D+00, BaH  
     7       3*0.0D+00/                                                 BaH  
      DATA TQ_YbH/                                                      071215
     1  0.699999789529, 0.725699933338, 0.765500101279, 0.867199946667, YbH  
     2  0.997799953734, 1.146600022282, 1.339400150455, 1.545600091341, YbH  
     3  1.787700029470, 2.032999799525, 2.241299817440, 2.433600032107, YbH  
     4  2.568699793336, 2.699299792087, 2.900099608551, 3.282899810654, YbH  
     5  3.419600082400, 3.559400392714, 3.638300293145, 3.715400151868, YbH  
     6  3.902899684777, 3.961900089316, 4.000000000000,      4*0.0D+00/ YbH  
      DATA  Q_YbH/                                                      071215
     1  4.20672580D-01, 4.35577302D-01, 4.60105835D-01, 5.29498445D-01, YbH  
     2  6.29299377D-01, 7.53188873D-01, 9.24481729D-01, 1.11633071D+00, YbH  
     3  1.34850862D+00, 1.58829079D+00, 1.79412060D+00, 1.98606643D+00, YbH  
     4  2.12436758D+00, 2.26523295D+00, 2.50700708D+00, 3.08465557D+00, YbH  
     5  3.33135672D+00, 3.61243098D+00, 3.78813860D+00, 3.97275412D+00, YbH  
     6  4.45722163D+00, 4.61229615D+00, 4.71148632D+00,      4*0.0D+00/ YbH  
      DATA TQ_PtH/                                                      071215
     1  0.699999789529, 0.816800093391, 0.951999846135, 1.262700153368, PtH  
     2  1.384899979696, 1.512400164694, 1.991499806834, 2.269100433774, PtH  
     3  2.569699722164, 2.695499707283, 2.810800357997, 3.087200084785, PtH  
     4  3.253000067652, 3.520300434930, 3.645199960849, 3.805200301695, PtH  
     5  3.860499675897, 3.913699939289, 3.967000203611, 4.000000000000, PtH  
     6       7*0.0D+00/                                                 PtH  
      DATA  Q_PtH/                                                      071215
     1  1.30103301D+00, 1.30108369D+00, 1.30168186D+00, 1.32229433D+00, PtH  
     2  1.34921707D+00, 1.39335934D+00, 1.68725668D+00, 1.91761684D+00, PtH  
     3  2.19174329D+00, 2.31169554D+00, 2.42526840D+00, 2.72790940D+00, PtH  
     4  2.94317608D+00, 3.35606494D+00, 3.57949111D+00, 3.90931915D+00, PtH  
     5  4.03984808D+00, 4.17430542D+00, 4.31653325D+00, 4.40732846D+00, PtH  
     6       7*0.0D+00/                                                 PtH  
      DATA TQ_AuH/                                                      071215
     1  0.699999789529, 0.732199975045, 0.782900147024, 0.904000000058, AuH  
     2  1.056700026279, 1.205699922186, 1.454099897989, 1.699999814320, AuH  
     3  1.968199841958, 2.252400058053, 2.465000138991, 2.639700333227, AuH  
     4  2.806900349835, 2.987399707526, 3.203699881125, 3.448299758209, AuH  
     5  3.548300155879, 3.650099624097, 3.754400091461, 3.850200355406, AuH  
     6  3.941199629687, 3.977199753496, 4.000000000000,      4*0.0D+00/ AuH  
      DATA  Q_AuH/                                                      071215
     1  2.11808243D-02, 2.81626737D-02, 4.21708696D-02, 9.13773060D-02, AuH  
     2  1.81251946D-01, 2.89390533D-01, 4.95818368D-01, 7.18732567D-01, AuH  
     3  9.73263251D-01, 1.25002787D+00, 1.45978728D+00, 1.63357377D+00, AuH  
     4  1.80318386D+00, 1.99813913D+00, 2.26546574D+00, 2.62827512D+00, AuH  
     5  2.79683463D+00, 2.98230001D+00, 3.19054906D+00, 3.40230252D+00, AuH  
     6  3.62333976D+00, 3.71569602D+00, 3.77539475D+00,      4*0.0D+00/ AuH  
      DATA TQ_HgH/                                                      071215
     1  0.699999789529, 0.727499976445, 0.769999982586, 0.875700010893, HgH  
     2  1.011500177598, 1.160499852852, 1.369799886005, 1.590200162712, HgH  
     3  1.831099938897, 2.086200061175, 2.270600410867, 2.439299634288, HgH  
     4  2.807200356154, 3.031099771257, 3.290999995386, 3.493999850041, HgH  
     5  3.589100138648, 3.686099872402, 3.850000369506, 4.000000000000, HgH  
     6       7*0.0D+00/                                                 HgH  
      DATA  Q_HgH/                                                      071215
     1  3.56598533D-01, 3.67464105D-01, 3.86481655D-01, 4.44368422D-01, HgH  
     2  5.35954237D-01, 6.51333796D-01, 8.30000370D-01, 1.03108620D+00, HgH  
     3  1.25980992D+00, 1.50795808D+00, 1.68974098D+00, 1.85800315D+00, HgH  
     4  2.25617614D+00, 2.54968159D+00, 2.95193742D+00, 3.30812297D+00, HgH  
     5  3.48726719D+00, 3.67797954D+00, 4.01421422D+00, 4.32863906D+00, HgH  
     6       7*0.0D+00/                                                 HgH  
      DATA TQ_TlH/                                                      071215
     1  0.699999789529, 0.727199969261, 0.769300001049, 0.875199998784, TlH  
     2  1.014700101148, 1.169400054392, 1.369799886005, 1.583900014824, TlH  
     3  1.832199964787, 2.087000077594, 2.310600385538, 2.513100312559, TlH  
     4  2.784399836172, 2.896399873931, 3.020700402838, 3.199899795865, TlH  
     5  3.385600246284, 3.455299921758, 3.528499838828, 3.592900223814, TlH  
     6  3.665499967371, 3.727399738317, 3.794000061915, 3.918200054911, TlH  
     7  3.967000203611, 4.000000000000,      1*0.0D+00/                 TlH  
      DATA  Q_TlH/                                                      071215
     1  7.91720506D-02, 9.21886530D-02, 1.14324411D-01, 1.79206567D-01, TlH  
     2  2.79707360D-01, 4.04496810D-01, 5.79571852D-01, 7.77128425D-01, TlH  
     3  1.01433915D+00, 1.26304417D+00, 1.48399562D+00, 1.68669518D+00, TlH  
     4  1.97657607D+00, 2.11094161D+00, 2.27471268D+00, 2.54081625D+00, TlH  
     5  2.85507211D+00, 2.98414994D+00, 3.12877523D+00, 3.26771101D+00, TlH  
     6  3.44549905D+00, 3.62239956D+00, 3.84309319D+00, 4.31957014D+00, TlH  
     7  4.51623424D+00, 4.64855968D+00,      1*0.0D+00/                 TlH  
      DATA TQ_PbH/                                                      071215
     1  0.699999789529, 0.738299805155, 0.800100175966, 0.959000029290, PbH  
     2  1.115699992889, 1.277400159711, 1.655799943658, 2.037099899025, PbH  
     3  2.258700204732, 2.450699818627, 2.624399975033, 2.804800305606, PbH  
     4  2.924100164439, 3.042000021314, 3.326099830207, 3.419900090294, PbH  
     5  3.515700356320, 3.613699746762, 3.715600156182, 3.794500072398, PbH  
     6  3.895199973359, 3.957599994330, 4.000000000000,      4*0.0D+00/ PbH  
      DATA  Q_PbH/                                                      071215
     1  1.08178835D+00, 1.08336594D+00, 1.08743606D+00, 1.11110805D+00, PbH  
     2  1.15948108D+00, 1.23657244D+00, 1.49845906D+00, 1.82882003D+00, PbH  
     3  2.03600297D+00, 2.22074910D+00, 2.39284709D+00, 2.58385723D+00, PbH  
     4  2.72304581D+00, 2.87430324D+00, 3.30264700D+00, 3.46481959D+00, PbH  
     5  3.64250326D+00, 3.83967823D+00, 4.06465862D+00, 4.25250583D+00, PbH  
     6  4.50300569D+00, 4.65934948D+00, 4.76429709D+00,      4*0.0D+00/ PbH  
      DATA TQ_BiH/                                                      071215
     1  0.699999789529, 0.727399974050, 0.769799987861, 0.875800013315, BiH  
     2  1.013000141762, 1.165899975135, 1.383100023504, 1.611300107614, BiH  
     3  1.861000111240, 2.115799795849, 2.336900090728, 2.534999835748, BiH  
     4  2.676300224550, 2.975699854948, 3.262300275660, 3.381100132331, BiH  
     5  3.528299853367, 3.621499932975, 3.717400195008, 3.798300152066, BiH  
     6  3.905999754394, 3.962000091557, 3.985399673540, 4.000000000000, BiH  
     7       3*0.0D+00/                                                 BiH  
      DATA  Q_BiH/                                                      071215
     1  5.43368332D-01, 5.55308360D-01, 5.75918344D-01, 6.37427744D-01, BiH  
     2  7.33132929D-01, 8.54028605D-01, 1.04204157D+00, 1.25215843D+00, BiH  
     3  1.49059933D+00, 1.73916940D+00, 1.95742751D+00, 2.15489000D+00, BiH  
     4  2.29913406D+00, 2.63939375D+00, 3.05365831D+00, 3.26192945D+00, BiH  
     5  3.55492160D+00, 3.76139850D+00, 3.99068303D+00, 4.19530622D+00, BiH  
     6  4.47581673D+00, 4.62110818D+00, 4.68103928D+00, 4.71810923D+00, BiH  
     7       3*0.0D+00/                                                 BiH  
      DATA TQ_HeHp/                                                     071215
     1  0.699999789529, 0.876000018159, 1.056600028516, 1.201499832375, HeHp 
     2  1.373799963743, 1.490299788711, 1.598899921776, 1.738199823705, HeHp 
     3  1.868299926941, 2.081599966766, 2.272100304450, 2.580799944189, HeHp 
     4  2.709100033710, 2.835500066226, 2.981199568174, 3.158299792496, HeHp 
     5  3.382500167783, 3.488399721202, 3.589900157392, 3.819199742725, HeHp 
     6  3.928800296672, 3.972100119951, 4.000000000000,      4*0.0D+00/ HeHp 
      DATA  Q_HeHp/                                                     071215
     1  5.65708114D-09, 3.45824192D-06, 2.72771632D-04, 3.00695921D-03, HeHp 
     2  2.14820907D-02, 5.41631197D-02, 1.02954439D-01, 1.87535445D-01, HeHp 
     3  2.82109497D-01, 4.57394909D-01, 6.27378071D-01, 9.18527112D-01, HeHp 
     4  1.04340017D+00, 1.16882668D+00, 1.31916913D+00, 1.51987157D+00, HeHp 
     5  1.82628123D+00, 2.00114540D+00, 2.19061699D+00, 2.66617442D+00, HeHp 
     6  2.88783168D+00, 2.97099087D+00, 3.02305108D+00,      4*0.0D+00/ HeHp 
      DATA TQ_BeHp/                                                     071215
     1  0.699999789529, 0.734199919344, 0.789899955268, 0.938500010103, BeHp 
     2  1.064400054251, 1.210999992731, 1.375800007358, 1.578299960377, BeHp 
     3  1.775400104098, 1.990899793199, 2.258700204732, 2.526000011179, BeHp 
     4  2.712900121409, 2.909699824432, 3.127700063155, 3.379800100371, BeHp 
     5  3.551400225758, 3.674100156820, 3.778399858074, 3.858899742027, BeHp 
     6  3.936599847644, 4.000000000000,      5*0.0D+00/                 BeHp 
      DATA  Q_BeHp/                                                     071215
     1  2.87061247D-03, 4.55227900D-03, 8.94323924D-03, 3.66204297D-02, BeHp 
     2  8.47300469D-02, 1.68202694D-01, 2.86331106D-01, 4.52532277D-01, BeHp 
     3  6.27829726D-01, 8.28757281D-01, 1.08648598D+00, 1.34901216D+00, BeHp 
     4  1.53567278D+00, 1.74161485D+00, 1.99940977D+00, 2.35895858D+00, BeHp 
     5  2.64726411D+00, 2.87948298D+00, 3.09838266D+00, 3.28063723D+00, BeHp 
     6  3.46399301D+00, 3.61570656D+00,      5*0.0D+00/                 BeHp 
      DATA TQ_CHp/                                                      071215
     1  0.699999789529, 0.848299948110, 1.009500202762, 1.144799982602, CHp  
     2  1.286000056881, 1.457799978405, 1.638399901975, 1.847699962586, CHp  
     3  2.071399718333, 2.330499939019, 2.604700030912, 2.801700240316, CHp  
     4  3.046200104941, 3.170800091849, 3.287499916397, 3.428500263167, CHp  
     5  3.551800234105, 3.629200104244, 3.702899873736, 3.887600299577, CHp  
     6  3.954899935374, 4.000000000000,      5*0.0D+00/                 CHp  
      DATA  Q_CHp/                                                      071215
     1  4.37362567D-04, 4.40323284D-03, 2.50736663D-02, 6.84143450D-02, CHp  
     2  1.41944363D-01, 2.59877252D-01, 4.04296105D-01, 5.87332544D-01, CHp  
     3  7.94168118D-01, 1.04210421D+00, 1.30996623D+00, 1.50548667D+00, CHp  
     4  1.76256559D+00, 1.90876646D+00, 2.06126637D+00, 2.27703895D+00, CHp  
     5  2.50960850D+00, 2.68145896D+00, 2.86288601D+00, 3.37064440D+00, CHp  
     6  3.56489043D+00, 3.69549114D+00,      5*0.0D+00/                 CHp  
      DATA TQ_NHp/                                                      071215
     1  0.699999789529, 0.875600008471, 1.055100062075, 1.196099895392, NHp  
     2  1.345400029923, 1.501400025088, 1.794200059021, 2.130400100205, NHp  
     3  2.306300342441, 2.509400231817, 2.639400326119, 2.770100467443, NHp  
     4  2.933200067623, 3.176100205463, 3.342000209198, 3.579599915890, NHp  
     5  3.830299933263, 3.930300301710, 3.972800069653, 4.000000000000, NHp  
     6       7*0.0D+00/                                                 NHp  
      DATA  Q_NHp/                                                      071215
     1  1.07918127D+00, 1.07918852D+00, 1.07953891D+00, 1.08212469D+00, NHp  
     2  1.09362920D+00, 1.12526613D+00, 1.25439932D+00, 1.51019919D+00, NHp  
     3  1.69086571D+00, 1.92674923D+00, 2.08354942D+00, 2.24106143D+00, NHp  
     4  2.43709604D+00, 2.74285216D+00, 2.97490816D+00, 3.35115729D+00, NHp  
     5  3.81470954D+00, 4.02514600D+00, 4.11890453D+00, 4.17998927D+00, NHp  
     6       7*0.0D+00/                                                 NHp  
      DATA TQ_OHp/                                                      071215
     1  0.699999789529, 0.878200071439, 1.073300096805, 1.210999992731, OHp  
     2  1.351699953250, 1.476099923055, 1.607200070454, 1.930399833611, OHp  
     3  2.382900188906, 2.604400052814, 2.795300098076, 2.968700240954, OHp  
     4  3.185999853936, 3.546000100867, 3.807200345537, 3.919600090882, OHp  
     5  3.968900246192, 4.000000000000,      9*0.0D+00/                 OHp  
      DATA  Q_OHp/                                                      071215
     1  3.01134936D-01, 3.03528073D-01, 3.24483967D-01, 3.67295815D-01, OHp  
     2  4.39576411D-01, 5.22016544D-01, 6.21418265D-01, 8.98889034D-01, OHp  
     3  1.32572550D+00, 1.54224822D+00, 1.73157864D+00, 1.90878296D+00, OHp  
     4  2.15245306D+00, 2.66206167D+00, 3.15282107D+00, 3.40888674D+00, OHp  
     5  3.52889801D+00, 3.60628402D+00,      9*0.0D+00/                 OHp  
      DATA TQ_HFp/                                                      071215
     1  0.699999789529, 0.856000044729, 1.030200216135, 1.159099864543, HFp  
     2  1.284300103718, 1.470999803633, 1.674599961865, 1.932599880754, HFp  
     3  2.188499697080, 2.567299892978, 2.710400065075, 2.848700356279, HFp  
     4  3.173700154015, 3.479599563587, 3.787999930901, 3.911999895610, HFp  
     5  3.965600172236, 4.000000000000,      9*0.0D+00/                 HFp  
      DATA  Q_HFp/                                                      071215
     1  1.07918145D+00, 1.07920080D+00, 1.07981285D+00, 1.08302373D+00, HFp  
     2  1.09331555D+00, 1.13372049D+00, 1.21999558D+00, 1.38319730D+00, HFp  
     3  1.58519036D+00, 1.92530697D+00, 2.06142200D+00, 2.19628431D+00, HFp  
     4  2.53935329D+00, 2.93412638D+00, 3.44176876D+00, 3.68326193D+00, HFp  
     5  3.79366201D+00, 3.86594247D+00,      9*0.0D+00/                 HFp  
      DATA TQ_NeHp/                                                     071215
     1  0.699999789529, 0.885099989271, 1.093000089364, 1.221199825944, NeHp 
     2  1.348099962778, 1.481599977263, 1.618799920079, 1.779400204528, NeHp 
     3  1.955699946942, 2.184399991073, 2.421600123395, 2.605899943305, NeHp 
     4  2.781299757630, 2.939899597437, 3.131000143404, 3.465800160895, NeHp 
     5  3.571199715211, 3.678400259170, 3.807100343345, 3.903699702742, NeHp 
     6  3.963000113968, 4.000000000000,      5*0.0D+00/                 NeHp 
      DATA  Q_NeHp/                                                     071215
     1  6.20758530D-05, 1.95785361D-03, 2.26570614D-02, 6.08590429D-02, NeHp 
     2  1.22619663D-01, 2.08207559D-01, 3.10941283D-01, 4.43632355D-01, NeHp 
     3  5.99533730D-01, 8.11938486D-01, 1.03971625D+00, 1.22006215D+00, NeHp 
     4  1.39435461D+00, 1.55745765D+00, 1.77250826D+00, 2.24510398D+00, NeHp 
     5  2.42840932D+00, 2.63480419D+00, 2.90504909D+00, 3.11450861D+00, NeHp 
     6  3.24131590D+00, 3.31868760D+00,      5*0.0D+00/                 NeHp 
      DATA TQ_MgHp/                                                     071215
     1  0.699999789529, 0.729900033922, 0.776500141108, 0.889699875840, MgHp 
     2  1.032100165370, 1.185299946265, 1.415399974569, 1.650700079360, MgHp 
     3  1.905099954000, 2.174500169189, 2.367399824950, 2.537599898984, MgHp 
     4  2.706299964171, 2.905899738979, 3.117699822019, 3.373799977235, MgHp 
     5  3.554300286279, 3.659699848099, 3.761200253191, 3.865999801566, MgHp 
     6  3.948599796769, 4.000000000000,      5*0.0D+00/                 MgHp 
      DATA  Q_MgHp/                                                     071215
     1  3.37704362D-02, 4.25457175D-02, 5.89432868D-02, 1.12206015D-01, MgHp 
     2  2.01413871D-01, 3.15853037D-01, 5.09374390D-01, 7.23162115D-01, MgHp 
     3  9.64422251D-01, 1.22623562D+00, 1.41599617D+00, 1.58483751D+00, MgHp 
     4  1.75623869D+00, 1.97545644D+00, 2.24475064D+00, 2.63363564D+00, MgHp 
     5  2.95142050D+00, 3.15592449D+00, 3.36644812D+00, 3.59461119D+00, MgHp 
     6  3.77901439D+00, 3.89524166D+00,      5*0.0D+00/                 MgHp 
      DATA TQ_AlHp/                                                     071215
     1  0.699999789529, 0.729900033922, 0.776500141108, 0.889699875840, AlHp 
     2  1.033900117277, 1.178899861440, 1.429399831450, 1.644799975341, AlHp 
     3  1.901000056768, 2.175500194169, 2.369399873405, 2.539899954924, AlHp 
     4  2.705199936853, 2.882100156552, 3.096400293374, 3.286099884214, AlHp 
     5  3.450599808106, 3.533099800332, 3.614099756467, 3.799000166742, AlHp 
     6  3.889300339928, 4.000000000000,      5*0.0D+00/                 AlHp 
      DATA  Q_AlHp/                                                     071215
     1  3.30148582D-01, 3.38120506D-01, 3.53230018D-01, 4.03557957D-01, AlHp 
     2  4.91294958D-01, 5.97701629D-01, 8.06951068D-01, 1.00215187D+00, AlHp 
     3  1.24461872D+00, 1.51148926D+00, 1.70267579D+00, 1.87256761D+00, AlHp 
     4  2.04185256D+00, 2.23724373D+00, 2.51041490D+00, 2.79462855D+00, AlHp 
     5  3.07569543D+00, 3.22934793D+00, 3.38819279D+00, 3.77027883D+00, AlHp 
     6  3.95925905D+00, 4.18986102D+00,      5*0.0D+00/                 AlHp 
      DATA TQ_SiHp/                                                     071215
     1  0.699999789529, 0.732399969475, 0.783600127849, 0.905399967786, SiHp 
     2  1.057300012856, 1.197599858818, 1.461799981939, 1.736099874204, SiHp 
     3  1.989999772746, 2.256900162824, 2.456399945861, 2.627000034631, SiHp 
     4  2.788399937516, 2.949099790724, 3.180700238300, 3.433300059644, SiHp 
     5  3.528899809750, 3.627600068656, 3.733999640789, 3.832699985538, SiHp 
     6  3.930100316125, 4.000000000000,      5*0.0D+00/                 SiHp 
      DATA  Q_SiHp/                                                     071215
     1  1.66962513D-02, 2.26550007D-02, 3.49729100D-02, 8.00651985D-02, SiHp 
     2  1.65186708D-01, 2.64053712D-01, 4.80689055D-01, 7.29090666D-01, SiHp 
     3  9.70232804D-01, 1.23003514D+00, 1.42676195D+00, 1.59649956D+00, SiHp 
     4  1.76053867D+00, 1.93379581D+00, 2.21897711D+00, 2.59441144D+00, SiHp 
     5  2.75501045D+00, 2.93213322D+00, 3.13840034D+00, 3.34754602D+00, SiHp 
     6  3.57155930D+00, 3.74004276D+00,      5*0.0D+00/                 SiHp 
      DATA TQ_PHp/                                                      071215
     1  0.699999789529, 0.896500013925, 1.141299905448, 1.265100096992, PHp  
     2  1.387699911551, 1.777800164356, 2.043200042658, 2.334300029097, PHp  
     3  2.524100152273, 2.700999832545, 2.841200193564, 3.059200419395, PHp  
     4  3.213800111600, 3.381400139928, 3.612199710369, 3.709900033111, PHp  
     5  3.801400218396, 3.923100170033, 3.969800266362, 4.000000000000, PHp  
     6       7*0.0D+00/                                                 PHp  
      DATA  Q_PHp/                                                      071215
     1  1.07922907D+00, 1.08076685D+00, 1.10100029D+00, 1.13022430D+00, PHp  
     2  1.17550055D+00, 1.41222621D+00, 1.62659246D+00, 1.88780068D+00, PHp  
     3  2.06672209D+00, 2.23786860D+00, 2.37805496D+00, 2.61498410D+00, PHp  
     4  2.80636857D+00, 3.04057699D+00, 3.40699985D+00, 3.57619374D+00, PHp  
     5  3.74281551D+00, 3.98039448D+00, 4.07768926D+00, 4.14262838D+00, PHp  
     6       7*0.0D+00/                                                 PHp  
      DATA TQ_SHp/                                                      071215
     1  0.699999789529, 0.725299923758, 0.764800119742, 0.870399882535, SHp  
     2  0.976799932696, 1.092600099595, 1.264100120482, 1.469599790058, SHp  
     3  1.746899923086, 2.040399978561, 2.522100300793, 2.699699801014, SHp  
     4  2.870599905696, 3.018100408788, 3.176000203319, 3.339400149246, SHp  
     5  3.488399721202, 3.667200002074, 3.830099928907, 3.930900258466, SHp  
     6  3.973000055282, 4.000000000000,      5*0.0D+00/                 SHp  
      DATA  Q_SHp/                                                      071215
     1  3.07856728D-01, 3.10191052D-01, 3.15036560D-01, 3.37258521D-01, SHp  
     2  3.76133334D-01, 4.36944202D-01, 5.53210850D-01, 7.17681609D-01, SHp  
     3  9.63819602D-01, 1.24055863D+00, 1.71139502D+00, 1.88824306D+00, SHp  
     4  2.06326091D+00, 2.22462887D+00, 2.41577974D+00, 2.63861442D+00, SHp  
     5  2.86446577D+00, 3.16141632D+00, 3.45404946D+00, 3.64552218D+00, SHp  
     6  3.72790257D+00, 3.78146516D+00,      5*0.0D+00/                 SHp  
      DATA TQ_HClp/                                                     071215
     1  0.699999789529, 0.902600032329, 1.024800093064, 1.204699900803, HClp 
     2  1.319600154428, 1.439600033502, 1.856600058492, 2.115099781611, HClp 
     3  2.383200195354, 2.559300396371, 2.724999919775, 2.908099788452, HClp 
     4  3.140500312349, 3.288799946280, 3.442599639555, 3.600900315566, HClp 
     5  3.783399822449, 3.845300263772, 3.906399763377, 3.963800131897, HClp 
     6  4.000000000000,      6*0.0D+00/                                 HClp 
      DATA  Q_HClp/                                                     071215
     1  1.07919075D+00, 1.07980883D+00, 1.08270743D+00, 1.10031964D+00, HClp 
     2  1.12641899D+00, 1.16889018D+00, 1.42078845D+00, 1.63050790D+00, HClp 
     3  1.87074694D+00, 2.03616275D+00, 2.19578032D+00, 2.37889135D+00, HClp 
     4  2.63697713D+00, 2.82750323D+00, 3.05130061D+00, 3.31111418D+00, HClp 
     5  3.65277315D+00, 3.78136729D+00, 3.91536352D+00, 4.04715801D+00, HClp 
     6  4.13264512D+00,      6*0.0D+00/                                 HClp 
      DATA TQ_ZnHp/                                                     071215
     1  0.699999789529, 0.731599991755, 0.781400188115, 0.900800073821, ZnHp 
     2  1.050800158276, 1.192899973416, 1.429799820859, 1.671300037812, ZnHp 
     3  1.933899908611, 2.213800128833, 2.407599791865, 2.573399777188, ZnHp 
     4  2.741999816949, 2.936099864109, 3.130400129654, 3.342000209198, ZnHp 
     5  3.462400087335, 3.563800132133, 3.648799706535, 3.731499589142, ZnHp 
     6  3.910099846792, 3.964000136379, 4.000000000000,      4*0.0D+00/ ZnHp 
      DATA  Q_ZnHp/                                                     071215
     1  1.94232576D-02, 2.58654465D-02, 3.88667632D-02, 8.52972096D-02, ZnHp 
     2  1.71196401D-01, 2.72427126D-01, 4.66572408D-01, 6.83551740D-01, ZnHp 
     3  9.31501644D-01, 1.20337806D+00, 1.39434901D+00, 1.55929306D+00, ZnHp 
     4  1.73096239D+00, 1.94378688D+00, 2.18842866D+00, 2.50255010D+00, ZnHp 
     5  2.70578117D+00, 2.89289018D+00, 3.06320399D+00, 3.24169956D+00, ZnHp 
     6  3.65628334D+00, 3.78214677D+00, 3.86481462D+00,      4*0.0D+00/ ZnHp 
      DATA TQ_HBrp/                                                     071215
     1  0.699999789529, 0.891599904253, 1.112200083739, 1.261800174509, HBrp 
     2  1.419100059636, 1.653400007518, 1.951199833382, 2.249399989435, HBrp 
     3  2.583200000418, 2.708900028743, 2.824099799989, 3.112599712923, HBrp 
     4  3.281499778472, 3.463200104644, 3.612899727353, 3.794700076591, HBrp 
     5  3.910599859639, 3.965400167754, 4.000000000000,      8*0.0D+00/ HBrp 
      DATA  Q_HBrp/                                                     071215
     1  1.07925966D+00, 1.08121347D+00, 1.09994908D+00, 1.13636340D+00, HBrp 
     2  1.20062924D+00, 1.33956129D+00, 1.56800220D+00, 1.83053066D+00, HBrp 
     3  2.14529948D+00, 2.26765646D+00, 2.38261972D+00, 2.70040405D+00, HBrp 
     4  2.92099463D+00, 3.19409835D+00, 3.44826059D+00, 3.79837158D+00, HBrp 
     5  4.05192534D+00, 4.18005994D+00, 4.26300920D+00,      8*0.0D+00/ HBrp 
      DATA TQ_CdHp/                                                     071215
     1  0.699999789529, 0.727999988419, 0.771400016729, 0.878300073861, CdHp 
     2  1.014100115482, 1.154599979609, 1.408099896796, 1.631000086260, CdHp 
     3  1.887799922134, 2.162299904275, 2.365299774072, 2.539999957356, CdHp 
     4  2.708900028743, 2.939499625507, 3.122299928551, 3.323999784675, CdHp 
     5  3.433200066721, 3.555000300888, 3.638900305514, 3.731499589142, CdHp 
     6  3.921800141150, 3.968300232745, 4.000000000000,      4*0.0D+00/ CdHp 
      DATA  Q_CdHp/                                                     071215
     1  4.03084955D-02, 4.94520087D-02, 6.60014358D-02, 1.18480943D-01, CdHp 
     2  2.04691619D-01, 3.09550854D-01, 5.23044018D-01, 7.26109383D-01, CdHp 
     3  9.69936133D-01, 1.23717111D+00, 1.43739545D+00, 1.61125178D+00, CdHp 
     4  1.78320061D+00, 2.03889223D+00, 2.27429550D+00, 2.57724193D+00, CdHp 
     5  2.76169375D+00, 2.98697439D+00, 3.15627653D+00, 3.35738657D+00, CdHp 
     6  3.80277570D+00, 3.91270315D+00, 3.98678068D+00,      4*0.0D+00/ CdHp 
      DATA TQ_HgHp/                                                     071215
     1  0.699999789529, 0.730800014036, 0.779100204517, 0.895599993781, HgHp 
     2  1.041499987577, 1.197799853941, 1.435799947235, 1.679199855999, HgHp 
     3  1.937599987897, 2.210100056939, 2.412599913470, 2.584400028533, HgHp 
     4  2.758300191537, 2.954999925729, 3.166899999219, 3.427000233125, HgHp 
     5  3.591500192879, 3.684299997768, 3.784599850741, 3.862799728450, HgHp 
     6  3.924900210024, 3.971300177434, 4.000000000000,      4*0.0D+00/ HgHp 
      DATA  Q_HgHp/                                                     071215
     1  2.99999919D-02, 3.83943062D-02, 5.44477636D-02, 1.07553292D-01, HgHp 
     2  1.97909570D-01, 3.14310186D-01, 5.14662514D-01, 7.36387302D-01, HgHp 
     3  9.82080929D-01, 1.24755221D+00, 1.44733113D+00, 1.61822915D+00, HgHp 
     4  1.79485448D+00, 2.00939430D+00, 2.27613868D+00, 2.67100741D+00, HgHp 
     5  2.96293573D+00, 3.14458679D+00, 3.35662231D+00, 3.53230110D+00, HgHp 
     6  3.67605839D+00, 3.78439908D+00, 3.85129440D+00,      4*0.0D+00/ HgHp 
      DATA TQ_CHm/                                                      071215
     1  0.699999789529, 0.846799984302, 1.017500034254, 1.153000020522, CHm  
     2  1.300200166544, 1.451899850175, 1.612000090111, 1.797700130973, CHm  
     3  1.988499807752, 2.284599837226, 2.554100276283, 2.726099840218, CHm  
     4  2.891700223334, 3.208599990934, 3.359999634096, 3.516900384344, CHm  
     5  3.660899873471, 3.780899763508, 3.871699931419, 3.954999937558, CHm  
     6  4.000000000000,      6*0.0D+00/                                 CHm  
      DATA  Q_CHm/                                                      071215
     1  4.77498652D-01, 4.81009830D-01, 5.02208791D-01, 5.45648026D-01, CHm  
     2  6.22905282D-01, 7.26692093D-01, 8.52865689D-01, 1.01269946D+00, CHm  
     3  1.18657089D+00, 1.46791589D+00, 1.73101709D+00, 1.90135454D+00, CHm  
     4  2.06894319D+00, 2.42755686D+00, 2.63433293D+00, 2.88303388D+00, CHm  
     5  3.14584612D+00, 3.39127339D+00, 3.59177067D+00, 3.78331676D+00, CHm  
     6  3.88799314D+00,      6*0.0D+00/                                 CHm  
      DATA TQ_OHm/                                                      071215
     1  0.699999789529, 0.796200092479, 0.898300054213, 1.141499909856, OHm  
     2  1.232299954887, 1.323900073811, 1.507900177925, 1.705699948431, OHm  
     3  1.963799956518, 2.221900115926, 2.605499972507, 2.762500293781, OHm  
     4  2.918400036547, 3.238799745438, 3.360199639399, 3.487899710114, OHm  
     5  3.643500080941, 3.735499671778, 3.817099894333, 3.923800185585, OHm  
     6  3.970300249288, 3.987799727422, 4.000000000000,      4*0.0D+00/ OHm  
      DATA  Q_OHm/                                                      071215
     1  3.09589940D-05, 2.56710642D-04, 1.53139568D-03, 2.68304330D-02, OHm  
     2  5.39030298D-02, 9.38400042D-02, 2.06062425D-01, 3.57178686D-01, OHm  
     3  5.80212262D-01, 8.19441927D-01, 1.19015314D+00, 1.34490598D+00, OHm  
     4  1.50051911D+00, 1.84502639D+00, 1.99553984D+00, 2.17122434D+00, OHm  
     5  2.41220477D+00, 2.57058328D+00, 2.72376181D+00, 2.94769018D+00, OHm  
     6  3.05482823D+00, 3.09662622D+00, 3.12622536D+00,      4*0.0D+00/ OHm  
      DATA TQ_SiHm/                                                     071215
     1  0.699999789529, 0.732199975045, 0.782900147024, 0.904099997753, SiHm 
     2  1.057300012856, 1.201299828098, 1.439200024421, 1.679599846793, SiHm 
     3  1.944999921213, 2.225299877252, 2.436199850646, 2.615099761380, SiHm 
     4  2.789199957785, 2.967700217962, 3.107399845236, 3.245499899759, SiHm 
     5  3.389200337447, 3.533699813986, 3.664199940834, 3.769500453854, SiHm 
     6  3.856599904184, 3.926400243350, 3.971700148692, 4.000000000000, SiHm 
     7       3*0.0D+00/                                                 SiHm 
      DATA  Q_SiHm/                                                     071215
     1  4.95406092D-01, 5.01713282D-01, 5.14559359D-01, 5.60964793D-01, SiHm 
     2  6.48410862D-01, 7.51102598D-01, 9.46318456D-01, 1.16243734D+00, SiHm 
     3  1.41310476D+00, 1.68539008D+00, 1.89322994D+00, 2.07128837D+00, SiHm 
     4  2.24910848D+00, 2.44566165D+00, 2.61854216D+00, 2.81272561D+00, SiHm 
     5  3.04411115D+00, 3.30864736D+00, 3.57401779D+00, 3.80528862D+00, SiHm 
     6  4.00645992D+00, 4.17178291D+00, 4.27968089D+00, 4.34689134D+00, SiHm 
     7       3*0.0D+00/                                                 SiHm 
      DATA TQ_HSm/                                                      071215
     1  0.699999789529, 0.734099922129, 0.788899982662, 0.922000037400, HSm  
     2  1.039899956967, 1.204299892249, 1.383300018637, 1.594400046398, HSm  
     3  1.787300039346, 1.996599922734, 2.265400352954, 2.563700149198, HSm  
     4  2.761200262379, 2.954199907298, 3.095700276508, 3.236299684720, HSm  
     5  3.562500225563, 3.820999709145, 3.923900187807, 3.970500234917, HSm  
     6  4.000000000000,      6*0.0D+00/                                 HSm  
      DATA  Q_HSm/                                                      071215
     1  6.11870692D-03, 9.13589177D-03, 1.62954007D-02, 4.97006511D-02, HSm  
     2  1.01297902D-01, 2.02542577D-01, 3.37796443D-01, 5.16870983D-01, HSm  
     3  6.91744598D-01, 8.88894895D-01, 1.14888005D+00, 1.44255163D+00, HSm  
     4  1.63950727D+00, 1.83888713D+00, 1.99707948D+00, 2.17125249D+00, HSm  
     5  2.65924743D+00, 3.14362135D+00, 3.36710767D+00, 3.47347569D+00, HSm  
     6  3.54193175D+00,      6*0.0D+00/                                 HSm  
      DATA TQ_CN/                                                       071215
     1  0.699999789529, 0.730800014036, 0.778500189884, 0.901200064601, CN   
     2  1.229099996184, 1.654699972927, 2.085700050913, 2.343000231691, CN   
     3  2.565899992619, 2.757600174302, 2.913299911859, 3.051200211750, CN   
     4  3.325499817198, 3.496199901386, 3.613199734631, 3.855499981738, CN   
     5  3.944099695165, 3.978299674457, 4.000000000000,      8*0.0D+00/ CN   
      DATA  Q_CN/                                                       071215
     1  6.46594758D-01, 6.71774343D-01, 7.11557067D-01, 8.17604575D-01, CN   
     2  1.11881632D+00, 1.52981618D+00, 1.95539416D+00, 2.21135443D+00, CN   
     3  2.43384470D+00, 2.62782462D+00, 2.79330935D+00, 2.95263238D+00, CN   
     4  3.32222357D+00, 3.59489726D+00, 3.80483067D+00, 4.30513988D+00, CN   
     5  4.51027706D+00, 4.59300642D+00, 4.64665528D+00,      8*0.0D+00/ CN   
      DATA TQ_CO/                                                       071215
     1  0.699999789529, 0.731200002896, 0.779600216711, 0.904199995448, CO   
     2  1.237399820311, 1.633700019021, 2.088500108379, 2.354900024839, CO   
     3  2.586600080076, 2.783199805769, 3.014900335448, 3.148799693004, CO   
     4  3.283899833642, 3.558500373931, 3.722200098645, 3.866299808421, CO   
     5  3.948399792254, 3.980299559040, 4.000000000000,      8*0.0D+00/ CO   
      DATA  Q_CO/                                                       071215
     1  3.39713497D-01, 3.65121675D-01, 4.05361205D-01, 5.12837927D-01, CO   
     2  8.18794237D-01, 1.20113568D+00, 1.64995527D+00, 1.91496020D+00, CO   
     3  2.14623676D+00, 2.34526807D+00, 2.59730264D+00, 2.76086831D+00, CO   
     4  2.94322845D+00, 3.36759419D+00, 3.65050185D+00, 3.91541984D+00, CO   
     5  4.07392170D+00, 4.13775878D+00, 4.17808165D+00,      8*0.0D+00/ CO   
      DATA TQ_CF/                                                       071215
     1  0.699999789529, 0.729800031527, 0.776000128914, 0.894099960208, CF   
     2  1.206099930740, 1.604700008252, 2.078899906352, 2.283299810877, CF   
     3  2.465400148013, 2.840800184886, 3.015300344616, 3.179300274060, CF   
     4  3.408799805524, 3.599700374070, 3.729999558153, 3.850700320154, CF   
     5  3.940699618398, 3.976899775052, 4.000000000000,      8*0.0D+00/ CF   
      DATA  Q_CF/                                                       071215
     1  1.21732623D+00, 1.23328082D+00, 1.25958248D+00, 1.33462890D+00, CF   
     2  1.57371316D+00, 1.92949283D+00, 2.38465127D+00, 2.58556536D+00, CF   
     3  2.76652935D+00, 3.17135319D+00, 3.39564056D+00, 3.63407823D+00, CF   
     4  4.00983685D+00, 4.35476869D+00, 4.60627029D+00, 4.85325426D+00, CF   
     5  5.04885703D+00, 5.13087051D+00, 5.18429996D+00,      8*0.0D+00/ CF   
      DATA TQ_SiC/                                                      071215
     1  0.699999789529, 0.727999988419, 0.771200011851, 0.882100063248, SiC  
     2  1.174799956587, 1.508400189682, 1.942799973183, 2.142400153184, SiC  
     3  2.332199979317, 2.478399642236, 2.618199832741, 2.832199991938, SiC  
     4  3.112999721479, 3.297700139571, 3.441499616657, 3.621899941872, SiC  
     5  3.807800358689, 3.923900187807, 4.000000000000,      8*0.0D+00/ SiC  
      DATA  Q_SiC/                                                      071215
     1  1.15915930D+00, 1.18113479D+00, 1.21573765D+00, 1.30799018D+00, SiC  
     2  1.56889251D+00, 1.88457729D+00, 2.30911991D+00, 2.50661824D+00, SiC  
     3  2.69571758D+00, 2.84461024D+00, 2.99462267D+00, 3.25138721D+00, SiC  
     4  3.65277792D+00, 3.95941454D+00, 4.22558788D+00, 4.59443211D+00, SiC  
     5  5.00236219D+00, 5.26247945D+00, 5.43285658D+00,      8*0.0D+00/ SiC  
      DATA TQ_CP/                                                       071215
     1  0.699999789529, 0.783200138806, 0.911799902048, 1.277200154905, CP   
     2  1.612900067607, 1.970299801424, 2.196899727746, 2.392900128340, CP   
     3  2.564300106495, 2.737799726818, 2.859999642519, 2.968500236356, CP   
     4  3.183900006231, 3.369499885991, 3.520900391313, 3.632600175638, CP   
     5  3.732699613933, 3.903699702742, 3.961200073629, 4.000000000000, CP   
     6       7*0.0D+00/                                                 CP   
      DATA  Q_CP/                                                       071215
     1  9.75539232D-01, 1.05288716D+00, 1.17440842D+00, 1.52818797D+00, CP   
     2  1.85918787D+00, 2.21438396D+00, 2.44034898D+00, 2.63642227D+00, CP   
     3  2.81098672D+00, 2.99874469D+00, 3.14357974D+00, 3.28356655D+00, CP   
     4  3.59688458D+00, 3.91282762D+00, 4.21553339D+00, 4.46764502D+00, CP   
     5  4.71060081D+00, 5.15015330D+00, 5.30326182D+00, 5.40808880D+00, CP   
     6       7*0.0D+00/                                                 CP   
      DATA TQ_CS/                                                       071215
     1  0.699999789529, 0.761200214697, 0.857100069962, 1.100499921844, CS   
     2  1.468999804818, 1.954699921706, 2.200499803104, 2.410199855268, CS   
     3  2.570099703058, 2.738899750604, 2.924600175266, 3.162799899390, CS   
     4  3.448299758209, 3.610699673977, 3.688199726142, 3.765000345061, CS   
     5  3.839200127116, 3.908299806045, 3.963600127415, 3.985899684766, CS   
     6  4.000000000000,      6*0.0D+00/                                 CS   
      DATA  Q_CS/                                                       071215
     1  6.63987746D-01, 7.20654575D-01, 8.10635559D-01, 1.04380422D+00, CS   
     2  1.40456786D+00, 1.88644857D+00, 2.13152462D+00, 2.34129203D+00, CS   
     3  2.50397100D+00, 2.68548517D+00, 2.90845224D+00, 3.24195223D+00, CS   
     4  3.70873669D+00, 4.00090836D+00, 4.14704257D+00, 4.29864441D+00, CS   
     5  4.45746707D+00, 4.62617810D+00, 4.78403411D+00, 4.85481380D+00, CS   
     6  4.90184717D+00,      6*0.0D+00/                                 CS   
      DATA TQ_CCl/                                                      071215
     1  0.699999789529, 0.735599880353, 0.791199979616, 0.933799918537, CCl  
     2  1.286400045861, 1.522700038242, 1.741199805231, 1.910999855627, CCl  
     3  2.080999954452, 2.222500073807, 2.361999694122, 2.494499872091, CCl  
     4  2.625399997955, 2.906099743477, 3.046600112905, 3.194199662280, CCl  
     5  3.501700029706, 3.640000328190, 3.817099894333, 3.924600203359, CCl  
     6  3.970600227731, 4.000000000000,      5*0.0D+00/                 CCl  
      DATA  Q_CCl/                                                      071215
     1  1.11563638D+00, 1.14270832D+00, 1.18628134D+00, 1.30421836D+00, CCl  
     2  1.62284740D+00, 1.85034245D+00, 2.06801284D+00, 2.24248953D+00, CCl  
     3  2.42461410D+00, 2.58524927D+00, 2.75129673D+00, 2.91115255D+00, CCl  
     4  3.06740357D+00, 3.40798659D+00, 3.59478806D+00, 3.81063981D+00, CCl  
     5  4.32773422D+00, 4.58609399D+00, 4.93713316D+00, 5.16073181D+00, CCl  
     6  5.25850881D+00, 5.32152038D+00,      5*0.0D+00/                 CCl  
      DATA TQ_CSe/                                                      071215
     1  0.699999789529, 0.760100243710, 0.854300005732, 1.098799941011, CSe  
     2  1.482899947109, 1.890399871855, 2.123499954682, 2.337600107322, CSe  
     3  2.483099605583, 2.631900148415, 2.840800184886, 2.994599873621, CSe  
     4  3.139900347367, 3.396999846599, 3.556000321758, 3.700999830477, CSe  
     5  3.792800036757, 3.884700230742, 3.954299922273, 3.982699612923, CSe  
     6  4.000000000000,      6*0.0D+00/                                 CSe  
      DATA  Q_CSe/                                                      071215
     1  8.07734963D-01, 8.64724736D-01, 9.54845683D-01, 1.19212062D+00, CSe  
     2  1.57064470D+00, 1.97581337D+00, 2.20835370D+00, 2.42270425D+00, CSe  
     3  2.57109687D+00, 2.73072014D+00, 2.98121439D+00, 3.19156282D+00, CSe  
     4  3.41098497D+00, 3.84203715D+00, 4.13112171D+00, 4.40900705D+00, CSe  
     5  4.59511183D+00, 4.79488696D+00, 4.95907041D+00, 5.02986203D+00, CSe  
     6  5.07412392D+00,      6*0.0D+00/                                 CSe  
      DATA TQ_CBr/                                                      071215
     1  0.699999789529, 0.729800031527, 0.775900126475, 0.893599949017, CBr  
     2  1.211299986309, 1.552600118261, 1.927499887696, 2.162599910217, CBr  
     3  2.354500054578, 2.489899779856, 2.640400311355, 2.790699993906, CBr  
     4  2.924900181761, 3.056400346719, 3.328899890918, 3.607099866868, CBr  
     5  3.836600070485, 3.936599847644, 4.000000000000,      8*0.0D+00/ CBr  
      DATA  Q_CBr/                                                      071215
     1  1.53248992D+00, 1.55704915D+00, 1.59569533D+00, 1.69750782D+00, CBr  
     2  1.98847387D+00, 2.31616161D+00, 2.68444347D+00, 2.91753724D+00, CBr  
     3  3.10901583D+00, 3.24684237D+00, 3.40804148D+00, 3.58424804D+00, CBr  
     4  3.75850645D+00, 3.94571262D+00, 4.37989215D+00, 4.87079559D+00, CBr  
     5  5.29853923D+00, 5.48914851D+00, 5.61099391D+00,      8*0.0D+00/ CBr  
      DATA TQ_RhC/                                                      071215
     1  0.699999789529, 0.762700175132, 0.861000110122, 1.115000011059, RhC  
     2  1.502700055655, 1.920100073158, 2.136700245006, 2.344600267680, RhC  
     3  2.503500074823, 2.664399958139, 2.780899747496, 2.901499640034, RhC  
     4  3.125500008316, 3.299600180459, 3.459600025738, 3.616999826826, RhC  
     5  3.736099684173, 3.844000234527, 3.939499638630, 3.976399810979, RhC  
     6  4.000000000000,      6*0.0D+00/                                 RhC  
      DATA  Q_RhC/                                                      071215
     1  1.08948844D+00, 1.14879531D+00, 1.24267892D+00, 1.48900699D+00, RhC  
     2  1.87104436D+00, 2.28607998D+00, 2.50216715D+00, 2.71032473D+00, RhC  
     3  2.87269653D+00, 3.04712124D+00, 3.18467271D+00, 3.33997089D+00, RhC  
     4  3.66563015D+00, 3.95019306D+00, 4.23727376D+00, 4.55421557D+00, RhC  
     5  4.82541448D+00, 5.09737662D+00, 5.35791847D+00, 5.46291951D+00, RhC  
     6  5.53114337D+00,      6*0.0D+00/                                 RhC  
      DATA TQ_IrC/                                                      071215
     1  0.699999789529, 0.728600002789, 0.772900053311, 0.886399957214, IrC  
     2  1.181499867144, 1.514600107301, 1.986999842759, 2.196999729749, IrC  
     3  2.389100322164, 2.625099991079, 2.842600223937, 3.063700164272, IrC  
     4  3.349300373396, 3.510800241890, 3.671700099694, 3.772800263319, IrC  
     5  3.863499744444, 3.947399769675, 3.979799566676, 4.000000000000, IrC  
     6       7*0.0D+00/                                                 IrC  
      DATA  Q_IrC/                                                      071215
     1  1.56981119D+00, 1.58951803D+00, 1.62104820D+00, 1.70680869D+00, IrC  
     2  1.95496774D+00, 2.26164636D+00, 2.71873120D+00, 2.92581488D+00, IrC  
     3  3.11716288D+00, 3.36344590D+00, 3.62122752D+00, 3.92837164D+00, IrC  
     4  4.39011475D+00, 4.67879802D+00, 4.99180251D+00, 5.21095266D+00, IrC  
     5  5.42905634D+00, 5.65099189D+00, 5.74166267D+00, 5.79946239D+00, IrC  
     6       7*0.0D+00/                                                 IrC  
      DATA TQ_PtC/                                                      071215
     1  0.699999789529, 0.763000167219, 0.861700091668, 1.117699940975, PtC  
     2  1.506300140304, 1.922300018021, 2.138300281781, 2.346500310416, PtC  
     3  2.506900165294, 2.666800011179, 2.787199907113, 2.910399840958, PtC  
     4  3.340900184456, 3.502200041363, 3.662399904091, 3.761300255609, PtC  
     5  3.855799960587, 3.945199720002, 3.978799638530, 4.000000000000, PtC  
     6       7*0.0D+00/                                                 PtC  
      DATA  Q_PtC/                                                      071215
     1  8.40812489D-01, 9.00816677D-01, 9.95622250D-01, 1.24484633D+00, PtC  
     2  1.62848848D+00, 2.04243106D+00, 2.25798642D+00, 2.46648641D+00, PtC  
     3  2.63048046D+00, 2.80401543D+00, 2.94653544D+00, 3.10597104D+00, PtC  
     4  3.77375553D+00, 4.06467610D+00, 4.38964140D+00, 4.62166921D+00, PtC  
     5  4.87260087D+00, 5.13598061D+00, 5.24057153D+00, 5.30785947D+00, PtC  
     6       7*0.0D+00/                                                 PtC  
      DATA TQ_CNp/                                                      071215
     1  0.699999789529, 0.731399997326, 0.780200220987, 0.905299970091, CNp  
     2  1.243099835125, 1.650100095325, 2.094100242272, 2.358099786927, CNp  
     3  2.578299887258, 2.746099904193, 2.931900158853, 3.092200192175, CNp  
     4  3.245099890655, 3.444599681188, 3.589100138648, 3.764300328137, CNp  
     5  3.862699726165, 3.945499726775, 3.978999624159, 4.000000000000, CNp  
     6       7*0.0D+00/                                                 CNp  
      DATA  Q_CNp/                                                      071215
     1  3.46327491D-01, 3.72014803D-01, 4.12755620D-01, 5.21030004D-01, CNp  
     2  8.31933974D-01, 1.22520479D+00, 1.66365146D+00, 1.92633293D+00, CNp  
     3  2.14624463D+00, 2.31622011D+00, 2.51504615D+00, 2.70553893D+00, CNp  
     4  2.91077434D+00, 3.22622975D+00, 3.50261214D+00, 3.89751145D+00, CNp  
     5  4.14245760D+00, 4.35891183D+00, 4.44917373D+00, 4.50659458D+00, CNp  
     6       7*0.0D+00/                                                 CNp  
      DATA TQ_COp/                                                      071215
     1  0.699999789529, 0.730200030747, 0.776900150863, 0.897400034069, COp  
     2  1.216399877143, 1.648300057890, 2.097700332805, 2.369899885518, COp  
     3  2.603700103918, 2.792200027875, 3.047900138790, 3.252500056740, COp  
     4  3.485799663541, 3.585100044933, 3.687999740071, 3.874900003810, COp  
     5  3.950999850215, 3.981299581491, 4.000000000000,      8*0.0D+00/ COp  
      DATA  Q_COp/                                                      071215
     1  6.32725224D-01, 6.57174262D-01, 6.95772290D-01, 7.99123721D-01, COp  
     2  1.09059420D+00, 1.50670479D+00, 1.95023004D+00, 2.22101291D+00, COp  
     3  2.45442543D+00, 2.64537261D+00, 2.92571007D+00, 3.18730928D+00, COp  
     4  3.53534033D+00, 3.69899208D+00, 3.87913813D+00, 4.24723751D+00, COp  
     5  4.41981203D+00, 4.49313629D+00, 4.53978032D+00,      8*0.0D+00/ COp  
      DATA TQ_CNm/                                                      071215
     1  0.699999789529, 0.731100005681, 0.779400211833, 0.903500011583, CNm  
     2  1.236799836143, 1.625800018306, 2.076799853707, 2.344400263181, CNm  
     3  2.571699739000, 2.765900375909, 3.000800017852, 3.127300053185, CNm  
     4  3.257600168051, 3.541399990842, 3.689799614705, 3.843300218779, CNm  
     5  3.938399717911, 3.975899846906, 4.000000000000,      8*0.0D+00/ CNm  
      DATA  Q_CNm/                                                      071215
     1  3.45564763D-01, 3.70991632D-01, 4.11286285D-01, 5.18606520D-01, CNm  
     2  8.25073437D-01, 1.20048586D+00, 1.64552835D+00, 1.91171951D+00, CNm  
     3  2.13862231D+00, 2.33535912D+00, 2.59161577D+00, 2.74650412D+00, CNm  
     4  2.92203731D+00, 3.36057547D+00, 3.61686606D+00, 3.89847178D+00, CNm  
     5  4.08096807D+00, 4.15470594D+00, 4.20266678D+00,      8*0.0D+00/ CNm  
      DATA TQ_CSm/                                                      071215
     1  0.699999789529, 0.730300027961, 0.777300160619, 0.897600038545, CSm  
     2  1.215599894267, 1.571200142665, 2.015600361429, 2.238699758309, CSm  
     3  2.439499620329, 2.850400354805, 3.006400142815, 3.163599918869, CSm  
     4  3.517800405362, 3.678800268691, 3.834300020388, 3.934699984585, CSm  
     5  3.974399954687, 4.000000000000,      9*0.0D+00/                 CSm  
      DATA  Q_CSm/                                                      071215
     1  1.36367057D+00, 1.38534951D+00, 1.42004294D+00, 1.51403048D+00, CSm  
     2  1.78850834D+00, 2.12125597D+00, 2.55409677D+00, 2.77466505D+00, CSm  
     3  2.97469780D+00, 3.41853907D+00, 3.61931343D+00, 3.84557677D+00, CSm  
     4  4.43381785D+00, 4.73005393D+00, 5.03067330D+00, 5.23276305D+00, CSm  
     5  5.31469412D+00, 5.36821889D+00,      9*0.0D+00/                 CSm  
      DATA TQ_BN/                                                       071215
     1  0.699999789529, 0.730800014036, 0.778600192323, 0.900900071516, BN   
     2  1.217399855738, 1.431699854157, 1.643499944679, 1.884800002625, BN   
     3  2.134400192142, 2.334100024356, 2.517400402410, 2.657899809434, BN   
     4  2.869899890457, 3.039199962594, 3.212100071654, 3.421000112953, BN   
     5  3.594700263588, 3.724499939269, 3.853600115694, 3.941199629687, BN   
     6  3.977199753496, 4.000000000000,      5*0.0D+00/                 BN   
      DATA  Q_BN/                                                       071215
     1  1.35923914D+00, 1.37366317D+00, 1.39789127D+00, 1.46918588D+00, BN   
     2  1.70123286D+00, 1.88395301D+00, 2.07644285D+00, 2.30459161D+00, BN   
     3  2.54649233D+00, 2.74264250D+00, 2.92458748D+00, 3.06748475D+00, BN   
     4  3.30009357D+00, 3.51174105D+00, 3.75725891D+00, 4.09238601D+00, BN   
     5  4.39955614D+00, 4.64507127D+00, 4.90522961D+00, 5.09405361D+00, BN   
     6  5.17526338D+00, 5.22788272D+00,      5*0.0D+00/                 BN   
      DATA TQ_NO/                                                       071215
     1  0.699999789529, 0.734399913773, 0.787800012795, 0.928399883083, NO   
     2  1.209299999167, 1.363400050199, 1.509100206141, 1.701099840201, NO   
     3  1.874099985452, 2.122699937903, 2.338300123915, 2.484699646588, NO   
     4  2.625399997955, 2.910499843402, 3.045400089012, 3.191099589629, NO   
     5  3.494099852375, 3.705899942040, 3.851100291953, 3.944599706455, NO   
     6  3.978399667271, 4.000000000000,      5*0.0D+00/                 NO   
      DATA  Q_NO/                                                       071215
     1  8.77169171D-01, 8.92955999D-01, 9.19780683D-01, 1.00256381D+00, NO   
     2  1.20856809D+00, 1.33765924D+00, 1.46812008D+00, 1.65638077D+00, NO   
     3  1.84755315D+00, 2.15091735D+00, 2.41974690D+00, 2.59737428D+00, NO   
     4  2.76329064D+00, 3.09867170D+00, 3.26887289D+00, 3.46859007D+00, NO   
     5  3.94302839D+00, 4.31766362D+00, 4.59328432D+00, 4.78167222D+00, NO   
     6  4.85307424D+00, 4.89993748D+00,      5*0.0D+00/                 NO   
      DATA TQ_NF/                                                       071215
     1  0.699999789529, 0.775700121598, 0.893099937826, 1.214999907110, NF   
     2  1.557100008637, 1.934299917183, 2.167800013206, 2.375800026349, NF   
     3  2.522300285941, 2.672500138507, 2.872199939655, 3.068899776537, NF   
     4  3.300200193469, 3.462300085172, 3.627000055310, 3.840300151289, NF   
     5  3.936499854852, 3.975099904389, 4.000000000000,      8*0.0D+00/ NF   
      DATA  Q_NF/                                                       071215
     1  9.91023640D-01, 1.05856838D+00, 1.16587690D+00, 1.47089204D+00, NF   
     2  1.80470063D+00, 2.17798808D+00, 2.41042604D+00, 2.61846724D+00, NF   
     3  2.76783583D+00, 2.92911709D+00, 3.16849909D+00, 3.44141799D+00, NF   
     4  3.80991568D+00, 4.09568500D+00, 4.40976849D+00, 4.86347621D+00, NF   
     5  5.08989788D+00, 5.18447265D+00, 5.24644048D+00,      8*0.0D+00/ NF   
      DATA TQ_AlN/                                                      071215
     1  0.699999789529, 0.726599954891, 0.767600045889, 0.872699938238, AlN  
     2  1.150800076776, 1.487899831130, 1.865200005205, 2.067399877035, AlN  
     3  2.248299966078, 2.485899677342, 2.710100058315, 2.968700240954, AlN  
     4  3.263200297270, 3.434999939322, 3.610899678829, 3.729899565083, AlN  
     5  3.846300286269, 3.937999746741, 3.975799854091, 4.000000000000, AlN  
     6       7*0.0D+00/                                                 AlN  
      DATA  Q_AlN/                                                      071215
     1  1.65703961D+00, 1.67823213D+00, 1.71149451D+00, 1.79969718D+00, AlN  
     2  2.04830930D+00, 2.36751103D+00, 2.73599163D+00, 2.93593894D+00, AlN  
     3  3.11658202D+00, 3.36567138D+00, 3.63411583D+00, 4.00209993D+00, AlN  
     4  4.49102386D+00, 4.80264653D+00, 5.13774239D+00, 5.37324359D+00, AlN  
     5  5.61151731D+00, 5.80608993D+00, 5.88826320D+00, 5.94145704D+00, AlN  
     6       7*0.0D+00/                                                 AlN  
      DATA TQ_SiN/                                                      071215
     1  0.699999789529, 0.759000220431, 0.851699946089, 1.085000043123, SiN  
     2  1.480699998140, 1.914299936307, 2.150299629321, 2.350300366837, SiN  
     3  2.512300295842, 2.706599971622, 2.860299650032, 3.098300339155, SiN  
     4  3.387700299463, 3.635600237483, 3.749099967163, 3.861299694176, SiN  
     5  3.944399701939, 3.978599652901, 4.000000000000,      8*0.0D+00/ SiN  
      DATA  Q_SiN/                                                      071215
     1  1.01134472D+00, 1.06644821D+00, 1.15401319D+00, 1.37837982D+00, SiN  
     2  1.76661806D+00, 2.19711772D+00, 2.43242763D+00, 2.63243415D+00, SiN  
     3  2.79692913D+00, 3.00635914D+00, 3.19082104D+00, 3.52043918D+00, SiN  
     4  3.99066151D+00, 4.44206675D+00, 4.66245356D+00, 4.89173960D+00, SiN  
     5  5.07219178D+00, 5.14979873D+00, 5.19944755D+00,      8*0.0D+00/ SiN  
      DATA TQ_PN/                                                       071215
     1  0.699999789529, 0.759500233390, 0.852899973616, 1.086700084934, PN   
     2  1.459600017526, 1.944799925937, 2.193599661638, 2.406099755270, PN   
     3  2.584100021504, 2.745599893554, 2.893900059784, 3.059500427182, PN   
     4  3.210800041106, 3.542600019544, 3.704199903334, 3.850500334255, PN   
     5  3.940999625171, 3.977099760681, 4.000000000000,      8*0.0D+00/ PN   
      DATA  Q_PN/                                                       071215
     1  6.80684077D-01, 7.35951479D-01, 8.23787399D-01, 1.04797450D+00, PN   
     2  1.41315329D+00, 1.89459853D+00, 2.14266473D+00, 2.35511258D+00, PN   
     3  2.53586659D+00, 2.70891304D+00, 2.88321013D+00, 3.10149075D+00, PN   
     4  3.32410751D+00, 3.88005754D+00, 4.17804397D+00, 4.46145394D+00, PN   
     5  4.64474993D+00, 4.72028000D+00, 4.76910347D+00,      8*0.0D+00/ PN   
      DATA TQ_NS/                                                       071215
     1  0.699999789529, 0.750800007916, 0.830399928302, 1.031700176058, NS   
     2  1.355900046448, 1.717199861924, 1.933699904326, 2.165599969634, NS   
     3  2.332099976947, 2.515800368977, 2.854900020921, 3.028899823049, NS   
     4  3.399599656760, 3.651899666098, 3.742699824802, 3.828999902553, NS   
     5  3.935899898096, 3.974699933131, 4.000000000000,      8*0.0D+00/ NS   
      DATA  Q_NS/                                                       071215
     1  1.08120207D+00, 1.11862953D+00, 1.18003228D+00, 1.34734223D+00, NS   
     2  1.64076457D+00, 1.98626783D+00, 2.19918210D+00, 2.44062600D+00, NS   
     3  2.63251018D+00, 2.86519192D+00, 3.34667423D+00, 3.62077012D+00, NS   
     4  4.26302147D+00, 4.73597149D+00, 4.91213619D+00, 5.08331277D+00, NS   
     5  5.30406072D+00, 5.38814873D+00, 5.44459157D+00,      8*0.0D+00/ NS   
      DATA TQ_NCl/                                                      071215
     1  0.699999789529, 0.759300228206, 0.852199957559, 1.093200084248, NCl  
     2  1.450299815400, 1.845600013255, 2.058300389466, 2.257400174465, NCl  
     3  2.406899774788, 2.548300143478, 2.662999927199, 2.782299782966, NCl  
     4  2.996999929789, 3.316199968200, 3.481299563743, 3.639500317883, NCl  
     5  3.751300017992, 3.847000302017, 3.938399717911, 3.976099832535, NCl  
     6  4.000000000000,      6*0.0D+00/                                 NCl  
      DATA  Q_NCl/                                                      071215
     1  1.23557823D+00, 1.29140849D+00, 1.37975096D+00, 1.61264667D+00, NCl  
     2  1.96367370D+00, 2.35622400D+00, 2.56830889D+00, 2.76779148D+00, NCl  
     3  2.92100719D+00, 3.07417881D+00, 3.20860490D+00, 3.36095640D+00, NCl  
     4  3.67024932D+00, 4.20578813D+00, 4.51233943D+00, 4.82477074D+00, NCl  
     5  5.05990087D+00, 5.27379777D+00, 5.48990008D+00, 5.58190125D+00, NCl  
     6  5.64084817D+00,      6*0.0D+00/                                 NCl  
      DATA TQ_TiN/                                                      071215
     1  0.699999789529, 0.758500207473, 0.849899909505, 1.089600156259, TiN  
     2  1.445199921057, 1.830299920067, 2.058000382557, 2.265100346401, TiN  
     3  2.404199708916, 2.529099780973, 2.745199885042, 2.899299658341, TiN  
     4  3.036299894090, 3.386100258946, 3.726699786823, 3.895599943884, TiN  
     5  3.959200029267, 4.000000000000,      9*0.0D+00/                 TiN  
      DATA  Q_TiN/                                                      071215
     1  1.07586004D+00, 1.13107029D+00, 1.21814882D+00, 1.45005355D+00, TiN  
     2  1.79978026D+00, 2.18222821D+00, 2.40925590D+00, 2.61672921D+00, TiN  
     3  2.75908499D+00, 2.89305903D+00, 3.15069445D+00, 3.36034966D+00, TiN  
     4  3.56524618D+00, 4.15321361D+00, 4.78653714D+00, 5.11541178D+00, TiN  
     5  5.24152629D+00, 5.32313419D+00,      9*0.0D+00/                 TiN  
      DATA TQ_AsN/                                                      071215
     1  0.699999789529, 0.758500207473, 0.850299913974, 1.084200023447, AsN  
     2  1.471599817683, 1.888199911401, 2.126100009214, 2.331399960353, AsN  
     3  2.489599772167, 2.656999787378, 2.831999987435, 3.029299794767, AsN  
     4  3.205699925945, 3.394700014534, 3.568899765599, 3.720700202585, AsN  
     5  3.858799749077, 3.943699686134, 3.978299674457, 4.000000000000, AsN  
     6       7*0.0D+00/                                                 AsN  
      DATA  Q_AsN/                                                      071215
     1  8.29276751D-01, 8.84900324D-01, 9.72905570D-01, 1.20011290D+00, AsN  
     2  1.58200789D+00, 1.99629800D+00, 2.23364582D+00, 2.43907133D+00, AsN  
     3  2.60002187D+00, 2.77987396D+00, 2.98875820D+00, 3.25841631D+00, AsN  
     4  3.53046473D+00, 3.84982925D+00, 4.16497224D+00, 4.45321068D+00, AsN  
     5  4.72656964D+00, 4.90233613D+00, 4.97659972D+00, 5.02417823D+00, AsN  
     6       7*0.0D+00/                                                 AsN  
      DATA TQ_SeN/                                                      071215
     1  0.699999789529, 0.729100014763, 0.774200085016, 0.888899895567, SeN  
     2  1.203599877281, 1.510500214260, 1.848399945697, 2.082399983185, SeN  
     3  2.272600268978, 2.404599718674, 2.545800087417, 2.703399892149, SeN  
     4  2.880500117430, 3.090300146394, 3.323299769497, 3.449999793597, SeN  
     5  3.566399945272, 3.667900016363, 3.753000058282, 3.833199996429, SeN  
     6  3.900299626388, 3.961200073629, 3.985099666805, 4.000000000000, SeN  
     7       3*0.0D+00/                                                 SeN  
      DATA  Q_SeN/                                                      071215
     1  1.51552356D+00, 1.53926035D+00, 1.57670273D+00, 1.67506893D+00, SeN  
     2  1.96163912D+00, 2.25503950D+00, 2.58584600D+00, 2.81740191D+00, SeN  
     3  3.00719372D+00, 3.14275977D+00, 3.29917595D+00, 3.50001957D+00, SeN  
     4  3.76739260D+00, 4.13292787D+00, 4.57729642D+00, 4.82914451D+00, SeN  
     5  5.06582364D+00, 5.27810473D+00, 5.46377827D+00, 5.64941154D+00, SeN  
     6  5.81538300D+00, 5.97525968D+00, 6.04028176D+00, 6.08141217D+00, SeN  
     7       3*0.0D+00/                                                 SeN  
      DATA TQ_ZrN/                                                      071215
     1  0.699999789529, 0.758600210065, 0.850299913974, 1.089900163638, ZrN  
     2  1.456999961018, 1.856400053904, 2.106599920221, 2.313000213222, ZrN  
     3  2.453199874432, 2.609399687784, 2.785699869109, 2.921400105978, ZrN  
     4  3.054200289617, 3.316399953976, 3.597700329877, 3.833900011675, ZrN  
     5  3.935399934133, 4.000000000000,      9*0.0D+00/                 ZrN  
      DATA  Q_ZrN/                                                      071215
     1  1.17908965D+00, 1.23513224D+00, 1.32346665D+00, 1.55701098D+00, ZrN  
     2  1.91942943D+00, 2.31671410D+00, 2.56630935D+00, 2.77293967D+00, ZrN  
     3  2.91588534D+00, 3.08354518D+00, 3.29290235D+00, 3.47307744D+00, ZrN  
     4  3.66604959D+00, 4.08925836D+00, 4.58882870D+00, 5.03099881D+00, ZrN  
     5  5.22513812D+00, 5.34966087D+00,      9*0.0D+00/                 ZrN  
      DATA TQ_NOp/                                                      071215
     1  0.699999789529, 0.731200002896, 0.779800221589, 0.904799981617, NOp  
     2  1.238299796562, 1.647500039022, 2.110099679915, 2.614199740663, NOp  
     3  2.811400315150, 3.039099960232, 3.168500038177, 3.300700204473, NOp  
     4  3.568199815908, 3.728899634376, 3.869499881537, 3.949899826122, NOp  
     5  3.980899572511, 4.000000000000,      9*0.0D+00/                 NOp  
      DATA  Q_NOp/                                                      071215
     1  3.28116913D-01, 3.53317493D-01, 3.93433964D-01, 5.00668897D-01, NOp  
     2  8.06006630D-01, 1.20051724D+00, 1.65706250D+00, 2.15931841D+00, NOp  
     3  2.35864557D+00, 2.60469249D+00, 2.76081505D+00, 2.93660997D+00, NOp  
     4  3.34399372D+00, 3.61835218D+00, 3.87416429D+00, 4.02787368D+00, NOp  
     5  4.08934486D+00, 4.12810823D+00,      9*0.0D+00/                 NOp  
      DATA TQ_NSp/                                                      071215
     1  0.699999789529, 0.729500024342, 0.775300111842, 0.892699928873, NSp  
     2  1.208599984199, 1.564200036894, 1.999599990910, 2.238899763159, NSp  
     3  2.453799887825, 2.593800241171, 2.755900132446, 2.893300104388, NSp  
     4  3.039599972043, 3.184399969970, 3.534999843571, 3.666899995950, NSp  
     5  3.818599786041, 3.924600203359, 3.970700220546, 4.000000000000, NSp  
     6       7*0.0D+00/                                                 NSp  
      DATA  Q_NSp/                                                      071215
     1  3.26509502D-01, 3.50297160D-01, 3.88002703D-01, 4.88307124D-01, NSp  
     2  7.76120942D-01, 1.11720763D+00, 1.54539890D+00, 1.78293672D+00, NSp  
     3  1.99728780D+00, 2.13921079D+00, 2.31166955D+00, 2.47103818D+00, NSp  
     4  2.65913500D+00, 2.86576503D+00, 3.44241770D+00, 3.68277291D+00, NSp  
     5  3.97387211D+00, 4.18693554D+00, 4.28182812D+00, 4.34268386D+00, NSp  
     6       7*0.0D+00/                                                 NSp  
      DATA TQ_LiO/                                                      071215
     1  0.699999789529, 0.727699981235, 0.770399992341, 0.879600105345, LiO  
     2  1.171300037811, 1.488899807934, 1.951399838429, 2.148599720034, LiO  
     3  2.318499818330, 2.546700107599, 2.750499999493, 2.918000026767, LiO  
     4  3.063400186641, 3.299800184763, 3.433200066721, 3.676000202044, LiO  
     5  3.802200235933, 3.899699641757, 3.961700084834, 4.000000000000, LiO  
     6       7*0.0D+00/                                                 LiO  
      DATA  Q_LiO/                                                      071215
     1  1.25472781D+00, 1.27117370D+00, 1.29773376D+00, 1.37167114D+00, LiO  
     2  1.60144377D+00, 1.88470819D+00, 2.32607952D+00, 2.51932492D+00, LiO  
     3  2.68856647D+00, 2.93027001D+00, 3.17799129D+00, 3.41407547D+00, LiO  
     4  3.64569154D+00, 4.07308314D+00, 4.33781514D+00, 4.85689127D+00, LiO  
     5  5.14348094D+00, 5.36812007D+00, 5.50997280D+00, 5.59650658D+00, LiO  
     6       7*0.0D+00/                                                 LiO  
      DATA TQ_BeO/                                                      071215
     1  0.699999789529, 0.731699988970, 0.780900201812, 0.906799935515, BeO  
     2  1.254700112920, 1.609500127680, 2.060100421162, 2.266400374797, BeO  
     3  2.453999892289, 2.634700214758, 2.845800293363, 3.093000211451, BeO  
     4  3.169900072265, 3.254500100391, 3.397499810092, 3.453399875814, BeO  
     5  3.514100318955, 3.601300286617, 3.706699960254, 3.793900059819, BeO  
     6  3.915399982968, 3.967700219299, 4.000000000000,      4*0.0D+00/ BeO  
      DATA  Q_BeO/                                                      071215
     1  3.96222926D-01, 4.22943087D-01, 4.65126665D-01, 5.76403739D-01, BeO  
     2  9.00606046D-01, 1.24476103D+00, 1.69005820D+00, 1.89536476D+00, BeO  
     3  2.08276141D+00, 2.26639700D+00, 2.49677637D+00, 2.81194962D+00, BeO  
     4  2.92300069D+00, 3.05441449D+00, 3.31112834D+00, 3.42928920D+00, BeO  
     5  3.57185782D+00, 3.80235850D+00, 4.11169830D+00, 4.38032223D+00, BeO  
     6  4.75679431D+00, 4.91696010D+00, 5.01503110D+00,      4*0.0D+00/ BeO  
      DATA TQ_BO/                                                       071215
     1  0.699999789529, 0.730600019606, 0.778100180129, 0.900000092262, BO   
     2  1.226499940156, 1.623999978509, 2.055800331892, 2.324899816734, BO   
     3  2.555200301686, 2.747999944624, 2.994499871280, 3.142500163109, BO   
     4  3.294100062098, 3.481999579267, 3.666199981661, 3.770500429759, BO   
     5  3.873399969876, 3.949899826122, 3.980899572511, 4.000000000000, BO   
     6       7*0.0D+00/                                                 BO   
      DATA  Q_BO/                                                       071215
     1  6.69500767D-01, 6.94882878D-01, 7.35008693D-01, 8.41406225D-01, BO   
     2  1.14292098D+00, 1.52721180D+00, 1.95355422D+00, 2.22127210D+00, BO   
     3  2.45128381D+00, 2.64738448D+00, 2.92062743D+00, 3.10799996D+00, BO   
     4  3.32211423D+00, 3.61760073D+00, 3.93634008D+00, 4.13061900D+00, BO   
     5  4.33681211D+00, 4.50377800D+00, 4.57567021D+00, 4.62133112D+00, BO   
     6       7*0.0D+00/                                                 BO   
      DATA TQ_FO/                                                       071215
     1  0.699999789529, 0.727599978840, 0.770299989902, 0.879900112611, FO   
     2  1.161099866439, 1.519299984690, 1.995399895463, 2.208500015224, FO   
     3  2.388500309268, 2.635800240821, 2.809900413020, 2.945699715885, FO   
     4  3.074399804708, 3.315300032208, 3.459000011230, 3.585700058990, FO   
     5  3.679900294873, 3.779399785709, 3.920300107825, 3.967600217058, FO   
     6  4.000000000000,      6*0.0D+00/                                 FO   
      DATA  Q_FO/                                                       071215
     1  1.29021589D+00, 1.30787835D+00, 1.33629365D+00, 1.41465946D+00, FO   
     2  1.64285667D+00, 1.96769848D+00, 2.42627170D+00, 2.63599107D+00, FO   
     3  2.81512250D+00, 3.07426131D+00, 3.27993548D+00, 3.45997355D+00, FO   
     4  3.64753799D+00, 4.03978025D+00, 4.29659663D+00, 4.53754357D+00, FO   
     5  4.72770520D+00, 4.94074178D+00, 5.26041351D+00, 5.36982858D+00, FO   
     6  5.44459616D+00,      6*0.0D+00/                                 FO   
      DATA TQ_NaO/                                                      071215
     1  0.699999789529, 0.766100085453, 0.868799904485, 1.149600088415, NaO  
     2  1.438400006259, 1.765700084890, 1.978699978690, 2.164699951808, NaO  
     3  2.387200281327, 2.619499862667, 2.806600343517, 2.985999676060, NaO  
     4  3.233099607001, 3.504800101978, 3.794700076591, 3.928000278898, NaO  
     5  4.000000000000,     10*0.0D+00/                                 NaO  
      DATA  Q_NaO/                                                      071215
     1  1.58528954D+00, 1.64169175D+00, 1.73190438D+00, 1.99024263D+00, NaO  
     2  2.26695038D+00, 2.58744739D+00, 2.79825020D+00, 2.98525557D+00, NaO  
     3  3.22431727D+00, 3.51405598D+00, 3.78843333D+00, 4.08759120D+00, NaO  
     4  4.54597188D+00, 5.08350329D+00, 5.66852860D+00, 5.93788826D+00, NaO  
     5  6.08312998D+00,     10*0.0D+00/                                 NaO  
      DATA TQ_MgO/                                                      071215
     1  0.699999789529, 0.768300027425, 0.874899991518, 1.157699900341, MgO  
     2  1.864300027927, 2.140200306882, 2.264400331110, 2.393000120990, MgO  
     3  2.519900454649, 2.622399929189, 2.716500202529, 2.808400381427, MgO  
     4  2.924200166605, 3.035299870469, 3.171500106854, 3.271500353067, MgO  
     5  3.385200236155, 3.498499955066, 3.596000292313, 3.691299627677, MgO  
     6  3.847700317764, 3.938899681874, 4.000000000000,      4*0.0D+00/ MgO  
      DATA  Q_MgO/                                                      071215
     1  8.08401812D-01, 8.73206242D-01, 9.75348780D-01, 1.25051298D+00, MgO  
     2  1.95056925D+00, 2.22598493D+00, 2.35095483D+00, 2.48336126D+00, MgO  
     3  2.62099602D+00, 2.74101901D+00, 2.86239221D+00, 2.99756781D+00, MgO  
     4  3.20486076D+00, 3.45367835D+00, 3.81668569D+00, 4.10392123D+00, MgO  
     5  4.43248790D+00, 4.75164772D+00, 5.01687604D+00, 5.26814663D+00, MgO  
     6  5.66947703D+00, 5.90209964D+00, 6.05902067D+00,      4*0.0D+00/ MgO  
      DATA TQ_AlO/                                                      071215
     1  0.699999789529, 0.784200111412, 0.914399960255, 1.292299998225, AlO  
     2  1.931299852897, 2.145599929623, 2.346200303668, 2.660399869739, AlO  
     3  2.792400032404, 2.934599969375, 3.095000259641, 3.273100238417, AlO  
     4  3.354900019169, 3.435499903933, 3.673400140158, 3.767800412754, AlO  
     5  3.883000190390, 3.952399880785, 4.000000000000,      8*0.0D+00/ AlO  
      DATA  Q_AlO/                                                      071215
     1  1.06459376D+00, 1.14406355D+00, 1.26855260D+00, 1.63700765D+00, AlO  
     2  2.27087606D+00, 2.48467755D+00, 2.68587569D+00, 3.02093073D+00, AlO  
     3  3.18105039D+00, 3.37218568D+00, 3.61384269D+00, 3.92507039D+00, AlO  
     4  4.08885656D+00, 4.26477326D+00, 4.85179097D+00, 5.10169956D+00, AlO  
     5  5.41431896D+00, 5.60652559D+00, 5.74003569D+00,      8*0.0D+00/ AlO  
      DATA TQ_SiO/                                                      071215
     1  0.699999789529, 0.760700227885, 0.856100047022, 1.098199956358, SiO  
     2  1.501300022737, 1.940100036964, 2.186299854832, 2.394200032783, SiO  
     3  2.553300257808, 2.724999919775, 2.904999718740, 3.146599857168, SiO  
     4  3.452799861305, 3.615999802564, 3.694499693896, 3.771800335685, SiO  
     5  3.843800230028, 3.913299929012, 3.965400167754, 3.986499698236, SiO  
     6  4.000000000000,      6*0.0D+00/                                 SiO  
      DATA  Q_SiO/                                                      071215
     1  7.12454195D-01, 7.69174496D-01, 8.59354622D-01, 1.09241058D+00, SiO  
     2  1.48822753D+00, 1.92406116D+00, 2.16959066D+00, 2.37756696D+00, SiO  
     3  2.53937224D+00, 2.72393906D+00, 2.93966700D+00, 3.27699433D+00, SiO  
     4  3.77907271D+00, 4.07388858D+00, 4.22173968D+00, 4.37270612D+00, SiO  
     5  4.52259404D+00, 4.68496674D+00, 4.82647408D+00, 4.89023688D+00, SiO  
     6  4.93324857D+00,      6*0.0D+00/                                 SiO  
      DATA TQ_PO/                                                       071215
     1  0.699999789529, 0.728900009973, 0.773500067944, 0.888099915294, PO   
     2  1.187599994154, 1.534600114673, 1.996299915916, 2.234499656452, PO   
     3  2.434599962314, 2.708700023776, 2.850500347385, 3.041800017332, PO   
     4  3.237299709008, 3.427000233125, 3.598600349764, 3.724599932340, PO   
     5  3.798700160452, 3.868199851834, 3.947299767417, 3.979799566676, PO   
     6  4.000000000000,      6*0.0D+00/                                 PO   
      DATA  Q_PO/                                                       071215
     1  1.39874693D+00, 1.42023256D+00, 1.45426642D+00, 1.54601662D+00, PO   
     2  1.80727012D+00, 2.13281403D+00, 2.58303237D+00, 2.81876845D+00, PO   
     3  3.01838286D+00, 3.30547819D+00, 3.47136870D+00, 3.72433757D+00, PO   
     4  4.01926891D+00, 4.33663268D+00, 4.64570240D+00, 4.88444898D+00, PO   
     5  5.02989317D+00, 5.17086824D+00, 5.33955723D+00, 5.41257554D+00, PO   
     6  5.45936748D+00,      6*0.0D+00/                                 PO   
      DATA TQ_SO/                                                       071215
     1  0.699999789529, 0.781900174418, 0.908099905548, 1.272000029951, SO   
     2  1.600499903753, 1.931799863611, 2.166999997362, 2.382800186757, SO   
     3  2.534699828451, 2.684999946435, 2.902299658024, 3.156499751500, SO   
     4  3.390300335801, 3.552700252888, 3.730599570549, 3.858699756127, SO   
     5  3.946799756128, 3.979399595417, 4.000000000000,      8*0.0D+00/ SO   
      DATA  Q_SO/                                                       071215
     1  1.19309773D+00, 1.26980206D+00, 1.38970786D+00, 1.74307177D+00, SO   
     2  2.06734628D+00, 2.39670209D+00, 2.63125817D+00, 2.84727615D+00, SO   
     3  3.00232365D+00, 3.16415527D+00, 3.42716781D+00, 3.79215025D+00, SO   
     4  4.18647505D+00, 4.49791303D+00, 4.87501072D+00, 5.16716773D+00, SO   
     5  5.37991693D+00, 5.46212402D+00, 5.51528257D+00,      8*0.0D+00/ SO   
      DATA TQ_ClO/                                                      071215
     1  0.699999789529, 0.749399973420, 0.826600001417, 1.025200102942, ClO  
     2  1.310699972318, 1.652200039448, 1.878100084408, 2.153999716824, ClO  
     3  2.296700113529, 2.551900225477, 2.713200128169, 2.871099916308, ClO  
     4  3.187599737901, 3.395499956122, 3.584900040247, 3.701199835031, ClO  
     5  3.835700050882, 3.933200092696, 3.973899990614, 4.000000000000, ClO  
     6       7*0.0D+00/                                                 ClO  
      DATA  Q_ClO/                                                      071215
     1  1.15120842D+00, 1.18998735D+00, 1.25274020D+00, 1.42381579D+00, ClO  
     2  1.68619400D+00, 2.01410098D+00, 2.23640935D+00, 2.52549079D+00, ClO  
     3  2.69204164D+00, 3.02542145D+00, 3.25982041D+00, 3.50740723D+00, ClO  
     4  4.05353452D+00, 4.44286923D+00, 4.81694084D+00, 5.05836746D+00, ClO  
     5  5.35577634D+00, 5.58839614D+00, 5.69006216D+00, 5.75653085D+00, ClO  
     6       7*0.0D+00/                                                 ClO  
      DATA TQ_KO/                                                       071215
     1  0.699999789529, 0.747699934427, 0.822000112744, 1.018000022309, KO   
     2  1.550400171855, 1.716999867137, 1.874499995348, 2.009900227328, KO   
     3  2.132300143876, 2.217300196840, 2.308800400668, 2.432800087942, KO   
     4  2.572199750232, 2.758300191537, 2.968000224860, 3.167600016263, KO   
     5  3.623199970788, 3.815999973747, 3.933100099903, 4.000000000000, KO   
     6       7*0.0D+00/                                                 KO   
      DATA  Q_KO/                                                       071215
     1  1.07460388D+00, 1.11957838D+00, 1.19019271D+00, 1.37903737D+00, KO   
     2  1.90260898D+00, 2.06807996D+00, 2.22568259D+00, 2.36600377D+00, KO   
     3  2.50482755D+00, 2.61306779D+00, 2.74333741D+00, 2.94282168D+00, KO   
     4  3.19245207D+00, 3.55048951D+00, 3.96819406D+00, 4.36968317D+00, KO   
     5  5.28610892D+00, 5.67297271D+00, 5.90740735D+00, 6.04081850D+00, KO   
     6       7*0.0D+00/                                                 KO   
      DATA TQ_CaO/                                                      071215
     1  0.699999789529, 0.782500157982, 0.909999861751, 1.280700202902, CaO  
     2  1.919400060995, 2.082999995499, 2.239699782561, 2.592700217597, CaO  
     3  2.794300075431, 2.968700240954, 3.086000061480, 3.201499831824, CaO  
     4  3.283399822148, 3.356699883261, 3.412699900839, 3.469900249600, CaO  
     5  3.554600292540, 3.679000273451, 3.788299937974, 3.917700042064, CaO  
     6  3.967600217058, 4.000000000000,      5*0.0D+00/                 CaO  
      DATA  Q_CaO/                                                      071215
     1  9.12672764D-01, 9.91951029D-01, 1.11553931D+00, 1.47967522D+00, CaO  
     2  2.11463501D+00, 2.27796437D+00, 2.43544490D+00, 2.81879034D+00, CaO  
     3  3.07975476D+00, 3.33663901D+00, 3.52465305D+00, 3.72257206D+00, CaO  
     4  3.87476734D+00, 4.02550265D+00, 4.15413282D+00, 4.30030804D+00, CaO  
     5  4.54509417D+00, 4.95093155D+00, 5.32638269D+00, 5.76596400D+00, CaO  
     6  5.93027238D+00, 6.03484780D+00,      5*0.0D+00/                 CaO  
      DATA TQ_ScO/                                                      071215
     1  0.699999789529, 0.761700201508, 0.858300097489, 1.109500129319, ScO  
     2  1.488399819532, 1.891499895915, 2.109199734902, 2.315000069625, ScO  
     3  2.470600208225, 2.622899940650, 2.740499785030, 2.861299675076, ScO  
     4  3.073399779665, 3.342700224943, 3.502900057683, 3.591000181831, ScO  
     5  3.676000202044, 3.771300371867, 3.862299717025, 3.946899758386, ScO  
     6  3.979599581047, 4.000000000000,      5*0.0D+00/                 ScO  
      DATA  Q_ScO/                                                      071215
     1  1.15533939D+00, 1.21419683D+00, 1.30708940D+00, 1.55177352D+00, ScO  
     2  1.92586227D+00, 2.32692997D+00, 2.54416831D+00, 2.75030973D+00, ScO  
     3  2.90948881D+00, 3.07447490D+00, 3.21286562D+00, 3.36793891D+00, ScO  
     4  3.67422160D+00, 4.11923077D+00, 4.40932823D+00, 4.57875964D+00, ScO  
     5  4.75290766D+00, 4.96665013D+00, 5.19418306D+00, 5.42744269D+00, ScO  
     6  5.52270428D+00, 5.58338618D+00,      5*0.0D+00/                 ScO  
      DATA TQ_TiO/                                                      071215
     1  0.699999789529, 0.730000036317, 0.776700145986, 0.894899978114, TiO  
     2  1.071900133719, 1.247799961605, 1.428099865869, 1.612500077608, TiO  
     3  1.761100209381, 1.908299873790, 2.264100324557, 2.458900001666, TiO  
     4  2.665499982449, 2.816399958096, 2.968300231757, 3.314000124664, TiO  
     5  3.473599993697, 3.620099901835, 3.744999875963, 3.845800275021, TiO  
     6  3.939999602593, 3.976599796608, 4.000000000000,      4*0.0D+00/ TiO  
      DATA  Q_TiO/                                                      071215
     1  1.27749662D+00, 1.29842465D+00, 1.33208880D+00, 1.42253683D+00, TiO  
     2  1.56926590D+00, 1.72497185D+00, 1.89355029D+00, 2.08055929D+00, TiO  
     3  2.24775130D+00, 2.42916933D+00, 2.90988399D+00, 3.17869902D+00, TiO  
     4  3.46714990D+00, 3.68778799D+00, 3.92313213D+00, 4.51017970D+00, TiO  
     5  4.80226416D+00, 5.08512780D+00, 5.34696031D+00, 5.58019248D+00, TiO  
     6  5.81958892D+00, 5.91825529D+00, 5.98290173D+00,      4*0.0D+00/ TiO  
      DATA TQ_VO/                                                       071215
     1  0.699999789529, 0.761900196233, 0.858700106665, 1.110000140845, VO   
     2  1.492799841149, 1.899700075272, 2.124599977753, 2.334200026726, VO   
     3  2.485699672217, 2.637100271623, 2.750700004417, 2.868099845377, VO   
     4  3.056700354506, 3.324799802020, 3.485899665759, 3.645999904334, VO   
     5  3.752200039322, 3.851600256701, 3.943499681618, 3.978099688827, VO   
     6  4.000000000000,      6*0.0D+00/                                 VO   
      DATA  Q_VO/                                                       071215
     1  1.42932781D+00, 1.48818100D+00, 1.58101412D+00, 1.82535678D+00, VO   
     2  2.20298746D+00, 2.60771800D+00, 2.83211217D+00, 3.04202860D+00, VO   
     3  3.19688845D+00, 3.36044058D+00, 3.49338914D+00, 3.64281912D+00, VO   
     4  3.91073979D+00, 4.34615762D+00, 4.63528803D+00, 4.95006477D+00, VO   
     5  5.18237000D+00, 5.42124514D+00, 5.66045950D+00, 5.75456396D+00, VO   
     6  5.81511108D+00,      6*0.0D+00/                                 VO   
      DATA TQ_CrO/                                                      071215
     1  0.699999789529, 0.726899962076, 0.768500022150, 0.875099996362, CrO  
     2  1.156099941254, 1.510400216869, 1.905099954000, 2.124299971461, CrO  
     3  2.316699947567, 2.566899921447, 2.672200131714, 2.792500034668, CrO  
     4  3.016900381285, 3.322099743478, 3.485599659106, 3.647899770113, CrO  
     5  3.747399929348, 3.842700205281, 3.936999818815, 4.000000000000, CrO  
     6       7*0.0D+00/                                                 CrO  
      DATA  Q_CrO/                                                      071215
     1  1.90390201D+00, 1.92584562D+00, 1.96033826D+00, 2.05146331D+00, CrO  
     2  2.30569538D+00, 2.64342480D+00, 3.03023723D+00, 3.24735335D+00, CrO  
     3  3.43951493D+00, 3.70201864D+00, 3.82326886D+00, 3.97358611D+00, CrO  
     4  4.29151840D+00, 4.79611083D+00, 5.09470450D+00, 5.41123949D+00, CrO  
     5  5.61972202D+00, 5.83376122D+00, 6.05973588D+00, 6.21637864D+00, CrO  
     6       7*0.0D+00/                                                 CrO  
      DATA TQ_MnO/                                                      071215
     1  0.699999789529, 0.759900243756, 0.853499987380, 1.099299928222, MnO  
     2  1.459300011006, 1.835000030689, 2.054100292742, 2.259200216374, MnO  
     3  2.404399713795, 2.541899999962, 2.656499775124, 2.771900332397, MnO  
     4  2.929300277032, 3.068799783994, 3.349600380144, 3.518000410032, MnO  
     5  3.670000059230, 3.776300010041, 3.863099735304, 3.945899735807, MnO  
     6  3.979199609788, 4.000000000000,      5*0.0D+00/                 MnO  
      DATA  Q_MnO/                                                      071215
     1  1.69982693D+00, 1.75737932D+00, 1.84790495D+00, 2.08819221D+00, MnO  
     2  2.44412062D+00, 2.81799693D+00, 3.03658351D+00, 3.24197720D+00, MnO  
     3  3.39037362D+00, 3.53816316D+00, 3.67081374D+00, 3.81583001D+00, MnO  
     4  4.03406214D+00, 4.24664710D+00, 4.72037343D+00, 5.02765873D+00, MnO  
     5  5.31876074D+00, 5.53284125D+00, 5.71672260D+00, 5.90088125D+00, MnO  
     6  5.97721329D+00, 6.02545274D+00,      5*0.0D+00/                 MnO  
      DATA TQ_FeO/                                                      071215
     1  0.699999789529, 0.727999988419, 0.771100009413, 0.881900068180, FeO  
     2  1.169400054392, 1.511500188172, 1.962399992968, 2.165699971614, FeO  
     3  2.353200151229, 2.599700367612, 2.805100311925, 3.185299904701, FeO  
     4  3.320999719627, 3.452699858887, 3.660599867347, 3.827599868707, FeO  
     5  3.933300085488, 3.973799997799, 4.000000000000,      8*0.0D+00/ FeO  
      DATA  Q_FeO/                                                      071215
     1  1.98605315D+00, 2.00580996D+00, 2.03713127D+00, 2.12217639D+00, FeO  
     2  2.36577683D+00, 2.68182660D+00, 3.11851951D+00, 3.31892706D+00, FeO  
     3  3.50559541D+00, 3.76324627D+00, 4.00692977D+00, 4.56464848D+00, FeO  
     4  4.79977504D+00, 5.04505598D+00, 5.46011019D+00, 5.81387012D+00, FeO  
     5  6.04775720D+00, 6.13965547D+00, 6.19978502D+00,      8*0.0D+00/ FeO  
      DATA TQ_NiO/                                                      071215
     1  0.699999789529, 0.758000194515, 0.848399945697, 1.088200121827, NiO  
     2  1.413199923989, 1.769599979344, 1.964999925274, 2.151299652970, NiO  
     3  2.441999628947, 2.559600403299, 2.689699606253, 2.933300060605, NiO  
     4  3.262600282863, 3.455299921758, 3.662099897967, 3.756100131751, NiO  
     5  3.850700320154, 3.942799665813, 3.977899703198, 4.000000000000, NiO  
     6       7*0.0D+00/                                                 NiO  
      DATA  Q_NiO/                                                      071215
     1  8.68912890D-01, 9.24322760D-01, 1.01132049D+00, 1.24489745D+00, NiO  
     2  1.56545991D+00, 1.91963254D+00, 2.11442793D+00, 2.30115178D+00, NiO  
     3  2.60876668D+00, 2.74780166D+00, 2.91598843D+00, 3.27308034D+00, NiO  
     4  3.82744034D+00, 4.17869914D+00, 4.57981498D+00, 4.77610597D+00, NiO  
     5  4.98539470D+00, 5.20012930D+00, 5.28431831D+00, 5.33784502D+00, NiO  
     6       7*0.0D+00/                                                 NiO  
      DATA TQ_CuO/                                                      071215
     1  0.699999789529, 0.747699934427, 0.822100110324, 1.013300134595, CuO  
     2  1.279100200561, 1.567200107539, 1.706999979017, 1.831399945958, CuO  
     3  2.141200237019, 2.323799792915, 2.507700186581, 2.712200105635, CuO  
     4  3.115999785654, 3.354400056921, 3.466200169549, 3.576299837052, CuO  
     5  3.689699621670, 3.805100299503, 3.917400034356, 4.000000000000, CuO  
     6       7*0.0D+00/                                                 CuO  
      DATA  Q_CuO/                                                      071215
     1  1.57079621D+00, 1.61107642D+00, 1.67537327D+00, 1.84719334D+00, CuO  
     2  2.09684402D+00, 2.37566022D+00, 2.51292396D+00, 2.63661530D+00, CuO  
     3  2.96698610D+00, 3.19153866D+00, 3.44581592D+00, 3.76112280D+00, CuO  
     4  4.46827151D+00, 4.92696888D+00, 5.15208572D+00, 5.38185441D+00, CuO  
     5  5.63026047D+00, 5.89923124D+00, 6.17708946D+00, 6.38826286D+00, CuO  
     6       7*0.0D+00/                                                 CuO  
      DATA TQ_GaO/                                                      071215
     1  0.699999789529, 0.758200199698, 0.849299923982, 1.086000067718, GaO  
     2  1.438600010800, 1.813400027824, 2.006900158441, 2.203099872043, GaO  
     3  2.358199779492, 2.516400381514, 2.744899878658, 3.005800129426, GaO  
     4  3.254300096026, 3.433600038410, 3.596700307781, 3.711400065586, GaO  
     5  3.819499721066, 3.931800193599, 4.000000000000,      8*0.0D+00/ GaO  
      DATA  Q_GaO/                                                      071215
     1  1.23034010D+00, 1.28629742D+00, 1.37444164D+00, 1.60583495D+00, GaO  
     2  1.95438623D+00, 2.32735713D+00, 2.52042015D+00, 2.71684427D+00, GaO  
     3  2.87511227D+00, 3.04581036D+00, 3.32497267D+00, 3.70385645D+00, GaO  
     4  4.11903758D+00, 4.44480648D+00, 4.75801979D+00, 4.98882881D+00, GaO  
     5  5.21524483D+00, 5.45867538D+00, 5.60930468D+00,      8*0.0D+00/ GaO  
      DATA TQ_GeO/                                                      071215
     1  0.699999789529, 0.761000219972, 0.856500056198, 1.104400011750, GeO  
     2  1.484299914635, 1.885299989210, 2.108199806179, 2.318599811150, GeO  
     3  2.470400222738, 2.621199901683, 2.846200302041, 3.080199948836, GeO  
     4  3.399199685966, 3.575699822718, 3.656699778099, 3.741899807007, GeO  
     5  3.817999829358, 3.897299818612, 3.959300031451, 3.984499653334, GeO  
     6  4.000000000000,      6*0.0D+00/                                 GeO  
      DATA  Q_GeO/                                                      071215
     1  8.77245139D-01, 9.35585835D-01, 1.02760994D+00, 1.26938152D+00, GeO  
     2  1.64466778D+00, 2.04372898D+00, 2.26617166D+00, 2.47687210D+00, GeO  
     3  2.63189139D+00, 2.79440578D+00, 3.06780956D+00, 3.40265849D+00, GeO  
     4  3.93467307D+00, 4.25843633D+00, 4.41383678D+00, 4.58442210D+00, GeO  
     5  4.74768757D+00, 4.93706129D+00, 5.10515615D+00, 5.17939698D+00, GeO  
     6  5.22682198D+00,      6*0.0D+00/                                 GeO  
      DATA TQ_AsO/                                                      071215
     1  0.699999789529, 0.730000036317, 0.776400138669, 0.894599971399, AsO  
     2  1.219999800085, 1.531200201690, 1.860000136486, 2.099600380586, AsO  
     3  2.299500183258, 2.434299983252, 2.590900179022, 2.719000258863, AsO  
     4  2.881900151661, 3.148299730314, 3.378400071639, 3.511000246561, AsO  
     5  3.631700157084, 3.728799641306, 3.819499721066, 3.930600280088, AsO  
     6  3.972600084024, 3.988499743138, 4.000000000000,      4*0.0D+00/ AsO  
      DATA  Q_AsO/                                                      071215
     1  1.23791483D+00, 1.26273177D+00, 1.30176618D+00, 1.40431824D+00, AsO  
     2  1.70309702D+00, 2.00200237D+00, 2.32458283D+00, 2.56190990D+00, AsO  
     3  2.76141283D+00, 2.89973276D+00, 3.07409599D+00, 3.23732280D+00, AsO  
     4  3.47950955D+00, 3.94783728D+00, 4.39493604D+00, 4.66204226D+00, AsO  
     5  4.90980845D+00, 5.11318894D+00, 5.30915523D+00, 5.56584425D+00, AsO  
     6  5.67083835D+00, 5.71206076D+00, 5.74242838D+00,      4*0.0D+00/ AsO  
      DATA TQ_SeO/                                                      071215
     1  0.699999789529, 0.808299988314, 0.989099808396, 1.333500015170, SeO  
     2  1.467299846638, 1.594900032552, 1.852899973616, 2.019500453221, SeO  
     3  2.202499856134, 2.360799665049, 2.490499792382, 2.755900132446, SeO  
     4  2.940599603626, 3.221700133940, 3.537799907291, 3.679700290113, SeO  
     5  3.820899706727, 3.928100281119, 3.971800141507, 4.000000000000, SeO  
     6       7*0.0D+00/                                                 SeO  
      DATA  Q_SeO/                                                      071215
     1  1.37209386D+00, 1.47610598D+00, 1.65176432D+00, 1.99074464D+00, SeO  
     2  2.12348055D+00, 2.25110802D+00, 2.52170077D+00, 2.71460111D+00, SeO  
     3  2.94220751D+00, 3.14496459D+00, 3.31290480D+00, 3.67298434D+00, SeO  
     4  3.94874248D+00, 4.41678845D+00, 5.01606774D+00, 5.31389524D+00, SeO  
     5  5.63385093D+00, 5.89693794D+00, 6.00987054D+00, 6.08450805D+00, SeO  
     6       7*0.0D+00/                                                 SeO  
      DATA TQ_BrO/                                                      071215
     1  0.699999789529, 0.773100058189, 0.886199962146, 1.198099846626, BrO  
     2  1.489799787058, 1.798200141251, 2.021600346866, 2.212400101630, BrO  
     3  2.349000366648, 2.504700106754, 2.640400311355, 2.799300188660, BrO  
     4  3.053600274043, 3.276000030615, 3.403799704405, 3.516200367997, BrO  
     5  3.723100036280, 3.880200123929, 4.000000000000,      8*0.0D+00/ BrO  
      DATA  Q_BrO/                                                      071215
     1  1.58291667D+00, 1.64530740D+00, 1.74493906D+00, 2.03346025D+00, BrO  
     2  2.31421561D+00, 2.61682553D+00, 2.83803955D+00, 3.02852082D+00, BrO  
     3  3.16911045D+00, 3.34399790D+00, 3.52036816D+00, 3.76326308D+00, BrO  
     4  4.22140916D+00, 4.66132738D+00, 4.92183035D+00, 5.15421168D+00, BrO  
     5  5.59278725D+00, 5.93965947D+00, 6.21154061D+00,      8*0.0D+00/ BrO  
      DATA TQ_RbO/                                                      071215
     1  0.699999789529, 0.822900090963, 1.020399984405, 1.613100062606, RbO  
     2  1.802900110029, 2.000700016074, 2.132400146174, 2.291199976560, RbO  
     3  2.427000226231, 2.520900389905, 2.612999713039, 2.740199778646, RbO  
     4  2.882000154107, 3.065999992774, 3.234699645861, 3.483699616969, RbO  
     5  3.787899928543, 3.920500112268, 4.000000000000,      8*0.0D+00/ RbO  
      DATA  Q_RbO/                                                      071215
     1  1.22751688D+00, 1.34597959D+00, 1.53853350D+00, 2.12485537D+00, RbO  
     2  2.31390727D+00, 2.51222714D+00, 2.64882281D+00, 2.83137324D+00, RbO  
     3  3.01820710D+00, 3.16914377D+00, 3.33398765D+00, 3.58307046D+00, RbO  
     4  3.87827322D+00, 4.27092026D+00, 4.63029295D+00, 5.15305915D+00, RbO  
     5  5.77979490D+00, 6.04966966D+00, 6.21018431D+00,      8*0.0D+00/ RbO  
      DATA TQ_SrO/                                                      071215
     1  0.699999789529, 0.782400160721, 0.909299877887, 1.286800034840, SrO  
     2  1.889399879205, 2.048100154827, 2.199699783837, 2.549600172629, SrO  
     3  2.755400120136, 2.934699962357, 3.080399952720, 3.213100095152, SrO  
     4  3.308800382734, 3.418000040299, 3.499599980738, 3.594500259169, SrO  
     5  3.713800117355, 3.810700356378, 3.924700205581, 3.970700220546, SrO  
     6  4.000000000000,      6*0.0D+00/                                 SrO  
      DATA  Q_SrO/                                                      071215
     1  1.02858613D+00, 1.10855423D+00, 1.23250955D+00, 1.60499566D+00, SrO  
     2  2.20498686D+00, 2.36357994D+00, 2.51621979D+00, 2.89814495D+00, SrO  
     3  3.16695525D+00, 3.43437031D+00, 3.67254629D+00, 3.90515424D+00, SrO  
     4  4.08652066D+00, 4.32158492D+00, 4.52908796D+00, 4.81105871D+00, SrO  
     5  5.21341228D+00, 5.55763540D+00, 5.96146752D+00, 6.12093554D+00, SrO  
     6  6.22097187D+00,      6*0.0D+00/                                 SrO  
      DATA TQ_YO/                                                       071215
     1  0.699999789529, 0.761700201508, 0.858000090607, 1.114500024038, YO   
     2  1.470999803633, 1.849499919156, 2.071599723347, 2.278499850407, YO   
     3  2.566099978385, 2.675800213229, 2.801600238210, 3.034099842122, YO   
     4  3.341800204700, 3.516400372667, 3.607899808972, 3.687799754001, YO   
     5  3.872599951779, 3.950499839297, 3.981099577001, 4.000000000000, YO   
     6       7*0.0D+00/                                                 YO   
      DATA  Q_YO/                                                       071215
     1  1.27126296D+00, 1.33081788D+00, 1.42432667D+00, 1.67582553D+00, YO   
     2  2.02887016D+00, 2.40583600D+00, 2.62755277D+00, 2.83494402D+00, YO   
     3  3.13763434D+00, 3.26557609D+00, 3.42533103D+00, 3.76033008D+00, YO   
     4  4.27346584D+00, 4.59230076D+00, 4.76749270D+00, 4.92700174D+00, YO   
     5  5.33456359D+00, 5.52852200D+00, 5.60865942D+00, 5.65923555D+00, YO   
     6       7*0.0D+00/                                                 YO   
      DATA TQ_ZrO/                                                      071215
     1  0.699999789529, 0.843200071162, 1.066000091269, 1.797200120694, ZrO  
     2  2.051000221351, 2.260600248105, 2.367999839486, 2.467500195377, ZrO  
     3  2.588100115220, 2.681300214237, 2.766800397648, 2.844800271667, ZrO  
     4  2.970300248954, 3.075799839769, 3.164799948087, 3.262500280462, ZrO  
     5  3.370299905405, 3.478399649609, 3.725299883834, 3.806300325808, ZrO  
     6  3.882800185643, 3.954899935374, 3.982899617413, 4.000000000000, ZrO  
     7       3*0.0D+00/                                                 ZrO  
      DATA  Q_ZrO/                                                      071215
     1  9.34826487D-01, 1.07307345D+00, 1.29079953D+00, 2.01590141D+00, ZrO  
     2  2.26913837D+00, 2.47892375D+00, 2.58835408D+00, 2.69490053D+00, ZrO  
     3  2.84077564D+00, 2.97588828D+00, 3.12270265D+00, 3.27621971D+00, ZrO  
     4  3.55446661D+00, 3.80598186D+00, 4.02293544D+00, 4.26173682D+00, ZrO  
     5  4.52388712D+00, 4.78569846D+00, 5.39822378D+00, 5.61179241D+00, ZrO  
     6  5.82197161D+00, 6.02769967D+00, 6.10945693D+00, 6.15986384D+00, ZrO  
     7       3*0.0D+00/                                                 ZrO  
      DATA TQ_NbO/                                                      071215
     1  0.699999789529, 0.761100217334, 0.856600058492, 1.105700041718, NbO  
     2  1.484099919274, 1.878200086882, 2.106999891710, 2.320499721459, NbO  
     3  2.469000229209, 2.615899779796, 2.834000032459, 3.061800305944, NbO  
     4  3.356299913463, 3.512500281590, 3.668800034735, 3.867199828985, NbO  
     5  3.948399792254, 3.980199556795, 4.000000000000,      8*0.0D+00/ NbO  
      DATA  Q_NbO/                                                      071215
     1  1.52753513D+00, 1.58626615D+00, 1.67867437D+00, 1.92230113D+00, NbO  
     2  2.29661060D+00, 2.68898479D+00, 2.91734219D+00, 3.13113182D+00, NbO  
     3  3.28268285D+00, 3.44049831D+00, 3.70325040D+00, 4.02499473D+00, NbO  
     4  4.50758455D+00, 4.78817425D+00, 5.08681200D+00, 5.50777787D+00, NbO  
     5  5.69919645D+00, 5.77746738D+00, 5.82712300D+00,      8*0.0D+00/ NbO  
      DATA TQ_InO/                                                      071215
     1  0.699999789529, 0.758700212656, 0.850399916268, 1.093800068901, InO  
     2  1.429199836745, 1.790499982959, 1.998799972729, 2.195499699700, InO  
     3  2.474499925231, 2.582199976990, 2.712100103382, 2.860299650032, InO  
     4  2.997499941491, 3.335500051963, 3.525100085992, 3.693399671133, InO  
     5  3.875700021907, 3.950999850215, 3.981299581491, 4.000000000000, InO  
     6       7*0.0D+00/                                                 InO  
      DATA  Q_InO/                                                      071215
     1  1.22558543D+00, 1.28199771D+00, 1.37069191D+00, 1.60861827D+00, InO  
     2  1.94010757D+00, 2.29951407D+00, 2.50727765D+00, 2.70437179D+00, InO  
     3  2.99786399D+00, 3.12306687D+00, 3.28759810D+00, 3.49506921D+00, InO  
     4  3.70531606D+00, 4.28348045D+00, 4.63535695D+00, 4.96157020D+00, InO  
     5  5.33551697D+00, 5.49911125D+00, 5.56660642D+00, 5.60871680D+00, InO  
     6       7*0.0D+00/                                                 InO  
      DATA TQ_SnO/                                                      071215
     1  0.699999789529, 0.856100047022, 1.109100120098, 1.831899957726, SnO  
     2  2.043500049525, 2.245599908746, 2.397099819615, 2.543900044811, SnO  
     3  2.657699804533, 2.777599904754, 3.009700216454, 3.365299774627, SnO  
     4  3.550600209062, 3.637200270468, 3.723200029351, 3.810900341940, SnO  
     5  3.892000209165, 3.957199985596, 3.983699635374, 4.000000000000, SnO  
     6       7*0.0D+00/                                                 SnO  
      DATA  Q_SnO/                                                      071215
     1  1.00700858D+00, 1.15863512D+00, 1.40708263D+00, 2.12528741D+00, SnO  
     2  2.33653424D+00, 2.53905394D+00, 2.69411009D+00, 2.85289455D+00, SnO  
     3  2.98611282D+00, 3.13890540D+00, 3.47390206D+00, 4.07451870D+00, SnO  
     4  4.41951899D+00, 4.58743099D+00, 4.75955122D+00, 4.94374612D+00, SnO  
     5  5.12710955D+00, 5.28798516D+00, 5.35752894D+00, 5.40158395D+00, SnO  
     6       7*0.0D+00/                                                 SnO  
      DATA TQ_SbO/                                                      071215
     1  0.699999789529, 0.758900217840, 0.851199934619, 1.088600131665, SbO  
     2  1.469399794978, 1.880900107265, 2.119399869070, 2.331599965094, SbO  
     3  2.566899921447, 2.810300393702, 2.976299811170, 3.133400198405, SbO  
     4  3.320099700113, 3.457599977376, 3.607499837920, 3.747999942695, SbO  
     5  3.898199752291, 3.960200051218, 3.984699657824, 4.000000000000, SbO  
     6       7*0.0D+00/                                                 SbO  
      DATA  Q_SbO/                                                      071215
     1  1.34909864D+00, 1.40065072D+00, 1.48322725D+00, 1.70323668D+00, SbO  
     2  2.06999977D+00, 2.47539969D+00, 2.71238392D+00, 2.92570349D+00, SbO  
     3  3.17790924D+00, 3.48431050D+00, 3.73466178D+00, 4.00836029D+00, SbO  
     4  4.37471604D+00, 4.66381551D+00, 4.99002515D+00, 5.30701533D+00, SbO  
     5  5.67398113D+00, 5.84014557D+00, 5.90859307D+00, 5.95213592D+00, SbO  
     6       7*0.0D+00/                                                 SbO  
      DATA TQ_TeO/                                                      071215
     1  0.699999789529, 0.847899957761, 1.087300099691, 1.741399809366, TeO  
     2  1.958200010031, 2.151899667160, 2.300000195710, 2.441799624595, TeO  
     3  2.532999787104, 2.628200062137, 2.805300316137, 2.975999833059, TeO  
     4  3.121499908610, 3.283999835941, 3.440899604168, 3.591900201718, TeO  
     5  3.707099969361, 3.816499937650, 3.927000256680, 3.971400170248, TeO  
     6  4.000000000000,      6*0.0D+00/                                 TeO  
      DATA  Q_TeO/                                                      071215
     1  1.00752051D+00, 1.15114916D+00, 1.38610121D+00, 2.03555263D+00, TeO  
     2  2.25188387D+00, 2.44593343D+00, 2.59806653D+00, 2.75489369D+00, TeO  
     3  2.86639669D+00, 2.99436932D+00, 3.26519159D+00, 3.55989301D+00, TeO  
     4  3.82898867D+00, 4.14219645D+00, 4.45357312D+00, 4.76042027D+00, TeO  
     5  5.00018270D+00, 5.23530914D+00, 5.48565869D+00, 5.59133706D+00, TeO  
     6  5.66115548D+00,      6*0.0D+00/                                 TeO  
      DATA TQ_IO/                                                       071215
     1  0.699999789529, 0.758200199698, 0.849199926395, 1.087500104610, IO   
     2  1.443399963124, 1.838800120128, 2.076299841172, 2.284099827092, IO   
     3  2.507300175938, 2.761100259963, 2.930300271136, 3.076699862308, IO   
     4  3.205099912499, 3.443999668698, 3.632600175638, 3.804100277582, IO   
     5  3.919100078035, 3.968500237227, 4.000000000000,      8*0.0D+00/ IO   
      DATA  Q_IO/                                                       071215
     1  1.66967385D+00, 1.72096577D+00, 1.80282653D+00, 2.02448147D+00, IO   
     2  2.36758676D+00, 2.75700370D+00, 2.99300666D+00, 3.20277989D+00, IO   
     3  3.44495954D+00, 3.76940046D+00, 4.02789565D+00, 4.28522561D+00, IO   
     4  4.53628225D+00, 5.04718532D+00, 5.47009232D+00, 5.86643892D+00, IO   
     5  6.14294362D+00, 6.26462171D+00, 6.34291053D+00,      8*0.0D+00/ IO   
      DATA TQ_BaO/                                                      071215
     1  0.699999789529, 0.845400018081, 1.080399929986, 1.789299989964, BaO  
     2  1.987899821755, 2.173400141711, 2.525200070587, 2.677200244929, BaO  
     3  2.830799960421, 3.135000235073, 3.308700380533, 3.397299824695, BaO  
     4  3.487299696807, 3.572199739102, 3.655099740765, 3.739099746150, BaO  
     5  3.825499817937, 3.920500112268, 4.000000000000,      8*0.0D+00/ BaO  
      DATA  Q_BaO/                                                      071215
     1  1.06098107D+00, 1.20267713D+00, 1.43379639D+00, 2.13843038D+00, BaO  
     2  2.33672613D+00, 2.52276661D+00, 2.89974596D+00, 3.08875581D+00, BaO  
     3  3.30204930D+00, 3.78649969D+00, 4.09298906D+00, 4.25617504D+00, BaO  
     4  4.42752079D+00, 4.59814652D+00, 4.78160804D+00, 4.99690746D+00, BaO  
     5  5.25892196D+00, 5.58988767D+00, 5.88588688D+00,      8*0.0D+00/ BaO  
      DATA TQ_LaO/                                                      071215
     1  0.699999789529, 0.764700122380, 0.865799983576, 1.133900026194, LaO  
     2  1.871699926078, 2.054100292742, 2.236599707381, 2.380500137323, LaO  
     3  2.586300073048, 2.717000213796, 2.859999642519, 3.049800176621, LaO  
     4  3.235199658004, 3.358999709601, 3.579999925446, 3.698499776669, LaO  
     5  3.821399718815, 3.935099955755, 3.974399954687, 4.000000000000, LaO  
     6       7*0.0D+00/                                                 LaO  
      DATA  Q_LaO/                                                      071215
     1  1.31130115D+00, 1.37396517D+00, 1.47243614D+00, 1.73592909D+00, LaO  
     2  2.46939969D+00, 2.65152655D+00, 2.83437742D+00, 2.98118412D+00, LaO  
     3  3.20555322D+00, 3.36404719D+00, 3.55568556D+00, 3.83981562D+00, LaO  
     4  4.14740535D+00, 4.37058429D+00, 4.82524699D+00, 5.11085541D+00, LaO  
     5  5.43944593D+00, 5.76701610D+00, 5.88414276D+00, 5.96128800D+00, LaO  
     6       7*0.0D+00/                                                 LaO  
      DATA TQ_TbO/                                                      071215
     1  0.699999789529, 0.852799971322, 1.100399919539, 1.817500130117, TbO  
     2  2.046600120490, 2.256500153511, 2.397799768161, 2.528999788399, TbO  
     3  2.742399825460, 2.906899761467, 3.054800305190, 3.418800061349, TbO  
     4  3.655299745432, 3.829899924311, 3.935299941340, 4.000000000000, TbO  
     5      11*0.0D+00/                                                 TbO  
      DATA  Q_TbO/                                                      071215
     1  1.00895561D+00, 1.15737872D+00, 1.40047767D+00, 2.11285468D+00, TbO  
     2  2.34150592D+00, 2.55175164D+00, 2.69612007D+00, 2.83665944D+00, TbO  
     3  3.09077001D+00, 3.31475524D+00, 3.53722850D+00, 4.15227605D+00, TbO  
     4  4.58648579D+00, 4.91751363D+00, 5.12031490D+00, 5.24566267D+00, TbO  
     5      11*0.0D+00/                                                 TbO  
      DATA TQ_LuO/                                                      071215
     1  0.699999789529, 0.761300212059, 0.857000067668, 1.111100112292, LuO  
     2  1.470699796609, 1.843200071162, 2.057600373345, 2.261400265580, LuO  
     3  2.413399932870, 2.559200394061, 2.674300179264, 2.794600082224, LuO  
     4  3.029899752343, 3.355199996518, 3.535899864052, 3.625700026395, LuO  
     5  3.709200017174, 3.881200147665, 3.953499904804, 3.982299603942, LuO  
     6  4.000000000000,      6*0.0D+00/                                 LuO  
      DATA  Q_LuO/                                                      071215
     1  1.30496005D+00, 1.36429396D+00, 1.45743173D+00, 1.70694005D+00, LuO  
     2  2.06330935D+00, 2.43440723D+00, 2.64845666D+00, 2.85268275D+00, LuO  
     3  3.00827079D+00, 3.16597796D+00, 3.30072463D+00, 3.45409127D+00, LuO  
     4  3.79399832D+00, 4.33971324D+00, 4.67191651D+00, 4.84447500D+00, LuO  
     5  5.01153913D+00, 5.39284126D+00, 5.57563829D+00, 5.65274050D+00, LuO  
     6  5.70131071D+00,      6*0.0D+00/                                 LuO  
      DATA TQ_HfO/                                                      071215
     1  0.699999789529, 0.762200188320, 0.859500125016, 1.115000011059, HfO  
     2  1.492999845344, 1.885399986527, 2.107099884583, 2.317299904488, HfO  
     3  2.471200164688, 2.621599910852, 2.736199692220, 2.855100006081, HfO  
     4  3.080499954662, 3.368799867430, 3.540799976491, 3.629900119814, HfO  
     5  3.707799985299, 3.796900122715, 3.879200101084, 3.953099896070, HfO  
     6  3.982099599452, 4.000000000000,      5*0.0D+00/                 HfO  
      DATA  Q_HfO/                                                      071215
     1  9.71881514D-01, 1.03192861D+00, 1.12642559D+00, 1.37697395D+00, HfO  
     2  1.75140170D+00, 2.14231009D+00, 2.36365053D+00, 2.57421494D+00, HfO  
     3  2.73152235D+00, 2.89395296D+00, 3.02799371D+00, 3.17931520D+00, HfO  
     4  3.50319949D+00, 3.98020249D+00, 4.29181442D+00, 4.46233535D+00, HfO  
     5  4.61984809D+00, 4.81551455D+00, 5.01737945D+00, 5.21958767D+00, HfO  
     6  5.30447345D+00, 5.35835758D+00,      5*0.0D+00/                 HfO  
      DATA TQ_TaO/                                                      071215
     1  0.699999789529, 0.728399997999, 0.772300038678, 0.884799996669, TaO  
     2  1.179999835913, 1.556400025690, 1.976599934373, 2.202599858785, TaO  
     3  2.397299804914, 2.648699710015, 2.908699801944, 3.164099931043, TaO  
     4  3.396199905011, 3.499099969069, 3.609599685942, 3.723200029351, TaO  
     5  3.832199974647, 3.935499926925, 3.974799925945, 4.000000000000, TaO  
     6       7*0.0D+00/                                                 TaO  
      DATA  Q_TaO/                                                      071215
     1  1.65294388D+00, 1.67442794D+00, 1.70843196D+00, 1.79947878D+00, TaO  
     2  2.05823191D+00, 2.41255467D+00, 2.82278184D+00, 3.04638203D+00, TaO  
     3  3.24100672D+00, 3.50675112D+00, 3.82883331D+00, 4.21875512D+00, TaO  
     4  4.64991129D+00, 4.86426290D+00, 5.11149591D+00, 5.38778840D+00, TaO  
     5  5.67767861D+00, 5.97540604D+00, 6.09395991D+00, 6.17129714D+00, TaO  
     6       7*0.0D+00/                                                 TaO  
      DATA TQ_WO/                                                       071215
     1  0.699999789529, 0.758600210065, 0.850599920856, 1.086300075096, WO   
     2  1.471399813000, 1.881100101898, 2.118799856866, 2.325299825396, WO   
     3  2.484999654277, 2.659299843744, 2.829099920686, 2.986999698536, WO   
     4  3.141800215343, 3.394899999931, 3.626600046413, 3.742999831475, WO   
     5  3.854900024040, 3.942499659039, 3.977699717569, 4.000000000000, WO   
     6       7*0.0D+00/                                                 WO   
      DATA  Q_WO/                                                       071215
     1  9.53518855D-01, 1.00998934D+00, 1.09917774D+00, 1.32989688D+00, WO   
     2  1.71094137D+00, 2.11890877D+00, 2.35614737D+00, 2.56273628D+00, WO   
     3  2.72502824D+00, 2.91215880D+00, 3.11453355D+00, 3.32628170D+00, WO   
     4  3.55644393D+00, 3.97397829D+00, 4.39027724D+00, 4.60993875D+00, WO   
     5  4.82950082D+00, 5.00909961D+00, 5.08356443D+00, 5.13147476D+00, WO   
     6       7*0.0D+00/                                                 WO   
      DATA TQ_PtO/                                                      071215
     1  0.699999789529, 0.761800198871, 0.858300097489, 1.115799990293, PtO  
     2  1.469699787598, 1.846100001191, 2.066999906850, 2.272200297356, PtO  
     3  2.561800284426, 2.671900124921, 2.797800154691, 3.022100303850, PtO  
     4  3.342400218195, 3.518400419374, 3.606099939239, 3.685499914191, PtO  
     5  3.871699931419, 3.949499817090, 3.980699568021, 4.000000000000, PtO  
     6       7*0.0D+00/                                                 PtO  
      DATA  Q_PtO/                                                      071215
     1  9.77175813D-01, 1.03686313D+00, 1.13061371D+00, 1.38318464D+00, PtO  
     2  1.73371249D+00, 2.10860398D+00, 2.32912580D+00, 2.53482704D+00, PtO  
     3  2.83991444D+00, 2.96867437D+00, 3.12913874D+00, 3.45298674D+00, PtO  
     4  3.99050564D+00, 4.31607536D+00, 4.48614612D+00, 4.64602910D+00, PtO  
     5  5.05702616D+00, 5.25011259D+00, 5.33143359D+00, 5.38277749D+00, PtO  
     6       7*0.0D+00/                                                 PtO  
      DATA TQ_PbO/                                                      071215
     1  0.699999789529, 0.846799984302, 1.083600008690, 1.802400121792, PbO  
     2  2.007700176810, 2.199199773820, 2.545400078447, 2.691899626942, PbO  
     3  2.841500200072, 3.147799767624, 3.312600224232, 3.398899707870, PbO  
     4  3.486299674630, 3.569899693730, 3.651499656764, 3.735999682107, PbO  
     5  3.832299976826, 3.933100099903, 3.974699933131, 4.000000000000, PbO  
     6       7*0.0D+00/                                                 PbO  
      DATA  Q_PbO/                                                      071215
     1  1.06857052D+00, 1.21170317D+00, 1.44468144D+00, 2.15929132D+00, PbO  
     2  2.36428554D+00, 2.55628683D+00, 2.92601414D+00, 3.10656721D+00, PbO  
     3  3.31200867D+00, 3.79640958D+00, 4.08647404D+00, 4.24538239D+00, PbO  
     4  4.41215036D+00, 4.58050082D+00, 4.75988433D+00, 4.97138890D+00, PbO  
     5  5.25477069D+00, 5.59691393D+00, 5.74710577D+00, 5.83976234D+00, PbO  
     6       7*0.0D+00/                                                 PbO  
      DATA TQ_BiO/                                                      071215
     1  0.699999789529, 0.771300014290, 0.881800070646, 1.181299862980, BiO  
     2  1.491299809686, 1.843000075988, 2.016900392026, 2.182300141655, BiO  
     3  2.514700345992, 2.642600151964, 2.780699742428, 3.025200084661, BiO  
     4  3.229799546223, 3.417200019248, 3.600200366225, 3.721300161009, BiO  
     5  3.831599961579, 3.933800049451, 3.974099976243, 4.000000000000, BiO  
     6       7*0.0D+00/                                                 BiO  
      DATA  Q_BiO/                                                      071215
     1  1.71336996D+00, 1.77721018D+00, 1.87828314D+00, 2.16148830D+00, BiO  
     2  2.46308921D+00, 2.81034686D+00, 2.98313956D+00, 3.14849099D+00, BiO  
     3  3.50180883D+00, 3.65689298D+00, 3.84168949D+00, 4.21351582D+00, BiO  
     4  4.56339817D+00, 4.91458651D+00, 5.29868033D+00, 5.58478844D+00, BiO  
     5  5.87140756D+00, 6.15834404D+00, 6.27615423D+00, 6.35291660D+00, BiO  
     6       7*0.0D+00/                                                 BiO  
      DATA TQ_ThO/                                                      071215
     1  0.699999789529, 0.788499993619, 0.924599974709, 1.319900160566, ThO  
     2  1.949499814910, 2.108199806179, 2.265500355138, 2.417200025023, ThO  
     3  2.660699876369, 2.805500320349, 2.996399915747, 3.130400129654, ThO  
     4  3.286699898007, 3.425100195070, 3.576899851386, 3.651399654431, ThO  
     5  3.724599932340, 3.889100335181, 3.954699931007, 4.000000000000, ThO  
     6       7*0.0D+00/                                                 ThO  
      DATA  Q_ThO/                                                      071215
     1  1.03472731D+00, 1.12066990D+00, 1.25373073D+00, 1.64411471D+00, ThO  
     2  2.27124700D+00, 2.42976834D+00, 2.58735711D+00, 2.74185044D+00, ThO  
     3  3.00957691D+00, 3.19025668D+00, 3.45925025D+00, 3.67038213D+00, ThO  
     4  3.94592722D+00, 4.22772914D+00, 4.59742667D+00, 4.80858278D+00, ThO  
     5  5.03592665D+00, 5.60027193D+00, 5.83520612D+00, 5.99773881D+00, ThO  
     6       7*0.0D+00/                                                 ThO  
      DATA TQ_BOp/                                                      071215
     1  0.699999789529, 0.730900011251, 0.778800197201, 0.901900048465, BOp  
     2  1.232999936416, 1.601499928634, 2.044000060971, 2.305700328467, BOp  
     3  2.528299840381, 2.718500247596, 2.957899992541, 3.092900209041, BOp  
     4  3.227499713105, 3.513100295602, 3.666699991867, 3.844000234527, BOp  
     5  3.937099811607, 3.975399882833, 4.000000000000,      8*0.0D+00/ BOp  
      DATA  Q_BOp/                                                      071215
     1  3.67053731D-01, 3.92664596D-01, 4.33106444D-01, 5.40534669D-01, BOp  
     2  8.46427998D-01, 1.20246684D+00, 1.63911466D+00, 1.89942156D+00, BOp  
     3  2.12168034D+00, 2.31472295D+00, 2.57749156D+00, 2.74470472D+00, BOp  
     4  2.92849512D+00, 3.37241615D+00, 3.63626841D+00, 3.95824755D+00, BOp  
     5  4.13491977D+00, 4.20935694D+00, 4.25776090D+00,      8*0.0D+00/ BOp  
      DATA TQ_SiOp/                                                     071215
     1  0.699999789529, 0.757100171190, 0.845700010842, 1.090500153308, SiOp 
     2  1.478299974570, 2.001300029851, 2.495599894011, 2.698299769770, SiOp 
     3  2.971200183287, 3.094800254822, 3.215800158597, 3.495099875714, SiOp 
     4  3.835200039991, 3.934499998999, 4.000000000000,     12*0.0D+00/ SiOp 
      DATA  Q_SiOp/                                                     071215
     1  1.01708940D+00, 1.07046438D+00, 1.15418287D+00, 1.38972618D+00, SiOp 
     2  1.77035379D+00, 2.28996804D+00, 2.78367230D+00, 2.98886557D+00, SiOp 
     3  3.28849092D+00, 3.44250033D+00, 3.60722393D+00, 4.03806927D+00, SiOp 
     4  4.63661015D+00, 4.82232221D+00, 4.94701459D+00,     12*0.0D+00/ SiOp 
      DATA TQ_POp/                                                      071215
     1  0.699999789529, 0.759800241165, 0.853799994262, 1.088400126746, POp  
     2  1.464399917979, 1.950199808146, 2.200699808407, 2.420600104351, POp  
     3  2.603100147722, 2.746399910577, 2.899899613737, 3.030499757083, POp  
     4  3.162399889651, 3.454999914504, 3.828699895300, 3.932200164770, POp  
     5  4.000000000000,     10*0.0D+00/                                 POp  
      DATA  Q_POp/                                                      071215
     1  1.01697027D+00, 1.07287903D+00, 1.16176020D+00, 1.38755092D+00, POp  
     2  1.75645438D+00, 2.23879289D+00, 2.48855874D+00, 2.70829628D+00, POp  
     3  2.89330857D+00, 3.04563797D+00, 3.22321444D+00, 3.39036747D+00, POp  
     4  3.57563282D+00, 4.04050477D+00, 4.71310446D+00, 4.91024317D+00, POp  
     5  5.04159290D+00,     10*0.0D+00/                                 POp  
      DATA TQ_SOp/                                                      071215
     1  0.699999789529, 0.748599955070, 0.824600049820, 1.016400060534, SOp  
     2  1.336000072494, 1.701999861377, 1.924599960377, 2.153799712094, SOp  
     3  2.329899925002, 2.537999908713, 2.666900013389, 2.864599757722, SOp  
     4  3.029799759414, 3.224899901756, 3.396099912313, 3.642400158648, SOp  
     5  3.737999723425, 3.826699846948, 3.933100099903, 3.973700004984, SOp  
     6  4.000000000000,      6*0.0D+00/                                 SOp  
      DATA  Q_SOp/                                                      071215
     1  1.10323261D+00, 1.13982025D+00, 1.19941430D+00, 1.36023061D+00, SOp  
     2  1.65035539D+00, 2.00079332D+00, 2.21990885D+00, 2.45871200D+00, SOp  
     3  2.66160661D+00, 2.92239677D+00, 3.09151507D+00, 3.36315814D+00, SOp  
     4  3.60653849D+00, 3.91671095D+00, 4.20783699D+00, 4.65134016D+00, SOp  
     5  4.83056534D+00, 5.00218217D+00, 5.21997524D+00, 5.30819910D+00, SOp  
     6  5.36722725D+00,      6*0.0D+00/                                 SOp  
      DATA TQ_AsOp/                                                     071215
     1  0.699999789529, 0.758700212656, 0.850699923150, 1.085100045582, AsOp 
     2  1.474199878564, 1.891699900290, 2.132000136980, 2.338900138138, AsOp 
     3  2.497399929880, 2.668300044329, 2.837700115752, 3.024100162438, AsOp 
     4  3.187599737901, 3.551600229932, 3.710900054801, 3.851800242600, AsOp 
     5  3.941699640976, 3.977399739125, 4.000000000000,      8*0.0D+00/ AsOp 
      DATA  Q_AsOp/                                                     071215
     1  8.49043237D-01, 9.04994123D-01, 9.93373067D-01, 1.22139069D+00, AsOp 
     2  1.60523850D+00, 2.02054045D+00, 2.26030669D+00, 2.46732302D+00, AsOp 
     3  2.62851478D+00, 2.81212426D+00, 3.01410700D+00, 3.26708479D+00, AsOp 
     4  3.51600738D+00, 4.14465905D+00, 4.44502262D+00, 4.72225404D+00, AsOp 
     5  4.90650471D+00, 4.98203174D+00, 5.03074907D+00,      8*0.0D+00/ AsOp 
      DATA TQ_TaOp/                                                     071215
     1  0.699999789529, 0.759100223023, 0.851699946089, 1.089200146421, TaOp 
     2  1.472699843440, 1.883000050920, 2.114799775509, 2.317399897308, TaOp 
     3  2.477699693030, 2.658399821688, 2.834600045966, 2.968600238655, TaOp 
     4  3.096600298193, 3.370599911562, 3.509000199894, 3.636100247791, TaOp 
     5  3.850600327204, 3.941899645492, 3.977399739125, 4.000000000000, TaOp 
     6       7*0.0D+00/                                                 TaOp 
      DATA  Q_TaOp/                                                     071215
     1  9.44343509D-01, 1.00124956D+00, 1.09096279D+00, 1.32335912D+00, TaOp 
     2  1.70276705D+00, 2.11134140D+00, 2.34272894D+00, 2.54546339D+00, TaOp 
     3  2.70834913D+00, 2.90256773D+00, 3.11388768D+00, 3.29373865D+00, TaOp 
     4  3.48153197D+00, 3.92985545D+00, 4.17595267D+00, 4.41206268D+00, TaOp 
     5  4.83734755D+00, 5.03206646D+00, 5.11036627D+00, 5.16098908D+00, TaOp 
     6       7*0.0D+00/                                                 TaOp 
      DATA TQ_FeOm/                                                     071215
     1  0.699999789529, 0.758400204881, 0.849499919156, 1.090600150751, FeOm 
     2  1.436299958586, 1.794600067244, 2.012500288465, 2.214900150206, FeOm 
     3  2.357699816666, 2.488499743976, 2.709100033710, 2.878700077613, FeOm 
     4  3.028699837190, 3.398699722473, 3.642500151584, 3.824799801014, FeOm 
     5  3.933100099903, 4.000000000000,      9*0.0D+00/                 FeOm 
      DATA  Q_FeOm/                                                     071215
     1  9.43292827D-01, 9.99516985D-01, 1.08775916D+00, 1.32365355D+00, FeOm 
     2  1.66551186D+00, 2.02204058D+00, 2.23940988D+00, 2.44214356D+00, FeOm 
     3  2.58821745D+00, 2.72874739D+00, 2.99326767D+00, 3.22675174D+00, FeOm 
     4  3.45490102D+00, 4.08532175D+00, 4.53545332D+00, 4.88259854D+00, FeOm 
     5  5.09173694D+00, 5.22175125D+00,      9*0.0D+00/                 FeOm 
      DATA TQ_LiF/                                                      071215
     1  0.699999789529, 0.771300014290, 0.881800070646, 1.183599910869, LiF  
     2  1.519399982082, 1.890999884979, 2.109099742030, 2.305100314492, LiF  
     3  2.583200000418, 2.690099586771, 2.812500236599, 2.966100181175, LiF  
     4  3.273600202589, 3.425800209090, 3.573599772548, 3.656999785099, LiF  
     5  3.741399795885, 3.821599723650, 3.911699887902, 3.964100138620, LiF  
     6  4.000000000000,      6*0.0D+00/                                 LiF  
      DATA  Q_LiF/                                                      071215
     1  4.72783440D-01, 5.35455230D-01, 6.35154873D-01, 9.18450043D-01, LiF  
     2  1.24444750D+00, 1.61137794D+00, 1.82828838D+00, 2.02452302D+00, LiF  
     3  2.31748977D+00, 2.44219398D+00, 2.59759079D+00, 2.81362917D+00, LiF  
     4  3.31384385D+00, 3.59232723D+00, 3.88364654D+00, 4.05920808D+00, LiF  
     5  4.24658472D+00, 4.43319309D+00, 4.64889806D+00, 4.77457608D+00, LiF  
     6  4.85975237D+00,      6*0.0D+00/                                 LiF  
      DATA TQ_BeF/                                                      071215
     1  0.699999789529, 0.729800031527, 0.776100131353, 0.894299964685, BeF  
     2  1.216399877143, 1.571800127260, 1.957499992366, 2.214100134662, BeF  
     3  2.420200096734, 2.560700362716, 2.724999919775, 2.879400092470, BeF  
     4  3.123599960956, 3.380900127267, 3.623799984133, 3.743099833699, BeF  
     5  3.860499675897, 3.943999692907, 3.978399667271, 4.000000000000, BeF  
     6       7*0.0D+00/                                                 BeF  
      DATA  Q_BeF/                                                      071215
     1  7.35244645D-01, 7.60832510D-01, 8.01154360D-01, 9.06763203D-01, BeF  
     2  1.20798994D+00, 1.55290278D+00, 1.93384784D+00, 2.18910551D+00, BeF  
     3  2.39524324D+00, 2.53853276D+00, 2.71557721D+00, 2.89928707D+00, BeF  
     4  3.23610784D+00, 3.65228556D+00, 4.09425403D+00, 4.32798463D+00, BeF  
     5  4.57184904D+00, 4.75728022D+00, 4.83740137D+00, 4.88894429D+00, BeF  
     6       7*0.0D+00/                                                 BeF  
      DATA TQ_BF/                                                       071215
     1  0.699999789529, 0.729800031527, 0.776100131353, 0.894399966923, BF   
     2  1.214499917813, 1.582199972065, 1.985199884767, 2.248999980941, BF   
     3  2.460400035242, 2.604700030912, 2.787199907113, 2.919900073220, BF   
     4  3.045900098967, 3.167900023568, 3.465400152241, 3.574599796438, BF   
     5  3.680800241535, 3.780499754078, 3.878600087511, 3.951699865500, BF   
     6  3.981699590472, 4.000000000000,      5*0.0D+00/                 BF   
      DATA  Q_BF/                                                       071215
     1  4.27501528D-01, 4.53009140D-01, 4.93217560D-01, 5.98679646D-01, BF   
     2  8.97609321D-01, 1.25434006D+00, 1.65247753D+00, 1.91496700D+00, BF   
     3  2.12636084D+00, 2.27332431D+00, 2.47044959D+00, 2.62850507D+00, BF   
     4  2.79379317D+00, 2.96896179D+00, 3.45451398D+00, 3.65113527D+00, BF   
     5  3.85139731D+00, 4.04933784D+00, 4.25887511D+00, 4.43024764D+00, BF   
     6  4.50563433D+00, 4.55328129D+00,      5*0.0D+00/                 BF   
      DATA TQ_NaF/                                                      071215
     1  0.699999789529, 0.757100171190, 0.846199998778, 1.083099996392, NaF  
     2  1.735599886227, 1.908699863763, 2.081899972923, 2.395499937225, NaF  
     3  2.521800323071, 2.656799782477, 2.820799720329, 3.171800113285, NaF  
     4  3.362499700384, 3.566299952459, 3.740399773641, 3.839200127116, NaF  
     5  3.917200029217, 3.968500237227, 4.000000000000,      8*0.0D+00/ NaF  
      DATA  Q_NaF/                                                      071215
     1  9.22049693D-01, 9.76909193D-01, 1.06305616D+00, 1.29452547D+00, NaF  
     2  1.94135851D+00, 2.11404676D+00, 2.28779894D+00, 2.62182837D+00, NaF  
     3  2.77441054D+00, 2.95415336D+00, 3.19653941D+00, 3.79297828D+00, NaF  
     4  4.15336676D+00, 4.56491816D+00, 4.94148434D+00, 5.16566946D+00, NaF  
     5  5.34574194D+00, 5.46402672D+00, 5.53603049D+00,      8*0.0D+00/ NaF  
      DATA TQ_MgF/                                                      071215
     1  0.699999789529, 0.759400230798, 0.852199957559, 1.096799992167, MgF  
     2  1.436999974477, 1.808899968874, 2.009900227328, 2.201699834922, MgF  
     3  2.497399929880, 2.612099692321, 2.744299865891, 2.938299709720, MgF  
     4  3.280699760082, 3.463800117625, 3.636000245730, 3.748199947143, MgF  
     5  3.847800320014, 3.941399634203, 3.977299746311, 4.000000000000, MgF  
     6       7*0.0D+00/                                                 MgF  
      DATA  Q_MgF/                                                      071215
     1  1.15126032D+00, 1.20788982D+00, 1.29705913D+00, 1.53512022D+00, MgF  
     2  1.87066394D+00, 2.24038739D+00, 2.44088325D+00, 2.63319284D+00, MgF  
     3  2.94622384D+00, 3.08204907D+00, 3.25356026D+00, 3.53621843D+00, MgF  
     4  4.11510541D+00, 4.45929849D+00, 4.80530178D+00, 5.04631415D+00, MgF  
     5  5.27368371D+00, 5.49953947D+00, 5.58900960D+00, 5.64628688D+00, MgF  
     6       7*0.0D+00/                                                 MgF  
      DATA TQ_AlF/                                                      071215
     1  0.699999789529, 0.759500233390, 0.852699969029, 1.095300030534, AlF  
     2  1.453599887122, 1.836900075409, 2.043600051814, 2.240599802577, AlF  
     3  2.392700143042, 2.537599898984, 2.654899735913, 2.776899957271, AlF  
     4  2.995399892344, 3.349000366648, 3.532099777575, 3.617799846235, AlF  
     5  3.703199880566, 3.796900122715, 3.884300221248, 3.954299922273, AlF  
     6  3.982699612923, 4.000000000000,      5*0.0D+00/                 AlF  
      DATA  Q_AlF/                                                      071215
     1  8.24647925D-01, 8.81192929D-01, 9.70512808D-01, 1.20620003D+00, AlF  
     2  1.55934990D+00, 1.94036796D+00, 2.14653818D+00, 2.34390343D+00, AlF  
     3  2.49969153D+00, 2.65668681D+00, 2.79446516D+00, 2.95090396D+00, AlF  
     4  3.26730315D+00, 3.86646663D+00, 4.21055722D+00, 4.37895099D+00, AlF  
     5  4.55246234D+00, 4.75248163D+00, 4.95345796D+00, 5.12874786D+00, AlF  
     6  5.20419114D+00, 5.25145581D+00,      5*0.0D+00/                 AlF  
      DATA TQ_SiF/                                                      071215
     1  0.699999789529, 0.726999964471, 0.768800014237, 0.875700010893, SiF  
     2  1.157699900341, 1.507500168520, 1.898400046838, 2.110399686017, SiF  
     3  2.296400106058, 2.537699901416, 2.761200262379, 3.022200296779, SiF  
     4  3.297400133115, 3.476399792979, 3.649099685342, 3.757400162560, SiF  
     5  3.863199737589, 3.945799733549, 3.979099616974, 4.000000000000, SiF  
     6       7*0.0D+00/                                                 SiF  
      DATA  Q_SiF/                                                      071215
     1  1.47585120D+00, 1.49728609D+00, 1.53109784D+00, 1.62064151D+00, SiF  
     2  1.87262086D+00, 2.20409167D+00, 2.58622660D+00, 2.79600814D+00, SiF  
     3  2.98172988D+00, 3.23432052D+00, 3.50138257D+00, 3.87288721D+00, SiF  
     4  4.33015249D+00, 4.65653260D+00, 4.99014869D+00, 5.20979028D+00, SiF  
     5  5.43538278D+00, 5.62242806D+00, 5.70126998D+00, 5.75187264D+00, SiF  
     6       7*0.0D+00/                                                 SiF  
      DATA TQ_PF/                                                       071215
     1  0.699999789529, 0.760800225247, 0.856000044729, 1.102899977171, PF   
     2  1.464799908139, 1.856500056198, 2.050700214442, 2.247899957584, PF   
     3  2.407199782107, 2.569399743516, 2.685399917483, 2.807200356154, PF   
     4  2.962600100702, 3.151899646734, 3.329699908264, 3.603500127402, PF   
     5  3.726799779893, 3.848100326763, 3.938999674667, 3.976199825350, PF   
     6  4.000000000000,      6*0.0D+00/                                 PF   
      DATA  Q_PF/                                                       071215
     1  1.29108565D+00, 1.34878866D+00, 1.43993843D+00, 1.67969243D+00, PF   
     2  2.03633110D+00, 2.42563243D+00, 2.61925194D+00, 2.81654244D+00, PF   
     3  2.97907841D+00, 3.15457256D+00, 3.29124567D+00, 3.44778568D+00, PF   
     4  3.66833439D+00, 3.96691851D+00, 4.27458014D+00, 4.80350244D+00, PF   
     5  5.06699118D+00, 5.34441091D+00, 5.56740223D+00, 5.66316743D+00, PF   
     6  5.72590160D+00,      6*0.0D+00/                                 PF   
      DATA TQ_SF/                                                       071215
     1  0.699999789529, 0.747399927545, 0.821600122425, 1.010400203878, SF   
     2  1.304000085678, 1.624499989564, 1.820100189696, 2.056800354922, SF   
     3  2.476499780105, 2.653299696702, 2.835900075231, 3.258900196424, SF   
     4  3.556400330105, 3.779699763999, 3.911399880194, 3.965400167754, SF   
     5  4.000000000000,     10*0.0D+00/                                 SF   
      DATA  Q_SF/                                                       071215
     1  1.49152800D+00, 1.52977920D+00, 1.59145405D+00, 1.75635029D+00, SF   
     2  2.02815920D+00, 2.33684472D+00, 2.52966655D+00, 2.77638235D+00, SF   
     3  3.29687927D+00, 3.55377102D+00, 3.84426824D+00, 4.59584611D+00, SF   
     4  5.16624874D+00, 5.61015165D+00, 5.87944691D+00, 5.99175239D+00, SF   
     5  6.06424857D+00,     10*0.0D+00/                                 SF   
      DATA TQ_KF/                                                       071215
     1  0.699999789529, 0.771000006974, 0.879800110189, 1.192699978292, KF   
     2  1.681799879715, 1.962799982554, 2.321899751774, 2.465000138991, KF   
     3  2.610199648584, 3.108799744043, 3.325899825871, 3.543600043463, KF   
     4  3.625400019722, 3.715400151868, 3.801900229357, 3.895199973359, KF   
     5  3.958200007432, 3.984099644354, 4.000000000000,      8*0.0D+00/ KF   
      DATA  Q_KF/                                                       071215
     1  1.10839860D+00, 1.17764602D+00, 1.28426333D+00, 1.59325071D+00, KF   
     2  2.07995370D+00, 2.36111802D+00, 2.74453446D+00, 2.92132434D+00, KF   
     3  3.12078034D+00, 3.94722389D+00, 4.36117051D+00, 4.80958977D+00, KF   
     4  4.98997455D+00, 5.19721092D+00, 5.40347758D+00, 5.62880113D+00, KF   
     5  5.77897265D+00, 5.83959616D+00, 5.87639343D+00,      8*0.0D+00/ KF   
      DATA TQ_CaF/                                                      071215
     1  0.699999789529, 0.854400008025, 1.108500106266, 1.767300041589, CaF  
     2  1.943699951922, 2.115099781611, 2.437299773874, 2.573399777188, CaF  
     3  2.712200105635, 3.279399786984, 3.478199663946, 3.587400098819, CaF  
     4  3.681900164923, 3.775900038987, 3.864199760438, 3.946499749354, CaF  
     5  3.979299602603, 4.000000000000,      9*0.0D+00/                 CaF  
      DATA  Q_CaF/                                                      071215
     1  1.32922322D+00, 1.47941212D+00, 1.72916017D+00, 2.38379254D+00, CaF  
     2  2.55989673D+00, 2.73186484D+00, 3.07517295D+00, 3.24057216D+00, CaF  
     3  3.42722960D+00, 4.36489183D+00, 4.74495250D+00, 4.96535463D+00, CaF  
     4  5.16680678D+00, 5.38288716D+00, 5.60579985D+00, 5.83348488D+00, CaF  
     5  5.92925520D+00, 5.99095388D+00,      9*0.0D+00/                 CaF  
      DATA TQ_ScF/                                                      071215
     1  0.699999789529, 0.781100196333, 0.906299947040, 1.272100032354, ScF  
     2  1.845600013255, 2.028899807929, 2.207199980754, 2.541699995478, ScF  
     3  2.677300247193, 2.823099775849, 3.060100432703, 3.275100095105, ScF  
     4  3.445999710331, 3.635200229237, 3.739199748216, 3.840300151289, ScF  
     5  3.937799761155, 3.975799854091, 4.000000000000,      8*0.0D+00/ ScF  
      DATA  Q_ScF/                                                      071215
     1  9.63300754D-01, 1.04159436D+00, 1.16337547D+00, 1.52336401D+00, ScF  
     2  2.09364786D+00, 2.27656874D+00, 2.45524413D+00, 2.81097968D+00, ScF  
     3  2.97577061D+00, 3.17222008D+00, 3.53410637D+00, 3.90151389D+00, ScF  
     4  4.21888618D+00, 4.61609182D+00, 4.87003402D+00, 5.14797784D+00, ScF  
     5  5.44376450D+00, 5.56486531D+00, 5.64319428D+00,      8*0.0D+00/ ScF  
      DATA TQ_MnF/                                                      071215
     1  0.699999789529, 0.758400204881, 0.849499919156, 1.093800068901, MnF  
     2  1.759500226552, 1.957399989843, 2.146099894691, 2.433800018149, MnF  
     3  2.547600127781, 2.676000217757, 2.908499797447, 3.109399700674, MnF  
     4  3.309700402541, 3.511100248896, 3.704499910165, 3.882700183270, MnF  
     5  3.953999915722, 3.982499608432, 4.000000000000,      8*0.0D+00/ MnF  
      DATA  Q_MnF/                                                      071215
     1  1.80698973D+00, 1.86330946D+00, 1.95167697D+00, 2.19095000D+00, MnF  
     2  2.85144281D+00, 3.04884430D+00, 3.23799834D+00, 3.54186241D+00, MnF  
     3  3.67568355D+00, 3.84052439D+00, 4.17866724D+00, 4.50758086D+00, MnF  
     4  4.86186989D+00, 5.23843367D+00, 5.61688492D+00, 5.98558759D+00, MnF  
     5  6.14144019D+00, 6.20536192D+00, 6.24508992D+00,      8*0.0D+00/ MnF  
      DATA TQ_NiF/                                                      071215
     1  0.699999789529, 0.758500207473, 0.849899909505, 1.092100112384, NiF  
     2  1.437399983558, 1.794000054910, 2.012100279051, 2.214700146320, NiF  
     3  2.357599824101, 2.488599746539, 2.709300038677, 2.879000083980, NiF  
     4  3.028999815979, 3.398599729775, 3.642400158648, 3.824699798596, NiF  
     5  3.933100099903, 4.000000000000,      9*0.0D+00/                 NiF  
      DATA  Q_NiF/                                                      071215
     1  9.61891686D-01, 1.01830810D+00, 1.10696815D+00, 1.34417808D+00, NiF  
     2  1.68581973D+00, 2.04073515D+00, 2.25832699D+00, 2.46127191D+00, NiF  
     3  2.60744964D+00, 2.74819937D+00, 3.01287263D+00, 3.24653653D+00, NiF  
     4  3.47472535D+00, 4.10449299D+00, 4.55461707D+00, 4.90175870D+00, NiF  
     5  5.11108946D+00, 5.24110383D+00,      9*0.0D+00/                 NiF  
      DATA TQ_CuF/                                                      071215
     1  0.699999789529, 0.760400235797, 0.854600012613, 1.107300078603, CuF  
     2  1.785200091198, 1.964799930481, 2.144400013458, 2.461000048774, CuF  
     3  2.594000245457, 2.732399610050, 3.277599915965, 3.472200094056, CuF  
     4  3.568399801534, 3.665199961247, 3.758700193370, 3.852400200298, CuF  
     5  3.942299654524, 3.977699717569, 4.000000000000,      8*0.0D+00/ CuF  
      DATA  Q_CuF/                                                      071215
     1  9.80527907D-01, 1.03887729D+00, 1.13040083D+00, 1.37824313D+00, CuF  
     2  2.05139397D+00, 2.23066079D+00, 2.41082601D+00, 2.74780842D+00, CuF  
     3  2.90875514D+00, 3.09393491D+00, 3.99024825D+00, 4.36176248D+00, CuF  
     4  4.55577177D+00, 4.76163735D+00, 4.97579414D+00, 5.21041947D+00, CuF  
     5  5.45530471D+00, 5.55620730D+00, 5.62074518D+00,      8*0.0D+00/ CuF  
      DATA TQ_ZnF/                                                      071215
     1  0.699999789529, 0.846100001191, 1.079599930695, 1.761800190437, ZnF  
     2  1.939200022183, 2.124199969363, 2.281399772367, 2.439499620329, ZnF  
     3  2.556800338636, 2.679300292480, 2.914499941197, 3.244199870172, ZnF  
     4  3.441899624984, 3.651699661431, 3.752300041692, 3.850200355406, ZnF  
     5  3.941899645492, 4.000000000000,      9*0.0D+00/                 ZnF  
      DATA  Q_ZnF/                                                      071215
     1  1.28020950D+00, 1.42178127D+00, 1.65058241D+00, 2.32756733D+00, ZnF  
     2  2.50452366D+00, 2.68972193D+00, 2.85028875D+00, 3.02142777D+00, ZnF  
     3  3.15954516D+00, 3.31689110D+00, 3.65885931D+00, 4.21369858D+00, ZnF  
     4  4.57752977D+00, 4.98389465D+00, 5.18675462D+00, 5.39003976D+00, ZnF  
     5  5.58564160D+00, 5.71172404D+00,      9*0.0D+00/                 ZnF  
      DATA TQ_GaF/                                                      071215
     1  0.699999789529, 0.852799971322, 1.102799974865, 1.771900016223, GaF  
     2  1.957099982272, 2.140200306882, 2.452099849878, 2.571999745739, GaF  
     3  2.705099934369, 2.878900081858, 3.099200360841, 3.291700010450, GaF  
     4  3.493699843039, 3.689299649529, 3.787099909682, 3.881000142918, GaF  
     5  3.952699887336, 3.982099599452, 4.000000000000,      8*0.0D+00/ GaF  
      DATA  Q_GaF/                                                      071215
     1  1.00296659D+00, 1.15133121D+00, 1.39674546D+00, 2.06132223D+00, GaF  
     2  2.24618300D+00, 2.42980432D+00, 2.76068876D+00, 2.90403652D+00, GaF  
     3  3.07873730D+00, 3.33272563D+00, 3.69294869D+00, 4.03662844D+00, GaF  
     4  4.42149358D+00, 4.81869827D+00, 5.03014445D+00, 5.24641554D+00, GaF  
     5  5.42416209D+00, 5.50093139D+00, 5.54887572D+00,      8*0.0D+00/ GaF  
      DATA TQ_GeF/                                                      071215
     1  0.699999789529, 0.755200121949, 0.841600109767, 1.061999998723, GeF  
     2  1.404199991992, 1.786500059099, 2.003100071183, 2.195499699700, GeF  
     3  2.338200121544, 2.500700000318, 2.647299811446, 2.748799961647, GeF  
     4  2.842300217429, 3.096100286146, 3.305700314510, 3.450499805687, GeF  
     5  3.589000136306, 3.706199948870, 3.812500226428, 3.925100214468, GeF  
     6  3.970700220546, 4.000000000000,      5*0.0D+00/                 GeF  
      DATA  Q_GeF/                                                      071215
     1  1.34124243D+00, 1.38938604D+00, 1.46635869D+00, 1.66961702D+00, GeF  
     2  1.99760555D+00, 2.37299587D+00, 2.58776676D+00, 2.78028398D+00, GeF  
     3  2.92801417D+00, 3.11236283D+00, 3.30584958D+00, 3.45905412D+00, GeF  
     4  3.61402965D+00, 4.08621364D+00, 4.50907956D+00, 4.80991771D+00, GeF  
     5  5.10317316D+00, 5.35744211D+00, 5.59702338D+00, 5.86684661D+00, GeF  
     6  5.98230077D+00, 6.05849670D+00,      5*0.0D+00/                 GeF  
      DATA TQ_AsF/                                                      071215
     1  0.699999789529, 0.780400215509, 0.907899910158, 1.279000198158, AsF  
     2  1.398500058671, 1.511300193390, 1.784500108482, 1.942199987356, AsF  
     3  2.138100277184, 2.276999956823, 2.529699736417, 2.680200293854, AsF  
     4  2.880500117430, 3.087900098380, 3.256900152773, 3.448599764454, AsF  
     5  3.627800073104, 3.740899784763, 3.858299784329, 3.942099650008, AsF  
     6  3.977599724754, 4.000000000000,      5*0.0D+00/                 AsF  
      DATA  Q_AsF/                                                      071215
     1  1.47368130D+00, 1.55151322D+00, 1.67580598D+00, 2.04154525D+00, AsF  
     2  2.16023470D+00, 2.27312893D+00, 2.56031506D+00, 2.74368103D+00, AsF  
     3  2.98831167D+00, 3.16742957D+00, 3.50421954D+00, 3.71888449D+00, AsF  
     4  4.02849812D+00, 4.37830475D+00, 4.68298393D+00, 5.04915118D+00, AsF  
     5  5.41452907D+00, 5.66176384D+00, 5.94361194D+00, 6.16981561D+00, AsF  
     6  6.27318688D+00, 6.34068929D+00,      5*0.0D+00/                 AsF  
      DATA TQ_SeF/                                                      071215
     1  0.699999789529, 0.771700024046, 0.882600050919, 1.186599973333, SeF  
     2  1.496799925049, 1.837000077762, 2.054000290439, 2.249099983065, SeF  
     3  2.485999679905, 2.703399892149, 2.871499924798, 3.023500204862, SeF  
     4  3.414499948202, 3.654699731432, 3.829599917059, 3.934999962962, SeF  
     5  4.000000000000,     10*0.0D+00/                                 SeF  
      DATA  Q_SeF/                                                      071215
     1  1.64437156D+00, 1.70711749D+00, 1.80672444D+00, 2.09105024D+00, SeF  
     2  2.39128109D+00, 2.72627217D+00, 2.94157551D+00, 3.13662159D+00, SeF  
     3  3.38443482D+00, 3.64231174D+00, 3.87101570D+00, 4.10003348D+00, SeF  
     4  4.76479651D+00, 5.20858640D+00, 5.54158336D+00, 5.74500844D+00, SeF  
     5  5.87125416D+00,     10*0.0D+00/                                 SeF  
      DATA TQ_BrF/                                                      071215
     1  0.699999789529, 0.843900054273, 1.075000051982, 1.782400160334, BrF  
     2  1.983399926775, 2.171100084257, 2.514500341813, 2.660999882999, BrF  
     3  2.809400402489, 3.114699757845, 3.285599872721, 3.373199964921, BrF  
     4  3.461300063537, 3.542500017152, 3.623099968564, 3.697999766323, BrF  
     5  3.780299749363, 3.914699964983, 3.965600172236, 4.000000000000, BrF  
     6       7*0.0D+00/                                                 BrF  
      DATA  Q_BrF/                                                      071215
     1  1.00714241D+00, 1.14686585D+00, 1.37359295D+00, 2.07608391D+00, BrF  
     2  2.27672250D+00, 2.46493859D+00, 2.83212483D+00, 3.01298875D+00, BrF  
     3  3.21714387D+00, 3.70107825D+00, 4.00346634D+00, 4.16593704D+00, BrF  
     4  4.33506664D+00, 4.49875148D+00, 4.67456631D+00, 4.85721183D+00, BrF  
     5  5.08523573D+00, 5.50796958D+00, 5.67453126D+00, 5.78634031D+00, BrF  
     6       7*0.0D+00/                                                 BrF  
      DATA TQ_RbF/                                                      071215
     1  0.699999789529, 0.868199920303, 1.178399873043, 1.644399965906, RbF  
     2  1.801600140613, 1.950099805623, 2.181600191849, 2.321199736616, RbF  
     3  2.450699818627, 2.571899743493, 2.848200345432, 3.082900001274, RbF  
     4  3.312700217120, 3.547200129569, 3.628000077553, 3.718400216579, RbF  
     5  3.805000297311, 3.896599870194, 3.958800020533, 3.984399651089, RbF  
     6  4.000000000000,      6*0.0D+00/                                 RbF  
      DATA  Q_RbF/                                                      071215
     1  1.22874896D+00, 1.39414710D+00, 1.70134628D+00, 2.16554619D+00, RbF  
     2  2.32260678D+00, 2.47198690D+00, 2.71553179D+00, 2.87760907D+00, RbF  
     3  3.04331368D+00, 3.21314124D+00, 3.65049340D+00, 4.06716611D+00, RbF  
     4  4.50727122D+00, 4.99186429D+00, 5.17030479D+00, 5.37801243D+00, RbF  
     5  5.58348227D+00, 5.80334778D+00, 5.95085319D+00, 6.01053436D+00, RbF  
     6  6.04651428D+00,      6*0.0D+00/                                 RbF  
      DATA TQ_SrF/                                                      071215
     1  0.699999789529, 0.775100106965, 0.890899888586, 1.230899991830, SrF  
     2  1.732199967988, 1.883800029456, 2.040099971693, 2.174700174185, SrF  
     3  2.388400307119, 2.528599818103, 2.678100265308, 2.979199599574, SrF  
     4  3.252200050192, 3.461000057046, 3.566099966833, 3.672200111595, SrF  
     5  3.764800340226, 3.857299854832, 3.944099695165, 3.978499660086, SrF  
     6  4.000000000000,      6*0.0D+00/                                 SrF  
      DATA  Q_SrF/                                                      071215
     1  1.45593317D+00, 1.52938007D+00, 1.64313328D+00, 1.97951657D+00, SrF  
     2  2.47881141D+00, 2.63021867D+00, 2.78699287D+00, 2.92457243D+00, SrF  
     3  3.15877262D+00, 3.33110592D+00, 3.53542983D+00, 4.00919148D+00, SrF  
     4  4.49640761D+00, 4.89821600D+00, 5.11164532D+00, 5.34008673D+00, SrF  
     5  5.55642638D+00, 5.79413331D+00, 6.03860427D+00, 6.14060087D+00, SrF  
     6  6.20554024D+00,      6*0.0D+00/                                 SrF  
      DATA TQ_YF/                                                       071215
     1  0.699999789529, 0.858900111253, 1.124099979854, 1.787200041816, YF   
     2  1.970699809865, 2.151999669525, 2.472500070356, 2.600600330237, YF   
     3  2.738899750604, 2.938299709720, 3.302200237484, 3.493999850041, YF   
     4  3.601900243195, 3.696299731144, 3.787299914398, 3.872699954041, YF   
     5  3.949799823864, 3.980799570266, 4.000000000000,      8*0.0D+00/ YF   
      DATA  Q_YF/                                                       071215
     1  1.09228122D+00, 1.24747970D+00, 1.50886360D+00, 2.16849304D+00, YF   
     2  2.35173666D+00, 2.53363434D+00, 2.87494718D+00, 3.02988155D+00, YF   
     3  3.21433484D+00, 3.51227951D+00, 4.13705012D+00, 4.49896276D+00, YF   
     4  4.71379317D+00, 4.91337883D+00, 5.12286875D+00, 5.34155251D+00, YF   
     5  5.56103869D+00, 5.65516804D+00, 5.71500350D+00,      8*0.0D+00/ YF   
      DATA TQ_AgF/                                                      071215
     1  0.699999789529, 0.838600116278, 1.065600082015, 1.721299819092, AgF  
     2  1.893899948410, 2.066999906850, 2.381200152368, 2.512200293753, AgF  
     3  2.648799702770, 3.213700109251, 3.422000132982, 3.638100289022, AgF  
     4  3.743899851494, 3.847100304266, 3.940199607108, 3.976699789423, AgF  
     5  4.000000000000,     10*0.0D+00/                                 AgF  
      DATA  Q_AgF/                                                      071215
     1  1.13022809D+00, 1.26581346D+00, 1.48956028D+00, 2.14166216D+00, AgF  
     2  2.31401960D+00, 2.48777436D+00, 2.82268585D+00, 2.98132334D+00, AgF  
     3  3.16388750D+00, 4.09284570D+00, 4.48951988D+00, 4.92913093D+00, AgF  
     4  5.15827084D+00, 5.39295697D+00, 5.61436773D+00, 5.70342067D+00, AgF  
     5  5.76082205D+00,     10*0.0D+00/                                 AgF  
      DATA TQ_CdF/                                                      071215
     1  0.699999789529, 0.846000003604, 1.086700084934, 1.711000023535, CdF  
     2  1.890799880604, 2.071699725854, 2.224199954470, 2.372799954759, CdF  
     3  2.490799798360, 2.614099738361, 2.871199918431, 3.093800230727, CdF  
     4  3.318899776177, 3.601200293854, 3.826999854201, 3.933300085488, CdF  
     5  4.000000000000,     10*0.0D+00/                                 CdF  
      DATA  Q_CdF/                                                      071215
     1  1.42924096D+00, 1.57207097D+00, 1.80940063D+00, 2.43023629D+00, CdF  
     2  2.60969232D+00, 2.79096423D+00, 2.94710496D+00, 3.10815071D+00, CdF  
     3  3.24687251D+00, 3.40494852D+00, 3.77990907D+00, 4.14719892D+00, CdF  
     4  4.54737875D+00, 5.07542737D+00, 5.51075935D+00, 5.71911155D+00, CdF  
     5  5.85140105D+00,     10*0.0D+00/                                 CdF  
      DATA TQ_InF/                                                      071215
     1  0.699999789529, 0.839100127740, 1.063500033428, 1.738299821300, InF  
     2  1.906599916401, 2.073399768471, 2.403699696717, 2.534599826019, InF  
     3  2.675100197379, 3.011900266693, 3.272300295742, 3.490599770689, InF  
     4  3.595500281265, 3.697799762184, 3.793700055626, 3.882400176149, InF  
     5  3.953299900437, 3.982299603942, 4.000000000000,      8*0.0D+00/ InF  
      DATA  Q_InF/                                                      071215
     1  1.13563129D+00, 1.27174564D+00, 1.49296463D+00, 2.16414704D+00, InF  
     2  2.33221910D+00, 2.49955229D+00, 2.85135410D+00, 3.01034478D+00, InF  
     3  3.19893700D+00, 3.72555775D+00, 4.19119152D+00, 4.61237965D+00, InF  
     4  4.82523587D+00, 5.04174226D+00, 5.25662309D+00, 5.47133785D+00, InF  
     5  5.65880230D+00, 5.74038292D+00, 5.79166867D+00,      8*0.0D+00/ InF  
      DATA TQ_SnF/                                                      071215
     1  0.699999789529, 0.772500043556, 0.884100013930, 1.199899802738, SnF  
     2  1.483599930872, 1.799800174143, 1.973599871064, 2.139400307064, SnF  
     3  2.435799878563, 2.569499736398, 2.731899599238, 2.922500129796, SnF  
     4  3.078699912395, 3.445599702005, 3.580199930132, 3.708800008067, SnF  
     5  3.803000253469, 3.880500131050, 3.953099896070, 3.982299603942, SnF  
     6  4.000000000000,      6*0.0D+00/                                 SnF  
      DATA  Q_SnF/                                                      071215
     1  1.45327608D+00, 1.51894698D+00, 1.62200134D+00, 1.92266506D+00, SnF  
     2  2.19956708D+00, 2.51199800D+00, 2.68474883D+00, 2.85092369D+00, SnF  
     3  3.16760188D+00, 3.32985177D+00, 3.55037115D+00, 3.84704783D+00, SnF  
     4  4.12488935D+00, 4.88507521D+00, 5.18641411D+00, 5.48274613D+00, SnF  
     5  5.70843823D+00, 5.90255767D+00, 6.09263845D+00, 6.17128290D+00, SnF  
     6  6.21950301D+00,      6*0.0D+00/                                 SnF  
      DATA TQ_SbF/                                                      071215
     1  0.699999789529, 0.851399939207, 1.106600062466, 1.730400011273, SbF  
     2  1.933899908611, 2.140200306882, 2.286999885871, 2.431600171693, SbF  
     3  2.525200070587, 2.623599956696, 2.794500079960, 2.965200160482, SbF  
     4  3.091200168079, 3.223400010592, 3.365199771975, 3.539899955081, SbF  
     5  3.655299745432, 3.755600119901, 3.901799660074, 3.961600082593, SbF  
     6  4.000000000000,      6*0.0D+00/                                 SbF  
      DATA  Q_SbF/                                                      071215
     1  1.58637521D+00, 1.73436034D+00, 1.98593463D+00, 2.60634840D+00, SbF  
     2  2.80952211D+00, 3.01659203D+00, 3.16866028D+00, 3.33145694D+00, SbF  
     3  3.44915598D+00, 3.58656116D+00, 3.86021789D+00, 4.17037036D+00, SbF  
     4  4.41497531D+00, 4.68063760D+00, 4.97344128D+00, 5.34611866D+00, SbF  
     5  5.60179922D+00, 5.83268838D+00, 6.18967749D+00, 6.34397345D+00, SbF  
     6  6.44546314D+00,      6*0.0D+00/                                 SbF  
      DATA TQ_IF/                                                       071215
     1  0.699999789529, 0.778100180129, 0.898500058689, 1.238699786007, IF   
     2  1.767300041589, 2.071099710812, 2.240499800453, 2.499699975713, IF   
     3  2.662299911729, 2.834300039212, 3.089000119744, 3.287999927890, IF   
     4  3.376200026489, 3.464400130606, 3.546500112826, 3.630400130284, IF   
     5  3.699599799432, 3.776400002805, 3.839300129294, 3.912099898179, IF   
     6  3.964700152066, 3.986299693746, 4.000000000000,      4*0.0D+00/ IF   
      DATA  Q_IF/                                                       071215
     1  1.10838390D+00, 1.18456969D+00, 1.30262036D+00, 1.63884362D+00, IF   
     2  2.16515854D+00, 2.46884136D+00, 2.64088539D+00, 2.92618153D+00, IF   
     3  3.13188576D+00, 3.37719944D+00, 3.78896818D+00, 4.14365908D+00, IF   
     4  4.30891282D+00, 4.48039860D+00, 4.65002927D+00, 4.84326747D+00, IF   
     5  5.02767256D+00, 5.26571818D+00, 5.48460157D+00, 5.75448582D+00, IF   
     6  5.95290191D+00, 6.03389921D+00, 6.08492927D+00,      4*0.0D+00/ IF   
      DATA TQ_CsF/                                                      071215
     1  0.699999789529, 0.864200025758, 1.169100047599, 1.626600035993, CsF  
     2  1.780000219593, 1.926499912758, 2.150999645876, 2.285699859522, CsF  
     3  2.422200134821, 2.555000297068, 2.831699980682, 3.083400010984, CsF  
     4  3.321199723964, 3.450799812942, 3.561600290245, 3.642600144520, CsF  
     5  3.728499662094, 3.907699792571, 3.963000113968, 4.000000000000, CsF  
     6       7*0.0D+00/                                                 CsF  
      DATA  Q_CsF/                                                      071215
     1  1.28537680D+00, 1.44717356D+00, 1.74945587D+00, 2.20536802D+00, CsF  
     2  2.35866492D+00, 2.50610315D+00, 2.74224288D+00, 2.89799969D+00, CsF  
     3  3.07201789D+00, 3.25859518D+00, 3.69857361D+00, 4.14744312D+00, CsF  
     4  4.60497929D+00, 4.86832338D+00, 5.10390799D+00, 5.28388345D+00, CsF  
     5  5.48246146D+00, 5.91327495D+00, 6.04568677D+00, 6.13273939D+00, CsF  
     6       7*0.0D+00/                                                 CsF  
      DATA TQ_BaF/                                                      071215
     1  0.699999789529, 0.829199938494, 1.034000114605, 1.688000024819, BaF  
     2  1.840400138721, 1.993099843195, 2.139800316258, 2.373999983395, BaF  
     3  2.527299914641, 2.693499662649, 2.943699671862, 3.233399614288, BaF  
     4  3.331499952185, 3.427400241136, 3.526499984219, 3.621399930751, BaF  
     5  3.726699786823, 3.830099928907, 3.932900114318, 3.974099976243, BaF  
     6  4.000000000000,      6*0.0D+00/                                 BaF  
      DATA  Q_BaF/                                                      071215
     1  1.51907159D+00, 1.64596542D+00, 1.84827076D+00, 2.49912784D+00, BaF  
     2  2.65132966D+00, 2.80436213D+00, 2.95396254D+00, 3.21144253D+00, BaF  
     3  3.40282656D+00, 3.63566018D+00, 4.03309178D+00, 4.54895319D+00, BaF  
     4  4.73453235D+00, 4.92172929D+00, 5.12382675D+00, 5.33079957D+00, BaF  
     5  5.58483952D+00, 5.86798840D+00, 6.18388455D+00, 6.31802744D+00, BaF  
     6  6.40380171D+00,      6*0.0D+00/                                 BaF  
      DATA TQ_LaF/                                                      071215
     1  0.699999789529, 0.849899909505, 1.100799928760, 1.727599965223, LaF  
     2  1.918600041436, 2.105100027135, 2.396599856368, 2.508200199886, LaF  
     3  2.637000269254, 2.772000324895, 2.897699777287, 3.276599987621, LaF  
     4  3.545000076948, 3.778499850838, 3.915899995815, 4.000000000000, LaF  
     5      11*0.0D+00/                                                 LaF  
      DATA  Q_LaF/                                                      071215
     1  1.16198638D+00, 1.30889846D+00, 1.55662817D+00, 2.18039912D+00, LaF  
     2  2.37113273D+00, 2.55816976D+00, 2.86547536D+00, 2.99617756D+00, LaF  
     3  3.16063325D+00, 3.34982716D+00, 3.54082410D+00, 4.18533677D+00, LaF  
     4  4.68416819D+00, 5.13566055D+00, 5.40730302D+00, 5.57566810D+00, LaF  
     5      11*0.0D+00/                                                 LaF  
      DATA TQ_HoF/                                                      071215
     1  0.699999789529, 0.857000067668, 1.120899902913, 1.766900052414, HoF  
     2  1.955199934324, 2.139900318556, 2.445399702920, 2.567399885860, HoF  
     3  2.703899904567, 2.906499752472, 3.299200171851, 3.498799962067, HoF  
     4  3.605799960950, 3.701099832754, 3.873999983450, 3.950199832747, HoF  
     5  3.980999574756, 4.000000000000,      9*0.0D+00/                 HoF  
      DATA  Q_HoF/                                                      071215
     1  1.13426845D+00, 1.28795168D+00, 1.54839448D+00, 2.19126140D+00, HoF  
     2  2.37930797D+00, 2.56460326D+00, 2.88856928D+00, 3.03412968D+00, HoF  
     3  3.21318046D+00, 3.51177441D+00, 4.18408698D+00, 4.56135988D+00, HoF  
     4  4.77272977D+00, 4.96861834D+00, 5.35757322D+00, 5.54901071D+00, HoF  
     5  5.63005610D+00, 5.68099431D+00,      9*0.0D+00/                 HoF  
      DATA TQ_YbF/                                                      071215
     1  0.699999789529, 0.837600093354, 1.060699968646, 1.725599918833, YbF  
     2  1.888099914084, 2.049800193743, 2.381000148069, 2.518900433753, YbF  
     3  2.664899969189, 2.985599667069, 3.250400010905, 3.467200191184, YbF  
     4  3.569299736852, 3.674400163960, 3.762200277367, 3.853400129795, YbF  
     5  3.943099672587, 3.978099688827, 4.000000000000,      8*0.0D+00/ YbF  
      DATA  Q_YbF/                                                      071215
     1  1.47168911D+00, 1.60656685D+00, 1.82674262D+00, 2.48831806D+00, YbF  
     2  2.65061947D+00, 2.81289572D+00, 3.16629911D+00, 3.33477239D+00, YbF  
     3  3.53270845D+00, 4.03621522D+00, 4.50901088D+00, 4.92574551D+00, YbF  
     4  5.13199486D+00, 5.35517070D+00, 5.55507270D+00, 5.78037905D+00, YbF  
     5  6.02025180D+00, 6.11787924D+00, 6.17980455D+00,      8*0.0D+00/ YbF  
      DATA TQ_LuF/                                                      071215
     1  0.699999789529, 0.845800008430, 1.083700011149, 1.786400061569, LuF  
     2  1.963799956518, 2.140700271950, 2.468000206654, 2.603600111219, LuF  
     3  2.745499891426, 3.127300053185, 3.302200237484, 3.401199651822, LuF  
     4  3.499299973737, 3.596000292313, 3.694799700104, 3.780599756436, LuF  
     5  3.866699817560, 3.947999783222, 3.980199556795, 4.000000000000, LuF  
     6       7*0.0D+00/                                                 LuF  
      DATA  Q_LuF/                                                      071215
     1  1.12685349D+00, 1.26948033D+00, 1.50403579D+00, 2.20314130D+00, LuF  
     2  2.38032254D+00, 2.55785998D+00, 2.90743374D+00, 3.07302321D+00, LuF  
     3  3.26485561D+00, 3.87133412D+00, 4.18451743D+00, 4.36956560D+00, LuF  
     4  4.55861945D+00, 4.75266579D+00, 4.96493807D+00, 5.16912173D+00, LuF  
     5  5.40058851D+00, 5.64627662D+00, 5.75012879D+00, 5.81550224D+00, LuF  
     6       7*0.0D+00/                                                 LuF  
      DATA TQ_HgF/                                                      071215
     1  0.699999789529, 0.848299948110, 1.097599971705, 1.698299858063, HgF  
     2  1.871699926078, 2.044800079285, 2.352200225577, 2.475399859924, HgF  
     3  2.611299673906, 2.803000267696, 3.080099946894, 3.248499968036, HgF  
     4  3.402099670024, 3.509300206888, 3.611899703091, 3.818399800480, HgF  
     5  3.926400243350, 3.971300177434, 4.000000000000,      8*0.0D+00/ HgF  
      DATA  Q_HgF/                                                      071215
     1  1.42607136D+00, 1.57113517D+00, 1.81696799D+00, 2.41432357D+00, HgF  
     2  2.58739348D+00, 2.76103604D+00, 3.08770390D+00, 3.23542617D+00, HgF  
     3  3.41462520D+00, 3.69741378D+00, 4.16042869D+00, 4.46691366D+00, HgF  
     4  4.76100468D+00, 4.97504342D+00, 5.18674669D+00, 5.62280034D+00, HgF  
     5  5.84650188D+00, 5.93737578D+00, 5.99480025D+00,      8*0.0D+00/ HgF  
      DATA TQ_TlF/                                                      071215
     1  0.699999789529, 0.835100036045, 1.053800091159, 1.705099934314, TlF  
     2  1.862500073370, 2.021900324718, 2.359599675406, 2.492899840208, TlF  
     3  2.635400231343, 2.956199953375, 3.222000112173, 3.444499679107, TlF  
     4  3.665799973495, 3.868499858688, 3.947599774191, 3.979899559490, TlF  
     5  4.000000000000,     10*0.0D+00/                                 TlF  
      DATA  Q_TlF/                                                      071215
     1  1.20415072D+00, 1.33677509D+00, 1.55279995D+00, 2.20101147D+00, TlF  
     2  2.35822454D+00, 2.51814811D+00, 2.87825911D+00, 3.04095418D+00, TlF  
     3  3.23343743D+00, 3.73545589D+00, 4.20962742D+00, 4.63798920D+00, TlF  
     4  5.09471693D+00, 5.55615473D+00, 5.75354792D+00, 5.83755559D+00, TlF  
     5  5.89082627D+00,     10*0.0D+00/                                 TlF  
      DATA TQ_PbF/                                                      071215
     1  0.699999789529, 0.766500074903, 0.869499886030, 1.149000075188, PbF  
     2  1.427599879107, 1.746699918951, 1.892999928725, 2.041500003742, PbF  
     3  2.176000206659, 2.386300261983, 2.523600189403, 2.663099929409, PbF  
     4  2.998299960214, 3.172100119716, 3.365799787884, 3.592500214976, PbF  
     5  3.710100037545, 3.822499745409, 3.930200308918, 3.972700076838, PbF  
     6  4.000000000000,      6*0.0D+00/                                 PbF  
      DATA  Q_PbF/                                                      071215
     1  1.82372842D+00, 1.88492743D+00, 1.98114610D+00, 2.24858227D+00, PbF  
     2  2.52085318D+00, 2.83633089D+00, 2.98175049D+00, 3.13021806D+00, PbF  
     3  3.26738503D+00, 3.49737414D+00, 3.66549737D+00, 3.85462834D+00, PbF  
     4  4.38208561D+00, 4.68912774D+00, 5.05430213D+00, 5.52148472D+00, PbF  
     5  5.78782446D+00, 6.06045145D+00, 6.33706179D+00, 6.44942520D+00, PbF  
     6  6.52226016D+00,      6*0.0D+00/                                 PbF  
      DATA TQ_LiNa/                                                     071215
     1  0.699999789529, 0.844400042209, 1.093300081690, 1.532700163300, LiNa 
     2  1.675299945755, 1.809299959463, 2.009500218143, 2.147899768938, LiNa 
     3  2.304600302847, 2.445899713799, 2.865999792784, 3.045600092994, LiNa 
     4  3.221800126685, 3.307200347522, 3.394300043740, 3.540699974099, LiNa 
     5  3.685599907226, 3.869699886107, 3.949299812574, 3.980499563531, LiNa 
     6  4.000000000000,      6*0.0D+00/                                 LiNa 
      DATA  Q_LiNa/                                                     071215
     1  9.62738653D-01, 1.10247375D+00, 1.34622662D+00, 1.78141256D+00, LiNa 
     2  1.92352540D+00, 2.05839776D+00, 2.26914532D+00, 2.42888752D+00, LiNa 
     3  2.63019199D+00, 2.83217437D+00, 3.53406850D+00, 3.87188027D+00, LiNa 
     4  4.22825503D+00, 4.41348032D+00, 4.61107317D+00, 4.95343365D+00, LiNa 
     5  5.28391213D+00, 5.66489519D+00, 5.81360164D+00, 5.86932389D+00, LiNa 
     6  5.90344635D+00,      6*0.0D+00/                                 LiNa 
      DATA TQ_AsP/                                                      071215
     1  0.699999789529, 0.859600127310, 1.139799881688, 1.749099968574, AsP  
     2  1.947699857431, 2.137300258797, 2.430400255445, 2.546200096387, AsP  
     3  2.675600208700, 2.894899985442, 3.108299780183, 3.312100259791, AsP  
     4  3.519500445062, 3.721600140221, 3.810400378037, 3.895399958621, AsP  
     5  3.958600016166, 3.984199646599, 4.000000000000,      8*0.0D+00/ AsP  
      DATA  Q_AsP/                                                      071215
     1  1.26648657D+00, 1.42362686D+00, 1.70120155D+00, 2.30836703D+00, AsP  
     2  2.50677237D+00, 2.69705635D+00, 3.00745672D+00, 3.14448483D+00, AsP  
     3  3.31191304D+00, 3.63201484D+00, 3.98263123D+00, 4.34564606D+00, AsP  
     4  4.73688945D+00, 5.13685932D+00, 5.31980859D+00, 5.50182261D+00, AsP  
     5  5.64388251D+00, 5.70358553D+00, 5.74115384D+00,      8*0.0D+00/ AsP  
      DATA TQ_SbP/                                                      071215
     1  0.699999789529, 0.851399939207, 1.094600048439, 1.702199866082, SbP  
     2  1.884700005308, 2.060900361532, 2.364499754690, 2.486899702971, SbP  
     3  2.623699958988, 2.817099908109, 3.049300166665, 3.274000173927, SbP  
     4  3.518300417038, 3.631100144715, 3.735499671778, 3.893500098631, SbP  
     5  3.958200007432, 3.984099644354, 4.000000000000,      8*0.0D+00/ SbP  
      DATA  Q_SbP/                                                      071215
     1  1.40064371D+00, 1.55032570D+00, 1.79176485D+00, 2.39763955D+00, SbP  
     2  2.58001152D+00, 2.75696709D+00, 3.08006748D+00, 3.22688772D+00, SbP  
     3  3.40727419D+00, 3.69252128D+00, 4.07596312D+00, 4.47972936D+00, SbP  
     4  4.94594961D+00, 5.16964552D+00, 5.38225329D+00, 5.72108706D+00, SbP  
     5  5.86946951D+00, 5.93087513D+00, 5.96916977D+00,      8*0.0D+00/ SbP  
      DATA TQ_BeS/                                                      071215
     1  0.699999789529, 0.781300190854, 0.907099928599, 1.259900214831, BeS  
     2  1.590300159943, 1.946299890503, 2.138500286378, 2.323299782089, BeS  
     3  2.467500195377, 2.678800281158, 2.810400386561, 2.955599939552, BeS  
     4  3.110399665862, 3.265500352498, 3.368699864779, 3.479399577924, BeS  
     5  3.578999901556, 3.677900247269, 3.781199770581, 3.915199977830, BeS  
     6  3.966100183442, 4.000000000000,      5*0.0D+00/                 BeS  
      DATA  Q_BeS/                                                      071215
     1  6.78876628D-01, 7.54509713D-01, 8.73410578D-01, 1.21485581D+00, BeS  
     2  1.54048005D+00, 1.89421665D+00, 2.08586364D+00, 2.27089824D+00, BeS  
     3  2.41802412D+00, 2.64906378D+00, 2.80978034D+00, 3.00626391D+00, BeS  
     4  3.23846502D+00, 3.49487976D+00, 3.68200298D+00, 3.90479007D+00, BeS  
     5  4.13221534D+00, 4.38788722D+00, 4.68555606D+00, 5.10579865D+00, BeS  
     6  5.27072194D+00, 5.38088913D+00,      5*0.0D+00/                 BeS  
      DATA TQ_BS/                                                       071215
     1  0.699999789529, 0.761400209421, 0.857600081432, 1.104900023276, BS   
     2  1.499899990072, 1.938700011469, 2.175000181679, 2.390900275352, BS   
     3  2.541299986508, 2.692799647027, 2.907299770462, 3.112299706505, BS   
     4  3.374699995705, 3.527899882445, 3.682400130099, 3.772200306739, BS   
     5  3.872199942730, 3.950099830563, 3.980899572511, 4.000000000000, BS   
     6       7*0.0D+00/                                                 BS   
      DATA  Q_BS/                                                       071215
     1  9.77555310D-01, 1.03455098D+00, 1.12500379D+00, 1.36228550D+00, BS   
     2  1.74954024D+00, 2.18509171D+00, 2.42068487D+00, 2.63675144D+00, BS   
     3  2.79013003D+00, 2.95293761D+00, 3.21164446D+00, 3.49937749D+00, BS   
     4  3.92340765D+00, 4.19652024D+00, 4.49689210D+00, 4.69042964D+00, BS   
     5  4.92982698D+00, 5.13683252D+00, 5.22376976D+00, 5.27907124D+00, BS   
     6       7*0.0D+00/                                                 BS   
      DATA TQ_MgS/                                                      071215
     1  0.699999789529, 0.838900123156, 1.064200049624, 1.730700004059, MgS  
     2  1.900400071808, 2.072899755937, 2.392600150392, 2.521000382479, MgS  
     3  2.661499894049, 2.824699814473, 3.001600035704, 3.233799624002, MgS  
     4  3.445099691596, 3.549400182189, 3.653399701098, 3.751700027472, MgS  
     5  3.848900344760, 3.940299609366, 3.976899775052, 4.000000000000, MgS  
     6       7*0.0D+00/                                                 MgS  
      DATA  Q_MgS/                                                      071215
     1  1.12649192D+00, 1.26234462D+00, 1.48438282D+00, 2.14721234D+00, MgS  
     2  2.31667128D+00, 2.48975087D+00, 2.83002554D+00, 2.98519929D+00, MgS  
     3  3.17267745D+00, 3.41430515D+00, 3.70250072D+00, 4.11493596D+00, MgS  
     4  4.51815866D+00, 4.72696145D+00, 4.94335479D+00, 5.15804334D+00, MgS  
     5  5.38261710D+00, 5.60546344D+00, 5.69740925D+00, 5.75604548D+00, MgS  
     6       7*0.0D+00/                                                 MgS  
      DATA TQ_AlS/                                                      071215
     1  0.699999789529, 0.856900065374, 1.119099904635, 1.771500006180, AlS  
     2  1.956599969654, 2.139700313959, 2.450099805234, 2.575999835592, AlS  
     3  2.712000101129, 2.895699925969, 3.273600202589, 3.468700223637, AlS  
     4  3.566999902151, 3.665899975537, 3.762400282203, 3.856799890084, AlS  
     5  3.944499704197, 3.978599652901, 4.000000000000,      8*0.0D+00/ AlS  
      DATA  Q_AlS/                                                      071215
     1  1.40884925D+00, 1.56221938D+00, 1.82075505D+00, 2.46979370D+00, AlS  
     2  2.65463217D+00, 2.83831735D+00, 3.16786443D+00, 3.31883073D+00, AlS  
     3  3.49841451D+00, 3.76927417D+00, 4.41370364D+00, 4.78262371D+00, AlS  
     4  4.97768310D+00, 5.18258111D+00, 5.39463522D+00, 5.61728417D+00, AlS  
     5  5.83866090D+00, 5.92818850D+00, 5.98522204D+00,      8*0.0D+00/ AlS  
      DATA TQ_SiS/                                                      071215
     1  0.699999789529, 0.859100115841, 1.122099931766, 1.814500055268, SiS  
     2  2.021200376397, 2.218100212384, 2.518500425395, 2.631500138937, SiS  
     3  2.760400243055, 2.999099978937, 3.379000083953, 3.585800061333, SiS  
     4  3.673400140158, 3.767000393413, 3.837200083553, 3.907599790325, SiS  
     5  3.963400122932, 3.985799682521, 4.000000000000,      8*0.0D+00/ SiS  
      DATA  Q_SiS/                                                      071215
     1  1.07349325D+00, 1.22871726D+00, 1.48774988D+00, 2.17641957D+00, SiS  
     2  2.38282253D+00, 2.58023534D+00, 2.89734400D+00, 3.03044139D+00, SiS  
     3  3.19624452D+00, 3.54506075D+00, 4.19426326D+00, 4.58304781D+00, SiS  
     4  4.75428941D+00, 4.94392873D+00, 5.09526879D+00, 5.26269278D+00, SiS  
     5  5.41295905D+00, 5.47870824D+00, 5.52212582D+00,      8*0.0D+00/ SiS  
      DATA TQ_PS/                                                       071215
     1  0.699999789529, 0.752900062340, 0.835100036045, 1.049900173955, PS   
     2  1.325500037379, 1.625600013884, 1.760500225619, 1.885899973112, PS   
     3  2.195399697697, 2.343900251935, 2.504400098771, 2.797400145633, PS   
     4  3.230099534140, 3.533199802608, 3.645999904334, 3.751800029842, PS   
     5  3.905199736428, 3.962700107245, 3.985599678030, 4.000000000000, PS   
     6       7*0.0D+00/                                                 PS   
      DATA  Q_PS/                                                       071215
     1  1.42004327D+00, 1.46740470D+00, 1.54221062D+00, 1.74310379D+00, PS   
     2  2.00844808D+00, 2.30278429D+00, 2.43628390D+00, 2.56156490D+00, PS   
     3  2.89126103D+00, 3.07080631D+00, 3.28501784D+00, 3.72942800D+00, PS   
     4  4.48749400D+00, 5.06609949D+00, 5.28905093D+00, 5.50301010D+00, PS   
     5  5.82799367D+00, 5.95696026D+00, 6.00969997D+00, 6.04329332D+00, PS   
     6       7*0.0D+00/                                                 PS   
      DATA TQ_CaS/                                                      071215
     1  0.699999789529, 0.835000033752, 1.058999974823, 1.695399932681, CaS  
     2  1.852899973616, 2.011700269636, 2.346300305918, 2.481899574829, CaS  
     3  2.626000011709, 2.957699987933, 3.225099887244, 3.448099754046, CaS  
     4  3.568599787160, 3.671800102074, 3.772100313975, 3.857599833681, CaS  
     5  3.944499704197, 3.978599652901, 4.000000000000,      8*0.0D+00/ CaS  
      DATA  Q_CaS/                                                      071215
     1  1.30322360D+00, 1.43626469D+00, 1.65810648D+00, 2.29211613D+00, CaS  
     2  2.44947478D+00, 2.60884689D+00, 2.96572837D+00, 3.13122806D+00, CaS  
     3  3.32603761D+00, 3.84629261D+00, 4.32324986D+00, 4.75008416D+00, CaS  
     4  4.99194134D+00, 5.20827447D+00, 5.43066044D+00, 5.63246032D+00, CaS  
     5  5.84985189D+00, 5.93820425D+00, 5.99438416D+00,      8*0.0D+00/ CaS  
      DATA TQ_ScS/                                                      071215
     1  0.699999789529, 0.852099955265, 1.094700045881, 1.737499840538, ScS  
     2  1.928999850102, 2.111099700254, 2.407699794305, 2.528599818103, ScS  
     3  2.666800011179, 2.886800271471, 3.245599902035, 3.437399769455, ScS  
     4  3.661599887760, 3.762400282203, 3.857799819580, 3.944199697423, ScS  
     5  3.978399667271, 4.000000000000,      9*0.0D+00/                 ScS  
      DATA  Q_ScS/                                                      071215
     1  1.60243460D+00, 1.75235499D+00, 1.99273285D+00, 2.63324127D+00, ScS  
     2  2.82453249D+00, 3.00728559D+00, 3.32175118D+00, 3.46547967D+00, ScS  
     3  3.64596325D+00, 3.97014949D+00, 4.57867568D+00, 4.93311394D+00, ScS  
     4  5.37343694D+00, 5.58521672D+00, 5.79605440D+00, 5.99615366D+00, ScS  
     5  6.07758130D+00, 6.12957197D+00,      9*0.0D+00/                 ScS  
      DATA TQ_TiS/                                                      071215
     1  0.699999789529, 0.778500189884, 0.901900048465, 1.226999950930, TiS  
     2  1.416399997560, 1.637999911936, 1.749099968574, 1.858900111253, TiS  
     3  2.015700363782, 2.173000131719, 2.321399740947, 2.504600104093, TiS  
     4  2.664799966979, 2.836300084236, 3.223899974313, 3.411899879788, TiS  
     5  3.592300210556, 3.710400044016, 3.831599961579, 3.933100099903, TiS  
     6  3.973700004984, 4.000000000000,      5*0.0D+00/                 TiS  
      DATA  Q_TiS/                                                      071215
     1  1.59670497D+00, 1.66581588D+00, 1.77741912D+00, 2.08363814D+00, TiS  
     2  2.26977973D+00, 2.50624484D+00, 2.63797351D+00, 2.77807818D+00, TiS  
     3  2.99213843D+00, 3.21604788D+00, 3.43109300D+00, 3.70402378D+00, TiS  
     4  3.95475223D+00, 4.23889041D+00, 4.93698429D+00, 5.29887560D+00, TiS  
     5  5.66394050D+00, 5.91891141D+00, 6.20006879D+00, 6.45333162D+00, TiS  
     6  6.55900501D+00, 6.62859385D+00,      5*0.0D+00/                 TiS  
      DATA TQ_CrS/                                                      071215
     1  0.699999789529, 0.858700106665, 1.133300040890, 1.756200143403, CrS  
     2  1.953299886377, 2.142900118252, 2.437599752936, 2.554300280902, CrS  
     3  2.684799960910, 2.905799736730, 3.266900386114, 3.453599880650, CrS  
     4  3.643800059748, 3.743799849270, 3.842900209781, 3.939599631422, CrS  
     5  4.000000000000,     10*0.0D+00/                                 CrS  
      DATA  Q_CrS/                                                      071215
     1  1.24506930D+00, 1.40119146D+00, 1.67305591D+00, 2.29359627D+00, CrS  
     2  2.49045036D+00, 2.68062311D+00, 2.99215326D+00, 3.12993659D+00, CrS  
     3  3.29850588D+00, 3.62088721D+00, 4.23238175D+00, 4.57932387D+00, CrS  
     4  4.95144728D+00, 5.15623506D+00, 5.36664795D+00, 5.57875973D+00, CrS  
     5  5.71382724D+00,     10*0.0D+00/                                 CrS  
      DATA TQ_CuS/                                                      071215
     1  0.699999789529, 0.750299994957, 0.828799948174, 1.033700122621, CuS  
     2  1.307500011196, 1.604800010740, 1.777100146781, 1.940200034602, CuS  
     3  2.091300171857, 2.224299947450, 2.327299868703, 2.433400046066, CuS  
     4  2.573299774941, 2.758400193999, 2.916499990094, 3.092500199403, CuS  
     5  3.258700192059, 3.445099691596, 3.570299693710, 3.683500053486, CuS  
     6  3.852500193248, 4.000000000000,      5*0.0D+00/                 CuS  
      DATA  Q_CuS/                                                      071215
     1  1.59758035D+00, 1.64448476D+00, 1.71840402D+00, 1.91458342D+00, CuS  
     2  2.18159678D+00, 2.47505922D+00, 2.64611455D+00, 2.80919309D+00, CuS  
     3  2.96537730D+00, 3.11519744D+00, 3.24520250D+00, 3.39512858D+00, CuS  
     4  3.61770544D+00, 3.94682218D+00, 4.24784052D+00, 4.59471763D+00, CuS  
     5  4.92887678D+00, 5.30939736D+00, 5.56928001D+00, 5.80946253D+00, CuS  
     6  6.18288387D+00, 6.52255044D+00,      5*0.0D+00/                 CuS  
      DATA TQ_GeS/                                                      071215
     1  0.699999789529, 0.851499941501, 1.092700097037, 1.737599838133, GeS  
     2  1.933999910754, 2.118799856866, 2.412299906195, 2.527799877511, GeS  
     3  2.662099907309, 2.873799973614, 3.093700228317, 3.303800272696, GeS  
     4  3.521900318617, 3.728799641306, 3.814900053161, 3.896799855456, GeS  
     5  3.959100027084, 3.984399651089, 4.000000000000,      8*0.0D+00/ GeS  
      DATA  Q_GeS/                                                      071215
     1  1.27980878D+00, 1.42902604D+00, 1.66790243D+00, 2.31043511D+00, GeS  
     2  2.50664978D+00, 2.69215813D+00, 3.00328295D+00, 3.14018309D+00, GeS  
     3  3.31454370D+00, 3.62418563D+00, 3.98583050D+00, 4.36098128D+00, GeS  
     4  4.77380685D+00, 5.18498757D+00, 5.36324114D+00, 5.53949430D+00, GeS  
     5  5.68022614D+00, 5.73950552D+00, 5.77675523D+00,      8*0.0D+00/ GeS  
      DATA TQ_AsS/                                                      071215
     1  0.699999789529, 0.771200011851, 0.881400080510, 1.190800024619, AsS  
     2  1.756000138363, 1.937699990040, 2.115499789747, 2.411299881944, AsS  
     3  2.528599818103, 2.661799900679, 2.862199697616, 3.072699762134, AsS  
     4  3.286599895708, 3.503300067008, 3.605899953713, 3.705999944316, AsS  
     5  3.880500131050, 3.953199898253, 3.982199601697, 4.000000000000, AsS  
     6       7*0.0D+00/                                                 AsS  
      DATA  Q_AsS/                                                      071215
     1  1.90888638D+00, 1.97550797D+00, 2.07993622D+00, 2.37912670D+00, AsS  
     2  2.93719855D+00, 3.11805640D+00, 3.29615348D+00, 3.60985644D+00, AsS  
     3  3.74940990D+00, 3.92311347D+00, 4.21632925D+00, 4.56113977D+00, AsS  
     4  4.94194387D+00, 5.35198372D+00, 5.55401728D+00, 5.75737512D+00, AsS  
     5  6.13662481D+00, 6.30857585D+00, 6.37989016D+00, 6.42443632D+00, AsS  
     6       7*0.0D+00/                                                 AsS  
      DATA TQ_SeS/                                                      071215
     1  0.699999789529, 0.835800052091, 1.069800179188, 1.473599864515, SeS  
     2  1.692899997008, 1.849699914331, 1.991699811379, 2.261900276502, SeS  
     3  2.464900136736, 2.640200325845, 2.812400243740, 3.193999657593, SeS  
     4  3.479099599429, 3.604000091217, 3.719300235992, 3.885500249731, SeS  
     5  3.954899935374, 3.982899617413, 4.000000000000,      8*0.0D+00/ SeS  
      DATA  Q_SeS/                                                      071215
     1  1.77402132D+00, 1.90782146D+00, 2.13955748D+00, 2.54147599D+00, SeS  
     2  2.76135295D+00, 2.92343126D+00, 3.07954556D+00, 3.41249735D+00, SeS  
     3  3.69494715D+00, 3.96179559D+00, 4.24415789D+00, 4.92861083D+00, SeS  
     4  5.47884700D+00, 5.72874530D+00, 5.96538990D+00, 6.32500411D+00, SeS  
     5  6.48617983D+00, 6.55353902D+00, 6.59535253D+00,      8*0.0D+00/ SeS  
      DATA TQ_SrS/                                                      071215
     1  0.699999789529, 0.827399982056, 1.039199975669, 1.641299892791, SrS  
     2  1.778800189464, 1.926599910252, 2.061800294448, 2.278899822029, SrS  
     3  2.415399981372, 2.562000270192, 2.901999651278, 3.188299687136, SrS  
     4  3.442899645800, 3.573899779715, 3.690599613192, 3.867499835839, SrS  
     5  3.948399792254, 3.980099554550, 4.000000000000,      8*0.0D+00/ SrS  
      DATA  Q_SrS/                                                      071215
     1  1.46604579D+00, 1.59217361D+00, 1.80253855D+00, 2.40298363D+00, SrS  
     2  2.54040167D+00, 2.68869679D+00, 2.82682132D+00, 3.06462435D+00, SrS  
     3  3.23215780D+00, 3.43161269D+00, 3.96837656D+00, 4.48213266D+00, SrS  
     4  4.97250775D+00, 5.23624394D+00, 5.47966304D+00, 5.87178599D+00, SrS  
     5  6.06277569D+00, 6.13956452D+00, 6.18827216D+00,      8*0.0D+00/ SrS  
      DATA TQ_YS/                                                       071215
     1  0.699999789529, 0.853199980498, 1.100099912623, 1.696699899232, YS   
     2  1.875400017613, 2.049300182297, 2.352000240446, 2.479699547904, YS   
     3  2.620899894806, 2.853000161894, 3.226199807431, 3.440699600004, YS   
     4  3.673000130637, 3.760300231432, 3.851400270802, 3.943199674844, YS   
     5  4.000000000000,     10*0.0D+00/                                 YS   
      DATA  Q_YS/                                                       071215
     1  1.76528851D+00, 1.91699098D+00, 2.16235660D+00, 2.75746135D+00, YS   
     2  2.93603154D+00, 3.11065868D+00, 3.43259473D+00, 3.58578921D+00, YS   
     3  3.77244224D+00, 4.11897569D+00, 4.75729815D+00, 5.15402475D+00, YS   
     4  5.61012747D+00, 5.79465620D+00, 5.99744976D+00, 6.21167885D+00, YS   
     5  6.34797463D+00,     10*0.0D+00/                                 YS   
      DATA TQ_SnS/                                                      071215
     1  0.699999789529, 0.846799984302, 1.080399929986, 1.705399941372, SnS  
     2  1.875100010191, 2.039399954843, 2.362899715926, 2.495199886040, SnS  
     3  2.636800264515, 2.816099979519, 3.010400232315, 3.290799991082, SnS  
     4  3.518900431050, 3.633700198315, 3.741499798109, 3.892800150214, SnS  
     5  3.957699996514, 3.983899639864, 4.000000000000,      8*0.0D+00/ SnS  
      DATA  Q_SnS/                                                      071215
     1  1.41222299D+00, 1.55739393D+00, 1.78931385D+00, 2.41256248D+00, SnS  
     2  2.58215406D+00, 2.74709121D+00, 3.09163812D+00, 3.25198247D+00, SnS  
     3  3.44162065D+00, 3.70882076D+00, 4.02839129D+00, 4.53167511D+00, SnS  
     4  4.96850143D+00, 5.19715077D+00, 5.42076276D+00, 5.76804908D+00, SnS  
     5  5.93768685D+00, 6.01044796D+00, 6.05637840D+00,      8*0.0D+00/ SnS  
      DATA TQ_TeS/                                                      071215
     1  0.699999789529, 0.845800008430, 1.077699980792, 1.700199819026, TeS  
     2  1.864200030451, 2.024000169681, 2.359399690275, 2.493099844193, TeS  
     3  2.634600212388, 2.805500320349, 2.983999631107, 3.274500138099, TeS  
     4  3.510800241890, 3.621099924078, 3.722000112504, 3.806000319232, TeS  
     5  3.877000051316, 3.951599863317, 3.981699590472, 4.000000000000, TeS  
     6       7*0.0D+00/                                                 TeS  
      DATA  Q_TeS/                                                      071215
     1  1.42722754D+00, 1.57146380D+00, 1.80174707D+00, 2.42254998D+00, TeS  
     2  2.58645293D+00, 2.74688450D+00, 3.10519109D+00, 3.26869070D+00, TeS  
     3  3.45998558D+00, 3.71611643D+00, 4.00948750D+00, 4.53058251D+00, TeS  
     4  4.98475298D+00, 5.20596786D+00, 5.41743480D+00, 5.60612046D+00, TeS  
     5  5.77940135D+00, 5.97723726D+00, 6.06123559D+00, 6.11325933D+00, TeS  
     6       7*0.0D+00/                                                 TeS  
      DATA TQ_BaS/                                                      071215
     1  0.699999789529, 0.847799960174, 1.087000092313, 1.635199981666, BaS  
     2  1.798200141251, 1.948699833809, 2.188999661227, 2.329399914175, BaS  
     3  2.505200120059, 2.666800011179, 2.953599893474, 3.199499786491, BaS  
     4  3.303300261692, 3.404899726651, 3.498799962067, 3.599200363022, BaS  
     5  3.674000154440, 3.757800172040, 3.903999709479, 3.961500080352, BaS  
     6  4.000000000000,      6*0.0D+00/                                 BaS  
      DATA  Q_BaS/                                                      071215
     1  1.53283734D+00, 1.67940239D+00, 1.91731832D+00, 2.46431203D+00, BaS  
     2  2.62726562D+00, 2.77859778D+00, 3.03109185D+00, 3.19405618D+00, BaS  
     3  3.42266747D+00, 3.65862728D+00, 4.13044215D+00, 4.57547221D+00, BaS  
     4  4.77182947D+00, 4.96932703D+00, 5.16016033D+00, 5.38304840D+00, BaS  
     5  5.57049400D+00, 5.80752285D+00, 6.27462450D+00, 6.46573283D+00, BaS  
     6  6.59296939D+00,      6*0.0D+00/                                 BaS  
      DATA TQ_LaS/                                                      071215
     1  0.699999789529, 0.845800008430, 1.078699954425, 1.683499919501, LaS  
     2  1.848699938459, 2.010300236685, 2.335700062283, 2.469300235975, LaS  
     3  2.607899797293, 2.787799922314, 3.002200049093, 3.224499930779, LaS  
     4  3.439099649133, 3.568199815908, 3.686599837578, 3.872699954041, LaS  
     5  3.950399837114, 3.980999574756, 4.000000000000,      8*0.0D+00/ LaS  
      DATA  Q_LaS/                                                      071215
     1  1.78060976D+00, 1.92502620D+00, 2.15649122D+00, 2.75981423D+00, LaS  
     2  2.92492491D+00, 3.08715127D+00, 3.43364319D+00, 3.59563395D+00, LaS  
     3  3.78117363D+00, 4.04885109D+00, 4.40232367D+00, 4.79924630D+00, LaS  
     4  5.20441770D+00, 5.45828751D+00, 5.70072140D+00, 6.11051627D+00, LaS  
     5  6.29509219D+00, 6.37024323D+00, 6.41762581D+00,      8*0.0D+00/ LaS  
      DATA TQ_PbS/                                                      071215
     1  0.699999789529, 0.828599953015, 1.040099956514, 1.663499914443, PbS  
     2  1.804100081798, 1.959700047884, 2.101800262347, 2.332999998281, PbS  
     3  2.483599618397, 2.644100043288, 2.955599939552, 3.252800063287, PbS  
     4  3.360199639399, 3.464300128442, 3.573199762992, 3.666499987784, PbS  
     5  3.769300449019, 3.852800172097, 3.941499636461, 3.977499731940, PbS  
     6  4.000000000000,      6*0.0D+00/                                 PbS  
      DATA  Q_PbS/                                                      071215
     1  1.48201093D+00, 1.60937469D+00, 1.81949589D+00, 2.44126902D+00, PbS  
     2  2.58177856D+00, 2.73779484D+00, 2.88275905D+00, 3.13655642D+00, PbS  
     3  3.32373417D+00, 3.54685371D+00, 4.04460130D+00, 4.57966706D+00, PbS  
     4  4.78367853D+00, 4.98651452D+00, 5.20601550D+00, 5.40480217D+00, PbS  
     5  5.64565073D+00, 5.86728274D+00, 6.13230918D+00, 6.24742834D+00, PbS  
     6  6.32102791D+00,      6*0.0D+00/                                 PbS  
      DATA TQ_BiS/                                                      071215
     1  0.699999789529, 0.765100111829, 0.865899980940, 1.151300063991, BiS  
     2  1.679299853698, 1.836500065994, 1.986199861429, 2.217500200726, BiS  
     3  2.353400136360, 2.481799572266, 2.605399979808, 2.933200067623, BiS  
     4  3.204099890089, 3.321399728300, 3.434799953477, 3.545500088908, BiS  
     5  3.658199813099, 3.747699936021, 3.842300196283, 3.937999746741, BiS  
     6  3.976099832535, 4.000000000000,      5*0.0D+00/                 BiS  
      DATA  Q_BiS/                                                      071215
     1  2.11147749D+00, 2.17399382D+00, 2.27148070D+00, 2.55076685D+00, BiS  
     2  3.07417542D+00, 3.23089945D+00, 3.38121201D+00, 3.62404355D+00, BiS  
     3  3.78117383D+00, 3.94439084D+00, 4.11644071D+00, 4.63823794D+00, BiS  
     4  5.12505043D+00, 5.34752788D+00, 5.56899571D+00, 5.79296527D+00, BiS  
     5  6.03357199D+00, 6.23865734D+00, 6.47292081D+00, 6.72798025D+00, BiS  
     6  6.83335553D+00, 6.90015243D+00,      5*0.0D+00/                 BiS  
      DATA TQ_LiCl/                                                     071215
     1  0.699999789529, 0.773200060628, 0.886499954749, 1.207299956400, LiCl 
     2  1.789199992433, 1.981199978118, 2.166499987459, 2.460800044263, LiCl 
     3  2.578799898490, 2.709900053578, 2.888900322818, 3.223899974313, LiCl 
     4  3.398599729775, 3.574799801216, 3.742499820353, 3.839400131472, LiCl 
     5  3.914699964983, 3.967600217058, 4.000000000000,      8*0.0D+00/ LiCl 
      DATA  Q_LiCl/                                                     071215
     1  7.24768596D-01, 7.93383593D-01, 9.00966103D-01, 1.21181033D+00, LiCl 
     2  1.78718176D+00, 1.97850695D+00, 2.16433509D+00, 2.47700676D+00, LiCl 
     3  2.61781006D+00, 2.78946576D+00, 3.05126784D+00, 3.61758245D+00, LiCl 
     4  3.94607611D+00, 4.29981128D+00, 4.66258291D+00, 4.88479388D+00, LiCl 
     5  5.06109594D+00, 5.18476048D+00, 5.25971710D+00,      8*0.0D+00/ LiCl 
      DATA TQ_BeCl/                                                     071215
     1  0.699999789529, 0.758900217840, 0.850999930032, 1.091600125173, BeCl 
     2  1.452999874082, 1.852299959853, 2.067199891942, 2.267600401009, BeCl 
     3  2.414499959546, 2.553100253189, 2.785799871642, 3.064800082251, BeCl 
     4  3.340700179958, 3.511800265243, 3.666499987784, 3.869599883822, BeCl 
     5  3.947899780964, 3.980099554550, 4.000000000000,      8*0.0D+00/ BeCl 
      DATA  Q_BeCl/                                                     071215
     1  1.01303756D+00, 1.06806478D+00, 1.15508193D+00, 1.38657140D+00, BeCl 
     2  1.74108007D+00, 2.13729822D+00, 2.35151499D+00, 2.55227789D+00, BeCl 
     3  2.70277274D+00, 2.85259730D+00, 3.13652656D+00, 3.54506746D+00, BeCl 
     4  4.01564067D+00, 4.33466354D+00, 4.64046036D+00, 5.07869482D+00, BeCl 
     5  5.26546061D+00, 5.34574114D+00, 5.39634240D+00,      8*0.0D+00/ BeCl 
      DATA TQ_BCl/                                                      071215
     1  0.699999789529, 0.759700238573, 0.853199980498, 1.095100035650, BCl  
     2  1.455899937111, 1.858100092901, 2.066399951573, 2.263100302714, BCl  
     3  2.416800015323, 2.564700078026, 2.683000091193, 2.804600301394, BCl  
     4  3.022800254356, 3.344600267680, 3.525200078723, 3.615299785581, BCl  
     5  3.696299731144, 3.874699999285, 3.950899848032, 3.981199579246, BCl  
     6  4.000000000000,      6*0.0D+00/                                 BCl  
      DATA  Q_BCl/                                                      071215
     1  7.37441906D-01, 7.93462072D-01, 8.82135936D-01, 1.11549152D+00, BCl  
     2  1.46987335D+00, 1.86916669D+00, 2.07683017D+00, 2.27387117D+00, BCl  
     3  2.43142766D+00, 2.59219536D+00, 2.73196314D+00, 2.88891880D+00, BCl  
     4  3.20677904D+00, 3.75099877D+00, 4.08856235D+00, 4.26602586D+00, BCl  
     5  4.43323403D+00, 4.84721604D+00, 5.05170137D+00, 5.13804169D+00, BCl  
     6  5.19295286D+00,      6*0.0D+00/                                 BCl  
      DATA TQ_NaCl/                                                     071215
     1  0.699999789529, 0.866199973031, 1.171300037811, 1.632200056376, NaCl 
     2  1.789099994902, 1.937299981469, 2.307100361074, 2.554100276283, NaCl 
     3  2.843400241294, 3.077299877334, 3.320199702282, 3.455399924176, NaCl 
     4  3.565799988394, 3.734799657316, 3.908299806045, 3.963200118450, NaCl 
     5  4.000000000000,     10*0.0D+00/                                 NaCl 
      DATA  Q_NaCl/                                                     071215
     1  1.21411312D+00, 1.37744154D+00, 1.67945731D+00, 2.13847530D+00, NaCl 
     2  2.29522156D+00, 2.44429778D+00, 2.84840354D+00, 3.17831697D+00, NaCl 
     3  3.63659400D+00, 4.05287419D+00, 4.51907729D+00, 4.79304082D+00, NaCl 
     4  5.02583867D+00, 5.40156744D+00, 5.80618057D+00, 5.93345807D+00, NaCl 
     5  6.01730268D+00,     10*0.0D+00/                                 NaCl 
      DATA TQ_MgCl/                                                     071215
     1  0.699999789529, 0.833700003951, 1.049100156205, 1.697099888939, MgCl 
     2  1.852699969029, 2.010100231978, 2.346400308167, 2.481099554326, MgCl 
     3  2.624899986494, 2.942599647649, 3.204599901294, 3.429900291207, MgCl 
     4  3.539999957356, 3.650299628764, 3.747799938246, 3.845800275021, MgCl 
     5  3.939399645837, 3.976499803794, 4.000000000000,      8*0.0D+00/ MgCl 
      DATA  Q_MgCl/                                                     071215
     1  1.46542681D+00, 1.59643102D+00, 1.80891701D+00, 2.45349015D+00, MgCl 
     2  2.60887852D+00, 2.76679972D+00, 3.12556080D+00, 3.29010679D+00, MgCl 
     3  3.48468756D+00, 3.98219015D+00, 4.44898938D+00, 4.88190063D+00, MgCl 
     4  5.10405312D+00, 5.33569142D+00, 5.55066474D+00, 5.77880286D+00, MgCl 
     5  6.00824711D+00, 6.10178589D+00, 6.16159750D+00,      8*0.0D+00/ MgCl 
      DATA TQ_AlCl/                                                     071215
     1  0.699999789529, 0.835100036045, 1.052900111294, 1.710100046995, AlCl 
     2  1.868399924416, 2.027999874373, 2.364199747422, 2.498799957779, AlCl 
     3  2.642100188189, 2.986699691793, 3.259800216067, 3.495199878048, AlCl 
     4  3.604800033320, 3.703599889674, 3.877300058102, 3.951699865500, AlCl 
     5  3.981599588227, 4.000000000000,      9*0.0D+00/                 AlCl 
      DATA  Q_AlCl/                                                     071215
     1  1.16632093D+00, 1.29871353D+00, 1.51359746D+00, 2.16740008D+00, AlCl 
     2  2.32549531D+00, 2.48562556D+00, 2.84424191D+00, 3.00861617D+00, AlCl 
     3  3.20241502D+00, 3.74543787D+00, 4.23703049D+00, 4.69376611D+00, AlCl 
     4  4.91705855D+00, 5.12663802D+00, 5.53056875D+00, 5.72513759D+00, AlCl 
     5  5.80751116D+00, 5.85933540D+00,      9*0.0D+00/                 AlCl 
      DATA TQ_SiCl/                                                     071215
     1  0.699999789529, 0.745499883965, 0.816400084922, 1.000700014940, SiCl 
     2  1.486499863604, 1.704099910786, 1.871599923605, 2.015200352014, SiCl 
     3  2.455499925772, 2.624799984202, 2.793100048256, 3.165399962696, SiCl 
     4  3.413999935046, 3.528299853367, 3.639800324067, 3.750499999032, SiCl 
     5  3.860499675897, 3.943799688392, 3.978299674457, 4.000000000000, SiCl 
     6       7*0.0D+00/                                                 SiCl 
      DATA  Q_SiCl/                                                     071215
     1  1.47805109D+00, 1.51940163D+00, 1.58463141D+00, 1.75783720D+00, SiCl 
     2  2.22949671D+00, 2.44560655D+00, 2.61853715D+00, 2.77768337D+00, SiCl 
     3  3.35534318D+00, 3.61474028D+00, 3.89222021D+00, 4.56324532D+00, SiCl 
     4  5.04534641D+00, 5.27525778D+00, 5.50542211D+00, 5.74194749D+00, SiCl 
     5  5.98966491D+00, 6.18986546D+00, 6.27661231D+00, 6.33234771D+00, SiCl 
     6       7*0.0D+00/                                                 SiCl 
      DATA TQ_PCl/                                                      071215
     1  0.699999789529, 0.854200003438, 1.114900013655, 1.742699836246, PCl  
     2  1.930899844326, 2.114299765340, 2.413699940146, 2.531799757918, PCl  
     3  2.667900035489, 2.881200134546, 3.225999821942, 3.423300159019, PCl  
     4  3.620699915181, 3.707199971638, 3.803000253469, 3.932200164770, PCl  
     5  4.000000000000,     10*0.0D+00/                                 PCl  
      DATA  Q_PCl/                                                      071215
     1  1.44546491D+00, 1.59647195D+00, 1.85379618D+00, 2.47846421D+00, PCl  
     2  2.66635381D+00, 2.85030511D+00, 3.16749977D+00, 3.30774341D+00, PCl  
     3  3.48505725D+00, 3.79826909D+00, 4.38163433D+00, 4.74756760D+00, PCl  
     4  5.13558355D+00, 5.31457969D+00, 5.52025140D+00, 5.80766546D+00, PCl  
     5  5.96049072D+00,     10*0.0D+00/                                 PCl  
      DATA TQ_KCl/                                                      071215
     1  0.699999789529, 0.830799937472, 1.038899983685, 1.549700175445, KCl  
     2  1.693399984143, 1.830599927128, 2.052600258198, 2.187699754444, KCl  
     3  2.329499916340, 2.468900226953, 2.721600165677, 2.978899621463, KCl  
     4  3.215400149197, 3.433500045488, 3.534999843571, 3.634600216868, KCl  
     5  3.845800275021, 3.939799617007, 3.976499803794, 4.000000000000, KCl  
     6       7*0.0D+00/                                                 KCl  
      DATA  Q_KCl/                                                      071215
     1  1.43932125D+00, 1.56873247D+00, 1.77532901D+00, 2.28444451D+00, KCl  
     2  2.42805737D+00, 2.56615929D+00, 2.79959015D+00, 2.95566395D+00, KCl  
     3  3.13659876D+00, 3.33329598D+00, 3.73417908D+00, 4.18937739D+00, KCl  
     4  4.63955876D+00, 5.08108618D+00, 5.29759897D+00, 5.51700611D+00, KCl  
     5  5.98850873D+00, 6.19044326D+00, 6.26666816D+00, 6.31459969D+00, KCl  
     6       7*0.0D+00/                                                 KCl  
      DATA TQ_CaCl/                                                     071215
     1  0.699999789529, 0.817200101861, 0.999499989485, 1.546300105701, CaCl 
     2  1.717999841071, 1.879600121517, 2.031999775256, 2.266900385719, CaCl 
     3  2.417900041999, 2.581499960589, 2.886200256801, 3.187399752405, CaCl 
     4  3.415399971884, 3.524400136879, 3.634400212745, 3.738999744084, CaCl 
     5  3.842600203032, 3.938399717911, 3.975999839721, 4.000000000000, CaCl 
     6       7*0.0D+00/                                                 CaCl 
      DATA  Q_CaCl/                                                     071215
     1  1.66876689D+00, 1.78447365D+00, 1.96512205D+00, 2.50965714D+00, CaCl 
     2  2.68108246D+00, 2.84290434D+00, 2.99809172D+00, 3.25585844D+00, CaCl 
     3  3.44374933D+00, 3.67172435D+00, 4.15916064D+00, 4.70183692D+00, CaCl 
     4  5.14185183D+00, 5.36148931D+00, 5.59329265D+00, 5.82970233D+00, CaCl 
     5  6.08621286D+00, 6.34603967D+00, 6.45353276D+00, 6.52351564D+00, CaCl 
     6       7*0.0D+00/                                                 CaCl 
      DATA TQ_ScCl/                                                     071215
     1  0.699999789529, 0.828299960275, 1.034200109262, 1.673099996386, ScCl 
     2  1.816400102672, 1.974899898498, 2.118499850765, 2.351600270185, ScCl 
     3  2.502200040231, 2.664499960349, 2.941799630040, 3.235499665291, ScCl 
     4  3.338900136774, 3.437899734066, 3.535599857225, 3.631300148838, ScCl 
     5  3.733499630460, 3.834200018210, 3.935099955755, 3.974899918760, ScCl 
     6  4.000000000000,      6*0.0D+00/                                 ScCl 
      DATA  Q_ScCl/                                                     071215
     1  1.31365199D+00, 1.44012237D+00, 1.64401919D+00, 2.28044805D+00, ScCl 
     2  2.42360079D+00, 2.58247681D+00, 2.72894749D+00, 2.98504366D+00, ScCl 
     3  3.17262489D+00, 3.39900948D+00, 3.84062402D+00, 4.36749718D+00, ScCl 
     4  4.56442510D+00, 4.75865242D+00, 4.95761629D+00, 5.16380267D+00, ScCl 
     5  5.40442992D+00, 5.67086729D+00, 5.96899838D+00, 6.09343264D+00, ScCl 
     6  6.17325746D+00,      6*0.0D+00/                                 ScCl 
      DATA TQ_MnCl/                                                     071215
     1  0.699999789529, 0.828999943334, 1.047300116266, 1.621999934291, MnCl 
     2  1.782100167741, 1.944499933024, 2.258600202404, 2.383100193205, MnCl 
     3  2.518800431664, 2.718900256609, 2.965100158182, 3.206599946114, MnCl 
     4  3.557700357236, 3.684699969909, 3.809400393762, 3.925100214468, MnCl 
     5  3.970600227731, 4.000000000000,      9*0.0D+00/                 MnCl 
      DATA  Q_MnCl/                                                     071215
     1  2.15626301D+00, 2.28341017D+00, 2.49959565D+00, 3.07192301D+00, MnCl 
     2  3.23180616D+00, 3.39480872D+00, 3.72929608D+00, 3.87930912D+00, MnCl 
     3  4.05893847D+00, 4.35483925D+00, 4.76184712D+00, 5.19390663D+00, MnCl 
     4  5.85664132D+00, 6.10268246D+00, 6.34769839D+00, 6.58074052D+00, MnCl 
     5  6.67484595D+00, 6.73659810D+00,      9*0.0D+00/                 MnCl 
      DATA TQ_FeCl/                                                     071215
     1  0.699999789529, 0.761700201508, 0.857600081432, 1.113000062974, FeCl 
     2  1.397100025220, 1.716399882777, 1.858300097489, 1.995299893190, FeCl 
     3  2.322999775593, 2.468100208910, 2.616399791306, 2.893400096954, FeCl 
     4  3.098000331927, 3.251800041461, 3.409299815636, 3.629600113141, FeCl 
     5  3.862899730735, 3.946899758386, 4.000000000000,      8*0.0D+00/ FeCl 
      DATA  Q_FeCl/                                                     071215
     1  2.43275981D+00, 2.48797185D+00, 2.57541241D+00, 2.81548409D+00, FeCl 
     2  3.09032522D+00, 3.40438797D+00, 3.54503081D+00, 3.68191416D+00, FeCl 
     3  4.03205812D+00, 4.21065535D+00, 4.41440208D+00, 4.85873196D+00, FeCl 
     4  5.24437436D+00, 5.56470911D+00, 5.91124958D+00, 6.40660002D+00, FeCl 
     5  6.92431408D+00, 7.10693342D+00, 7.22107432D+00,      8*0.0D+00/ FeCl 
      DATA TQ_CuCl/                                                     071215
     1  0.699999789529, 0.841500112180, 1.066500102838, 1.658199879798, CuCl 
     2  1.825100049954, 1.978099966028, 2.233399629775, 2.382300176010, CuCl 
     3  2.522900241385, 2.648799702770, 2.971300175990, 3.243599856516, CuCl 
     4  3.353400132426, 3.457799982212, 3.557300348888, 3.656799780432, CuCl 
     5  3.753600072501, 3.848300331262, 3.939499638630, 3.976699789423, CuCl 
     6  4.000000000000,      6*0.0D+00/                                 CuCl 
      DATA  Q_CuCl/                                                     071215
     1  1.30434601D+00, 1.44381377D+00, 1.66668877D+00, 2.25613431D+00, CuCl 
     2  2.42289123D+00, 2.57659311D+00, 2.84495533D+00, 3.01920863D+00, CuCl 
     3  3.20193217D+00, 3.38172488D+00, 3.90503755D+00, 4.40129894D+00, CuCl 
     4  4.61246547D+00, 4.81923564D+00, 5.02356860D+00, 5.24004770D+00, CuCl 
     5  5.47173321D+00, 5.72957290D+00, 6.01159988D+00, 6.13479306D+00, CuCl 
     6  6.21370242D+00,      6*0.0D+00/                                 CuCl 
      DATA TQ_ZnCl/                                                     071215
     1  0.699999789529, 0.827099989317, 1.035800066512, 1.638699894504, ZnCl 
     2  1.784400110951, 1.935299938612, 2.275000098712, 2.407799796745, ZnCl 
     3  2.550800200074, 2.728099695569, 2.918700043881, 3.198399760711, ZnCl 
     4  3.451199822614, 3.589300143334, 3.720500216444, 3.885100240236, ZnCl 
     5  3.955099939741, 3.982899617413, 4.000000000000,      8*0.0D+00/ ZnCl 
      DATA  Q_ZnCl/                                                     071215
     1  1.60485350D+00, 1.73009236D+00, 1.93670546D+00, 2.53707059D+00, ZnCl 
     2  2.68257528D+00, 2.83392643D+00, 3.19588471D+00, 3.35772166D+00, ZnCl 
     3  3.55034200D+00, 3.81561999D+00, 4.12926280D+00, 4.62991492D+00, ZnCl 
     4  5.11231715D+00, 5.38634258D+00, 5.65543989D+00, 6.01102647D+00, ZnCl 
     5  6.17031198D+00, 6.23508859D+00, 6.27537285D+00,      8*0.0D+00/ ZnCl 
      DATA TQ_GaCl/                                                     071215
     1  0.699999789529, 0.866399967758, 1.185099942101, 1.635399976685, GaCl 
     2  1.791500003516, 1.939900037183, 2.167600009245, 2.301400228317, GaCl 
     3  2.435299913459, 2.565200042440, 2.918000026767, 3.211700062255, GaCl 
     4  3.472900043876, 3.595000270217, 3.705899942040, 3.878900094298, GaCl 
     5  3.951999872051, 3.981799592717, 4.000000000000,      8*0.0D+00/ GaCl 
      DATA  Q_GaCl/                                                     071215
     1  1.37358034D+00, 1.53800447D+00, 1.85452071D+00, 2.30358607D+00, GaCl 
     2  2.45961140D+00, 2.60892799D+00, 2.84829905D+00, 3.00290600D+00, GaCl 
     3  3.17319604D+00, 3.35483894D+00, 3.92280194D+00, 4.45875575D+00, GaCl 
     4  4.97092347D+00, 5.22239521D+00, 5.46057973D+00, 5.86676622D+00, GaCl 
     5  6.05895419D+00, 6.14145242D+00, 6.19297054D+00,      8*0.0D+00/ GaCl 
      DATA TQ_GeCl/                                                     071215
     1  0.699999789529, 0.754400101215, 0.839400134618, 1.059399965874, GeCl 
     2  1.671500033209, 1.883500037505, 2.078199888803, 2.354800032274, GeCl 
     3  2.554000273974, 2.689899591777, 2.824099799989, 3.063000216467, GeCl 
     4  3.303700270495, 3.581399958246, 3.777599915966, 3.909999844223, GeCl 
     5  3.965000158790, 4.000000000000,      9*0.0D+00/                 GeCl 
      DATA  Q_GeCl/                                                     071215
     1  1.69175625D+00, 1.74324912D+00, 1.82437068D+00, 2.03722381D+00, GeCl 
     2  2.64111141D+00, 2.85229576D+00, 3.04955725D+00, 3.35769760D+00, GeCl 
     3  3.62524297D+00, 3.83957194D+00, 4.07754166D+00, 4.54953835D+00, GeCl 
     4  5.05440381D+00, 5.64458944D+00, 6.06608571D+00, 6.35841084D+00, GeCl 
     5  6.48269014D+00, 6.56261873D+00,      9*0.0D+00/                 GeCl 
      DATA TQ_AsCl/                                                     071215
     1  0.699999789529, 0.841800104942, 1.067200119033, 1.667800015849, AsCl 
     2  1.834000007153, 2.000200004592, 2.314800083984, 2.442099631122, AsCl 
     3  2.580499937160, 2.779099792216, 3.175500192601, 3.437499762377, AsCl 
     4  3.666899995950, 3.793400049336, 3.930400294503, 4.000000000000, AsCl 
     5      11*0.0D+00/                                                 AsCl 
      DATA  Q_AsCl/                                                     071215
     1  1.37216097D+00, 1.51222047D+00, 1.73579659D+00, 2.33439047D+00, AsCl 
     2  2.50041785D+00, 2.66721854D+00, 3.00176572D+00, 3.15500183D+00, AsCl 
     3  3.33851976D+00, 3.63325880D+00, 4.31130095D+00, 4.80497045D+00, AsCl 
     4  5.26122075D+00, 5.52431673D+00, 5.81720546D+00, 5.96686839D+00, AsCl 
     5      11*0.0D+00/                                                 AsCl 
      DATA TQ_SeCl/                                                     071215
     1  0.699999789529, 0.755700134907, 0.842700083226, 1.069200165306, SeCl 
     2  1.390399865132, 1.748799962371, 1.939900037183, 2.124499975656, SeCl 
     3  2.270500417961, 2.409799845539, 2.641100260640, 2.935999871127, SeCl 
     4  3.139500338200, 3.362599703036, 3.622099946321, 3.826699846948, SeCl 
     5  3.931900186392, 4.000000000000,      9*0.0D+00/                 SeCl 
      DATA  Q_SeCl/                                                     071215
     1  1.99278625D+00, 2.04551399D+00, 2.12856471D+00, 2.34781187D+00, SeCl 
     2  2.66345942D+00, 3.01900596D+00, 3.20933658D+00, 3.39411942D+00, SeCl 
     3  3.54357186D+00, 3.69406191D+00, 3.97557270D+00, 4.40536058D+00, SeCl 
     4  4.74243123D+00, 5.13889154D+00, 5.62332279D+00, 6.01629911D+00, SeCl 
     5  6.22083224D+00, 6.35386397D+00,      9*0.0D+00/                 SeCl 
      DATA TQ_BrCl/                                                     071215
     1  0.699999789529, 0.846499991540, 1.080899942284, 1.676899908932, BrCl 
     2  1.847699962586, 2.004700107923, 2.263300307083, 2.413199928020, BrCl 
     3  2.565300035323, 2.709200036194, 2.963200114497, 3.223599996081, BrCl 
     4  3.330799934724, 3.432000151655, 3.523900173227, 3.602100228721, BrCl 
     5  3.665999977578, 3.731199582944, 3.800100189899, 3.901099644354, BrCl 
     6  3.960200051218, 3.984799660070, 4.000000000000,      4*0.0D+00/ BrCl 
      DATA  Q_BrCl/                                                     071215
     1  1.36626481D+00, 1.51095452D+00, 1.74347941D+00, 2.33758224D+00, BrCl 
     2  2.50826148D+00, 2.66595513D+00, 2.93773883D+00, 3.11330304D+00, BrCl 
     3  3.31204309D+00, 3.52052783D+00, 3.93234452D+00, 4.40041455D+00, BrCl 
     4  4.60399265D+00, 4.80258982D+00, 4.99231709D+00, 5.16784486D+00, BrCl 
     5  5.32689564D+00, 5.50719690D+00, 5.71659783D+00, 6.04515564D+00, BrCl 
     6  6.23970397D+00, 6.31964446D+00, 6.36853805D+00,      4*0.0D+00/ BrCl 
      DATA TQ_RbCl/                                                     071215
     1  0.699999789529, 0.827399982056, 1.027800167149, 1.494499876806, RbCl 
     2  1.630900088751, 1.759700231591, 2.084800032442, 2.265700359507, RbCl 
     3  2.435899871584, 2.690499595698, 2.920400084325, 3.191799606034, RbCl 
     4  3.458600001557, 3.549000172622, 3.649199678278, 3.744299860392, RbCl 
     5  3.835700050882, 3.936299869266, 3.974999911574, 4.000000000000, RbCl 
     6       7*0.0D+00/                                                 RbCl 
      DATA  Q_RbCl/                                                     071215
     1  1.60408491D+00, 1.73056684D+00, 1.92998004D+00, 2.39564790D+00, RbCl 
     2  2.53213890D+00, 2.66219704D+00, 3.01702712D+00, 3.24925269D+00, RbCl 
     3  3.49632007D+00, 3.91239302D+00, 4.32576435D+00, 4.84930021D+00, RbCl 
     4  5.40501343D+00, 5.60668023D+00, 5.83860110D+00, 6.06410144D+00, RbCl 
     5  6.28112401D+00, 6.51371519D+00, 6.59997099D+00, 6.65447410D+00, RbCl 
     6       7*0.0D+00/                                                 RbCl 
      DATA TQ_SrCl/                                                     071215
     1  0.699999789529, 0.840500136308, 1.067800132915, 1.566300086345, SrCl 
     2  1.716599877564, 1.857200072256, 2.086600069384, 2.222000108906, SrCl 
     3  2.375100009645, 2.520600412183, 2.864699760227, 3.184299977222, SrCl 
     4  3.420200096931, 3.542700021936, 3.650099624097, 3.756800148340, SrCl 
     5  3.847500313265, 3.940399611624, 3.976899775052, 4.000000000000, SrCl 
     6       7*0.0D+00/                                                 SrCl 
      DATA  Q_SrCl/                                                     071215
     1  1.90221702D+00, 2.04170507D+00, 2.26792570D+00, 2.76535821D+00, SrCl 
     2  2.91561685D+00, 3.05708467D+00, 3.29809359D+00, 3.45450002D+00, SrCl 
     3  3.65050992D+00, 3.85734699D+00, 4.41675015D+00, 4.99919204D+00, SrCl 
     4  5.45473064D+00, 5.70062729D+00, 5.92661229D+00, 6.16847815D+00, SrCl 
     5  6.39349898D+00, 6.64467581D+00, 6.74859366D+00, 6.81562508D+00, SrCl 
     6       7*0.0D+00/                                                 SrCl 
      DATA TQ_YCl/                                                      071215
     1  0.699999789529, 0.883600026260, 1.202699858035, 1.648900072042, YCl  
     2  1.807899992400, 1.956799974701, 2.185699897856, 2.318799796790, YCl  
     3  2.453499881128, 2.584400028533, 2.918700043881, 3.201299827342, YCl  
     4  3.327499860563, 3.446999731148, 3.567599859029, 3.671700099694, YCl  
     5  3.761200253191, 3.844800252524, 3.934899970170, 4.000000000000, YCl  
     6       7*0.0D+00/                                                 YCl  
      DATA  Q_YCl/                                                      071215
     1  1.48295340D+00, 1.66490139D+00, 1.98239146D+00, 2.42773788D+00, YCl  
     2  2.58674115D+00, 2.73661839D+00, 2.97732924D+00, 3.13115473D+00, YCl  
     3  3.30247546D+00, 3.48561539D+00, 4.02167716D+00, 4.53353791D+00, YCl  
     4  4.77510785D+00, 5.01077967D+00, 5.25708349D+00, 5.48041740D+00, YCl  
     5  5.68519434D+00, 5.89081085D+00, 6.12804239D+00, 6.30611403D+00, YCl  
     6       7*0.0D+00/                                                 YCl  
      DATA TQ_AgCl/                                                     071215
     1  0.699999789529, 0.876500030268, 1.181899875473, 1.612200085110, AgCl 
     2  1.768300014526, 1.914299936307, 2.140600278937, 2.274300148373, AgCl 
     3  2.407199782107, 2.536299867366, 2.864699760227, 3.140700297425, AgCl 
     4  3.402199672046, 3.537299895912, 3.653899712765, 3.749599978285, AgCl 
     5  3.840100146790, 3.935699912511, 4.000000000000,      8*0.0D+00/ AgCl 
      DATA  Q_AgCl/                                                     071215
     1  1.45833046D+00, 1.63313097D+00, 1.93683329D+00, 2.36614398D+00, AgCl 
     2  2.52219888D+00, 2.66913862D+00, 2.90706725D+00, 3.06151652D+00, AgCl 
     3  3.23040549D+00, 3.41067545D+00, 3.93557997D+00, 4.43356998D+00, AgCl 
     4  4.94030542D+00, 5.21519946D+00, 5.46207607D+00, 5.67367875D+00, AgCl 
     5  5.88255887D+00, 6.11196352D+00, 6.26976606D+00,      8*0.0D+00/ AgCl 
      DATA TQ_CdCl/                                                     071215
     1  0.699999789529, 0.817800114565, 1.006800145135, 1.588000117947, CdCl 
     2  1.846799984302, 1.990199777291, 2.213000113288, 2.348500355402, CdCl 
     3  2.495399890026, 2.845700291193, 3.155799735557, 3.489299741162, CdCl 
     4  3.754400091461, 3.900199624142, 3.961400078111, 4.000000000000, CdCl 
     5      11*0.0D+00/                                                 CdCl 
      DATA  Q_CdCl/                                                     071215
     1  1.75825208D+00, 1.87483548D+00, 2.06245829D+00, 2.64179849D+00, CdCl 
     2  2.90081530D+00, 3.04703945D+00, 3.29093260D+00, 3.45728669D+00, CdCl 
     3  3.65700593D+00, 4.20976881D+00, 4.76538449D+00, 5.40612212D+00, CdCl 
     4  5.94361425D+00, 6.25490684D+00, 6.38939787D+00, 6.47524972D+00, CdCl 
     5      11*0.0D+00/                                                 CdCl 
      DATA TQ_InCl/                                                     071215
     1  0.699999789529, 0.838700118571, 1.061899996410, 1.578199962945, InCl 
     2  1.730600006464, 1.875900029982, 2.108899756285, 2.244799891759, InCl 
     3  2.392100187145, 2.536199864934, 2.879000083980, 3.197199732588, InCl 
     4  3.438099719911, 3.547900146312, 3.655099740765, 3.756100131751, InCl 
     5  3.861199691891, 3.943699686134, 3.978299674457, 4.000000000000, InCl 
     6       7*0.0D+00/                                                 InCl 
      DATA  Q_InCl/                                                     071215
     1  1.50992258D+00, 1.64738947D+00, 1.86927963D+00, 2.38428403D+00, InCl 
     2  2.53664046D+00, 2.68282446D+00, 2.92776988D+00, 3.08512919D+00, InCl 
     3  3.27394640D+00, 3.47876777D+00, 4.03746586D+00, 4.62267196D+00, InCl 
     4  5.09726604D+00, 5.32319033D+00, 5.55209805D+00, 5.77965545D+00, InCl 
     5  6.03705946D+00, 6.26032462D+00, 6.36006672D+00, 6.42432029D+00, InCl 
     6       7*0.0D+00/                                                 InCl 
      DATA TQ_SnCl/                                                     071215
     1  0.699999789529, 0.764200135568, 0.863600041576, 1.144999987011, SnCl 
     2  1.654699972927, 1.800200173550, 1.940200034602, 2.147599789897, SnCl 
     3  2.287899904112, 2.442599642001, 2.598300337609, 2.773600204855, SnCl 
     4  2.973000051951, 3.077599884847, 3.184199984475, 3.336000064435, SnCl 
     5  3.491599794028, 3.629900119814, 3.805400306079, 3.910399854500, SnCl 
     6  4.000000000000,      6*0.0D+00/                                 SnCl 
      DATA  Q_SnCl/                                                     071215
     1  2.11525963D+00, 2.17693087D+00, 2.27307920D+00, 2.54843843D+00, SnCl 
     2  3.05355767D+00, 3.19861928D+00, 3.33943342D+00, 3.55745836D+00, SnCl 
     3  3.71930039D+00, 3.91765622D+00, 4.14037750D+00, 4.41802998D+00, SnCl 
     4  4.76962796D+00, 4.97106291D+00, 5.18836530D+00, 5.51501131D+00, SnCl 
     5  5.86346367D+00, 6.18024909D+00, 6.59452146D+00, 6.85242678D+00, SnCl 
     6  7.07758541D+00,      6*0.0D+00/                                 SnCl 
      DATA TQ_SbCl/                                                     071215
     1  0.699999789529, 0.865199999394, 1.146200013464, 1.627800062524, SbCl 
     2  1.839300131896, 2.031899772829, 2.326999862207, 2.511600281215, SbCl 
     3  2.610799662396, 2.702299864830, 3.005400120500, 3.135200239656, SbCl 
     4  3.274000173927, 3.414299942940, 3.561500297432, 3.765700361984, SbCl 
     5  3.910699862208, 3.964900156549, 4.000000000000,      8*0.0D+00/ SbCl 
      DATA  Q_SbCl/                                                     071215
     1  1.97457225D+00, 2.13829708D+00, 2.41777958D+00, 2.89826408D+00, SbCl 
     2  3.10973243D+00, 3.30499391D+00, 3.63478184D+00, 3.88538462D+00, SbCl 
     3  4.04029388D+00, 4.19618329D+00, 4.78058573D+00, 5.04835727D+00, SbCl 
     4  5.33821381D+00, 5.63202763D+00, 5.94029303D+00, 6.37045995D+00, SbCl 
     5  6.68250186D+00, 6.80176012D+00, 6.87995874D+00,      8*0.0D+00/ SbCl 
      DATA TQ_ICl/                                                      071215
     1  0.699999789529, 0.847099977063, 1.084700035744, 1.637499924388, ICl  
     2  1.803400098266, 1.955599944418, 2.193299655628, 2.338200121544, ICl  
     3  2.532899784672, 2.724999919775, 2.962400096104, 3.184599955466, ICl  
     4  3.282699806057, 3.384600220961, 3.462400087335, 3.552600250801, ICl  
     5  3.626500044189, 3.699099789085, 3.797600137390, 3.870099895224, ICl  
     6  3.949999828380, 3.980999574756, 4.000000000000,      4*0.0D+00/ ICl  
      DATA  Q_ICl/                                                      071215
     1  1.49026767D+00, 1.63600921D+00, 1.87219575D+00, 2.42364962D+00, ICl  
     2  2.58949094D+00, 2.74255572D+00, 2.99245264D+00, 3.16099433D+00, ICl  
     3  3.41704465D+00, 3.70481876D+00, 4.10206960D+00, 4.50616122D+00, ICl  
     4  4.69285377D+00, 4.89366237D+00, 5.05629454D+00, 5.26642487D+00, ICl  
     5  5.46621980D+00, 5.69018268D+00, 6.02538592D+00, 6.27908253D+00, ICl  
     6  6.55143698D+00, 6.65319998D+00, 6.71426483D+00,      4*0.0D+00/ ICl  
      DATA TQ_CsCl/                                                     071215
     1  0.699999789529, 1.039699962310, 1.474599887931, 1.610600125117, CsCl 
     2  1.738899806871, 2.057500371042, 2.245899915116, 2.421800127204, CsCl 
     3  2.689899591777, 2.943299663057, 3.101800250011, 3.253700082930, CsCl 
     4  3.497599934061, 3.588500124591, 3.685499914191, 3.853100150946, CsCl 
     5  3.942999670329, 3.977799710384, 4.000000000000,      8*0.0D+00/ CsCl 
      DATA  Q_CsCl/                                                     071215
     1  1.68815345D+00, 2.02625286D+00, 2.46037474D+00, 2.59650746D+00, CsCl 
     2  2.72617022D+00, 3.07422572D+00, 3.31655500D+00, 3.57331244D+00, CsCl 
     3  4.01476433D+00, 4.47454970D+00, 4.77833500D+00, 5.08035528D+00, CsCl 
     4  5.59459993D+00, 5.79960002D+00, 6.02697715D+00, 6.42999391D+00, CsCl 
     5  6.64152455D+00, 6.72086606D+00, 6.77053392D+00,      8*0.0D+00/ CsCl 
      DATA TQ_BaCl/                                                     071215
     1  0.699999789529, 0.842200095290, 1.073500091532, 1.544100060572, BaCl 
     2  1.693899971277, 1.834000007153, 2.183700041267, 2.361199674740, BaCl 
     3  2.539299940331, 2.810700365138, 3.131900164029, 3.247499945277, BaCl 
     4  3.360699652657, 3.472500072550, 3.586000066019, 3.685899886332, BaCl 
     5  3.792800036757, 3.900499630879, 4.000000000000,      8*0.0D+00/ BaCl 
      DATA  Q_BaCl/                                                     071215
     1  1.98870124D+00, 2.13006506D+00, 2.36047382D+00, 2.83024058D+00, BaCl 
     2  2.98004523D+00, 3.12117750D+00, 3.50252701D+00, 3.73118792D+00, BaCl 
     3  3.99040938D+00, 4.43547287D+00, 5.01811487D+00, 5.23840028D+00, BaCl 
     4  5.45957021D+00, 5.68641776D+00, 5.93257211D+00, 6.17020577D+00, BaCl 
     5  6.45114248D+00, 6.75893779D+00, 7.05574945D+00,      8*0.0D+00/ BaCl 
      DATA TQ_YbCl/                                                     071215
     1  0.699999789529, 0.844200047034, 1.077699980792, 1.569800168764, YbCl 
     2  1.714199940123, 1.850599920856, 2.201599832270, 2.354500054578, YbCl 
     3  2.502000034910, 2.830599955919, 3.142900133261, 3.278699837144, YbCl 
     4  3.406099750920, 3.544300060205, 3.666599989826, 3.762400282203, YbCl 
     5  3.849600360508, 3.940699618398, 3.977099760681, 4.000000000000, YbCl 
     6       7*0.0D+00/                                                 YbCl 
      DATA  Q_YbCl/                                                     071215
     1  1.98870124D+00, 2.13205509D+00, 2.36466152D+00, 2.85591672D+00, YbCl 
     2  3.00033007D+00, 3.13769986D+00, 3.51983944D+00, 3.71466303D+00, YbCl 
     3  3.92362236D+00, 4.45550815D+00, 5.02295936D+00, 5.28301959D+00, YbCl 
     4  5.53362987D+00, 5.81456575D+00, 6.07563503D+00, 6.29272032D+00, YbCl 
     5  6.50255545D+00, 6.73436722D+00, 6.82993447D+00, 6.89067011D+00, YbCl 
     6       7*0.0D+00/                                                 YbCl 
      DATA TQ_AuCl/                                                     071215
     1  0.699999789529, 1.047600122923, 1.634599996608, 1.912899902079, AuCl 
     2  2.274000169656, 2.409699843099, 2.555900317852, 2.734399653298, AuCl 
     3  2.913999928973, 3.216200167996, 3.365099769324, 3.499799985406, AuCl 
     4  3.614599768598, 3.714800138925, 3.803000253469, 3.881400152413, AuCl 
     5  3.954099917906, 3.982599610678, 4.000000000000,      8*0.0D+00/ AuCl 
      DATA  Q_AuCl/                                                     071215
     1  1.68767125D+00, 2.03362171D+00, 2.61962642D+00, 2.89830341D+00, AuCl 
     2  3.28309524D+00, 3.44951823D+00, 3.64805697D+00, 3.91722242D+00, AuCl 
     3  4.21374453D+00, 4.75602596D+00, 5.03790968D+00, 5.29979042D+00, AuCl 
     4  5.52904257D+00, 5.73654643D+00, 5.92845687D+00, 6.10867285D+00, AuCl 
     5  6.28426627D+00, 6.35505856D+00, 6.39872747D+00,      8*0.0D+00/ AuCl 
      DATA TQ_HgCl/                                                     071215
     1  0.699999789529, 0.843200071162, 1.077699980792, 1.554300076848, HgCl 
     2  1.702199866082, 1.841700107354, 2.072999758444, 2.214700146320, HgCl 
     3  2.347800339657, 2.469200233720, 2.935499906215, 3.193899655250, HgCl 
     4  3.414599950834, 3.571199715211, 3.762200277367, 3.856899883033, HgCl 
     5  3.943099672587, 4.000000000000,      9*0.0D+00/                 HgCl 
      DATA  Q_HgCl/                                                     071215
     1  1.98870124D+00, 2.13106007D+00, 2.36466152D+00, 2.84042972D+00, HgCl 
     2  2.98831648D+00, 3.12870658D+00, 3.37190378D+00, 3.53627596D+00, HgCl 
     3  3.70662378D+00, 3.87675140D+00, 4.64436483D+00, 5.12632288D+00, HgCl 
     4  5.56307943D+00, 5.88905714D+00, 6.29792041D+00, 6.49978972D+00, HgCl 
     5  6.68237976D+00, 6.80324245D+00,      9*0.0D+00/                 HgCl 
      DATA TQ_TlCl/                                                     071215
     1  0.699999789529, 0.837600093354, 1.058399988246, 1.557599996457, TlCl 
     2  1.701799856671, 1.837700094238, 2.187099797468, 2.338600131026, TlCl 
     3  2.485299661965, 2.801900244528, 3.104800033167, 3.387000281737, TlCl 
     4  3.526399991488, 3.653399701098, 3.743199835924, 3.832499981182, TlCl 
     5  3.920600114490, 4.000000000000,      9*0.0D+00/                 TlCl 
      DATA  Q_TlCl/                                                     071215
     1  1.58583029D+00, 1.72240494D+00, 1.94211581D+00, 2.44025378D+00, TlCl 
     2  2.58447516D+00, 2.72138624D+00, 3.10206548D+00, 3.29520589D+00, TlCl 
     3  3.50311406D+00, 4.01510750D+00, 4.56558486D+00, 5.11655910D+00, TlCl 
     4  5.40226848D+00, 5.67356159D+00, 5.87409119D+00, 6.08196750D+00, TlCl 
     5  6.29474042D+00, 6.49047534D+00,      9*0.0D+00/                 TlCl 
      DATA TQ_PbCl/                                                     071215
     1  0.699999789529, 0.845800008430, 1.094100061228, 1.585000042491, PbCl 
     2  1.729900018573, 1.871399918657, 2.210500064712, 2.354800032274, PbCl 
     3  2.499699975713, 2.775300077312, 3.074699812221, 3.199299781804, PbCl 
     4  3.313000195784, 3.598300343135, 3.699599799432, 3.792600032564, PbCl 
     5  4.000000000000,     10*0.0D+00/                                 PbCl 
      DATA  Q_PbCl/                                                     071215
     1  2.20010875D+00, 2.34161785D+00, 2.58522541D+00, 3.07200423D+00, PbCl 
     2  3.21647230D+00, 3.35870615D+00, 3.72696535D+00, 3.90912709D+00, PbCl 
     3  4.11186696D+00, 4.54797715D+00, 5.07833005D+00, 5.31137722D+00, PbCl 
     4  5.52966224D+00, 6.11097462D+00, 6.33381917D+00, 6.54682701D+00, PbCl 
     5  7.04138254D+00,     10*0.0D+00/                                 PbCl 
      DATA TQ_AlSe/                                                     071215
     1  0.699999789529, 0.855200026377, 1.107400080908, 1.681699877374, AlSe 
     2  1.855200026377, 2.026799962966, 2.337600107322, 2.464700132225, AlSe 
     3  2.601200286434, 2.790699993906, 3.010300230024, 3.225699843709, AlSe 
     4  3.463200104644, 3.587600103505, 3.699299793224, 3.867599838124, AlSe 
     5  3.948699799027, 4.000000000000,      9*0.0D+00/                 AlSe 
      DATA  Q_AlSe/                                                     071215
     1  1.88778319D+00, 2.04184538D+00, 2.29287901D+00, 2.86609309D+00, AlSe 
     2  3.03950235D+00, 3.21180013D+00, 3.54247109D+00, 3.69539105D+00, AlSe 
     3  3.87615337D+00, 4.15614134D+00, 4.51752562D+00, 4.90167159D+00, AlSe 
     4  5.35052343D+00, 5.59494300D+00, 5.82094856D+00, 6.17775272D+00, AlSe 
     5  6.35758086D+00, 6.47355958D+00,      9*0.0D+00/                 AlSe 
      DATA TQ_SiSe/                                                     071215
     1  0.699999789529, 0.844700034970, 1.088300124286, 1.750199992222, SiSe 
     2  1.934999932183, 2.116999820256, 2.427400233848, 2.550700197764, SiSe 
     3  2.686499837866, 2.880400114985, 3.098900353612, 3.318999769065, SiSe 
     4  3.531599766196, 3.639400315821, 3.740199769192, 3.892900142845, SiSe 
     5  3.957699996514, 3.983899639864, 4.000000000000,      8*0.0D+00/ SiSe 
      DATA  Q_SiSe/                                                     071215
     1  1.26754616D+00, 1.40998537D+00, 1.65113758D+00, 2.31054877D+00, SiSe 
     2  2.49517255D+00, 2.67782125D+00, 3.00760574D+00, 3.15541283D+00, SiSe 
     3  3.33437092D+00, 3.62020755D+00, 3.97999101D+00, 4.37382264D+00, SiSe 
     4  4.77718902D+00, 4.98999166D+00, 5.19708034D+00, 5.54472211D+00, SiSe 
     5  5.71361243D+00, 5.78633677D+00, 5.83229661D+00,      8*0.0D+00/ SiSe 
      DATA TQ_GeSe/                                                     071215
     1  0.699999789529, 0.839300132325, 1.060799970959, 1.652800023483, GeSe 
     2  1.800900157081, 1.952799873759, 2.296600111038, 2.430500248465, GeSe 
     3  2.574199795158, 2.751500024114, 2.939699611472, 3.224599923523, GeSe 
     4  3.472600065382, 3.602700185299, 3.722300091715, 3.886100263973, GeSe 
     5  3.955199941925, 3.982999619658, 4.000000000000,      8*0.0D+00/ GeSe 
      DATA  Q_GeSe/                                                     071215
     1  1.56286110D+00, 1.70106546D+00, 1.92141176D+00, 2.51216414D+00, GeSe 
     2  2.66020316D+00, 2.81265294D+00, 3.17927388D+00, 3.34287652D+00, GeSe 
     3  3.53712119D+00, 3.80343275D+00, 4.11426178D+00, 4.62696877D+00, GeSe 
     4  5.10421739D+00, 5.36509612D+00, 5.61326707D+00, 5.97549393D+00, GeSe 
     5  6.14018619D+00, 6.20883196D+00, 6.25149602D+00,      8*0.0D+00/ GeSe 
      DATA TQ_KBr/                                                      071215
     1  0.699999789529, 1.030600205448, 1.471599817683, 1.606900062990, KBr  
     2  1.734499912679, 2.054900311166, 2.241499821687, 2.417000020173, KBr  
     3  2.675100197379, 2.906599754720, 3.185499890196, 3.458299994303, KBr  
     4  3.548700165446, 3.649399664149, 3.743999853719, 3.835700050882, KBr  
     5  3.936199876474, 3.974999911574, 4.000000000000,      8*0.0D+00/ KBr  
      DATA  Q_KBr/                                                      071215
     1  1.63681329D+00, 1.96564083D+00, 2.40572971D+00, 2.54114356D+00, KBr  
     2  2.67006599D+00, 3.01998314D+00, 3.25989829D+00, 3.51586047D+00, KBr  
     3  3.93982970D+00, 4.35758977D+00, 4.89697073D+00, 5.46632307D+00, KBr  
     4  5.66829459D+00, 5.90160772D+00, 6.12605153D+00, 6.34372204D+00, KBr  
     5  6.57557821D+00, 6.66170959D+00, 6.71594235D+00,      8*0.0D+00/ KBr  
      DATA TQ_SiTe/                                                     071215
     1  0.699999789529, 0.862300075849, 1.129100100075, 1.685799973331, SiTe 
     2  1.869799889071, 2.047100131936, 2.347600335158, 2.467000184100, SiTe 
     3  2.602900162323, 2.809600406701, 3.070699712047, 3.296700118051, SiTe 
     4  3.582199976990, 3.686199865437, 3.784099838953, 3.916600013801, SiTe 
     5  3.967100205852, 3.986899707217, 4.000000000000,      8*0.0D+00/ SiTe 
      DATA  Q_SiTe/                                                     071215
     1  1.63689330D+00, 1.79814254D+00, 2.06387841D+00, 2.61967280D+00, SiTe 
     2  2.80359226D+00, 2.98167798D+00, 3.30138295D+00, 3.44419944D+00, SiTe 
     3  3.62252966D+00, 3.92686555D+00, 4.36010260D+00, 4.76725639D+00, SiTe 
     4  5.31046143D+00, 5.51494290D+00, 5.71178184D+00, 5.99139398D+00, SiTe 
     5  6.10460807D+00, 6.15026325D+00, 6.18088634D+00,      8*0.0D+00/ SiTe 
      DATA TQ_GeTe/                                                     071215
     1  0.699999789529, 0.850299913974, 1.098299953800, 1.589200148130, GeTe 
     2  1.740499790758, 1.883800029456, 2.118899858900, 2.256200146526, GeTe 
     3  2.385800251236, 2.509700239800, 2.894700000311, 3.201999843029, GeTe 
     4  3.368099848869, 3.521800325887, 3.638400295206, 3.744799871514, GeTe 
     5  3.898799708078, 3.960100048977, 3.984699657824, 4.000000000000, GeTe 
     6       7*0.0D+00/                                                 GeTe 
      DATA  Q_GeTe/                                                     071215
     1  1.73013766D+00, 1.87964710D+00, 2.12682105D+00, 2.61703403D+00, GeTe 
     2  2.76834991D+00, 2.91253798D+00, 3.15962908D+00, 3.31865157D+00, GeTe 
     3  3.48375028D+00, 3.65653699D+00, 4.27717892D+00, 4.83847660D+00, GeTe 
     4  5.15831854D+00, 5.46313138D+00, 5.70071068D+00, 5.92465729D+00, GeTe 
     5  6.27262195D+00, 6.42334404D+00, 6.48605363D+00, 6.52568038D+00, GeTe 
     6       7*0.0D+00/                                                 GeTe 
      DATA TQ_KI/                                                       071215
     1  0.699999789529, 1.028400181967, 1.426099918822, 1.560799956830, KI   
     2  1.687200006096, 1.999899997727, 2.198299755791, 2.384500223295, KI   
     3  2.650999640336, 2.909799826681, 3.075999844778, 3.233799624002, KI   
     4  3.489599747815, 3.581599962932, 3.679100275831, 3.847100304266, KI   
     5  3.940599616140, 3.976799782237, 4.000000000000,      8*0.0D+00/ KI   
      DATA  Q_KI/                                                       071215
     1  1.76107483D+00, 2.08815889D+00, 2.48523536D+00, 2.62011563D+00, KI   
     2  2.74801578D+00, 3.09040215D+00, 3.34689052D+00, 3.62153665D+00, KI   
     3  4.06430735D+00, 4.53639685D+00, 4.85626795D+00, 5.17118261D+00, KI   
     4  5.71227976D+00, 5.92017894D+00, 6.14864085D+00, 6.55147403D+00, KI   
     5  6.77076754D+00, 6.85305623D+00, 6.90482395D+00,      8*0.0D+00/ KI   
C
C Molecular equilibrium constants
C
      DATA TK_H2/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, H2   
     2  0.829499931233, 0.922100034988, 1.027700164680, 1.108600108571, H2   
     3  1.188500012893, 1.275700118861, 1.356300055324, 1.453899893643, H2   
     4  1.558999962352, 1.647200031946, 1.760300231031, 1.896700009654, H2   
     5  1.990799790926, 2.091200169343, 2.233299627350, 2.399399650552, H2   
     6  2.597300316178, 2.824499809645, 3.038399943696, 3.452699858887, H2   
     7  3.579299908723, 3.707899987575, 3.878600087511, 3.949499817090, H2   
     8  4.000000000000,     61*0.0D+00/                                 H2   
      DATA  K_H2/                                                       071215
     1 -9.16636527D-05, 1.51781634D-01, 3.05570688D-01, 7.76236213D-01, H2   
     2  1.85053975D+00, 2.92597870D+00, 3.95647613D+00, 4.63198308D+00, H2   
     3  5.22008271D+00, 5.78748199D+00, 6.25376023D+00, 6.75325121D+00, H2   
     4  7.21596819D+00, 7.54425819D+00, 7.89064668D+00, 8.22664256D+00, H2   
     5  8.42959955D+00, 8.63460352D+00, 8.91556173D+00, 9.23256049D+00, H2   
     6  9.58846069D+00, 9.97029263D+00, 1.03100736D+01, 1.08778722D+01, H2   
     7  1.10078948D+01, 1.11120526D+01, 1.12145557D+01, 1.12512419D+01, H2   
     8  1.12771036D+01,     61*0.0D+00/                                 H2   
      DATA TK_Li2/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720099799225, 0.752300046791, Li2  
     2  0.834300017705, 0.931499873728, 1.041399985358, 1.119199902039, Li2  
     3  1.197999849065, 1.347599975212, 1.509500215547, 1.673299991784, Li2  
     4  1.772400028776, 1.872099935974, 2.042800033501, 2.371699928509, Li2  
     5  2.770500437433, 2.920100077830, 3.065500030056, 3.253800085113, Li2  
     6  3.340400173210, 3.432200137499, 3.528099867906, 3.649999621764, Li2  
     7  3.740699780314, 3.817399872675, 3.897799781767, 3.931700200807, Li2  
     8  3.958700018349, 3.984299648844, 3.991699814756, 4.000000000000, Li2  
     9      58*0.0D+00/                                                 Li2  
      DATA  K_Li2/                                                      071215
     1 -5.24680115D-05, 1.50577231D-01, 3.08747885D-01, 7.77420515D-01, Li2  
     2  1.84020655D+00, 2.89285458D+00, 3.86597913D+00, 4.44304139D+00, Li2  
     3  4.95059842D+00, 5.74526068D+00, 6.41956567D+00, 6.96520972D+00, Li2  
     4  7.24727854D+00, 7.50314118D+00, 7.88913459D+00, 8.47754850D+00, Li2  
     5  8.95515479D+00, 9.08283843D+00, 9.18525674D+00, 9.28336719D+00, Li2  
     6  9.31096853D+00, 9.32688118D+00, 9.33474602D+00, 9.35246617D+00, Li2  
     7  9.39045996D+00, 9.45731869D+00, 9.59229491D+00, 9.67938624D+00, Li2  
     8  9.76437566D+00, 9.85830725D+00, 9.88792329D+00, 9.92230526D+00, Li2  
     9      58*0.0D+00/                                                 Li2  
      DATA TK_B2/                                                       071215
     1  0.699999789529, 0.710200044546, 0.722299851912, 0.757400178965, B2   
     2  0.793800038305, 0.846299996366, 0.901700053075, 0.961400019004, B2   
     3  1.084800038204, 1.229199998339, 1.393499939202, 1.552000132877, B2   
     4  1.710400039175, 1.917500014543, 2.141600209074, 2.322399762601, B2   
     5  2.483999628649, 2.691999629173, 2.870099895084, 3.044000061136, B2   
     6  3.234599643432, 3.445499699923, 3.659399841099, 3.729199613588, B2   
     7  3.802200235933, 3.864099758153, 3.921800141150, 3.969700264120, B2   
     8  3.987699725177, 4.000000000000,     60*0.0D+00/                 B2   
      DATA  K_B2/                                                       071215
     1 -4.42612809D-05, 1.62701996D-01, 3.51582662D-01, 8.75325522D-01, B2   
     2  1.38310859D+00, 2.05868924D+00, 2.70732775D+00, 3.34189656D+00, B2   
     3  4.47540149D+00, 5.55328541D+00, 6.52422265D+00, 7.26343611D+00, B2   
     4  7.85837558D+00, 8.48304396D+00, 9.02935847D+00, 9.40705574D+00, B2   
     5  9.71104144D+00, 1.00590308D+01, 1.03141491D+01, 1.05239547D+01, B2   
     6  1.07137400D+01, 1.08832242D+01, 1.10180291D+01, 1.10539804D+01, B2   
     7  1.10876687D+01, 1.11143310D+01, 1.11400856D+01, 1.11654608D+01, B2   
     8  1.11769022D+01, 1.11856090D+01,     60*0.0D+00/                 B2   
      DATA TK_C2/                                                       071215
     1  0.699999789529, 0.710600034435, 0.723899890230, 0.761200214697, C2   
     2  0.799600169226, 0.855700037847, 0.913199933390, 0.975499906389, C2   
     3  1.039299972998, 1.109200122403, 1.182899896294, 1.258900195233, C2   
     4  1.399600084954, 1.647800046098, 1.850699923150, 1.979800001903, C2   
     5  2.099700383101, 2.185199933709, 2.273800183845, 2.369999887941, C2   
     6  2.473400005050, 2.598400339752, 2.706899979073, 2.838500133762, C2   
     7  2.977999687131, 3.212900090452, 3.326599841049, 3.443399656209, C2   
     8  3.673000130637, 3.791600011598, 3.893900069155, 3.956699974678, C2   
     9  4.000000000000,     57*0.0D+00/                                 C2   
      DATA  K_C2/                                                       071215
     1  9.29252838D-05, 1.74382043D-01, 3.87897203D-01, 9.58097123D-01, C2   
     2  1.50450893D+00, 2.23778009D+00, 2.92048841D+00, 3.59385949D+00, C2   
     3  4.22293834D+00, 4.85239081D+00, 5.45792821D+00, 6.02797134D+00, C2   
     4  6.95692797D+00, 8.25731297D+00, 9.06237138D+00, 9.48300404D+00, C2   
     5  9.82279287D+00, 1.00371162D+01, 1.02324960D+01, 1.04097974D+01, C2   
     6  1.05590005D+01, 1.06941825D+01, 1.07890651D+01, 1.08940223D+01, C2   
     7  1.10016269D+01, 1.11731712D+01, 1.12477083D+01, 1.13160737D+01, C2   
     8  1.14288505D+01, 1.14808728D+01, 1.15237751D+01, 1.15487025D+01, C2   
     9  1.15651126D+01,     57*0.0D+00/                                 C2   
      DATA TK_N2/                                                       071215
     1  0.699999789529, 0.709300031396, 0.718799827163, 0.749299971126, N2   
     2  0.826400006258, 0.918300047565, 1.022000023918, 1.091800120057, N2   
     3  1.162099889084, 1.298600139423, 1.451399839308, 1.605600030645, N2   
     4  1.723099860844, 1.844700034970, 1.961400019004, 2.084200020128, N2   
     5  2.207899999315, 2.342400218195, 2.475799830899, 2.618699844251, N2   
     6  2.820599715501, 3.016100362951, 3.212000069304, 3.373299966973, N2   
     7  3.522300289539, 3.648799706535, 3.772300299502, 3.841000167037, N2   
     8  3.903999709479, 3.962100093798, 3.985399673540, 4.000000000000, N2   
     9      58*0.0D+00/                                                 N2   
      DATA  K_N2/                                                       071215
     1  3.60413895D-04, 1.78079341D-01, 3.56054520D-01, 9.03522469D-01, N2   
     2  2.13939765D+00, 3.37482944D+00, 4.51811950D+00, 5.16620116D+00, N2   
     3  5.73730888D+00, 6.65674529D+00, 7.45990595D+00, 8.09799411D+00, N2   
     4  8.50044645D+00, 8.86064110D+00, 9.16533794D+00, 9.45310709D+00, N2   
     5  9.71719620D+00, 9.98239864D+00, 1.02286875D+01, 1.04788800D+01, N2   
     6  1.08133044D+01, 1.11123639D+01, 1.13766794D+01, 1.15631676D+01, N2   
     7  1.17117863D+01, 1.18252794D+01, 1.19364821D+01, 1.20053527D+01, N2   
     8  1.20765591D+01, 1.21502115D+01, 1.21817729D+01, 1.22020452D+01, N2   
     9      58*0.0D+00/                                                 N2   
      DATA TK_O2/                                                       071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749799982595, O2   
     2  0.827599977216, 0.919600076668, 1.025400107881, 1.095300030534, O2   
     3  1.166099979664, 1.302400119727, 1.452399861042, 1.607900087871, O2   
     4  1.721499823732, 1.832299967141, 1.938400005040, 2.039499957270, O2   
     5  2.230199552170, 2.356399913318, 2.500399992335, 2.791800018816, O2   
     6  2.929500281363, 3.327699864899, 3.547900146312, 3.765300352314, O2   
     7  3.912999921304, 3.965300165513, 3.986399695991, 4.000000000000, O2   
     8      62*0.0D+00/                                                 O2   
      DATA  K_O2/                                                       071215
     1 -1.13492348D-04, 1.70648784D-01, 3.41699942D-01, 8.68123153D-01, O2   
     2  2.05871725D+00, 3.24373949D+00, 4.36268341D+00, 4.98512965D+00, O2   
     3  5.53720410D+00, 6.41995375D+00, 7.18242710D+00, 7.80943239D+00, O2   
     4  8.19464666D+00, 8.53091613D+00, 8.82911571D+00, 9.09922603D+00, O2   
     5  9.58238619D+00, 9.88401595D+00, 1.02086358D+01, 1.07926798D+01, O2   
     6  1.10285764D+01, 1.15509214D+01, 1.17436511D+01, 1.18724308D+01, O2   
     7  1.19308960D+01, 1.19450591D+01, 1.19495605D+01, 1.19520617D+01, O2   
     8      62*0.0D+00/                                                 O2   
      DATA TK_F2/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750399997549, F2   
     2  0.829199938494, 0.922800018110, 1.029200201723, 1.104000002529, F2   
     3  1.180299842159, 1.251500050205, 1.327399994117, 1.480000014377, F2   
     4  1.633500024002, 1.756500150962, 1.878400091830, 2.016500382612, F2   
     5  2.145799915650, 2.241199815317, 2.338800135767, 2.517800410768, F2   
     6  2.682200149096, 2.848500351940, 3.004200093722, 3.144799991483, F2   
     7  3.246999933897, 3.346900319413, 3.463100102480, 3.577599868109, F2   
     8  3.693699677341, 3.848300331262, 3.938399717911, 3.975999839721, F2   
     9  4.000000000000,     57*0.0D+00/                                 F2   
      DATA  K_F2/                                                       071215
     1  3.02356741D-05, 1.77388119D-01, 3.56684336D-01, 9.01702147D-01, F2   
     2  2.13391877D+00, 3.35997720D+00, 4.50042183D+00, 5.17161143D+00, F2   
     3  5.76523681D+00, 6.24928406D+00, 6.70267887D+00, 7.45835951D+00, F2   
     4  8.06112052D+00, 8.46232823D+00, 8.80732030D+00, 9.15199307D+00, F2   
     5  9.44465024D+00, 9.64923474D+00, 9.85263476D+00, 1.02132116D+01, F2   
     6  1.05223792D+01, 1.08006981D+01, 1.10235999D+01, 1.11934942D+01, F2   
     7  1.12993789D+01, 1.13888680D+01, 1.14738434D+01, 1.15339831D+01, F2   
     8  1.15708589D+01, 1.15974578D+01, 1.16126808D+01, 1.16208134D+01, F2   
     9  1.16267950D+01,     57*0.0D+00/                                 F2   
      DATA TK_Na2/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751300020874, Na2  
     2  0.831799960396, 0.926599926485, 1.034200109262, 1.110800120079, Na2  
     3  1.189500033715, 1.342800094581, 1.498799966999, 1.653500004857, Na2  
     4  1.815300075228, 1.986799847427, 2.339800159472, 2.612499701529, Na2  
     5  2.755100112750, 2.890100342280, 3.024300148297, 3.176700218325, Na2  
     6  3.284799854331, 3.399399671363, 3.518900431050, 3.607299852394, Na2  
     7  3.650399631097, 3.693099664925, 3.731799595339, 3.769600456272, Na2  
     8  3.836800074841, 3.864499767293, 3.893100128107, 3.927500267789, Na2  
     9  3.958000003064, 3.984399651089, 3.991699814756, 4.000000000000, Na2  
     A      54*0.0D+00/                                                 Na2  
      DATA  K_Na2/                                                      071215
     1  3.08584049D-06, 1.50645488D-01, 3.07590810D-01, 7.79316068D-01, Na2  
     2  1.84783888D+00, 2.90268682D+00, 3.88509881D+00, 4.47217897D+00, Na2  
     3  4.99594479D+00, 5.83320358D+00, 6.49976186D+00, 7.02934230D+00, Na2  
     4  7.47900979D+00, 7.86287681D+00, 8.42346440D+00, 8.71024423D+00, Na2  
     5  8.82648330D+00, 8.91964976D+00, 8.99624339D+00, 9.05888778D+00, Na2  
     6  9.08303644D+00, 9.09088853D+00, 9.09025753D+00, 9.09436855D+00, Na2  
     7  9.10100516D+00, 9.11340584D+00, 9.13343301D+00, 9.16707771D+00, Na2  
     8  9.29085320D+00, 9.37817786D+00, 9.49630065D+00, 9.67693143D+00, Na2  
     9  9.86860527D+00, 1.00530629D+01, 1.01064006D+01, 1.01679307D+01, Na2  
     A      54*0.0D+00/                                                 Na2  
      DATA TK_Mg2/                                                      071215
     1  0.699999789529, 0.709100026195, 0.717899849912, 0.746999918371, Mg2  
     2  0.820500149046, 0.906799935515, 1.006200132329, 1.075500038798, Mg2  
     3  1.146100011260, 1.279400207770, 1.412199900998, 1.554300076848, Mg2  
     4  1.664599940384, 1.772800038819, 1.911299862961, 2.040599983139, Mg2  
     5  2.275200084523, 2.401699647923, 2.579199907475, 2.774500137332, Mg2  
     6  2.906099743477, 3.031399778343, 3.199699791178, 3.377600055221, Mg2  
     7  3.480899554872, 3.580399934818, 3.684499983839, 3.787499919113, Mg2  
     8  3.856199932386, 3.918800070327, 3.967800221540, 3.978899631344, Mg2  
     9  3.987199713952, 4.000000000000,     56*0.0D+00/                 Mg2  
      DATA  K_Mg2/                                                      071215
     1 -2.34240358D-06, 1.27841637D-01, 2.49273156D-01, 6.35468636D-01, Mg2  
     2  1.51529983D+00, 2.39673766D+00, 3.24521869D+00, 3.74955017D+00, Mg2  
     3  4.20145244D+00, 4.91311310D+00, 5.47310841D+00, 5.94273998D+00, Mg2  
     4  6.23280146D+00, 6.46507934D+00, 6.69841864D+00, 6.86116410D+00, Mg2  
     5  7.06049120D+00, 7.14384527D+00, 7.26567858D+00, 7.43410023D+00, Mg2  
     6  7.57528559D+00, 7.72977478D+00, 7.96148451D+00, 8.22444966D+00, Mg2  
     7  8.38093934D+00, 8.53375154D+00, 8.70103413D+00, 8.89053132D+00, Mg2  
     8  9.04715920D+00, 9.23096076D+00, 9.42071952D+00, 9.47129550D+00, Mg2  
     9  9.51116471D+00, 9.57643742D+00,     56*0.0D+00/                 Mg2  
      DATA TK_Al2/                                                      071215
     1  0.699999789529, 0.710500036963, 0.723399878256, 0.760000246348, Al2  
     2  0.797800128595, 0.852499964441, 0.974999896272, 1.112200083739, Al2  
     3  1.200099802438, 1.281700175351, 1.506200137952, 1.616999965087, Al2  
     4  1.736199871799, 1.847699962586, 1.956799974701, 2.056800354922, Al2  
     5  2.154699733379, 2.256000141870, 2.355299995100, 2.566899921447, Al2  
     6  2.790499989377, 2.958199999453, 3.121299903624, 3.242299826929, Al2  
     7  3.351100306086, 3.442999647882, 3.544000053030, 3.658699824766, Al2  
     8  3.749399973836, 3.828599892883, 3.883500202259, 3.936299869266, Al2  
     9  3.974699933131, 3.988999754363, 4.000000000000,     55*0.0D+00/ Al2  
      DATA  K_Al2/                                                      071215
     1 -3.85222404D-05, 1.60349381D-01, 3.52633176D-01, 8.71205675D-01, Al2  
     2  1.36746058D+00, 2.02152623D+00, 3.25266536D+00, 4.32911252D+00, Al2  
     3  4.89041252D+00, 5.34099399D+00, 6.32893150D+00, 6.73143872D+00, Al2  
     4  7.13709221D+00, 7.50466705D+00, 7.85450920D+00, 8.16136743D+00, Al2  
     5  8.44292571D+00, 8.70983174D+00, 8.94490776D+00, 9.35898624D+00, Al2  
     6  9.68602640D+00, 9.87426201D+00, 1.00215361D+01, 1.01103721D+01, Al2  
     7  1.01725509D+01, 1.02079314D+01, 1.02266060D+01, 1.02278717D+01, Al2  
     8  1.02244026D+01, 1.02298343D+01, 1.02476393D+01, 1.02865453D+01, Al2  
     9  1.03359208D+01, 1.03601752D+01, 1.03813262D+01,     55*0.0D+00/ Al2  
      DATA TK_Si2/                                                      071215
     1  0.699999789529, 0.710500036963, 0.723599883045, 0.760400235797, Si2  
     2  0.798600146653, 0.854400008025, 0.979099979238, 1.106700064771, Si2  
     3  1.251000040406, 1.382300042974, 1.496499918757, 1.612800070107, Si2  
     4  1.764700111953, 1.927699882684, 2.043100040369, 2.182500127314, Si2  
     5  2.283899823038, 2.402699672320, 2.514700345992, 2.623399952111, Si2  
     6  3.033299823225, 3.324899804189, 3.430200279055, 3.543900050638, Si2  
     7  3.655699754765, 3.761700265279, 3.862399719310, 3.907299783588, Si2  
     8  3.950599841481, 3.980999574756, 3.990599790206, 4.000000000000, Si2  
     9      58*0.0D+00/                                                 Si2  
      DATA  K_Si2/                                                      071215
     1 -9.61348901D-05, 1.48107497D-01, 3.28542120D-01, 8.10381578D-01, Si2  
     2  1.27390954D+00, 1.89032941D+00, 3.04718130D+00, 3.98129512D+00, Si2  
     3  4.81687620D+00, 5.44122799D+00, 5.92600903D+00, 6.39581328D+00, Si2  
     4  6.99552151D+00, 7.62484629D+00, 8.05662520D+00, 8.55136669D+00, Si2  
     5  8.88412552D+00, 9.23731505D+00, 9.53067650D+00, 9.77870923D+00, Si2  
     6  1.04430148D+01, 1.07434632D+01, 1.08364474D+01, 1.09339017D+01, Si2  
     7  1.10269020D+01, 1.11077396D+01, 1.11705030D+01, 1.11937197D+01, Si2  
     8  1.12156685D+01, 1.12331486D+01, 1.12394596D+01, 1.12461484D+01, Si2  
     9      58*0.0D+00/                                                 Si2  
      DATA TK_P2/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750700005324, P2   
     2  0.830099921425, 0.924399979531, 1.031300186745, 1.105400034803, P2   
     3  1.180799852570, 1.326200021441, 1.480100012057, 1.635399976685, P2   
     4  1.760300231031, 1.891399893728, 2.016900392026, 2.154999740474, P2   
     5  2.275100091618, 2.389400328612, 2.638500304795, 2.871499924798, P2   
     6  3.148299730314, 3.332799984613, 3.499699983072, 3.587100091791, P2   
     7  3.686399851508, 3.772800263319, 3.863999755868, 3.905899752148, P2   
     8  3.943399679360, 3.978099688827, 4.000000000000,     59*0.0D+00/ P2   
      DATA  K_P2/                                                       071215
     1  1.65804377D-05, 1.74226100D-01, 3.52118941D-01, 8.90708136D-01, P2   
     2  2.10902600D+00, 3.31995690D+00, 4.44254272D+00, 5.09450887D+00, P2   
     3  5.67091734D+00, 6.58427247D+00, 7.33680494D+00, 7.93853761D+00, P2   
     4  8.34029018D+00, 8.70458442D+00, 9.01220925D+00, 9.31622069D+00, P2   
     5  9.55796966D+00, 9.77139256D+00, 1.01803767D+01, 1.04891168D+01, P2   
     6  1.07723756D+01, 1.09234560D+01, 1.10500738D+01, 1.11222118D+01, P2   
     7  1.12154237D+01, 1.13036908D+01, 1.13905653D+01, 1.14227530D+01, P2   
     8  1.14454969D+01, 1.14611225D+01, 1.14685255D+01,     59*0.0D+00/ P2   
      DATA TK_S2/                                                       071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749599978008, S2   
     2  0.827099989317, 0.919200067714, 1.024000073308, 1.097199981936, S2   
     3  1.172000021566, 1.315600072581, 1.463299945039, 1.612400080109, S2   
     4  1.736099874204, 1.856300051610, 2.142100174142, 2.235699685554, S2   
     5  2.333900019615, 2.429900281457, 2.526599966623, 2.658599826589, S2   
     6  2.797200141103, 2.967400211064, 3.130100122779, 3.327799867067, S2   
     7  3.520600413121, 3.611299688534, 3.702499864629, 3.787099909682, S2   
     8  3.861099689607, 3.944299699681, 3.978299674457, 4.000000000000, S2   
     9      58*0.0D+00/                                                 S2   
      DATA  K_S2/                                                       071215
     1 -1.48565768D-04, 1.66402786D-01, 3.31449899D-01, 8.42774088D-01, S2   
     2  1.99718409D+00, 3.15078509D+00, 4.22959856D+00, 4.86318896D+00, S2   
     3  5.42679837D+00, 6.31837619D+00, 7.03737181D+00, 7.61786257D+00, S2   
     4  8.01934235D+00, 8.35821849D+00, 9.03492390D+00, 9.23557592D+00, S2   
     5  9.44204631D+00, 9.64105592D+00, 9.83740633D+00, 1.00935395D+01, S2   
     6  1.03401089D+01, 1.06040443D+01, 1.08156189D+01, 1.10265274D+01, S2   
     7  1.11992219D+01, 1.12753246D+01, 1.13507817D+01, 1.14195955D+01, S2   
     8  1.14771738D+01, 1.15355861D+01, 1.15566861D+01, 1.15692431D+01, S2   
     9      58*0.0D+00/                                                 S2   
      DATA TK_Cl2/                                                      071215
     1  0.699999789529, 0.710000049601, 0.721499832753, 0.755400127132, Cl2  
     2  0.791199979616, 0.842100097703, 0.946599877282, 1.063900042683, Cl2  
     3  1.141199903243, 1.218099840754, 1.363300052765, 1.520599982388, Cl2  
     4  1.681699877374, 1.816300100178, 1.963499964328, 2.101500283730, Cl2  
     5  2.238099743758, 2.414099949846, 2.619399860365, 2.740599787158, Cl2  
     6  2.864099745200, 3.052800253279, 3.266400374108, 3.398399744378, Cl2  
     7  3.460800052719, 3.520900391313, 3.608299780024, 3.689399642564, Cl2  
     8  3.761000248356, 3.876200033218, 3.950799845848, 3.981399583736, Cl2  
     9  4.000000000000,     57*0.0D+00/                                 Cl2  
      DATA  K_Cl2/                                                      071215
     1  1.08748534D-04, 1.83245444D-01, 3.89108213D-01, 9.67810378D-01, Cl2  
     2  1.53596123D+00, 2.27446231D+00, 3.57002267D+00, 4.73681745D+00, Cl2  
     3  5.37376746D+00, 5.92310532D+00, 6.77758224D+00, 7.50080454D+00, Cl2  
     4  8.08842222D+00, 8.49560248D+00, 8.87873133D+00, 9.19457546D+00, Cl2  
     5  9.47514980D+00, 9.79689196D+00, 1.01262208D+01, 1.03035384D+01, Cl2  
     6  1.04733920D+01, 1.07102205D+01, 1.09397924D+01, 1.10587420D+01, Cl2  
     7  1.11064395D+01, 1.11431826D+01, 1.11681762D+01, 1.11457596D+01, Cl2  
     8  1.10868315D+01, 1.09479488D+01, 1.08564639D+01, 1.08224352D+01, Cl2  
     9  1.08031304D+01,     57*0.0D+00/                                 Cl2  
      DATA TK_K2/                                                       071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749699980301, K2   
     2  0.827499979636, 0.919100065475, 1.023600063430, 1.096899989609, K2   
     3  1.172600007642, 1.317500111458, 1.465499890919, 1.614300032600, K2   
     4  1.774100071459, 1.956399964607, 2.132200141577, 2.318499818330, K2   
     5  2.480699544075, 2.642300173699, 2.810300393702, 2.971300175990, K2   
     6  3.190199568536, 3.284999858928, 3.389800352641, 3.487699705678, K2   
     7  3.575899827496, 3.649299671214, 3.688599698283, 3.723200029351, K2   
     8  3.792100022081, 3.826599844531, 3.862099712455, 3.897499803874, K2   
     9  3.937899753948, 3.980599565776, 4.000000000000,     55*0.0D+00/ K2   
      DATA  K_K2/                                                       071215
     1  1.05457459D-05, 1.47169990D-01, 2.94554528D-01, 7.46602450D-01, K2   
     2  1.77239701D+00, 2.79093787D+00, 3.75019713D+00, 4.31827499D+00, K2   
     3  4.83034948D+00, 5.64075869D+00, 6.29158024D+00, 6.81062039D+00, K2   
     4  7.25073933D+00, 7.63759038D+00, 7.92274957D+00, 8.15656032D+00, K2   
     5  8.31863560D+00, 8.45171811D+00, 8.56609320D+00, 8.65433257D+00, K2   
     6  8.73400185D+00, 8.75001425D+00, 8.75753724D+00, 8.76594101D+00, K2   
     7  8.78860119D+00, 8.83201667D+00, 8.87180534D+00, 8.92193078D+00, K2   
     8  9.08573509D+00, 9.21043510D+00, 9.37257689D+00, 9.56656560D+00, K2   
     9  9.81950525D+00, 1.01098888D+01, 1.02456293D+01,     55*0.0D+00/ K2   
      DATA TK_Cu2/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719799801886, 0.751700031241, Cu2  
     2  0.832599978735, 0.928599878261, 1.037200029106, 1.113500049995, Cu2  
     3  1.191400009990, 1.341900116962, 1.499799987974, 1.659399847868, Cu2  
     4  1.767500036176, 1.875600022561, 2.075099811089, 2.333300005392, Cu2  
     5  2.567499878743, 2.791000000700, 3.021100374556, 3.252200050192, Cu2  
     6  3.335600054457, 3.416600003460, 3.606799888579, 3.692299648371, Cu2  
     7  3.778899821891, 3.874900003810, 3.951299856766, 3.981499585982, Cu2  
     8  3.990799794669, 4.000000000000,     60*0.0D+00/                 Cu2  
      DATA  K_Cu2/                                                      071215
     1  3.69935037D-05, 1.63384074D-01, 3.33391521D-01, 8.42256178D-01, Cu2  
     2  1.99162136D+00, 3.13131852D+00, 4.18471563D+00, 4.80417757D+00, Cu2  
     3  5.35278020D+00, 6.22213146D+00, 6.93263524D+00, 7.50411234D+00, Cu2  
     4  7.83047981D+00, 8.11891154D+00, 8.56972242D+00, 9.02063379D+00, Cu2  
     5  9.32465563D+00, 9.54640258D+00, 9.72570085D+00, 9.86925588D+00, Cu2  
     6  9.91413350D+00, 9.95644671D+00, 1.00760296D+01, 1.01540773D+01, Cu2  
     7  1.02538182D+01, 1.03889803D+01, 1.05216367D+01, 1.05842221D+01, Cu2  
     8  1.06051076D+01, 1.06265955D+01,     60*0.0D+00/                 Cu2  
      DATA TK_As2/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751400023466, As2  
     2  0.831799960396, 0.927299909607, 1.035400077199, 1.110200135654, As2  
     3  1.186199965004, 1.333000003705, 1.490299788711, 1.649800093269, As2  
     4  1.773800063927, 1.905499943973, 2.022600273039, 2.133800178352, As2  
     5  2.441999628947, 2.709800051095, 2.937399772879, 3.195799699778, As2  
     6  3.344500265430, 3.504900104309, 3.586400075391, 3.668100020445, As2  
     7  3.848300331262, 3.929600314445, 4.000000000000,     63*0.0D+00/ As2  
      DATA  K_As2/                                                      071215
     1  8.97922537D-05, 1.78327397D-01, 3.62003324D-01, 9.13707074D-01, As2  
     2  2.15951186D+00, 3.39491541D+00, 4.53554750D+00, 5.19551256D+00, As2  
     3  5.77737924D+00, 6.69927279D+00, 7.46538033D+00, 8.07797097D+00, As2  
     4  8.47291177D+00, 8.83567350D+00, 9.12061040D+00, 9.36444879D+00, As2  
     5  9.92172218D+00, 1.02790717D+01, 1.05105369D+01, 1.07179626D+01, As2  
     6  1.08201770D+01, 1.09299127D+01, 1.09922085D+01, 1.10616207D+01, As2  
     7  1.12249997D+01, 1.12907668D+01, 1.13418698D+01,     63*0.0D+00/ As2  
      DATA TK_Se2/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719899799358, 0.751800033832, Se2  
     2  0.832799983320, 0.928799873439, 1.037300026434, 1.110400130462, Se2  
     3  1.184099921280, 1.328299973624, 1.491999824369, 1.657099909067, Se2  
     4  1.854900019495, 2.027899881756, 2.149999622226, 2.280799760206, Se2  
     5  2.425800203378, 2.523100226533, 2.617999828137, 2.815000058071, Se2  
     6  2.918800046326, 3.032099794879, 3.144600006407, 3.260200225235, Se2  
     7  3.411199861369, 3.605799960950, 3.700899828200, 3.788799949762, Se2  
     8  3.908899819520, 3.961800087075, 4.000000000000,     59*0.0D+00/ Se2  
      DATA  K_Se2/                                                      071215
     1 -5.91028265D-06, 1.71610860D-01, 3.51938939D-01, 8.86288630D-01, Se2  
     2  2.09385124D+00, 3.28834872D+00, 4.38984675D+00, 5.01190859D+00, Se2  
     3  5.55874757D+00, 6.44153680D+00, 7.21857459D+00, 7.83335276D+00, Se2  
     4  8.42020873D+00, 8.84234762D+00, 9.09994926D+00, 9.33960264D+00, Se2  
     5  9.55910554D+00, 9.67969319D+00, 9.77927405D+00, 9.94741317D+00, Se2  
     6  1.00274749D+01, 1.01174082D+01, 1.02132392D+01, 1.03184003D+01, Se2  
     7  1.04615722D+01, 1.06490343D+01, 1.07399081D+01, 1.08206664D+01, Se2  
     8  1.09189710D+01, 1.09562303D+01, 1.09813106D+01,     59*0.0D+00/ Se2  
      DATA TK_Sb2/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751400023466, Sb2  
     2  0.831799960396, 0.926999916840, 1.034900090559, 1.109300124708, Sb2  
     3  1.184799935855, 1.330699950967, 1.487699835769, 1.649900095627, Sb2  
     4  1.781900172679, 1.912099882520, 2.048400161695, 2.184100012585, Sb2  
     5  2.560300391185, 2.771200384915, 3.186099846684, 3.340100166462, Sb2  
     6  3.525700042375, 3.597600327668, 3.668400026569, 3.844800252524, Sb2  
     7  3.907399785834, 3.962700107245, 3.985599678030, 4.000000000000, Sb2  
     8      62*0.0D+00/                                                 Sb2  
      DATA  K_Sb2/                                                      071215
     1  1.18622765D-04, 1.78621473D-01, 3.62570228D-01, 9.15090091D-01, Sb2  
     2  2.16271586D+00, 3.39636239D+00, 4.53735397D+00, 5.19554130D+00, Sb2  
     3  5.77555428D+00, 6.69585373D+00, 7.46432327D+00, 8.08896347D+00, Sb2  
     4  8.50707612D+00, 8.86025112D+00, 9.17865912D+00, 9.45026888D+00, Sb2  
     5  1.00085455D+01, 1.02287339D+01, 1.05499109D+01, 1.06496730D+01, Sb2  
     6  1.07868506D+01, 1.08541925D+01, 1.09304236D+01, 1.11493667D+01, Sb2  
     7  1.12289564D+01, 1.12985953D+01, 1.13277676D+01, 1.13464377D+01, Sb2  
     8      62*0.0D+00/                                                 Sb2  
      DATA TK_Te2/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720099799225, 0.752400049382, Te2  
     2  0.834500022290, 0.931599875676, 1.041499987577, 1.119099904635, Te2  
     3  1.198199844188, 1.272400039563, 1.351399946593, 1.512800154259, Te2  
     4  1.677599892822, 1.807799994752, 1.933099891469, 2.253800090648, Te2  
     5  2.527599892363, 2.639200321380, 2.751300019190, 2.894100044915, Te2  
     6  3.036899908264, 3.152999671786, 3.260600234840, 3.353600117325, Te2  
     7  3.446299716576, 3.531699768472, 3.609099722127, 3.777699908730, Te2  
     8  3.871499926895, 3.953299900437, 3.982199601697, 4.000000000000, Te2  
     9      58*0.0D+00/                                                 Te2  
      DATA  K_Te2/                                                      071215
     1  2.85526903D-06, 1.72453752D-01, 3.53486672D-01, 8.91179230D-01, Te2  
     2  2.10578489D+00, 3.30228077D+00, 4.40484833D+00, 5.05409629D+00, Te2  
     3  5.62625909D+00, 6.09408815D+00, 6.53084940D+00, 7.26810991D+00, Te2  
     4  7.86243309D+00, 8.25118325D+00, 8.57358419D+00, 9.21450539D+00, Te2  
     5  9.59540149D+00, 9.71766745D+00, 9.82444318D+00, 9.93789485D+00, Te2  
     6  1.00268099D+01, 1.00852196D+01, 1.01351395D+01, 1.01807388D+01, Te2  
     7  1.02326310D+01, 1.02881057D+01, 1.03447656D+01, 1.04802928D+01, Te2  
     8  1.05548727D+01, 1.06179156D+01, 1.06407363D+01, 1.06553188D+01, Te2  
     9      58*0.0D+00/                                                 Te2  
      DATA TK_I2/                                                       071215
     1  0.699999789529, 0.710200044546, 0.722199849517, 0.734699905418, I2   
     2  0.757100171190, 0.794100045076, 0.846599989127, 0.954299906314, I2   
     3  1.075000051982, 1.158499879885, 1.242699824361, 1.320900142120, I2   
     4  1.403100018843, 1.564000032184, 1.730800001654, 1.848299948110, I2   
     5  1.956899977225, 2.247399946967, 2.493999862128, 2.748099946752, I2   
     6  2.928300255380, 3.126700038229, 3.201499831824, 3.275700052112, I2   
     7  3.346900319413, 3.407999789345, 3.501400022712, 3.594100250330, I2   
     8  3.714300128140, 3.773700198191, 3.837000079197, 3.941499636461, I2   
     9  3.977199753496, 3.989599767834, 4.000000000000,     55*0.0D+00/ I2   
      DATA  K_I2/                                                       071215
     1 -4.94887188D-05, 1.87099446D-01, 4.02141457D-01, 6.20473107D-01, I2   
     2  9.97636982D-01, 1.58321671D+00, 2.34067464D+00, 3.66230217D+00, I2   
     3  4.84349763D+00, 5.51477101D+00, 6.09500978D+00, 6.56193158D+00, I2   
     4  6.99080489D+00, 7.68608028D+00, 8.25779836D+00, 8.59328474D+00, I2   
     5  8.86323829D+00, 9.42895748D+00, 9.76922071D+00, 1.00284275D+01, I2   
     6  1.01743073D+01, 1.03094317D+01, 1.03544028D+01, 1.03948025D+01, I2   
     7  1.04268501D+01, 1.04451746D+01, 1.04465482D+01, 1.04078240D+01, I2   
     8  1.03129486D+01, 1.02609402D+01, 1.02099875D+01, 1.01462972D+01, I2   
     9  1.01323165D+01, 1.01286218D+01, 1.01260530D+01,     55*0.0D+00/ I2   
      DATA TK_Cs2/                                                      071215
     1  0.699999789529, 0.709100026195, 0.717899849912, 0.746999918371, Cs2  
     2  0.820300153887, 0.907499919379, 1.007400157941, 1.139099898833, Cs2  
     3  1.262900148670, 1.400200089630, 1.543200042110, 1.664199930951, Cs2  
     4  1.794400063133, 1.923399990452, 2.057100361831, 2.290999971580, Cs2  
     5  2.533799806562, 2.662999927199, 2.798200163749, 2.912299887410, Cs2  
     6  3.022900247285, 3.091200168079, 3.159899828936, 3.240299781411, Cs2  
     7  3.313800138888, 3.432200137499, 3.486999690154, 3.540299964532, Cs2  
     8  3.586700082419, 3.629900119814, 3.719600242463, 3.768700434513, Cs2  
     9  3.817299879894, 3.860399673612, 3.925700227798, 3.970300249288, Cs2  
     A  3.987899729667, 4.000000000000,     52*0.0D+00/                 Cs2  
      DATA  K_Cs2/                                                      071215
     1 -7.40185442D-07, 1.39312680D-01, 2.71615033D-01, 6.92217013D-01, Cs2  
     2  1.64693141D+00, 2.61337017D+00, 3.53471188D+00, 4.50908879D+00, Cs2  
     3  5.23066444D+00, 5.86097166D+00, 6.37039589D+00, 6.71077721D+00, Cs2  
     4  7.00535186D+00, 7.24046065D+00, 7.43863607D+00, 7.70589003D+00, Cs2  
     5  7.91168078D+00, 8.00035467D+00, 8.07976427D+00, 8.13431040D+00, Cs2  
     6  8.17132273D+00, 8.18377195D+00, 8.18764122D+00, 8.18329038D+00, Cs2  
     7  8.17591917D+00, 8.17548283D+00, 8.18849902D+00, 8.21523760D+00, Cs2  
     8  8.25484412D+00, 8.31132243D+00, 8.52529747D+00, 8.72239054D+00, Cs2  
     9  8.98288886D+00, 9.26237401D+00, 9.74200640D+00, 1.00852543D+01, Cs2  
     A  1.02205789D+01, 1.03130464D+01,     52*0.0D+00/                 Cs2  
      DATA TK_H2p/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718599832218, 0.748899961951, H2p  
     2  0.824800044980, 0.920900063923, 1.028400181967, 1.103399988697, H2p  
     3  1.174299968191, 1.312400007103, 1.372699939755, 1.430199820105, H2p  
     4  1.560299945056, 1.645699996568, 1.737699835728, 1.835300037750, H2p  
     5  1.933499900040, 2.066599936665, 2.200399800452, 2.377200059758, H2p  
     6  2.530999738461, 2.651699657491, 2.762100284118, 3.070699712047, H2p  
     7  3.207799973006, 3.353900094674, 3.498299950398, 3.626700048637, H2p  
     8  3.764500332973, 3.873399969876, 3.951199854582, 4.000000000000, H2p  
     9      58*0.0D+00/                                                 H2p  
      DATA  K_H2p/                                                      071215
     1  4.26172519D-05, 1.35984554D-01, 2.69575890D-01, 6.88783740D-01, H2p  
     2  1.64068741D+00, 2.67305343D+00, 3.64199031D+00, 4.22214394D+00, H2p  
     3  4.70764688D+00, 5.49379219D+00, 5.77576646D+00, 6.01386970D+00, H2p  
     4  6.46405571D+00, 6.71293014D+00, 6.95613543D+00, 7.19617646D+00, H2p  
     5  7.42468459D+00, 7.71676181D+00, 7.99042578D+00, 8.32477889D+00, H2p  
     6  8.59656888D+00, 8.80072686D+00, 8.98118819D+00, 9.44150918D+00, H2p  
     7  9.61526707D+00, 9.77463441D+00, 9.90478672D+00, 9.99779275D+00, H2p  
     8  1.00770179D+01, 1.01303524D+01, 1.01674419D+01, 1.01915386D+01, H2p  
     9      58*0.0D+00/                                                 H2p  
      DATA TK_He2p/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719899799358, 0.752000039016, He2p 
     2  0.832299971858, 0.935899959449, 1.045400074109, 1.229099996184, He2p 
     3  1.308299994172, 1.417100013654, 1.519499979473, 1.630000111164, He2p 
     4  1.737599838133, 1.850499918562, 2.074299791034, 2.315500033725, He2p 
     5  2.545800087417, 2.722400107818, 2.966700194970, 3.379200088058, He2p 
     6  3.518800428715, 3.614799773450, 3.683100081345, 3.749999987182, He2p 
     7  3.869599883822, 3.949799823864, 3.980799570266, 4.000000000000, He2p 
     8      62*0.0D+00/                                                 He2p 
      DATA  K_He2p/                                                     071215
     1 -6.04135557D-05, 1.46710345D-01, 3.01212117D-01, 7.63650962D-01, He2p 
     2  1.80524721D+00, 2.94161815D+00, 3.93557007D+00, 5.22303038D+00, He2p 
     3  5.65713475D+00, 6.16029933D+00, 6.55765132D+00, 6.92604385D+00, He2p 
     4  7.24004898D+00, 7.53347672D+00, 8.03727136D+00, 8.50494313D+00, He2p 
     5  8.90789928D+00, 9.19508278D+00, 9.55068330D+00, 9.99295088D+00, He2p 
     6  1.00913859D+01, 1.01406838D+01, 1.01659406D+01, 1.01838895D+01, He2p 
     7  1.02078671D+01, 1.02264345D+01, 1.02355408D+01, 1.02418841D+01, He2p 
     8      62*0.0D+00/                                                 He2p 
      DATA TK_C2p/                                                      071215
     1  0.699999789529, 0.710500036963, 0.723199873466, 0.759500233390, C2p  
     2  0.798400142139, 0.852899973616, 0.964899927878, 1.089100143962, C2p  
     3  1.181699871309, 1.277700166920, 1.360500124600, 1.450399817574, C2p  
     4  1.540299982622, 1.629400097898, 1.818500155067, 2.029299778398, C2p  
     5  2.247999959707, 2.385400242639, 2.524000159699, 2.659299843744, C2p  
     6  2.846300304210, 2.978999614167, 3.131900164029, 3.360699652657, C2p  
     7  3.618699868071, 3.816899908772, 3.951099852399, 4.000000000000, C2p  
     8      62*0.0D+00/                                                 C2p  
      DATA  K_C2p/                                                      071215
     1  7.08488687D-05, 1.60565962D-01, 3.50271784D-01, 8.67455891D-01, C2p  
     2  1.38328658D+00, 2.04613730D+00, 3.22340033D+00, 4.29868548D+00, C2p  
     3  4.97910435D+00, 5.59918999D+00, 6.07791109D+00, 6.55012819D+00, C2p  
     4  6.98129900D+00, 7.37379071D+00, 8.10578585D+00, 8.78036698D+00, C2p  
     5  9.35534356D+00, 9.66813084D+00, 9.95489711D+00, 1.02106538D+01, C2p  
     6  1.05247706D+01, 1.07185591D+01, 1.09124091D+01, 1.11513876D+01, C2p  
     7  1.13706302D+01, 1.15236471D+01, 1.16251055D+01, 1.16621581D+01, C2p  
     8      62*0.0D+00/                                                 C2p  
      DATA TK_N2p/                                                      071215
     1  0.699999789529, 0.710300042018, 0.722399854307, 0.757500181557, N2p  
     2  0.794800060877, 0.847699962586, 0.951899843518, 1.071200152175, N2p  
     3  1.166899997780, 1.260500205046, 1.362100083551, 1.457799978405, N2p  
     4  1.542500027751, 1.630700093731, 1.784100118358, 1.984299905771, N2p  
     5  2.139700313959, 2.373599973850, 2.548600150205, 2.665499982449, N2p  
     6  2.777599904754, 2.959100020188, 3.110799674418, 3.293100040578, N2p  
     7  3.435199925166, 3.544900074557, 3.621099924078, 3.689699621670, N2p  
     8  3.840000144540, 3.933300085488, 3.974099976243, 4.000000000000, N2p  
     9      58*0.0D+00/                                                 N2p  
      DATA  K_N2p/                                                      071215
     1 -4.23698053D-04, 1.74668097D-01, 3.75492586D-01, 9.29889210D-01, N2p  
     2  1.47588790D+00, 2.18114319D+00, 3.36680372D+00, 4.45878979D+00, N2p  
     3  5.17896963D+00, 5.78305451D+00, 6.35511268D+00, 6.83447519D+00, N2p  
     4  7.22149451D+00, 7.59429876D+00, 8.18220118D+00, 8.85074127D+00, N2p  
     5  9.30061572D+00, 9.88192858D+00, 1.02585902D+01, 1.04893527D+01, N2p  
     6  1.06973706D+01, 1.10065792D+01, 1.12355279D+01, 1.14706672D+01, N2p  
     7  1.16193497D+01, 1.17086784D+01, 1.17564512D+01, 1.17903634D+01, N2p  
     8  1.18461839D+01, 1.18788227D+01, 1.18943129D+01, 1.19045863D+01, N2p  
     9      58*0.0D+00/                                                 N2p  
      DATA TK_O2p/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750299994957, O2p  
     2  0.829099940914, 0.922500025344, 1.028700189375, 1.100799928760, O2p  
     3  1.173799979794, 1.314900058258, 1.469799785138, 1.627300051469, O2p  
     4  1.743699856922, 1.859500125016, 2.118699854832, 2.368599854023, O2p  
     5  2.488099733725, 2.606199921403, 2.772100317392, 2.984699646841, O2p  
     6  3.145299954174, 3.334500027018, 3.488099714549, 3.639100309637, O2p  
     7  3.753200063022, 3.864499767293, 3.946399747096, 3.979399595417, O2p  
     8  4.000000000000,     61*0.0D+00/                                 O2p  
      DATA  K_O2p/                                                      071215
     1 -1.90292313D-04, 1.74121790D-01, 3.52261075D-01, 8.93876115D-01, O2p  
     2  2.11913345D+00, 3.33691024D+00, 4.47126053D+00, 5.11808214D+00, O2p  
     3  5.68923735D+00, 6.60077672D+00, 7.38035727D+00, 8.00566514D+00, O2p  
     4  8.39038992D+00, 8.72472010D+00, 9.34680508D+00, 9.83705092D+00, O2p  
     5  1.00487760D+01, 1.02486169D+01, 1.05156663D+01, 1.08287065D+01, O2p  
     6  1.10359223D+01, 1.12440995D+01, 1.13864085D+01, 1.15068833D+01, O2p  
     7  1.15882324D+01, 1.16604532D+01, 1.17053310D+01, 1.17196769D+01, O2p  
     8  1.17271288D+01,     61*0.0D+00/                                 O2p  
      DATA TK_Ne2p/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751400023466, Ne2p 
     2  0.831899962688, 0.927499904784, 1.035700069184, 1.110900117484, Ne2p 
     3  1.187399989990, 1.335300056443, 1.493599857929, 1.653000018161, Ne2p 
     4  1.776300126695, 1.909099853737, 2.025700044175, 2.137000251902, Ne2p 
     5  2.390100334157, 2.648699710015, 3.116499796349, 3.476899757136, Ne2p 
     6  3.818799771602, 3.925400221133, 4.000000000000,     67*0.0D+00/ Ne2p 
      DATA  K_Ne2p/                                                     071215
     1  5.39797620D-05, 1.55972170D-01, 3.16689858D-01, 7.99712596D-01, Ne2p 
     2  1.89347695D+00, 2.98090256D+00, 3.98867203D+00, 4.57659491D+00, Ne2p 
     3  5.09724215D+00, 5.92730831D+00, 6.62266103D+00, 7.18162450D+00, Ne2p 
     4  7.54485609D+00, 7.88706977D+00, 8.15577200D+00, 8.39027724D+00, Ne2p 
     5  8.85624733D+00, 9.24789017D+00, 9.78973480D+00, 1.00929190D+01, Ne2p 
     6  1.03223169D+01, 1.03867651D+01, 1.04306183D+01,     67*0.0D+00/ Ne2p 
      DATA TK_P2p/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720299804015, 0.752900062340, P2p  
     2  0.835400042922, 0.932699897106, 1.049500165080, 1.172200016925, P2p  
     3  1.298000125976, 1.434499917722, 1.621599925448, 1.736799857371, P2p  
     4  1.852499964441, 1.960000055455, 2.074699801061, 2.303500277227, P2p  
     5  2.511500279126, 2.649299666545, 2.793700061843, 3.022900247285, P2p  
     6  3.228999604269, 3.383400190574, 3.522800253191, 3.611499693386, P2p  
     7  3.703299882843, 3.787599921470, 3.882900188017, 3.951099852399, P2p  
     8  3.981499585982, 4.000000000000,     60*0.0D+00/                 P2p  
      DATA  K_P2p/                                                      071215
     1  5.26508938D-06, 1.50329210D-01, 3.11181173D-01, 7.84337843D-01, P2p  
     2  1.84968748D+00, 2.89905547D+00, 3.92261279D+00, 4.78332210D+00, P2p  
     3  5.49326611D+00, 6.11662981D+00, 6.79803609D+00, 7.15405836D+00, P2p  
     4  7.48572699D+00, 7.78094449D+00, 8.08532603D+00, 8.65447382D+00, P2p  
     5  9.11411407D+00, 9.37853151D+00, 9.61723828D+00, 9.91863123D+00, P2p  
     6  1.01213986D+01, 1.02437077D+01, 1.03456464D+01, 1.04134379D+01, P2p  
     7  1.04897174D+01, 1.05643540D+01, 1.06472708D+01, 1.07000006D+01, P2p  
     8  1.07209785D+01, 1.07330308D+01,     60*0.0D+00/                 P2p  
      DATA TK_S2p/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, S2p  
     2  0.830499930595, 0.925199960242, 1.032500154683, 1.107000071687, S2p  
     3  1.182799894212, 1.328799962240, 1.483199940150, 1.639099884543, S2p  
     4  1.762900160667, 1.895899992156, 2.016000370843, 2.142900118252, S2p  
     5  2.246999938473, 2.350600344533, 2.493399850171, 2.649599644810, S2p  
     6  2.892100193598, 3.187999708892, 3.382200160186, 3.565000045889, S2p  
     7  3.657399794432, 3.754200086721, 3.906499765622, 3.960100048977, S2p  
     8  4.000000000000,     61*0.0D+00/                                 S2p  
      DATA  K_S2p/                                                      071215
     1  1.09973595D-04, 1.69504100D-01, 3.44217058D-01, 8.69520389D-01, S2p  
     2  2.05702815D+00, 3.23938855D+00, 4.33502229D+00, 4.97249095D+00, S2p  
     3  5.53611384D+00, 6.42874084D+00, 7.16472270D+00, 7.75490600D+00, S2p  
     4  8.14527180D+00, 8.50854410D+00, 8.79940065D+00, 9.07814379D+00, S2p  
     5  9.28933708D+00, 9.48557859D+00, 9.73322462D+00, 9.97230786D+00, S2p  
     6  1.02789878D+01, 1.05664664D+01, 1.07200390D+01, 1.08543512D+01, S2p  
     7  1.09245813D+01, 1.10031949D+01, 1.11368301D+01, 1.11847799D+01, S2p  
     8  1.12201504D+01,     61*0.0D+00/                                 S2p  
      DATA TK_H2m/                                                      071215
     1  0.699999789529, 0.710400039490, 0.723199873466, 0.759300228206, H2m  
     2  0.852499964441, 0.963399966932, 1.087600107070, 1.174299968191, H2m  
     3  1.262600155717, 1.342500102041, 1.427899871165, 1.603299973419, H2m  
     4  1.789099994902, 1.914999953421, 2.067199891942, 2.165299963692, H2m  
     5  2.416200000772, 2.536999884391, 2.663799944879, 2.901499640034, H2m  
     6  3.109999657305, 3.502800055351, 3.608999729364, 3.718600220893, H2m  
     7  3.883200195138, 3.951899869867, 4.000000000000,     63*0.0D+00/ H2m  
      DATA  K_H2m/                                                      071215
     1  2.23458582D-06, 1.40601361D-01, 3.09899384D-01, 7.66329669D-01, H2m  
     2  1.81524337D+00, 2.86066423D+00, 3.82791098D+00, 4.40364944D+00, H2m  
     3  4.92328636D+00, 5.34568268D+00, 5.75552505D+00, 6.48968233D+00, H2m  
     4  7.14402720D+00, 7.52957491D+00, 7.94006625D+00, 8.17331366D+00, H2m  
     5  8.68013006D+00, 8.89601503D+00, 9.11243669D+00, 9.49945740D+00, H2m  
     6  9.82116028D+00, 1.03320806D+01, 1.04336571D+01, 1.05183344D+01, H2m  
     7  1.06156331D+01, 1.06511555D+01, 1.06758249D+01,     63*0.0D+00/ H2m  
      DATA TK_C2m/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720099799225, 0.752400049382, C2m  
     2  0.834200015413, 0.934299928278, 1.044000043046, 1.123899975046, C2m  
     3  1.202999864450, 1.345000039870, 1.520399977068, 1.706099957842, C2m  
     4  1.823500094671, 1.940100036964, 2.071699725854, 2.195299695693, C2m  
     5  2.433400046066, 2.652299672195, 2.925300190422, 3.061200350682, C2m  
     6  3.207299961801, 3.360999660611, 3.452299849214, 3.540999981275, C2m  
     7  3.617099829252, 3.688099733107, 3.832599983360, 3.932300157562, C2m  
     8  3.973600012170, 4.000000000000,     60*0.0D+00/                 C2m  
      DATA  K_C2m/                                                      071215
     1  1.55017511D-04, 1.64328801D-01, 3.36723841D-01, 8.49203624D-01, C2m  
     2  2.00775555D+00, 3.20062432D+00, 4.28764507D+00, 4.96529083D+00, C2m  
     3  5.55900692D+00, 6.46793657D+00, 7.37116651D+00, 8.12692692D+00, C2m  
     4  8.52358571D+00, 8.86946927D+00, 9.21447592D+00, 9.50449632D+00, C2m  
     5  9.99723596D+00, 1.03981717D+01, 1.08396106D+01, 1.10298930D+01, C2m  
     6  1.12101008D+01, 1.13754786D+01, 1.14663525D+01, 1.15536397D+01, C2m  
     7  1.16295929D+01, 1.17007300D+01, 1.18329791D+01, 1.19018333D+01, C2m  
     8  1.19243434D+01, 1.19372692D+01,     60*0.0D+00/                 C2m  
      DATA TK_LiH/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.751000013099, LiH  
     2  0.830799937472, 0.924099986765, 1.031400184073, 1.102399965644, LiH  
     3  1.174999951946, 1.320200158058, 1.480400005098, 1.641699902225, LiH  
     4  1.756000138363, 1.874599997822, 1.995799904553, 2.120599893858, LiH  
     5  2.356699891013, 2.469700244997, 2.571899743493, 2.854100080278, LiH  
     6  3.010500234607, 3.160499843388, 3.375200005967, 3.512400279255, LiH  
     7  3.655599752432, 3.790699992730, 3.915799993246, 3.966500192406, LiH  
     8  3.986799704971, 4.000000000000,     60*0.0D+00/                 LiH  
      DATA  K_LiH/                                                      071215
     1 -8.77919083D-05, 1.44509553D-01, 2.93751423D-01, 7.44474693D-01, LiH  
     2  1.76536738D+00, 2.76854865D+00, 3.71603542D+00, 4.24479674D+00, LiH  
     3  4.71919477D+00, 5.50711388D+00, 6.18989668D+00, 6.73929936D+00, LiH  
     4  7.06881351D+00, 7.37207894D+00, 7.65098489D+00, 7.91338032D+00, LiH  
     5  8.36163462D+00, 8.56012992D+00, 8.73215082D+00, 9.16245435D+00, LiH  
     6  9.36362313D+00, 9.52746843D+00, 9.71068947D+00, 9.79412216D+00, LiH  
     7  9.85954672D+00, 9.92565695D+00, 1.00274749D+01, 1.00912664D+01, LiH  
     8  1.01212468D+01, 1.01420812D+01,     60*0.0D+00/                 LiH  
      DATA TK_BeH/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751300020874, BeH  
     2  0.831899962688, 0.921800042222, 1.033800119949, 1.178999859119, BeH  
     3  1.324700055595, 1.481099988861, 1.639399877072, 1.753000062773, BeH  
     4  1.868599919367, 2.093600229698, 2.326299847049, 2.564200113612, BeH  
     5  2.681600192523, 2.786399886844, 3.042500031269, 3.208599990934, BeH  
     6  3.368199851521, 3.651199649764, 3.717500197165, 3.784999860172, BeH  
     7  3.897199825981, 3.961200073629, 3.984999664560, 4.000000000000, BeH  
     8      62*0.0D+00/                                                 BeH  
      DATA  K_BeH/                                                      071215
     1  2.06258583D-05, 1.41500452D-01, 2.87519670D-01, 7.26168520D-01, BeH  
     2  1.73006053D+00, 2.67924082D+00, 3.65323882D+00, 4.64667111D+00, BeH  
     3  5.41931103D+00, 6.07309398D+00, 6.60544335D+00, 6.93002445D+00, BeH  
     4  7.22384967D+00, 7.72078712D+00, 8.16688419D+00, 8.58065926D+00, BeH  
     5  8.77392285D+00, 8.94028793D+00, 9.31150454D+00, 9.51520825D+00, BeH  
     6  9.67897623D+00, 9.88580886D+00, 9.91713959D+00, 9.94389586D+00, BeH  
     7  9.98664005D+00, 1.00175923D+01, 1.00317165D+01, 1.00416097D+01, BeH  
     8      62*0.0D+00/                                                 BeH  
      DATA TK_BH/                                                       071215
     1  0.699999789529, 0.709700041799, 0.720299804015, 0.752800059749, BH   
     2  0.835500045214, 0.929699851738, 1.047600122923, 1.191200014866, BH   
     3  1.343200084633, 1.504500097980, 1.668800039432, 1.783800125766, BH   
     4  1.901300049249, 2.124599977753, 2.351800255316, 2.591500191880, BH   
     5  2.714800164222, 2.826499857924, 3.120699888668, 3.437499762377, BH   
     6  3.724699925410, 3.834800031278, 3.928800296672, 3.972600084024, BH   
     7  3.988499743138, 4.000000000000,     64*0.0D+00/                 BH   
      DATA  K_BH/                                                       071215
     1 -1.17453728D-04, 1.51357957D-01, 3.13756941D-01, 7.92076168D-01, BH   
     2  1.88783751D+00, 2.95310650D+00, 4.06055315D+00, 5.13821952D+00, BH   
     3  6.03153562D+00, 6.77611268D+00, 7.38010455D+00, 7.73523704D+00, BH   
     4  8.05506696D+00, 8.57669343D+00, 9.03085589D+00, 9.45986337D+00, BH   
     5  9.66705516D+00, 9.84746522D+00, 1.02756784D+01, 1.06263754D+01, BH   
     6  1.08217290D+01, 1.08587367D+01, 1.08790698D+01, 1.08886837D+01, BH   
     7  1.08927979D+01, 1.08961151D+01,     64*0.0D+00/                 BH   
      DATA TK_CH/                                                       071215
     1  0.699999789529, 0.710600034435, 0.723799887835, 0.760900222609, CH   
     2  0.855000021789, 0.976199920554, 1.109100120098, 1.217699849316, CH   
     3  1.320600148951, 1.546100101598, 1.796400104248, 1.962699985157, CH   
     4  2.149899629212, 2.342500220445, 2.540199961841, 2.726599804056, CH   
     5  2.902899671516, 3.201599834065, 3.307300349723, 3.420600104942, CH   
     6  3.593600239282, 3.763300303961, 3.889000332808, 4.000000000000, CH   
     7      66*0.0D+00/                                                 CH   
      DATA  K_CH/                                                       071215
     1 -2.06404210D-01,-6.86456104D-02, 9.93421846D-02, 5.51709351D-01, CH   
     2  1.58331512D+00, 2.71525183D+00, 3.76483817D+00, 4.50734424D+00, CH   
     3  5.13127073D+00, 6.26490139D+00, 7.21431541D+00, 7.71253253D+00, CH   
     4  8.18538446D+00, 8.60427889D+00, 8.98673665D+00, 9.31798049D+00, CH   
     5  9.61116187D+00, 1.00495458D+01, 1.01780677D+01, 1.02958479D+01, CH   
     6  1.04345796D+01, 1.05255272D+01, 1.05651613D+01, 1.05855448D+01, CH   
     7      66*0.0D+00/                                                 CH   
      DATA TK_NH/                                                       071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751400023466, NH   
     2  0.831299948934, 0.930699858142, 1.038499994372, 1.191300012428, NH   
     3  1.342700097068, 1.485099896078, 1.626900042626, 1.750600002301, NH   
     4  1.878900104199, 1.998499965912, 2.120599893858, 2.338300123915, NH   
     5  2.591800198310, 2.761600272041, 2.919700068330, 3.188399679884, NH   
     6  3.442899645800, 3.605200004372, 3.754400091461, 3.886400271093, NH   
     7  4.000000000000,     65*0.0D+00/                                 NH   
      DATA  K_NH/                                                       071215
     1  4.33593273D-05, 1.45343376D-01, 2.95357908D-01, 7.47796812D-01, NH   
     2  1.77510648D+00, 2.85878297D+00, 3.83338714D+00, 4.92849274D+00, NH   
     3  5.75877125D+00, 6.37343561D+00, 6.87263002D+00, 7.24162732D+00, NH   
     4  7.57624627D+00, 7.85509391D+00, 8.11502022D+00, 8.53457357D+00, NH   
     5  8.97752687D+00, 9.25717528D+00, 9.50796889D+00, 9.90193656D+00, NH   
     6  1.02094949D+01, 1.03615263D+01, 1.04673243D+01, 1.05345073D+01, NH   
     7  1.05788305D+01,     65*0.0D+00/                                 NH   
      DATA TK_OH/                                                       071215
     1  0.699999789529, 0.709800044400, 0.720399806410, 0.752900062340, OH   
     2  0.835600047507, 0.934099924381, 1.045700080766, 1.121499917340, OH   
     3  1.197899851503, 1.334100028928, 1.469599790058, 1.692700002154, OH   
     4  1.790799989126, 1.906299923921, 2.010500241392, 2.117499830425, OH   
     5  2.300900216672, 2.605799950606, 2.767700419388, 2.930800236047, OH   
     6  3.078099897369, 3.225899829198, 3.453599880650, 3.637100268406, OH   
     7  3.795500093363, 3.914699964983, 3.966800199129, 4.000000000000, OH   
     8      62*0.0D+00/                                                 OH   
      DATA  K_OH/                                                       071215
     1 -3.85549856D-05, 1.37146049D-01, 2.82700371D-01, 7.11385773D-01, OH   
     2  1.69338871D+00, 2.68961336D+00, 3.63612930D+00, 4.19033126D+00, OH   
     3  4.68938630D+00, 5.45759194D+00, 6.09679134D+00, 6.94117167D+00, OH   
     4  7.25090481D+00, 7.58179550D+00, 7.85678031D+00, 8.12145764D+00, OH   
     5  8.54199801D+00, 9.16562983D+00, 9.46585027D+00, 9.75100869D+00, OH   
     6  9.99336067D+00, 1.02178045D+01, 1.05147601D+01, 1.07055156D+01, OH   
     7  1.08352772D+01, 1.09099034D+01, 1.09359481D+01, 1.09506503D+01, OH   
     8      62*0.0D+00/                                                 OH   
      DATA TK_HF/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, HF   
     2  0.828999943334, 0.925899943363, 1.033000141324, 1.166799995515, HF   
     3  1.326200021441, 1.473699866856, 1.609300122704, 1.727699967543, HF   
     4  1.846100001191, 1.967299865391, 2.091200169343, 2.219400237644, HF   
     5  2.345700292422, 2.464900136736, 2.582899993390, 2.801300231891, HF   
     6  2.976299811170, 3.243099845136, 3.463100102480, 3.641500222226, HF   
     7  3.795800099653, 3.917200029217, 3.967600217058, 4.000000000000, HF   
     8      62*0.0D+00/                                                 HF   
      DATA  K_HF/                                                       071215
     1  1.69923625D-04, 1.54562320D-01, 3.10892585D-01, 7.89267882D-01, HF   
     2  1.87425083D+00, 3.01279113D+00, 4.06080152D+00, 5.12180365D+00, HF   
     3  6.10232050D+00, 6.80252087D+00, 7.32085493D+00, 7.70258674D+00, HF   
     4  8.03528053D+00, 8.33769978D+00, 8.61768074D+00, 8.88607382D+00, HF   
     5  9.13705420D+00, 9.36621953D+00, 9.58762801D+00, 9.98397733D+00, HF   
     6  1.02875457D+01, 1.07142154D+01, 1.10119684D+01, 1.12074024D+01, HF   
     7  1.13408662D+01, 1.14201874D+01, 1.14458840D+01, 1.14602246D+01, HF   
     8      62*0.0D+00/                                                 HF   
      DATA TK_NaH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750199992366, NaH  
     2  0.828799948174, 0.922300030166, 1.027900169619, 1.101699949507, NaH  
     3  1.176699912495, 1.320700146674, 1.472699843440, 1.626000022728, NaH  
     4  1.749099968574, 1.879600121517, 2.015800366136, 2.160199862683, NaH  
     5  2.393600076886, 2.616699798212, 2.861499680085, 3.153699687729, NaH  
     6  3.335100041985, 3.470300230258, 3.534499832192, 3.597600327668, NaH  
     7  3.650399631097, 3.700699823647, 3.749499976061, 3.796200108039, NaH  
     8  3.851900235550, 3.906199758885, 3.964800154308, 3.986399695991, NaH  
     9  4.000000000000,     57*0.0D+00/                                 NaH  
      DATA  K_NaH/                                                      071215
     1 -3.27547461D-01,-1.79925415D-01,-2.90507026D-02, 4.28321261D-01, NaH  
     2  1.46467543D+00, 2.50035440D+00, 3.46239651D+00, 4.02904205D+00, NaH  
     3  4.53181883D+00, 5.33187120D+00, 5.99837634D+00, 6.53941947D+00, NaH  
     4  6.90537504D+00, 7.24463405D+00, 7.55903101D+00, 7.85999468D+00, NaH  
     5  8.29788257D+00, 8.67328270D+00, 9.02795467D+00, 9.35391427D+00, NaH  
     6  9.50188463D+00, 9.58326430D+00, 9.61245818D+00, 9.63586099D+00, NaH  
     7  9.65294237D+00, 9.66975008D+00, 9.69113107D+00, 9.72359104D+00, NaH  
     8  9.79219482D+00, 9.90548577D+00, 1.00811779D+01, 1.01563550D+01, NaH  
     9  1.02055401D+01,     57*0.0D+00/                                 NaH  
      DATA TK_MgH/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720999820779, 0.754400101215, MgH  
     2  0.839300132325, 0.942499977946, 1.057500008381, 1.133100045788, MgH  
     3  1.207699964953, 1.347299982672, 1.496999929244, 1.658899861172, MgH  
     4  1.839300131896, 2.027499911287, 2.207799996663, 2.398099746110, MgH  
     5  2.517500404499, 2.629800098812, 2.902399660273, 3.035599877555, MgH  
     6  3.152399658121, 3.340300170960, 3.456599953194, 3.555300307149, MgH  
     7  3.635600237483, 3.725599863046, 3.790799994826, 3.857499840731, MgH  
     8  3.903099689268, 3.946599751612, 3.979199609788, 3.990199781278, MgH  
     9  4.000000000000,     57*0.0D+00/                                 MgH  
      DATA  K_MgH/                                                      071215
     1  5.37199455D-05, 1.36309384D-01, 2.85818709D-01, 7.15766763D-01, MgH  
     2  1.68598295D+00, 2.66464008D+00, 3.55010765D+00, 4.03876595D+00, MgH  
     3  4.46195209D+00, 5.12768077D+00, 5.70296742D+00, 6.21087470D+00, MgH  
     4  6.68077225D+00, 7.09901232D+00, 7.45440398D+00, 7.79796794D+00, MgH  
     5  8.00152025D+00, 8.18511218D+00, 8.58858984D+00, 8.75668515D+00, MgH  
     6  8.88597541D+00, 9.05713496D+00, 9.13783890D+00, 9.18771028D+00, MgH  
     7  9.21484120D+00, 9.23484578D+00, 9.24841017D+00, 9.26961987D+00, MgH  
     8  9.29444380D+00, 9.33212898D+00, 9.37339531D+00, 9.39042095D+00, MgH  
     9  9.40703946D+00,     57*0.0D+00/                                 MgH  
      DATA TK_AlH/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720999820779, 0.754300098624, AlH  
     2  0.839000125448, 0.936699975035, 1.056500030753, 1.193899949033, AlH  
     3  1.339200145869, 1.492899843246, 1.656299930354, 1.784600106013, AlH  
     4  1.908199876296, 2.180300285067, 2.440699600663, 2.550200186218, AlH  
     5  2.665299978029, 2.837400108999, 3.012800287320, 3.323399771665, AlH  
     6  3.431000222432, 3.551200221584, 3.642900123327, 3.737799719293, AlH  
     7  3.820599699474, 3.899199678602, 3.960100048977, 3.984799660070, AlH  
     8  4.000000000000,     61*0.0D+00/                                 AlH  
      DATA  K_AlH/                                                      071215
     1  1.04796459D-04, 1.50875166D-01, 3.16290108D-01, 7.90455880D-01, AlH  
     2  1.86025408D+00, 2.88710152D+00, 3.90599951D+00, 4.82732559D+00, AlH  
     3  5.59166538D+00, 6.23573273D+00, 6.80209046D+00, 7.19702496D+00, AlH  
     4  7.55143629D+00, 8.25453042D+00, 8.83096202D+00, 9.04917599D+00, AlH  
     5  9.26489638D+00, 9.56020797D+00, 9.82268709D+00, 1.01830598D+01, AlH  
     6  1.02776225D+01, 1.03648713D+01, 1.04173507D+01, 1.04588456D+01, AlH  
     7  1.04886690D+01, 1.05223535D+01, 1.05644674D+01, 1.05886801D+01, AlH  
     8  1.06062416D+01,     61*0.0D+00/                                 AlH  
      DATA TK_SiH/                                                      071215
     1  0.699999789529, 0.710700031907, 0.724099895020, 0.761600204146, SiH  
     2  0.858000090607, 0.974299882107, 1.100499921844, 1.354000004287, SiH  
     3  1.440400033234, 1.531100204249, 1.681399870353, 1.846599989127, SiH  
     4  2.021700339484, 2.149999622226, 2.288399914247, 2.544600060508, SiH  
     5  2.697399749685, 3.003800084796, 3.207299961801, 3.423500163025, SiH  
     6  3.557400350975, 3.686199865437, 3.811000334720, 3.887000285335, SiH  
     7  3.945599729033, 3.979799566676, 4.000000000000,     63*0.0D+00/ SiH  
      DATA  K_SiH/                                                      071215
     1  9.90292283D-05, 1.33128080D-01, 2.96036151D-01, 7.31398817D-01, SiH  
     2  1.72536553D+00, 2.72343022D+00, 3.60839738D+00, 4.93811812D+00, SiH  
     3  5.29898808D+00, 5.64734863D+00, 6.17862019D+00, 6.71674498D+00, SiH  
     4  7.24606705D+00, 7.60817151D+00, 7.97230574D+00, 8.57223255D+00, SiH  
     5  8.88828018D+00, 9.43233289D+00, 9.72051686D+00, 9.96656060D+00, SiH  
     6  1.00935549D+01, 1.01970219D+01, 1.02739111D+01, 1.03081015D+01, SiH  
     7  1.03307292D+01, 1.03447289D+01, 1.03541779D+01,     63*0.0D+00/ SiH  
      DATA TK_PH/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, PH   
     2  0.829799923973, 0.920300078390, 1.028400181967, 1.169800063450, PH   
     3  1.313700033703, 1.469199799898, 1.626900042626, 1.742599834178, PH   
     4  1.861700093567, 2.090800159283, 2.324199801577, 2.581499960589, PH   
     5  2.709700048611, 2.832700003193, 3.061600320857, 3.190699580254, PH   
     6  3.320899717459, 3.454999914504, 3.601000308329, 3.717400195008, PH   
     7  3.852500193248, 4.000000000000,     64*0.0D+00/                 PH   
      DATA  K_PH/                                                       071215
     1  4.00286740D-05, 1.43050839D-01, 2.87776308D-01, 7.30024558D-01, PH   
     2  1.73742051D+00, 2.70859692D+00, 3.66387754D+00, 4.65291084D+00, PH   
     3  5.43644184D+00, 6.10391503D+00, 6.64726901D+00, 6.98460857D+00, PH   
     4  7.29188498D+00, 7.80249665D+00, 8.25227211D+00, 8.70031218D+00, PH   
     5  8.91130692D+00, 9.10652749D+00, 9.44237348D+00, 9.60921493D+00, PH   
     6  9.75704188D+00, 9.88687073D+00, 1.00066652D+01, 1.00924955D+01, PH   
     7  1.01842930D+01, 1.02732779D+01,     64*0.0D+00/                 PH   
      DATA TK_HS/                                                       071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751500026057, HS   
     2  0.831699958103, 0.926099938541, 1.038200002388, 1.179399849837, HS   
     3  1.336300079373, 1.488199824171, 1.641999909301, 1.760600222913, HS   
     4  1.878400091830, 2.134300189844, 2.343800249685, 2.554700290139, HS   
     5  2.813800143764, 3.055100312977, 3.225799836454, 3.390900291992, HS   
     6  3.502800055351, 3.610899678829, 3.710500046173, 3.793400049336, HS   
     7  3.912699913596, 3.965600172236, 4.000000000000,     63*0.0D+00/ HS   
      DATA  K_HS/                                                       071215
     1  5.48586800D-05, 1.41463901D-01, 2.87480835D-01, 7.29350106D-01, HS   
     2  1.73399478D+00, 2.74380409D+00, 3.74598420D+00, 4.76805379D+00, HS   
     3  5.66180477D+00, 6.34529994D+00, 6.90446123D+00, 7.26837386D+00, HS   
     4  7.58647769D+00, 8.17693018D+00, 8.60500086D+00, 9.01780779D+00, HS   
     5  9.50562773D+00, 9.92105819D+00, 1.01759117D+01, 1.03852404D+01, HS   
     6  1.05071795D+01, 1.06111461D+01, 1.06947709D+01, 1.07537907D+01, HS   
     7  1.08177953D+01, 1.08382535D+01, 1.08496509D+01,     63*0.0D+00/ HS   
      DATA TK_HCl/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.749999987182, HCl  
     2  0.828499955435, 0.917100020700, 1.026000122698, 1.167000000044, HCl  
     3  1.306800026093, 1.456399947978, 1.608000090359, 1.721599826051, HCl  
     4  1.838500113067, 1.948899829084, 2.059300412495, 2.318599811150, HCl  
     5  2.589800155049, 2.711400087609, 2.831399973928, 3.149899610922, HCl  
     6  3.280499755484, 3.419800087662, 3.610099659420, 3.790899996923, HCl  
     7  3.913599936720, 3.966400190165, 4.000000000000,     63*0.0D+00/ HCl  
      DATA  K_HCl/                                                      071215
     1  1.73411817D-04, 1.52774103D-01, 3.07319671D-01, 7.80223676D-01, HCl  
     2  1.85955047D+00, 2.89340965D+00, 3.94155714D+00, 5.01128903D+00, HCl  
     3  5.83398979D+00, 6.52746407D+00, 7.09045752D+00, 7.44559427D+00, HCl  
     4  7.76668064D+00, 8.03813091D+00, 8.28577750D+00, 8.80245577D+00, HCl  
     5  9.28893959D+00, 9.50106581D+00, 9.70913458D+00, 1.02342924D+01, HCl  
     6  1.04233437D+01, 1.06009557D+01, 1.08013974D+01, 1.09462068D+01, HCl  
     7  1.10163367D+01, 1.10387664D+01, 1.10506946D+01,     63*0.0D+00/ HCl  
      DATA TK_KH/                                                       071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.750099989774, KH   
     2  0.828499955435, 0.921700044633, 1.027500159741, 1.100599924149, KH   
     3  1.174799956587, 1.317300107366, 1.467699836798, 1.619899892574, KH   
     4  1.746399912748, 1.880200126046, 2.016600384965, 2.172000106739, KH   
     5  2.296800116019, 2.417700037148, 2.672200131714, 2.884700220124, KH   
     6  3.131800161738, 3.315799996648, 3.438199712833, 3.497599934061, KH   
     7  3.556100323845, 3.607599830683, 3.656199766432, 3.704399907888, KH   
     8  3.743299838148, 3.809300391570, 3.843400221029, 3.877400060365, KH   
     9  3.954499926640, 3.983099621903, 4.000000000000,     55*0.0D+00/ KH   
      DATA  K_KH/                                                       071215
     1 -2.96911935D-06, 1.30759474D-01, 2.63093070D-01, 6.68712833D-01, KH   
     2  1.58758159D+00, 2.50914200D+00, 3.37366474D+00, 3.87958619D+00, KH   
     3  4.33033200D+00, 5.05360060D+00, 5.66312978D+00, 6.16561333D+00, KH   
     4  6.52037149D+00, 6.85038195D+00, 7.15146726D+00, 7.46299757D+00, KH   
     5  7.69516955D+00, 7.90749063D+00, 8.31120410D+00, 8.59208855D+00, KH   
     6  8.84778674D+00, 8.99058228D+00, 9.06230762D+00, 9.09025543D+00, KH   
     7  9.11445011D+00, 9.13507256D+00, 9.15727808D+00, 9.18709550D+00, KH   
     8  9.22189464D+00, 9.31737760D+00, 9.39026686D+00, 9.48027876D+00, KH   
     9  9.73845531D+00, 9.84708755D+00, 9.91305613D+00,     55*0.0D+00/ KH   
      DATA TK_CaH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749599978008, CaH  
     2  0.827299984476, 0.919700078907, 1.024400083186, 1.170500056376, CaH  
     3  1.311499988688, 1.460800006540, 1.611400105114, 1.733199943941, CaH  
     4  1.862900063272, 2.000900020666, 2.147299810856, 2.276200013579, CaH  
     5  2.401399640604, 2.523900167125, 2.637600283470, 2.770300452438, CaH  
     6  2.924600175266, 3.088700113917, 3.278499851475, 3.393400109454, CaH  
     7  3.505900127623, 3.626000033067, 3.739499754414, 3.827899875959, CaH  
     8  3.886200266346, 3.923400176698, 3.953799911355, 3.982599610678, CaH  
     9  4.000000000000,     57*0.0D+00/                                 CaH  
      DATA  K_CaH/                                                      071215
     1  1.45114700D-05, 1.26622276D-01, 2.52139725D-01, 6.41369451D-01, CaH  
     2  1.52473748D+00, 2.41257102D+00, 3.24599471D+00, 4.17341501D+00, CaH  
     3  4.87568189D+00, 5.47174074D+00, 5.96336809D+00, 6.30283387D+00, CaH  
     4  6.62239588D+00, 6.92717029D+00, 7.22156784D+00, 7.46280531D+00, CaH  
     5  7.68502941D+00, 7.89269732D+00, 8.07637159D+00, 8.27710222D+00, CaH  
     6  8.48765671D+00, 8.68147375D+00, 8.86759647D+00, 8.96216585D+00, CaH  
     7  9.04424018D+00, 9.13042832D+00, 9.23565785D+00, 9.36176042D+00, CaH  
     8  9.48188696D+00, 9.57792885D+00, 9.66766830D+00, 9.76093439D+00, CaH  
     9  9.82051560D+00,     57*0.0D+00/                                 CaH  
      DATA TK_TiH/                                                      071215
     1  0.699999789529, 0.710200044546, 0.722099847123, 0.756800163415, TiH  
     2  0.845100025319, 0.952399856601, 1.076100022978, 1.209700007721, TiH  
     3  1.346000015001, 1.500700008629, 1.691500033031, 1.871099911235, TiH  
     4  2.028799815312, 2.173600146707, 2.316099990646, 2.427900243370, TiH  
     5  2.533399796833, 2.736999709519, 2.876800037287, 3.022400282638, TiH  
     6  3.193799652906, 3.361299668566, 3.517200391350, 3.647799777178, TiH  
     7  3.746999920451, 3.851100291953, 3.942199652266, 3.977399739125, TiH  
     8  3.989699770079, 4.000000000000,     60*0.0D+00/                 TiH  
      DATA  K_TiH/                                                      071215
     1 -6.19873763D-05, 1.41730854D-01, 3.03725236D-01, 7.56086149D-01, TiH  
     2  1.78544099D+00, 2.83837757D+00, 3.83941462D+00, 4.72205985D+00, TiH  
     3  5.45721094D+00, 6.13297588D+00, 6.78681770D+00, 7.27096440D+00, TiH  
     4  7.62941207D+00, 7.92786961D+00, 8.20625750D+00, 8.41635059D+00, TiH  
     5  8.60590114D+00, 8.93257742D+00, 9.11291405D+00, 9.25516726D+00, TiH  
     6  9.36642642D+00, 9.43255607D+00, 9.47954105D+00, 9.52909735D+00, TiH  
     7  9.58897246D+00, 9.68692502D+00, 9.80934707D+00, 9.86627529D+00, TiH  
     8  9.88739956D+00, 9.90557202D+00,     60*0.0D+00/                 TiH  
      DATA TK_CrH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749799982595, CrH  
     2  0.827799972376, 0.920200080801, 1.024900095533, 1.097199981936, CrH  
     3  1.170900047093, 1.314100041888, 1.464399917979, 1.616599975089, CrH  
     4  1.741899819704, 1.874399992874, 2.009400215846, 2.145999901677, CrH  
     5  2.387900296372, 2.607399833796, 2.951899854308, 3.086000061480, CrH  
     6  3.212300076353, 3.418300048193, 3.504500094984, 3.598100338716, CrH  
     7  3.696399733213, 3.812500226428, 3.869099872398, 3.921200127820, CrH  
     8  3.969400257397, 3.987599722932, 4.000000000000,     59*0.0D+00/ CrH  
      DATA  K_CrH/                                                      071215
     1  5.35580450D-05, 1.38265403D-01, 2.76715324D-01, 7.02836877D-01, CrH  
     2  1.66922947D+00, 2.63377147D+00, 3.53543208D+00, 4.06210479D+00, CrH  
     3  4.53273748D+00, 5.29396804D+00, 5.92838630D+00, 6.44891206D+00, CrH  
     4  6.81127827D+00, 7.14708661D+00, 7.45214888D+00, 7.73286064D+00, CrH  
     5  8.18292025D+00, 8.55648747D+00, 9.07458187D+00, 9.24340203D+00, CrH  
     6  9.38362746D+00, 9.58682199D+00, 9.66920736D+00, 9.75690438D+00, CrH  
     7  9.84496524D+00, 9.95566767D+00, 1.00230337D+01, 1.00982758D+01, CrH  
     8  1.01807648D+01, 1.02151251D+01, 1.02395133D+01,     59*0.0D+00/ CrH  
      DATA TK_MnH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750299994957, MnH  
     2  0.828999943334, 0.922600022932, 1.028100174558, 1.101999956423, MnH  
     3  1.177099903212, 1.322200112519, 1.475199901980, 1.629500100109, MnH  
     4  1.751300019938, 1.879600121517, 2.011100255514, 2.144699992499, MnH  
     5  2.379900124190, 2.588800131620, 2.929700285693, 3.096700300603, MnH  
     6  3.287599918695, 3.443199652045, 3.573799777326, 3.631600155023, MnH  
     7  3.687399781860, 3.806400328000, 3.867499835839, 3.920800118933, MnH  
     8  3.969400257397, 3.987599722932, 4.000000000000,     59*0.0D+00/ MnH  
      DATA  K_MnH/                                                      071215
     1  1.73966005D-05, 1.34552432D-01, 2.72090700D-01, 6.90575078D-01, MnH  
     2  1.63803356D+00, 2.58629844D+00, 3.46757707D+00, 3.98937337D+00, MnH  
     3  4.45360217D+00, 5.19989098D+00, 5.82484532D+00, 6.33616271D+00, MnH  
     4  6.67923483D+00, 6.99832978D+00, 7.29176413D+00, 7.56408431D+00, MnH  
     5  7.99997660D+00, 8.35520862D+00, 8.86680730D+00, 9.07261006D+00, MnH  
     6  9.26589842D+00, 9.39013925D+00, 9.47148134D+00, 9.50254602D+00, MnH  
     7  9.53264473D+00, 9.61814766D+00, 9.68792558D+00, 9.77021903D+00, MnH  
     8  9.86455781D+00, 9.90462236D+00, 9.93336739D+00,     59*0.0D+00/ MnH  
      DATA TK_FeH/                                                      071215
     1  0.699999789529, 0.710300042018, 0.722599859097, 0.758100197107, FeH  
     2  0.848799936046, 0.958300010975, 1.084100020987, 1.237999804478, FeH  
     3  1.413499930887, 1.572000122125, 1.726799946667, 1.962100000779, FeH  
     4  2.094700257361, 2.227099750895, 2.331399960353, 2.434999934397, FeH  
     5  2.548200141235, 2.649599644810, 2.890900282807, 2.990899787028, FeH  
     6  3.091700180127, 3.194399666967, 3.285199863526, 3.416099990304, FeH  
     7  3.633000183884, 3.739799760611, 3.859399706775, 3.944099695165, FeH  
     8  3.978399667271, 3.989899774569, 4.000000000000,     59*0.0D+00/ FeH  
      DATA  K_FeH/                                                      071215
     1  9.92176001D-06, 1.52631280D-01, 3.30967249D-01, 8.23022513D-01, FeH  
     2  1.94204100D+00, 3.07108573D+00, 4.13184848D+00, 5.16589607D+00, FeH  
     3  6.08015766D+00, 6.72816022D+00, 7.24219569D+00, 7.85740776D+00, FeH  
     4  8.13569206D+00, 8.37803628D+00, 8.55163741D+00, 8.71404396D+00, FeH  
     5  8.88277689D+00, 9.02575890D+00, 9.31788219D+00, 9.41057780D+00, FeH  
     6  9.48369255D+00, 9.53754577D+00, 9.57049539D+00, 9.60255661D+00, FeH  
     7  9.65045078D+00, 9.68861514D+00, 9.76271550D+00, 9.84817543D+00, FeH  
     8  9.89290814D+00, 9.90934333D+00, 9.92439321D+00,     59*0.0D+00/ FeH  
      DATA TK_CoH/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718399837274, 0.748399950483, CoH  
     2  0.824100061921, 0.912499917719, 1.015000093981, 1.082999993933, CoH  
     3  1.152400035864, 1.288999974227, 1.432399870049, 1.578399957810, CoH  
     4  1.703499896669, 1.836100056579, 1.970799811976, 2.096000290053, CoH  
     5  2.405899750391, 2.509400231817, 2.620499885638, 2.777799889749, CoH  
     6  2.974299957098, 3.110599670140, 3.275700052112, 3.449599785270, CoH  
     7  3.588000112877, 3.710300041859, 3.829999926729, 3.932800121525, CoH  
     8  3.973700004984, 4.000000000000,     60*0.0D+00/                 CoH  
      DATA  K_CoH/                                                      071215
     1 -1.32577925D-05, 1.56080380D-01, 3.09288081D-01, 7.88991437D-01, CoH  
     2  1.87542316D+00, 2.94835951D+00, 3.97599292D+00, 4.55287356D+00, CoH  
     3  5.07026628D+00, 5.91815663D+00, 6.62200351D+00, 7.19744992D+00, CoH  
     4  7.60825473D+00, 7.98296379D+00, 8.31623972D+00, 8.59436945D+00, CoH  
     5  9.20009843D+00, 9.38881777D+00, 9.58983608D+00, 9.87431206D+00, CoH  
     6  1.02229474D+01, 1.04533751D+01, 1.07168373D+01, 1.09753667D+01, CoH  
     7  1.11664653D+01, 1.13263003D+01, 1.14813578D+01, 1.16234750D+01, CoH  
     8  1.16850062D+01, 1.17266112D+01,     60*0.0D+00/                 CoH  
      DATA TK_NiH/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720199801620, 0.752500051974, NiH  
     2  0.834500022290, 0.934299928278, 1.045600078547, 1.121699922148, NiH  
     3  1.195999897830, 1.335500061029, 1.500299999223, 1.671800026305, NiH  
     4  1.784900098605, 1.893099930912, 1.997299938641, 2.102500212454, NiH  
     5  2.210300060825, 2.311200342459, 2.573499779434, 2.781799770298, NiH  
     6  2.958300001757, 3.395099985328, 3.668500028611, 3.731699593274, NiH  
     7  3.799300173031, 3.909199826257, 3.964900156549, 3.986299693746, NiH  
     8  4.000000000000,     61*0.0D+00/                                 NiH  
      DATA  K_NiH/                                                      071215
     1  7.19439204D-05, 1.37236632D-01, 2.82909889D-01, 7.13433804D-01, NiH  
     2  1.69785344D+00, 2.71714304D+00, 3.66786760D+00, 4.22625668D+00, NiH  
     3  4.71138427D+00, 5.48844690D+00, 6.22707854D+00, 6.84259672D+00, NiH  
     4  7.18652767D+00, 7.48271130D+00, 7.74610865D+00, 7.99648167D+00, NiH  
     5  8.24075972D+00, 8.45974805D+00, 8.98564410D+00, 9.35623845D+00, NiH  
     6  9.63843860D+00, 1.02053330D+01, 1.04475437D+01, 1.04888068D+01, NiH  
     7  1.05280505D+01, 1.05908160D+01, 1.06295376D+01, 1.06469717D+01, NiH  
     8  1.06590873D+01,     61*0.0D+00/                                 NiH  
      DATA TK_CuH/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, CuH  
     2  0.829599928813, 0.921200056689, 1.028200177028, 1.097799966589, CuH  
     3  1.169300052128, 1.312400007103, 1.469099802358, 1.627900064735, CuH  
     4  1.744199867260, 1.865200005205, 1.984899891769, 2.107099884583, CuH  
     5  2.343500242938, 2.582299979332, 2.764400339676, 2.949399797328, CuH  
     6  3.092900209041, 3.239099752724, 3.381300137396, 3.517100389015, CuH  
     7  3.583300002761, 3.653399701098, 3.809400393762, 3.885600252105, CuH  
     8  3.953299900437, 3.982299603942, 4.000000000000,     59*0.0D+00/ CuH  
      DATA  K_CuH/                                                      071215
     1  7.38064273D-05, 1.46708058D-01, 2.95067914D-01, 7.48206979D-01, CuH  
     2  1.77673444D+00, 2.78037247D+00, 3.74464922D+00, 4.27398623D+00, CuH  
     3  4.75175631D+00, 5.54671296D+00, 6.23156577D+00, 6.78654398D+00, CuH  
     4  7.12947582D+00, 7.44432450D+00, 7.72378594D+00, 7.98422090D+00, CuH  
     5  8.43817648D+00, 8.85349315D+00, 9.14879365D+00, 9.42451615D+00, CuH  
     6  9.61472916D+00, 9.78343573D+00, 9.92363550D+00, 1.00417579D+01, CuH  
     7  1.00972108D+01, 1.01559229D+01, 1.02804261D+01, 1.03310016D+01, CuH  
     8  1.03708270D+01, 1.03884664D+01, 1.04001052D+01,     59*0.0D+00/ CuH  
      DATA TK_ZnH/                                                      071215
     1  0.699999789529, 0.710000049601, 0.721199825569, 0.754700108990, ZnH  
     2  0.840000148372, 0.943199960759, 1.058399988246, 1.133700031093, ZnH  
     3  1.208299977784, 1.350199919965, 1.500299999223, 1.662899900293, ZnH  
     4  1.845800008430, 2.038299928148, 2.229799561359, 2.423800165291, ZnH  
     5  2.622499931482, 2.910299838513, 3.183100064248, 3.381000129799, ZnH  
     6  3.524100158687, 3.661899893884, 3.834900033457, 3.887300292456, ZnH  
     7  3.938799689081, 3.975999839721, 3.989399763344, 4.000000000000, ZnH  
     8      62*0.0D+00/                                                 ZnH  
      DATA  K_ZnH/                                                      071215
     1  2.12832761D-05, 1.39268016D-01, 2.91885202D-01, 7.28217619D-01, ZnH  
     2  1.71451887D+00, 2.70416141D+00, 3.60006087D+00, 4.09124019D+00, ZnH  
     3  4.51815829D+00, 5.19912289D+00, 5.77796877D+00, 6.28923598D+00, ZnH  
     4  6.76579199D+00, 7.19260828D+00, 7.56787309D+00, 7.91557991D+00, ZnH  
     5  8.24754301D+00, 8.67673595D+00, 8.99683871D+00, 9.16299423D+00, ZnH  
     6  9.24240320D+00, 9.28886502D+00, 9.33523365D+00, 9.35446698D+00, ZnH  
     7  9.37931231D+00, 9.40309235D+00, 9.41336125D+00, 9.42229145D+00, ZnH  
     8      62*0.0D+00/                                                 ZnH  
      DATA TK_GaH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.749999987182, GaH  
     2  0.828199962695, 0.920500073567, 1.025200102942, 1.098599946127, GaH  
     3  1.173499986756, 1.318500131920, 1.468299822038, 1.623199960822, GaH  
     4  1.771099996137, 1.917800021877, 2.214200136605, 2.403799699157, GaH  
     5  2.521500345349, 2.620299881053, 2.717700229569, 2.846300304210, GaH  
     6  2.992199817452, 3.104000090992, 3.208199981970, 3.416399998198, GaH  
     7  3.542800024328, 3.680400269394, 3.742099811455, 3.804100277582, GaH  
     8  3.861499698746, 3.910499857069, 3.965400167754, 3.986399695991, GaH  
     9  4.000000000000,     57*0.0D+00/                                 GaH  
      DATA  K_GaH/                                                      071215
     1 -9.81381218D-05, 1.43071357D-01, 2.87951368D-01, 7.30521751D-01, GaH  
     2  1.73270041D+00, 2.72839545D+00, 3.65931682D+00, 4.21019473D+00, GaH  
     3  4.70146558D+00, 5.48998360D+00, 6.13467407D+00, 6.67277041D+00, GaH  
     4  7.10042921D+00, 7.46579896D+00, 8.08956546D+00, 8.44536815D+00, GaH  
     5  8.66476578D+00, 8.85315381D+00, 9.04236706D+00, 9.29093971D+00, GaH  
     6  9.55704513D+00, 9.74103185D+00, 9.89373789D+00, 1.01423079D+01, GaH  
     7  1.02547326D+01, 1.03333179D+01, 1.03494579D+01, 1.03539596D+01, GaH  
     8  1.03513040D+01, 1.03492007D+01, 1.03557690D+01, 1.03630905D+01, GaH  
     9  1.03697850D+01,     57*0.0D+00/                                 GaH  
      DATA TK_GeH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749599978008, GeH  
     2  0.827199986897, 0.915399982642, 1.025300105412, 1.166599990986, GeH  
     3  1.317700115550, 1.455599930590, 1.596299993780, 1.706199960195, GeH  
     4  1.810299950480, 2.036399882037, 2.166399985478, 2.293900043799, GeH  
     5  2.381300154517, 2.470100244507, 2.597800326893, 2.709800051095, GeH  
     6  2.955399934944, 3.135400244240, 3.316599939752, 3.417600029773, GeH  
     7  3.516800382009, 3.706399953424, 3.767400403084, 3.828899900135, GeH  
     8  3.881800161907, 3.928100281119, 3.972200112765, 3.988399740893, GeH  
     9  4.000000000000,     57*0.0D+00/                                 GeH  
      DATA  K_GeH/                                                      071215
     1  8.74044299D-05, 1.17334180D-01, 2.33802187D-01, 5.96536912D-01, GeH  
     2  1.42846077D+00, 2.24530213D+00, 3.10411118D+00, 4.00013701D+00, GeH  
     3  4.75893166D+00, 5.31760367D+00, 5.79094381D+00, 6.10969582D+00, GeH  
     4  6.38003829D+00, 6.89143890D+00, 7.15607555D+00, 7.41083318D+00, GeH  
     5  7.59012083D+00, 7.77995454D+00, 8.06777613D+00, 8.33076557D+00, GeH  
     6  8.90032931D+00, 9.26763087D+00, 9.57293895D+00, 9.71574488D+00, GeH  
     7  9.83822015D+00, 1.00185047D+01, 1.00591257D+01, 1.00919386D+01, GeH  
     8  1.01155442D+01, 1.01353112D+01, 1.01567521D+01, 1.01661475D+01, GeH  
     9  1.01736384D+01,     57*0.0D+00/                                 GeH  
      DATA TK_AsH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749699980301, AsH  
     2  0.827599977216, 0.918600054281, 1.023500060960, 1.162099889084, AsH  
     3  1.301700134623, 1.455299924070, 1.610400130118, 1.726999951306, AsH  
     4  1.847499967412, 1.964299943499, 2.086800073489, 2.344100256433, AsH  
     5  2.609299695084, 2.807200356154, 2.973999978987, 3.135400244240, AsH  
     6  3.267100390917, 3.385400241220, 3.652999691765, 3.785299867245, AsH  
     7  3.893300113369, 4.000000000000,     64*0.0D+00/                 AsH  
      DATA  K_AsH/                                                      071215
     1  4.87291527D-05, 1.40872897D-01, 2.81970806D-01, 7.15066179D-01, AsH  
     2  1.70024710D+00, 2.67160751D+00, 3.59556459D+00, 4.56793685D+00, AsH  
     3  5.33526669D+00, 6.00360786D+00, 6.54601777D+00, 6.89075908D+00, AsH  
     4  7.20507304D+00, 7.47871652D+00, 7.74081687D+00, 8.23401364D+00, AsH  
     5  8.69169440D+00, 9.01063757D+00, 9.26081955D+00, 9.47932426D+00, AsH  
     6  9.63751560D+00, 9.76521997D+00, 1.00339763D+01, 1.01782492D+01, AsH  
     7  1.03044023D+01, 1.04320537D+01,     64*0.0D+00/                 AsH  
      DATA TK_SeH/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719899799358, 0.751900036424, SeH  
     2  0.832899985612, 0.928299885495, 1.041399985358, 1.183699912951, SeH  
     3  1.341000139344, 1.493099847441, 1.647800046098, 1.766000076771, SeH  
     4  1.889799868472, 2.011700269636, 2.136600242708, 2.378600093167, SeH  
     5  2.611299673906, 2.779299777211, 2.951499845093, 3.089100121686, SeH  
     6  3.238199730866, 3.399799642157, 3.609099722127, 3.733299626328, SeH  
     7  3.837800096622, 3.965300165513, 4.000000000000,     63*0.0D+00/ SeH  
      DATA  K_SeH/                                                      071215
     1 -9.98068367D-05, 1.34973937D-01, 2.77204960D-01, 7.01900258D-01, SeH  
     2  1.67092403D+00, 2.64475632D+00, 3.60772027D+00, 4.58505391D+00, SeH  
     3  5.43262266D+00, 6.08063022D+00, 6.61525203D+00, 6.96187999D+00, SeH  
     4  7.28224992D+00, 7.56524581D+00, 7.82996469D+00, 8.29230646D+00, SeH  
     5  8.69611918D+00, 8.96916303D+00, 9.22929487D+00, 9.41852502D+00, SeH  
     6  9.60318628D+00, 9.78130810D+00, 9.98762453D+00, 1.01027699D+01, SeH  
     7  1.01970932D+01, 1.03095767D+01, 1.03401450D+01,     63*0.0D+00/ SeH  
      DATA TK_HBr/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.750099989774, HBr  
     2  0.828599953015, 0.918700056520, 1.025700115290, 1.164399941167, HBr  
     3  1.305200060142, 1.459500015353, 1.615999990092, 1.735099898251, HBr  
     4  1.858300097489, 1.977399951256, 2.099900388130, 2.326099842718, HBr  
     5  2.598300337609, 2.814300108059, 2.955999948767, 3.099100358432, HBr  
     6  3.289299957774, 3.475099886169, 3.642600144520, 3.799600179321, HBr  
     7  3.916500011231, 3.967400212576, 4.000000000000,     63*0.0D+00/ HBr  
      DATA  K_HBr/                                                      071215
     1 -1.10869222D-04, 1.51856512D-01, 3.05702873D-01, 7.77536436D-01, HBr  
     2  1.84812323D+00, 2.88544493D+00, 3.89903294D+00, 4.93904949D+00, HBr  
     3  5.76084545D+00, 6.46841255D+00, 7.04056139D+00, 7.40563861D+00, HBr  
     4  7.73630987D+00, 8.02157484D+00, 8.28832088D+00, 8.73086196D+00, HBr  
     5  9.20907670D+00, 9.56252094D+00, 9.78243319D+00, 9.99149945D+00, HBr  
     6  1.02447948D+01, 1.04624620D+01, 1.06298427D+01, 1.07554906D+01, HBr  
     7  1.08235178D+01, 1.08456201D+01, 1.08576083D+01,     63*0.0D+00/ HBr  
      DATA TK_RbH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749899984889, RbH  
     2  0.827899969956, 0.920600071156, 1.026000122698, 1.098399951243, RbH  
     3  1.171900023887, 1.312800015288, 1.461599986860, 1.612800070107, RbH  
     4  1.743399850719, 1.881000104582, 2.017500406148, 2.172500119229, RbH  
     5  2.297700138432, 2.415399981372, 2.700099810193, 2.894300030047, RbH  
     6  3.107299852464, 3.313100188672, 3.430100286132, 3.488599725638, RbH  
     7  3.545500088908, 3.594200252540, 3.644400017363, 3.689099663459, RbH  
     8  3.727299745247, 3.798000145777, 3.832199974647, 3.864899776432, RbH  
     9  3.954899935374, 3.983299626393, 4.000000000000,     55*0.0D+00/ RbH  
      DATA  K_RbH/                                                      071215
     1  5.53530731D-05, 1.36255387D-01, 2.72677634D-01, 6.93844961D-01, RbH  
     2  1.64618529D+00, 2.60102508D+00, 3.49771241D+00, 4.01912162D+00, RbH  
     3  4.48362117D+00, 5.22713877D+00, 5.85307381D+00, 6.36955778D+00, RbH  
     4  6.74604811D+00, 7.09214501D+00, 7.39784941D+00, 7.71214797D+00, RbH  
     5  7.94700420D+00, 8.15481737D+00, 8.60106886D+00, 8.84994741D+00, RbH  
     6  9.06806097D+00, 9.22709495D+00, 9.29451126D+00, 9.32218925D+00, RbH  
     7  9.34662041D+00, 9.36768428D+00, 9.39319508D+00, 9.42391347D+00, RbH  
     8  9.46085616D+00, 9.57112026D+00, 9.64955965D+00, 9.74075054D+00, RbH  
     9  1.00565110D+01, 1.01681545D+01, 1.02349619D+01,     55*0.0D+00/ RbH  
      DATA TK_SrH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749699980301, SrH  
     2  0.827499979636, 0.920000085623, 1.025100100473, 1.171200040131, SrH  
     3  1.312500009149, 1.461899979479, 1.612700072607, 1.733799929512, SrH  
     4  1.863000060747, 2.001600036740, 2.150499634051, 2.274000169656, SrH  
     5  2.404399713795, 2.524200144847, 2.649899623075, 2.773300227362, SrH  
     6  2.913099906969, 3.072099747108, 3.264700333288, 3.373599973130, SrH  
     7  3.484199628057, 3.543200033895, 3.600100373462, 3.710400044016, SrH  
     8  3.797100126908, 3.870899913321, 3.912699913596, 3.948599796769, SrH  
     9  3.980599565776, 4.000000000000,     56*0.0D+00/                 SrH  
      DATA  K_SrH/                                                      071215
     1 -1.62149101D-05, 1.25245880D-01, 2.50729694D-01, 6.35818530D-01, SrH  
     2  1.51118796D+00, 2.39120528D+00, 3.21970970D+00, 4.13848561D+00, SrH  
     3  4.83613770D+00, 5.42780177D+00, 5.91654472D+00, 6.25203544D+00, SrH  
     4  6.56887098D+00, 6.87371017D+00, 7.17177515D+00, 7.40205680D+00, SrH  
     5  7.63250156D+00, 7.83402429D+00, 8.03389649D+00, 8.21588191D+00, SrH  
     6  8.40167786D+00, 8.58445266D+00, 8.76579796D+00, 8.85009629D+00, SrH  
     7  8.92472069D+00, 8.96279497D+00, 9.00192689D+00, 9.10205796D+00, SrH  
     8  9.22561021D+00, 9.37741490D+00, 9.48500204D+00, 9.58925364D+00, SrH  
     9  9.69007804D+00, 9.75414107D+00,     56*0.0D+00/                 SrH  
      DATA TK_AgH/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.751000013099, AgH  
     2  0.830699935179, 0.925199960242, 1.031900170714, 1.105100027887, AgH  
     3  1.179399849837, 1.325500037379, 1.484799903037, 1.645900001285, AgH  
     4  1.759000213954, 1.875200012665, 2.105100027135, 2.334600036208, AgH  
     5  2.558000366349, 2.739199757091, 2.998399962554, 3.179100269772, AgH  
     6  3.363699732203, 3.467300193348, 3.578499889611, 3.662299902050, AgH  
     7  3.724599932340, 3.796600116425, 3.852300207349, 3.903499698251, AgH  
     8  3.963400122932, 3.985799682521, 4.000000000000,     59*0.0D+00/ AgH  
      DATA  K_AgH/                                                      071215
     1 -7.17040280D-05, 1.45110744D-01, 2.94914312D-01, 7.47074189D-01, AgH  
     2  1.76876843D+00, 2.78484008D+00, 3.72658657D+00, 4.27105837D+00, AgH  
     3  4.75450167D+00, 5.54272434D+00, 6.21828729D+00, 6.76514223D+00, AgH  
     4  7.09053321D+00, 7.38768096D+00, 7.89695337D+00, 8.33759710D+00, AgH  
     5  8.72734561D+00, 9.02170400D+00, 9.39673659D+00, 9.61396330D+00, AgH  
     6  9.79405967D+00, 9.87613653D+00, 9.94741543D+00, 9.98806326D+00, AgH  
     7  1.00108651D+01, 1.00311126D+01, 1.00451427D+01, 1.00603379D+01, AgH  
     8  1.00872808D+01, 1.01018278D+01, 1.01128089D+01,     59*0.0D+00/ AgH  
      DATA TK_CdH/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751200018282, CdH  
     2  0.831299948934, 0.926799921662, 1.034200109262, 1.108300101655, CdH  
     3  1.183099900459, 1.328899959963, 1.491399811784, 1.654799970266, CdH  
     4  1.767900025351, 1.885699978478, 2.005200119404, 2.128900067940, CdH  
     5  2.355000017404, 2.554200278593, 2.756000134908, 2.925100186092, CdH  
     6  3.118199832715, 3.316399953976, 3.483999623622, 3.662299902050, CdH  
     7  3.732299605669, 3.804100277582, 3.864499767293, 3.921300130042, CdH  
     8  3.969500259638, 3.987599722932, 4.000000000000,     59*0.0D+00/ CdH  
      DATA  K_CdH/                                                      071215
     1  4.61233026D-06, 1.30520849D-01, 2.66529629D-01, 6.74265867D-01, CdH  
     2  1.59765608D+00, 2.52223101D+00, 3.37810807D+00, 3.87729789D+00, CdH  
     3  4.31931764D+00, 5.03815795D+00, 5.67299983D+00, 6.18820594D+00, CdH  
     4  6.49335769D+00, 6.77820855D+00, 7.04087421D+00, 7.29159774D+00, CdH  
     5  7.71057375D+00, 8.04990825D+00, 8.36572703D+00, 8.60112225D+00, CdH  
     6  8.83121270D+00, 9.02574418D+00, 9.16137934D+00, 9.27835779D+00, CdH  
     7  9.31581466D+00, 9.34914075D+00, 9.37421590D+00, 9.39786775D+00, CdH  
     8  9.42177174D+00, 9.43276961D+00, 9.44126580D+00,     59*0.0D+00/ CdH  
      DATA TK_InH/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.751000013099, InH  
     2  0.830699935179, 0.925399955419, 1.032200162698, 1.105700041718, InH  
     3  1.179999835913, 1.324400062426, 1.484999898398, 1.646100006002, InH  
     4  1.758100191277, 1.874199987926, 1.992599831832, 2.116899818222, InH  
     5  2.362999718349, 2.563100191902, 2.747399931856, 2.944599691672, InH  
     6  3.034299846847, 3.126800040721, 3.379100086005, 3.475399864664, InH  
     7  3.573999782104, 3.657299792099, 3.702199857799, 3.747199924899, InH  
     8  3.798700160452, 3.849200351509, 3.926300241128, 3.971500163063, InH  
     9  3.988199736403, 4.000000000000,     56*0.0D+00/                 InH  
      DATA  K_InH/                                                      071215
     1  7.49535475D-05, 1.42833550D-01, 2.90107302D-01, 7.34482052D-01, InH  
     2  1.73822099D+00, 2.73892299D+00, 3.66642950D+00, 4.20487089D+00, InH  
     3  4.68130662D+00, 5.45099694D+00, 6.12516629D+00, 6.66680753D+00, InH  
     4  6.98659573D+00, 7.28172427D+00, 7.55345889D+00, 7.81458485D+00, InH  
     5  8.28083073D+00, 8.62609148D+00, 8.92218297D+00, 9.22043947D+00, InH  
     6  9.35307291D+00, 9.48832266D+00, 9.83124203D+00, 9.94097667D+00, InH  
     7  1.00307813D+01, 1.00786494D+01, 1.00901550D+01, 1.00906767D+01, InH  
     8  1.00789588D+01, 1.00587134D+01, 1.00275281D+01, 1.00206897D+01, InH  
     9  1.00221100D+01, 1.00246948D+01,     56*0.0D+00/                 InH  
      DATA TK_SnH/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718799827163, 0.749399973420, SnH  
     2  0.826699998997, 0.917300025178, 1.023300056021, 1.162899907200, SnH  
     3  1.297100105804, 1.440000042582, 1.584000017339, 1.696199912097, SnH  
     4  1.813600032814, 2.045000083863, 2.251500037099, 2.430500248465, SnH  
     5  2.557300350183, 2.680600264902, 2.859499679617, 2.944199682867, SnH  
     6  3.027699907896, 3.264700333288, 3.386900279204, 3.515300346979, SnH  
     7  3.753400067761, 3.852300207349, 3.934599991792, 3.974799925945, SnH  
     8  3.988999754363, 4.000000000000,     60*0.0D+00/                 SnH  
      DATA  K_SnH/                                                      071215
     1  8.69336027D-05, 1.15981259D-01, 2.32340642D-01, 5.93433521D-01, SnH  
     2  1.42019487D+00, 2.25316900D+00, 3.07261987D+00, 3.94786942D+00, SnH  
     3  4.62097474D+00, 5.20317397D+00, 5.68682861D+00, 6.01055371D+00, SnH  
     4  6.31174018D+00, 6.82549212D+00, 7.22480240D+00, 7.54278858D+00, SnH  
     5  7.75769118D+00, 7.96276603D+00, 8.26701246D+00, 8.41943108D+00, SnH  
     6  8.57628650D+00, 9.03775046D+00, 9.26547071D+00, 9.48603748D+00, SnH  
     7  9.83572242D+00, 9.95808760D+00, 1.00536911D+01, 1.01013769D+01, SnH  
     8  1.01189056D+01, 1.01328573D+01,     60*0.0D+00/                 SnH  
      DATA TK_SbH/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751200018282, SbH  
     2  0.831199946641, 0.926399931307, 1.033600125293, 1.108600108571, SbH  
     3  1.184699933773, 1.332099983068, 1.488399819532, 1.646300010719, SbH  
     4  1.776200124184, 1.915299960756, 2.039299952416, 2.167900015187, SbH  
     5  2.311200342459, 2.482799597895, 2.663899947089, 2.900999628790, SbH  
     6  3.156599753778, 3.263100294869, 3.364799761369, 3.560800347740, SbH  
     7  3.635100227176, 3.708800008067, 3.864999778717, 3.926300241128, SbH  
     8  4.000000000000,     61*0.0D+00/                                 SbH  
      DATA  K_SbH/                                                      071215
     1  1.15569325D-04, 1.40493319D-01, 2.86762571D-01, 7.25142243D-01, SbH  
     2  1.71589629D+00, 2.70452129D+00, 3.61873668D+00, 4.15754389D+00, SbH  
     3  4.63497130D+00, 5.40166152D+00, 6.04394732D+00, 6.56765331D+00, SbH  
     4  6.93156407D+00, 7.27303642D+00, 7.54628579D+00, 7.80619330D+00, SbH  
     5  8.07256240D+00, 8.36050077D+00, 8.62615772D+00, 8.91671722D+00, SbH  
     6  9.16716800D+00, 9.25634191D+00, 9.33640824D+00, 9.49261808D+00, SbH  
     7  9.55817922D+00, 9.62806404D+00, 9.78585869D+00, 9.84739518D+00, SbH  
     8  9.91991119D+00,     61*0.0D+00/                                 SbH  
      DATA TK_TeH/                                                      071215
     1  0.699999789529, 0.709700041799, 0.719999796830, 0.752000039016, TeH  
     2  0.833399997074, 0.929399858972, 1.040599967607, 1.193499958786, TeH  
     3  1.330399944088, 1.499099973292, 1.652100042108, 1.764500117366, TeH  
     4  1.883700032139, 2.001800041332, 2.123899963071, 2.364299749845, TeH  
     5  2.584400028533, 2.865599782767, 3.002000044630, 3.145299954174, TeH  
     6  3.325999828039, 3.479899542081, 3.622399952994, 3.746099900431, TeH  
     7  3.884700230742, 3.956999981229, 3.983599633128, 4.000000000000, TeH  
     8      62*0.0D+00/                                                 TeH  
      DATA  K_TeH/                                                      071215
     1 -1.02752124D-05, 1.36338779D-01, 2.78389596D-01, 7.02376501D-01, TeH  
     2  1.67263949D+00, 2.64423782D+00, 3.57770446D+00, 4.59731743D+00, TeH  
     3  5.31408673D+00, 6.01425046D+00, 6.52466956D+00, 6.84531965D+00, TeH  
     4  7.14785475D+00, 7.41839191D+00, 7.67484371D+00, 8.13146124D+00, TeH  
     5  8.51162031D+00, 8.94755463D+00, 9.13084811D+00, 9.29861135D+00, TeH  
     6  9.47426238D+00, 9.60094049D+00, 9.71088389D+00, 9.80677762D+00, TeH  
     7  9.91689710D+00, 9.97593939D+00, 9.99828319D+00, 1.00123479D+01, TeH  
     8      62*0.0D+00/                                                 TeH  
      DATA TK_HI/                                                       071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750199992366, HI   
     2  0.828799948174, 0.921800042222, 1.027200152332, 1.098499948685, HI   
     3  1.170900047093, 1.313000019380, 1.468799809738, 1.626500033782, HI   
     4  1.742899840381, 1.862600070845, 1.975599913270, 2.093400224668, HI   
     5  2.332899995910, 2.583600009790, 2.707599996457, 2.820599715501, HI   
     6  3.075599834760, 3.234899650718, 3.385300238688, 3.490199761354, HI   
     7  3.603000163587, 3.787499919113, 3.909299828503, 3.963200118450, HI   
     8  4.000000000000,     61*0.0D+00/                                 HI   
      DATA  K_HI/                                                       071215
     1 -5.85057696D-05, 1.50062739D-01, 3.03520589D-01, 7.68873550D-01, HI   
     2  1.82367507D+00, 2.87207794D+00, 3.84858227D+00, 4.40596290D+00, HI   
     3  4.90190360D+00, 5.71044201D+00, 6.40749747D+00, 6.97084322D+00, HI   
     4  7.32087467D+00, 7.63787360D+00, 7.90656046D+00, 8.16250893D+00, HI   
     5  8.62939880D+00, 9.06914050D+00, 9.27446661D+00, 9.45497462D+00, HI   
     6  9.82900907D+00, 1.00299502D+01, 1.01941034D+01, 1.02952782D+01, HI   
     7  1.03931156D+01, 1.05256327D+01, 1.05876957D+01, 1.06082957D+01, HI   
     8  1.06206913D+01,     61*0.0D+00/                                 HI   
      DATA TK_CsH/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751200018282, CsH  
     2  0.831399951226, 0.926799921662, 1.034700095902, 1.110200135654, CsH  
     3  1.187099983744, 1.335400058736, 1.492199828564, 1.649700090910, CsH  
     4  1.764600114660, 1.886299962379, 2.010700246100, 2.152799688445, CsH  
     5  2.271800325734, 2.388900317865, 2.613099715341, 2.729199616013, CsH  
     6  2.852000236090, 3.134500223614, 3.258400185511, 3.378400071639, CsH  
     7  3.474999893338, 3.632600175638, 3.710600048330, 3.805600310463, CsH  
     8  3.932700128733, 3.972600084024, 4.000000000000,     59*0.0D+00/ CsH  
      DATA  K_CsH/                                                      071215
     1 -2.54712938D-05, 1.36628867D-01, 2.79025090D-01, 7.05879737D-01, CsH  
     2  1.67373229D+00, 2.64074113D+00, 3.54024495D+00, 4.07079179D+00, CsH  
     3  4.54276085D+00, 5.29771026D+00, 5.92935175D+00, 6.44290692D+00, CsH  
     4  6.76271444D+00, 7.06424581D+00, 7.34246945D+00, 7.63263850D+00, CsH  
     5  7.85837153D+00, 8.06762075D+00, 8.43131395D+00, 8.59743301D+00, CsH  
     6  8.75491704D+00, 9.04838917D+00, 9.15008247D+00, 9.23562784D+00, CsH  
     7  9.29876981D+00, 9.40971255D+00, 9.47709739D+00, 9.56903910D+00, CsH  
     8  9.69786816D+00, 9.74001162D+00, 9.76965756D+00,     59*0.0D+00/ CsH  
      DATA TK_BaH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749799982595, BaH  
     2  0.827699974796, 0.920600071156, 1.026000122698, 1.098799941011, BaH  
     3  1.172700005321, 1.314600052119, 1.464499915519, 1.615899992592, BaH  
     4  1.738099826109, 1.868099931990, 2.005800133182, 2.155599754664, BaH  
     5  2.277099949729, 2.405299735752, 2.527399907215, 2.644300028798, BaH  
     6  2.848000341093, 3.064200126989, 3.218000210292, 3.310300387807, BaH  
     7  3.394999992630, 3.450299800851, 3.506400139279, 3.597600327668, BaH  
     8  3.683300067416, 3.760800243521, 3.861199691891, 3.905199736428, BaH  
     9  3.942999670329, 3.978199681642, 3.989899774569, 4.000000000000, BaH  
     A      54*0.0D+00/                                                 BaH  
      DATA  K_BaH/                                                      071215
     1  2.55093028D-07, 1.30952400D-01, 2.62125974D-01, 6.65874935D-01, BaH  
     2  1.58119327D+00, 2.50278311D+00, 3.36731223D+00, 3.87332252D+00, BaH  
     3  4.32433769D+00, 5.04822061D+00, 5.65898470D+00, 6.16157162D+00, BaH  
     4  6.50680843D+00, 6.83066758D+00, 7.13747183D+00, 7.44043056D+00, BaH  
     5  7.66884640D+00, 7.89695460D+00, 8.10315911D+00, 8.28947816D+00, BaH  
     6  8.58151090D+00, 8.83824686D+00, 8.98866508D+00, 9.07093901D+00, BaH  
     7  9.14822612D+00, 9.20332415D+00, 9.26457272D+00, 9.37301320D+00, BaH  
     8  9.47600755D+00, 9.56766827D+00, 9.70271316D+00, 9.77563088D+00, BaH  
     9  9.84680434D+00, 9.92000302D+00, 9.94568191D+00, 9.96834326D+00, BaH  
     A      54*0.0D+00/                                                 BaH  
      DATA TK_YbH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749899984889, YbH  
     2  0.828099965115, 0.921000061511, 1.026300130107, 1.099499923107, YbH  
     3  1.173799979794, 1.316300086904, 1.466599863858, 1.618299932581, YbH  
     4  1.742199825907, 1.874499995348, 2.014700340246, 2.164399945867, YbH  
     5  2.407799796745, 2.632100153154, 2.787299909646, 2.940299597022, YbH  
     6  3.170500085418, 3.311200323799, 3.440999606249, 3.500700006393, YbH  
     7  3.560100398049, 3.609999656994, 3.660599867347, 3.709700028558, YbH  
     8  3.756800148340, 3.821399718815, 3.885300244984, 3.961400078111, YbH  
     9  3.985099666805, 4.000000000000,     56*0.0D+00/                 YbH  
      DATA  K_YbH/                                                      071215
     1  1.39433316D-05, 1.26180881D-01, 2.52566430D-01, 6.42837010D-01, YbH  
     2  1.52803571D+00, 2.41620978D+00, 3.24964166D+00, 3.74112677D+00, YbH  
     3  4.17941826D+00, 4.88300905D+00, 5.47733219D+00, 5.96779283D+00, YbH  
     4  6.30966991D+00, 6.63216341D+00, 6.93847381D+00, 7.23635485D+00, YbH  
     5  7.67698997D+00, 8.04618873D+00, 8.27675703D+00, 8.47809516D+00, YbH  
     6  8.72741240D+00, 8.84830312D+00, 8.93787488D+00, 8.97118224D+00, YbH  
     7  8.99943727D+00, 9.02070683D+00, 9.04320239D+00, 9.07159234D+00, YbH  
     8  9.11261420D+00, 9.20610185D+00, 9.35269433D+00, 9.59011551D+00, YbH  
     9  9.67349866D+00, 9.72734693D+00,     56*0.0D+00/                 YbH  
      DATA TK_PtH/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751500026057, PtH  
     2  0.831899962688, 0.928999868616, 1.037600018419, 1.114000037016, PtH  
     3  1.189700037879, 1.333300010584, 1.486699858965, 1.649800093269, PtH  
     4  1.794400063133, 1.937599987897, 2.201399826967, 2.307400368061, PtH  
     5  2.408399811383, 2.562600227488, 2.662499916149, 2.760200238223, PtH  
     6  2.880900127211, 3.019100431706, 3.196299711496, 3.370499909510, PtH  
     7  3.479499570755, 3.585600056648, 3.683900025627, 3.785299867245, PtH  
     8  3.886500273467, 3.961300075870, 4.000000000000,     59*0.0D+00/ PtH  
      DATA  K_PtH/                                                      071215
     1 -5.91345970D-05, 1.33487332D-01, 2.71417547D-01, 6.89029506D-01, PtH  
     2  1.64214204D+00, 2.62551268D+00, 3.54959188D+00, 4.10887527D+00, PtH  
     3  4.60100736D+00, 5.39337288D+00, 6.07674612D+00, 6.66500343D+00, PtH  
     4  7.09898216D+00, 7.46970858D+00, 8.04925940D+00, 8.25823828D+00, PtH  
     5  8.45199256D+00, 8.74962485D+00, 8.94843288D+00, 9.14676367D+00, PtH  
     6  9.39105274D+00, 9.65790332D+00, 9.96249434D+00, 1.02132053D+01, PtH  
     7  1.03474425D+01, 1.04637083D+01, 1.05584867D+01, 1.06388675D+01, PtH  
     8  1.06959150D+01, 1.07252119D+01, 1.07388075D+01,     59*0.0D+00/ PtH  
      DATA TK_AuH/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749899984889, AuH  
     2  0.828099965115, 0.919600076668, 1.024800093064, 1.164099934374, AuH  
     3  1.304200081422, 1.458800000139, 1.615699997593, 1.735499888632, AuH  
     4  1.859300120429, 1.978699978690, 2.101500283730, 2.341200191204, AuH  
     5  2.600100366740, 2.772400294885, 3.013500303363, 3.109799671761, AuH  
     6  3.215000139798, 3.348300350903, 3.536899886809, 3.670200063991, AuH  
     7  3.754700098571, 3.841600180535, 3.910299851931, 3.971200184619, AuH  
     8  4.000000000000,     61*0.0D+00/                                 AuH  
      DATA  K_AuH/                                                      071215
     1 -1.19255964D-04, 1.44527278D-01, 2.89442955D-01, 7.36985029D-01, AuH  
     2  1.75135761D+00, 2.75165732D+00, 3.69918224D+00, 4.69660748D+00, AuH  
     3  5.48063227D+00, 6.16357196D+00, 6.71856044D+00, 7.07525468D+00, AuH  
     4  7.39910430D+00, 7.67883302D+00, 7.94130501D+00, 8.40226344D+00, AuH  
     5  8.85187254D+00, 9.13267318D+00, 9.49570643D+00, 9.62634206D+00, AuH  
     6  9.75760202D+00, 9.90764543D+00, 1.01024341D+01, 1.02391822D+01, AuH  
     7  1.03239147D+01, 1.04045211D+01, 1.04615481D+01, 1.05080941D+01, AuH  
     8  1.05297057D+01,     61*0.0D+00/                                 AuH  
      DATA TK_HgH/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720999820779, 0.754400101215, HgH  
     2  0.839300132325, 0.936899978931, 1.057700003907, 1.199499812491, HgH  
     3  1.348899942883, 1.499399979584, 1.659299850529, 1.829599924186, HgH  
     4  2.009200211254, 2.180400277896, 2.353600121490, 2.462600084861, HgH  
     5  2.561900277309, 2.810700365138, 2.954499914209, 3.094700252413, HgH  
     6  3.324099786843, 3.555500311323, 3.644300024427, 3.739999764743, HgH  
     7  3.824099784091, 3.895699936515, 3.961400078111, 3.985099666805, HgH  
     8  4.000000000000,     61*0.0D+00/                                 HgH  
      DATA  K_HgH/                                                      071215
     1 -3.71443565D-06, 1.36145959D-01, 2.85529422D-01, 7.15069383D-01, HgH  
     2  1.68424815D+00, 2.61381123D+00, 3.54827736D+00, 4.41444583D+00, HgH  
     3  5.13076125D+00, 5.70740046D+00, 6.20806281D+00, 6.65345070D+00, HgH  
     4  7.05666042D+00, 7.39833354D+00, 7.71535993D+00, 7.90369657D+00, HgH  
     5  8.06817988D+00, 8.44375662D+00, 8.63112260D+00, 8.79101194D+00, HgH  
     6  9.00750884D+00, 9.17556396D+00, 9.22709630D+00, 9.27541863D+00, HgH  
     7  9.31368775D+00, 9.34567732D+00, 9.37795581D+00, 9.39147099D+00, HgH  
     8  9.40086701D+00,     61*0.0D+00/                                 HgH  
      DATA TK_TlH/                                                      071215
     1  0.699999789529, 0.710000049601, 0.721299827964, 0.755000116765, TlH  
     2  0.840800129070, 0.940200034416, 1.062100001037, 1.203699879419, TlH  
     3  1.355400035353, 1.514700104693, 1.676599915836, 1.790399980903, TlH  
     4  1.904999956506, 2.124599977753, 2.373699976236, 2.565799999736, TlH  
     5  2.799500193189, 2.996199911067, 3.174400169021, 3.336800084391, TlH  
     6  3.405999748897, 3.482399588138, 3.583099998075, 3.654699731432, TlH  
     7  3.719600242463, 3.770300444233, 3.820999709145, 3.868499858688, TlH  
     8  3.911699887902, 3.966000181200, 3.986599700481, 4.000000000000, TlH  
     9      58*0.0D+00/                                                 TlH  
      DATA  K_TlH/                                                      071215
     1 -3.14526753D-05, 1.49805537D-01, 3.15391595D-01, 7.86807346D-01, TlH  
     2  1.84955752D+00, 2.87217524D+00, 3.88529905D+00, 4.81004513D+00, TlH  
     3  5.58212029D+00, 6.22082865D+00, 6.74430865D+00, 7.05813435D+00, TlH  
     4  7.34078331D+00, 7.81553720D+00, 8.28499467D+00, 8.61460967D+00, TlH  
     5  8.97756520D+00, 9.23872400D+00, 9.43474892D+00, 9.58334457D+00, TlH  
     6  9.64028405D+00, 9.69941815D+00, 9.76695772D+00, 9.79853416D+00, TlH  
     7  9.80561593D+00, 9.79331262D+00, 9.76657878D+00, 9.73340061D+00, TlH  
     8  9.70285314D+00, 9.67493186D+00, 9.67000909D+00, 9.66901722D+00, TlH  
     9      58*0.0D+00/                                                 TlH  
      DATA TK_PbH/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751400023466, PbH  
     2  0.831999964980, 0.926399931307, 1.036500047809, 1.185299946265, PbH  
     3  1.331099960139, 1.483299937831, 1.636199956763, 1.752300045135, PbH  
     4  1.875300015139, 2.004600105627, 2.138500286378, 2.394500010731, PbH  
     5  2.611799685415, 2.846400306380, 3.057700380461, 3.151099628513, PbH  
     6  3.243699858792, 3.426400221107, 3.492699819701, 3.559000384366, PbH  
     7  3.752600048802, 3.854900024040, 3.935199948548, 3.974899918760, PbH  
     8  3.989099756609, 4.000000000000,     60*0.0D+00/                 PbH  
      DATA  K_PbH/                                                      071215
     1 -4.85919998D-05, 1.19480780D-01, 2.42947152D-01, 6.15671742D-01, PbH  
     2  1.47064918D+00, 2.32496439D+00, 3.15596170D+00, 4.05637600D+00, PbH  
     3  4.75183227D+00, 5.33538571D+00, 5.81683121D+00, 6.13140157D+00, PbH  
     4  6.42902873D+00, 6.71199159D+00, 6.98091421D+00, 7.44852360D+00, PbH  
     5  7.81278254D+00, 8.17083292D+00, 8.44661025D+00, 8.55197512D+00, PbH  
     6  8.64742615D+00, 8.82390919D+00, 8.89124133D+00, 8.96320645D+00, PbH  
     7  9.19051965D+00, 9.30997768D+00, 9.40691616D+00, 9.45951959D+00, PbH  
     8  9.47961499D+00, 9.49561275D+00,     60*0.0D+00/                 PbH  
      DATA TK_BiH/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720999820779, 0.754400101215, BiH  
     2  0.839200130033, 0.936899978931, 1.057700003907, 1.199299817367, BiH  
     3  1.348299957804, 1.496999929244, 1.658299877137, 1.843000075988, BiH  
     4  2.034099826220, 2.234699661302, 2.437299773874, 2.611399676208, BiH  
     5  2.765100356584, 3.206299939391, 3.365799787884, 3.536799884533, BiH  
     6  3.681500192782, 3.763800316049, 3.839700138006, 3.935899898096, BiH  
     7  3.974999911574, 3.989099756609, 4.000000000000,     63*0.0D+00/ BiH  
      DATA  K_BiH/                                                      071215
     1 -9.22134132D-05, 1.45015576D-01, 3.04196046D-01, 7.61721860D-01, BiH  
     2  1.79190989D+00, 2.78064580D+00, 3.77129520D+00, 4.68528368D+00, BiH  
     3  5.43703321D+00, 6.03412319D+00, 6.56016546D+00, 7.05592948D+00, BiH  
     4  7.48935655D+00, 7.88831157D+00, 8.25409211D+00, 8.54768762D+00, BiH  
     5  8.79088653D+00, 9.35413920D+00, 9.48844573D+00, 9.58478028D+00, BiH  
     6  9.63470248D+00, 9.65616634D+00, 9.67677913D+00, 9.71279419D+00, BiH  
     7  9.73332925D+00, 9.74186981D+00, 9.74894165D+00,     63*0.0D+00/ BiH  
      DATA TK_HeHp/                                                     071215
     1  0.699999789529, 0.709700041799, 0.720199801620, 0.752500051974, HeHp 
     2  0.834500022290, 0.933399910744, 1.044400051921, 1.176299921777, HeHp 
     3  1.292199995984, 1.463499940119, 1.675799934248, 1.763600141723, HeHp 
     4  1.861900088518, 2.080899952399, 2.329299912010, 2.593100226169, HeHp 
     5  2.811500308009, 3.026300006885, 3.250400010905, 3.398199758981, HeHp 
     6  3.504500094984, 3.602500199773, 3.823499769585, 3.930800265673, HeHp 
     7  3.972800069653, 4.000000000000,     64*0.0D+00/                 HeHp 
      DATA  K_HeHp/                                                     071215
     1 -2.98845257D-05, 1.35161393D-01, 2.78747625D-01, 7.03160520D-01, HeHp 
     2  1.67399018D+00, 2.67173668D+00, 3.61091974D+00, 4.53408968D+00, HeHp 
     3  5.20947893D+00, 6.02094254D+00, 6.78427913D+00, 7.04303315D+00, HeHp 
     4  7.30449124D+00, 7.80942210D+00, 8.29811493D+00, 8.76101126D+00, HeHp 
     5  9.11805124D+00, 9.44708718D+00, 9.74781259D+00, 9.90757139D+00, HeHp 
     6  9.99651069D+00, 1.00573612D+01, 1.01527777D+01, 1.02052340D+01, HeHp 
     7  1.02299280D+01, 1.02473512D+01,     64*0.0D+00/                 HeHp 
      DATA TK_BeHp/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751300020874, BeHp 
     2  0.831799960396, 0.922300030166, 1.034000114605, 1.180399844241, BeHp 
     3  1.325500037379, 1.480500002779, 1.637599921898, 1.749299972709, BeHp 
     4  1.863100058222, 2.081999974975, 2.317699875769, 2.566499949916, BeHp 
     5  2.691999629173, 2.807300358260, 3.084800038174, 3.238299733295, BeHp 
     6  3.431100215355, 3.536499877706, 3.723600001634, 3.824099784091, BeHp 
     7  3.918400060049, 3.968400234986, 4.000000000000,     63*0.0D+00/ BeHp 
      DATA  K_BeHp/                                                     071215
     1 -1.06344925D-04, 1.48081705D-01, 3.01012636D-01, 7.60356707D-01, BeHp 
     2  1.81008968D+00, 2.80967775D+00, 3.82491151D+00, 4.86894957D+00, BeHp 
     3  5.66774826D+00, 6.33865326D+00, 6.88405840D+00, 7.21263597D+00, BeHp 
     4  7.50980137D+00, 8.00508208D+00, 8.46497075D+00, 8.90181944D+00, BeHp 
     5  9.10961727D+00, 9.29350008D+00, 9.69550851D+00, 9.88326642D+00, BeHp 
     6  1.00781634D+01, 1.01644955D+01, 1.02752852D+01, 1.03089090D+01, BeHp 
     7  1.03292167D+01, 1.03393538D+01, 1.03466859D+01,     63*0.0D+00/ BeHp 
      DATA TK_CHp/                                                      071215
     1  0.699999789529, 0.710300042018, 0.722699861492, 0.758300202290, CHp  
     2  0.796600101508, 0.849699914331, 0.956299958644, 1.080099922608, CHp  
     3  1.188700017058, 1.285200078922, 1.373299952839, 1.462099974559, CHp  
     4  1.533500142825, 1.606400050550, 1.749399974777, 1.973199862623, CHp  
     5  2.302100244620, 2.509600237139, 2.710600069582, 2.890400319978, CHp  
     6  3.082099985737, 3.222600068638, 3.352800177728, 3.501100015718, CHp  
     7  3.575499817940, 3.649899628828, 3.848500335761, 3.940199607108, CHp  
     8  4.000000000000,     61*0.0D+00/                                 CHp  
      DATA  K_CHp/                                                      071215
     1  9.39339470D-05, 1.59786855D-01, 3.47824867D-01, 8.63425311D-01, CHp  
     2  1.38027172D+00, 2.03754866D+00, 3.17507353D+00, 4.24204048D+00, CHp  
     3  4.99855872D+00, 5.56177998D+00, 6.00781875D+00, 6.40778267D+00, CHp  
     4  6.70139193D+00, 6.98029789D+00, 7.47644807D+00, 8.13964481D+00, CHp  
     5  8.92228363D+00, 9.33692220D+00, 9.70217613D+00, 1.00059772D+01, CHp  
     6  1.03025946D+01, 1.04953437D+01, 1.06477014D+01, 1.07732789D+01, CHp  
     7  1.08099926D+01, 1.08275930D+01, 1.08004805D+01, 1.07684977D+01, CHp  
     8  1.07467624D+01,     61*0.0D+00/                                 CHp  
      DATA TK_NHp/                                                      071215
     1  0.699999789529, 0.710300042018, 0.722699861492, 0.758200199698, NHp  
     2  0.848699938459, 0.964599935689, 1.097199981936, 1.184199923362, NHp  
     3  1.269200000684, 1.364900011716, 1.458299989272, 1.582499979611, NHp  
     4  1.706799974311, 1.817200122632, 1.927699882684, 2.079099911366, NHp  
     5  2.232299603099, 2.477099736567, 2.597500320464, 2.730899577615, NHp  
     6  2.839600158525, 2.942399643247, 3.172700132578, 3.352800177728, NHp  
     7  3.505700122960, 3.641200243419, 3.806700334576, 3.915999998385, NHp  
     8  3.967500214817, 4.000000000000,     60*0.0D+00/                 NHp  
      DATA  K_NHp/                                                      071215
     1 -2.69418706D-05, 1.27587777D-01, 2.78034732D-01, 6.90534843D-01, NHp  
     2  1.63178321D+00, 2.64485126D+00, 3.60158126D+00, 4.14171047D+00, NHp  
     3  4.61922508D+00, 5.11053687D+00, 5.55197431D+00, 6.09062217D+00, NHp  
     4  6.58219070D+00, 6.98276319D+00, 7.35032400D+00, 7.79774525D+00, NHp  
     5  8.18439185D+00, 8.68843423D+00, 8.90177800D+00, 9.12266220D+00, NHp  
     6  9.29419847D+00, 9.45013588D+00, 9.77119612D+00, 9.98629758D+00, NHp  
     7  1.01413801D+01, 1.02584573D+01, 1.03744892D+01, 1.04325590D+01, NHp  
     8  1.04552336D+01, 1.04684370D+01,     60*0.0D+00/                 NHp  
      DATA TK_OHp/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751500026057, OHp  
     2  0.831799960396, 0.931099865935, 1.039199975669, 1.192099992922, OHp  
     3  1.344400054791, 1.488099826491, 1.630500098712, 1.754700105607, OHp  
     4  1.883300042871, 2.003300075776, 2.125599998727, 2.340400173210, OHp  
     5  2.586300073048, 2.744299865891, 2.891800215900, 3.199999798209, OHp  
     6  3.496299903720, 3.770800408050, 3.860699680467, 3.953599906988, OHp  
     7  4.000000000000,     65*0.0D+00/                                 OHp  
      DATA  K_OHp/                                                      071215
     1  1.73752975D-04, 1.49441096D-01, 3.03536700D-01, 7.69607667D-01, OHp  
     2  1.82872518D+00, 2.93789105D+00, 3.93794573D+00, 5.05742332D+00, OHp  
     3  5.90860829D+00, 6.53855217D+00, 7.04610106D+00, 7.42035465D+00, OHp  
     4  7.75853460D+00, 8.04017773D+00, 8.30190582D+00, 8.71776090D+00, OHp  
     5  9.14937901D+00, 9.41085242D+00, 9.64605553D+00, 1.00949163D+01, OHp  
     6  1.04321809D+01, 1.06302472D+01, 1.06645509D+01, 1.06876333D+01, OHp  
     7  1.06975657D+01,     65*0.0D+00/                                 OHp  
      DATA TK_HFp/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720099799225, 0.752400049382, HFp  
     2  0.834200015413, 0.932999902951, 1.043600034171, 1.119299899443, HFp  
     3  1.193799951471, 1.335400058736, 1.505500121493, 1.672700005592, HFp  
     4  1.896700009654, 2.142800125238, 2.251600039427, 2.361299677163, HFp  
     5  2.581299955904, 2.859499679617, 3.162099882346, 3.462200083008, HFp  
     6  3.778699836364, 3.907999799308, 3.964100138620, 4.000000000000, HFp  
     7      66*0.0D+00/                                                 HFp  
      DATA  K_HFp/                                                      071215
     1  9.91785864D-05, 1.27237382D-01, 2.61032684D-01, 6.60571933D-01, HFp  
     2  1.57389085D+00, 2.51668208D+00, 3.40496508D+00, 3.93045446D+00, HFp  
     3  4.39293070D+00, 5.14667973D+00, 5.87712825D+00, 6.45432096D+00, HFp  
     4  7.07574085D+00, 7.63288558D+00, 7.85493528D+00, 8.07033085D+00, HFp  
     5  8.48604909D+00, 8.98486817D+00, 9.47625439D+00, 9.87493091D+00, HFp  
     6  1.01684828D+01, 1.02455220D+01, 1.02719464D+01, 1.02871843D+01, HFp  
     7      66*0.0D+00/                                                 HFp  
      DATA TK_NeHp/                                                     071215
     1  0.699999789529, 0.709700041799, 0.719999796830, 0.752100041607, NeHp 
     2  0.833199992489, 0.935799957501, 1.043700036390, 1.195099919774, NeHp 
     3  1.349899918014, 1.494199870514, 1.641999909301, 1.752000037576, NeHp 
     4  1.862900063272, 2.072799753430, 2.315500033725, 2.572399754724, NeHp 
     5  2.793100048256, 2.993699852558, 3.177800241905, 3.436399840233, NeHp 
     6  3.532699791229, 3.635600237483, 3.766800388578, 3.877300058102, NeHp 
     7  3.952799889519, 3.982099599452, 4.000000000000,     63*0.0D+00/ NeHp 
      DATA  K_NeHp/                                                     071215
     1 -9.55453697D-05, 1.37610976D-01, 2.81095921D-01, 7.10816364D-01, NeHp 
     2  1.68909236D+00, 2.73792813D+00, 3.65426191D+00, 4.68110968D+00, NeHp 
     3  5.48817547D+00, 6.08068545D+00, 6.57593048D+00, 6.89145354D+00, NeHp 
     4  7.17528263D+00, 7.64457538D+00, 8.11459680D+00, 8.56246600D+00, NeHp 
     5  8.92253966D+00, 9.23031913D+00, 9.48477238D+00, 9.77177428D+00, NeHp 
     6  9.85245546D+00, 9.92051769D+00, 9.98155551D+00, 1.00199730D+01, NeHp 
     7  1.00468982D+01, 1.00585952D+01, 1.00662249D+01,     63*0.0D+00/ NeHp 
      DATA TK_MgHp/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, MgHp 
     2  0.830599932887, 0.925099962653, 1.031700176058, 1.104800020971, MgHp 
     3  1.178899861440, 1.324600057872, 1.483899923913, 1.644999980058, MgHp 
     4  1.758400198836, 1.875100010191, 1.989899775079, 2.107399863200, MgHp 
     5  2.338700133397, 2.558600380205, 2.713500134929, 2.968900245552, MgHp 
     6  3.139100329033, 3.380600119670, 3.555200305062, 3.680500262430, MgHp 
     7  3.772300299502, 3.883200195138, 3.950799845848, 4.000000000000, MgHp 
     8      62*0.0D+00/                                                 MgHp 
      DATA  K_MgHp/                                                     071215
     1 -5.18624672D-05, 1.44617678D-01, 2.93891435D-01, 7.43064873D-01, MgHp 
     2  1.76135217D+00, 2.77409657D+00, 3.71209458D+00, 4.25434715D+00, MgHp 
     3  4.73536554D+00, 5.52011293D+00, 6.19490354D+00, 6.74125705D+00, MgHp 
     4  7.06723552D+00, 7.36536971D+00, 7.63030285D+00, 7.87919859D+00, MgHp 
     5  8.32268039D+00, 8.70639632D+00, 8.95946833D+00, 9.33538251D+00, MgHp 
     6  9.54726628D+00, 9.78948399D+00, 9.92175497D+00, 9.99210308D+00, MgHp 
     7  1.00307859D+01, 1.00678902D+01, 1.00887968D+01, 1.01044433D+01, MgHp 
     8      62*0.0D+00/                                                 MgHp 
      DATA TK_AlHp/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751200018282, AlHp 
     2  0.831399951226, 0.926099938541, 1.033300133308, 1.106700064771, AlHp 
     3  1.181299862980, 1.328899959963, 1.489699789377, 1.652300036787, AlHp 
     4  1.764900106541, 1.881300096532, 2.112199722628, 2.336600083617, AlHp 
     5  2.550100183908, 2.722300115050, 2.944699693873, 3.111099680836, AlHp 
     6  3.308300371730, 3.459800030575, 3.584400028533, 3.706499955700, AlHp 
     7  3.833700007319, 3.935499926925, 3.974699933131, 4.000000000000, AlHp 
     8      62*0.0D+00/                                                 AlHp 
      DATA  K_AlHp/                                                     071215
     1 -5.10432899D-05, 1.32000851D-01, 2.69643108D-01, 6.82458473D-01, AlHp 
     2  1.61904254D+00, 2.54825045D+00, 3.41379398D+00, 3.91455335D+00, AlHp 
     3  4.36095684D+00, 5.09613615D+00, 5.72977311D+00, 6.24698546D+00, AlHp 
     4  6.55325827D+00, 6.83690989D+00, 7.32934369D+00, 7.74953978D+00, AlHp 
     5  8.11654927D+00, 8.39361898D+00, 8.71676719D+00, 8.92331658D+00, AlHp 
     6  9.12387662D+00, 9.24457698D+00, 9.32160589D+00, 9.38065234D+00, AlHp 
     7  9.43529991D+00, 9.48405903D+00, 9.50603883D+00, 9.52159351D+00, AlHp 
     8      62*0.0D+00/                                                 AlHp 
      DATA TK_SiHp/                                                     071215
     1  0.699999789529, 0.709300031396, 0.718799827163, 0.749199968833, SiHp 
     2  0.826400006258, 0.915899993836, 1.020699991814, 1.091300132846, SiHp 
     3  1.163999932109, 1.307500011196, 1.454899915377, 1.601299923658, SiHp 
     4  1.867799939564, 1.972699852072, 2.074299791034, 2.167300003303, SiHp 
     5  2.254000095305, 2.451499836485, 2.602300206127, 2.952599870435, SiHp 
     6  3.147099819858, 3.357299837958, 3.465400152241, 3.575499817940, SiHp 
     7  3.760200229015, 3.819399728286, 3.889300339928, 3.951599863317, SiHp 
     8  4.000000000000,     61*0.0D+00/                                 SiHp 
      DATA  K_SiHp/                                                     071215
     1  2.12329692D-05, 1.43113102D-01, 2.86556135D-01, 7.27320379D-01, SiHp 
     2  1.73182462D+00, 2.71669035D+00, 3.66897252D+00, 4.21091836D+00, SiHp 
     3  4.70014420D+00, 5.50160769D+00, 6.15309130D+00, 6.67726524D+00, SiHp 
     4  7.42833313D+00, 7.68113589D+00, 7.91613712D+00, 8.12817690D+00, SiHp 
     5  8.32527582D+00, 8.76912327D+00, 9.09407432D+00, 9.76442960D+00, SiHp 
     6  1.00674570D+01, 1.03311360D+01, 1.04419342D+01, 1.05382381D+01, SiHp 
     7  1.06583343D+01, 1.06834634D+01, 1.07042250D+01, 1.07160661D+01, SiHp 
     8  1.07228088D+01,     61*0.0D+00/                                 SiHp 
      DATA TK_PHp/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749899984889, PHp  
     2  0.827299984476, 0.922400027755, 1.030200216135, 1.146900028895, PHp  
     3  1.270499993906, 1.408299891914, 1.595300021474, 1.691100043323, PHp  
     4  1.785700078852, 1.877000057195, 1.966699881013, 2.191599621573, PHp  
     5  2.370599902259, 2.564500092260, 2.780799744962, 3.016200365242, PHp  
     6  3.218100212642, 3.360699652657, 3.549100175014, 3.678200254409, PHp  
     7  3.831199952866, 3.930900258466, 3.972900062468, 4.000000000000, PHp  
     8      62*0.0D+00/                                                 PHp  
      DATA  K_PHp/                                                      071215
     1 -4.49287263D-05, 1.17009399D-01, 2.34508011D-01, 5.98979951D-01, PHp  
     2  1.42803176D+00, 2.30516196D+00, 3.14458091D+00, 3.90295535D+00, PHp  
     3  4.56888209D+00, 5.18013259D+00, 5.84813886D+00, 6.14160114D+00, PHp  
     4  6.41275791D+00, 6.66522366D+00, 6.90950012D+00, 7.51508685D+00, PHp  
     5  7.98650271D+00, 8.47261337D+00, 8.96656754D+00, 9.43112609D+00, PHp  
     6  9.75976143D+00, 9.95365499D+00, 1.01709297D+01, 1.03013073D+01, PHp  
     7  1.04385364D+01, 1.05137531D+01, 1.05408686D+01, 1.05567892D+01, PHp  
     8      62*0.0D+00/                                                 PHp  
      DATA TK_SHp/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.750099989774, SHp  
     2  0.828799948174, 0.917600031894, 1.025900120229, 1.164499943432, SHp  
     3  1.305000064398, 1.457799978405, 1.612900067607, 1.729400006975, SHp  
     4  1.849499919156, 1.965199920067, 2.085800052965, 2.329899925002, SHp  
     5  2.603600111219, 2.801000225573, 3.022700261426, 3.229599560734, SHp  
     6  3.359699656748, 3.491699796362, 3.613099732205, 3.730199562285, SHp  
     7  3.807600354305, 3.881900164281, 3.955999959393, 4.000000000000, SHp  
     8      62*0.0D+00/                                                 SHp  
      DATA  K_SHp/                                                      071215
     1  5.36079906D-05, 1.45731736D-01, 2.93253461D-01, 7.45977677D-01, SHp  
     2  1.77740792D+00, 2.76336327D+00, 3.75488583D+00, 4.75962202D+00, SHp  
     3  5.55408726D+00, 6.23563648D+00, 6.79008095D+00, 7.14110082D+00, SHp  
     4  7.45959070D+00, 7.73464716D+00, 7.99611357D+00, 8.47028272D+00, SHp  
     5  8.94719139D+00, 9.26817331D+00, 9.60169557D+00, 9.87323927D+00, SHp  
     6  1.00207918D+01, 1.01535381D+01, 1.02650887D+01, 1.03707945D+01, SHp  
     7  1.04438506D+01, 1.05181588D+01, 1.05963861D+01, 1.06443610D+01, SHp  
     8      62*0.0D+00/                                                 SHp  
      DATA TK_HClp/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719899799358, 0.751800033832, HClp 
     2  0.832599978735, 0.928599878261, 1.040399963170, 1.173499986756, HClp 
     3  1.325100046487, 1.486899854326, 1.649800093269, 1.762200179611, HClp 
     4  1.878600096778, 2.104100098412, 2.366799810413, 2.471400150175, HClp 
     5  2.581499960589, 2.705899954237, 2.817399886685, 3.105899953658, HClp 
     6  3.261900266055, 3.409399817659, 3.623199970788, 3.795900101749, HClp 
     7  3.915899995815, 3.967400212576, 4.000000000000,     63*0.0D+00/ HClp 
      DATA  K_HClp/                                                     071215
     1 -7.67661713D-05, 1.34797659D-01, 2.76823613D-01, 6.99658838D-01, HClp 
     2  1.66573223D+00, 2.64617353D+00, 3.60239036D+00, 4.53444827D+00, HClp 
     3  5.38325493D+00, 6.09991122D+00, 6.67763021D+00, 7.01462038D+00, HClp 
     4  7.32383449D+00, 7.84024947D+00, 8.35683647D+00, 8.55004710D+00, HClp 
     5  8.75119210D+00, 8.97813434D+00, 9.18077212D+00, 9.67674647D+00, HClp 
     6  9.90850473D+00, 1.00963628D+01, 1.03157587D+01, 1.04487748D+01, HClp 
     7  1.05136809D+01, 1.05345647D+01, 1.05461213D+01,     63*0.0D+00/ HClp 
      DATA TK_ZnHp/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751100015691, ZnHp 
     2  0.831199946641, 0.924799969886, 1.032300160026, 1.105000025581, ZnHp 
     3  1.179299852157, 1.327199998671, 1.486099872882, 1.646700020153, ZnHp 
     4  1.757400173639, 1.869999884022, 2.090900161798, 2.314900076804, ZnHp 
     5  2.554800292449, 2.752900058583, 3.018600420247, 3.417800035036, ZnHp 
     6  3.568399801534, 3.676700218706, 3.793000040950, 3.891200268117, ZnHp 
     7  3.958600016166, 3.984199646599, 4.000000000000,     63*0.0D+00/ ZnHp 
      DATA  K_ZnHp/                                                     071215
     1  8.68607121D-05, 1.46265730D-01, 2.98623271D-01, 7.54100478D-01, ZnHp 
     2  1.78901029D+00, 2.80440445D+00, 3.76137823D+00, 4.30614416D+00, ZnHp 
     3  4.79309766D+00, 5.59509496D+00, 6.27113285D+00, 6.81807204D+00, ZnHp 
     4  7.13793654D+00, 7.42777390D+00, 7.92185175D+00, 8.35643894D+00, ZnHp 
     5  8.77743980D+00, 9.10039936D+00, 9.48725563D+00, 9.91655822D+00, ZnHp 
     6  1.00227863D+01, 1.00749897D+01, 1.01076347D+01, 1.01240075D+01, ZnHp 
     7  1.01356895D+01, 1.01411941D+01, 1.01450292D+01,     63*0.0D+00/ ZnHp 
      DATA TK_HBrp/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719899799358, 0.751900036424, HBrp 
     2  0.832799983320, 0.928399883083, 1.041199980920, 1.182799894212, HBrp 
     3  1.340000164213, 1.492199828564, 1.647000027229, 1.764700111953, HBrp 
     4  1.887299935549, 2.006000137774, 2.130200095609, 2.380100128725, HBrp 
     5  2.613399722247, 2.774500137332, 2.985799671564, 3.218700226741, HBrp 
     6  3.501300020381, 3.640000328190, 3.802800249085, 3.915099975260, HBrp 
     7  3.967100205852, 4.000000000000,     64*0.0D+00/                 HBrp 
      DATA  K_HBrp/                                                     071215
     1 -5.98864081D-05, 1.35018742D-01, 2.77255543D-01, 7.01973270D-01, HBrp 
     2  1.66999132D+00, 2.64614823D+00, 3.60727598D+00, 4.58193454D+00, HBrp 
     3  5.43184943D+00, 6.08257004D+00, 6.61922768D+00, 6.96542267D+00, HBrp 
     4  7.28375221D+00, 7.56055929D+00, 7.82509100D+00, 8.30330085D+00, HBrp 
     5  8.70805325D+00, 8.97199367D+00, 9.29860297D+00, 9.62676361D+00, HBrp 
     6  9.97338696D+00, 1.01180410D+01, 1.02582760D+01, 1.03296540D+01, HBrp 
     7  1.03554735D+01, 1.03699728D+01,     64*0.0D+00/                 HBrp 
      DATA TK_CdHp/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751100015691, CdHp 
     2  0.831099944349, 0.926199936130, 1.033200135980, 1.107500083213, CdHp 
     3  1.182899896294, 1.330199939502, 1.489999782419, 1.651400060734, CdHp 
     4  1.761700193143, 1.874800002769, 2.097400325260, 2.325599831892, CdHp 
     5  2.550200186218, 2.723700013796, 2.983699624364, 3.126400030751, CdHp 
     6  3.268700429335, 3.413799929783, 3.558700378105, 3.678400259170, CdHp 
     7  3.798800162549, 3.893900069155, 3.959800042369, 3.984599655579, CdHp 
     8  4.000000000000,     61*0.0D+00/                                 CdHp 
      DATA  K_CdHp/                                                     071215
     1 -1.98611724D-06, 1.44646646D-01, 2.95364362D-01, 7.45673668D-01, CdHp 
     2  1.76647578D+00, 2.78337556D+00, 3.72211849D+00, 4.27100929D+00, CdHp 
     3  4.75749678D+00, 5.54453538D+00, 6.21555421D+00, 6.75851330D+00, CdHp 
     4  7.07382061D+00, 7.36225172D+00, 7.85613414D+00, 8.29580267D+00, CdHp 
     5  8.68867493D+00, 8.97185560D+00, 9.35227134D+00, 9.53005030D+00, CdHp 
     6  9.68233144D+00, 9.81112024D+00, 9.91012595D+00, 9.96518553D+00, CdHp 
     7  9.99628331D+00, 1.00104322D+01, 1.00201412D+01, 1.00246253D+01, CdHp 
     8  1.00277770D+01,     61*0.0D+00/                                 CdHp 
      DATA TK_HgHp/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, HgHp 
     2  0.830299926010, 0.924399979531, 1.030900197432, 1.103799997918, HgHp 
     3  1.177999882326, 1.324000071534, 1.482099965666, 1.641799904584, HgHp 
     4  1.754300095529, 1.868599919367, 2.090900161798, 2.324099799412, HgHp 
     5  2.561600298661, 2.680000308330, 2.785799871642, 3.038399943696, HgHp 
     6  3.226999749384, 3.409599821704, 3.527999875175, 3.716100166967, HgHp 
     7  3.822799752662, 3.921200127820, 3.969500259638, 4.000000000000, HgHp 
     8      62*0.0D+00/                                                 HgHp 
      DATA  K_HgHp/                                                     071215
     1 -1.18951500D-04, 1.45338309D-01, 2.93945845D-01, 7.45706898D-01, HgHp 
     2  1.76761340D+00, 2.78276317D+00, 3.72626439D+00, 4.27065731D+00, HgHp 
     3  4.75544604D+00, 5.54628072D+00, 6.21976030D+00, 6.76487049D+00, HgHp 
     4  7.09045226D+00, 7.38458922D+00, 7.88116379D+00, 8.33211395D+00, HgHp 
     5  8.74767974D+00, 8.94343852D+00, 9.11187392D+00, 9.47865937D+00, HgHp 
     6  9.70801564D+00, 9.88871055D+00, 9.98348549D+00, 1.00927995D+01, HgHp 
     7  1.01298569D+01, 1.01532245D+01, 1.01645055D+01, 1.01727331D+01, HgHp 
     8      62*0.0D+00/                                                 HgHp 
      DATA TK_CHm/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.751000013099, CHm  
     2  0.830599932887, 0.923300006054, 1.034800093230, 1.182499887966, CHm  
     3  1.322400107965, 1.470099782559, 1.619199910077, 1.737699835728, CHm  
     4  1.859800131898, 1.973299864733, 2.091500176887, 2.328199888191, CHm  
     5  2.599100354753, 2.738899750604, 2.875800016063, 3.108099794639, CHm  
     6  3.208899997657, 3.314500089104, 3.429800289204, 3.585800061333, CHm  
     7  3.687799754001, 3.847100304266, 4.000000000000,     63*0.0D+00/ CHm  
      DATA  K_CHm/                                                      071215
     1  2.18508180D-05, 1.43609709D-01, 2.91961186D-01, 7.41054333D-01, CHm  
     2  1.76319175D+00, 2.77679949D+00, 3.78676380D+00, 4.84259983D+00, CHm  
     3  5.61683538D+00, 6.26376153D+00, 6.79144035D+00, 7.14610334D+00, CHm  
     4  7.46730100D+00, 7.73524196D+00, 7.99038950D+00, 8.44974004D+00, CHm  
     5  8.92227577D+00, 9.15227711D+00, 9.36980507D+00, 9.71332477D+00, CHm  
     6  9.84713413D+00, 9.97466703D+00, 1.00990472D+01, 1.02485573D+01, CHm  
     7  1.03393547D+01, 1.04706312D+01, 1.05832169D+01,     63*0.0D+00/ CHm  
      DATA TK_OHm/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, OHm  
     2  0.829799923973, 0.928799873439, 1.035300079871, 1.172800003001, OHm  
     3  1.331899978482, 1.474099876223, 1.613200060105, 1.723199863164, OHm  
     4  1.839100127189, 1.939300024326, 2.046400115911, 2.166699991420, OHm  
     5  2.287299891951, 2.576199840085, 2.735599679246, 2.891500238202, OHm  
     6  3.063700164272, 3.257200159321, 3.386600271607, 3.511600260573, OHm  
     7  3.652499680098, 3.740199769192, 3.816399944869, 3.923400176698, OHm  
     8  3.970100263658, 3.987799727422, 4.000000000000,     59*0.0D+00/ OHm  
      DATA  K_OHm/                                                      071215
     1 -2.33219858D-04, 1.54029753D-01, 3.11805769D-01, 7.92649622D-01, OHm  
     2  1.88236063D+00, 3.03994777D+00, 4.07387632D+00, 5.14874252D+00, OHm  
     3  6.10634288D+00, 6.77005032D+00, 7.29557679D+00, 7.64883750D+00, OHm  
     4  7.97796383D+00, 8.23633145D+00, 8.49281722D+00, 8.76264336D+00, OHm  
     5  9.01780778D+00, 9.57978345D+00, 9.86594017D+00, 1.01327278D+01, OHm  
     6  1.04120444D+01, 1.06985377D+01, 1.08677616D+01, 1.10112973D+01, OHm  
     7  1.11479583D+01, 1.12179650D+01, 1.12670623D+01, 1.13123670D+01, OHm  
     8  1.13223736D+01, 1.13246448D+01, 1.13257460D+01,     59*0.0D+00/ OHm  
      DATA TK_SiHm/                                                     071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750199992366, SiHm 
     2  0.828699950594, 0.920300078390, 1.026200127637, 1.165999977399, SiHm 
     3  1.307300015452, 1.462899954879, 1.620799907760, 1.739799785229, SiHm 
     4  1.863100058222, 1.983199931443, 2.106099955859, 2.349000366648, SiHm 
     5  2.592900221883, 2.759600223544, 2.978499650649, 3.166399987045, SiHm 
     6  3.304600290302, 3.488399721202, 3.565100038702, 3.633000183884, SiHm 
     7  3.797100126908, 3.910499857069, 4.000000000000,     63*0.0D+00/ SiHm 
      DATA  K_SiHm/                                                     071215
     1  8.45595955D-05, 1.40788615D-01, 2.84679774D-01, 7.21400445D-01, SiHm 
     2  1.71216374D+00, 2.68706500D+00, 3.61579040D+00, 4.59049843D+00, SiHm 
     3  5.36078083D+00, 6.03092020D+00, 6.57650368D+00, 6.92367092D+00, SiHm 
     4  7.24080864D+00, 7.51817184D+00, 7.77766760D+00, 8.24018280D+00, SiHm 
     5  8.66148521D+00, 8.93169801D+00, 9.25807191D+00, 9.49922663D+00, SiHm 
     6  9.65077854D+00, 9.83281783D+00, 9.90798445D+00, 9.97493622D+00, SiHm 
     7  1.01313303D+01, 1.02279266D+01, 1.02981374D+01,     63*0.0D+00/ SiHm 
      DATA TK_HSm/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.750099989774, HSm  
     2  0.828699950594, 0.917400027417, 1.025600112820, 1.165199959283, HSm  
     3  1.306300036733, 1.459100006659, 1.613500052604, 1.730400011273, HSm  
     4  1.848299948110, 1.970599807755, 2.106399934476, 2.334000021985, HSm  
     5  2.556300327090, 2.706199961688, 2.839000145017, 3.097600322289, HSm  
     6  3.265100342893, 3.439599613744, 3.595400279055, 3.785599874318, HSm  
     7  3.905099734182, 3.963300120691, 4.000000000000,     63*0.0D+00/ HSm  
      DATA  K_HSm/                                                      071215
     1  3.32361642D-05, 1.52238632D-01, 3.06356842D-01, 7.79238396D-01, HSm  
     2  1.85482790D+00, 2.88259746D+00, 3.91546509D+00, 4.96774292D+00, HSm  
     3  5.79386772D+00, 6.49687380D+00, 7.06426770D+00, 7.42510452D+00, HSm  
     4  7.74466218D+00, 8.04045386D+00, 8.33754585D+00, 8.78925138D+00, HSm  
     5  9.20388321D+00, 9.47576848D+00, 9.71010475D+00, 1.01321851D+01, HSm  
     6  1.03669916D+01, 1.05718510D+01, 1.07194435D+01, 1.08504383D+01, HSm  
     7  1.09000686D+01, 1.09156532D+01, 1.09234868D+01,     63*0.0D+00/ HSm  
      DATA TK_CN/                                                       071215
     1  0.699999789529, 0.709800044400, 0.720399806410, 0.753100067524, CN   
     2  0.836100058969, 0.936699975035, 1.048300138454, 1.128600088053, CN   
     3  1.208399979922, 1.352799977659, 1.528800200485, 1.715399908843, CN   
     4  1.950999828335, 2.184100012585, 2.416099998347, 2.640200325845, CN   
     5  2.811800286586, 3.013700307946, 3.158099787941, 3.354500049371, CN   
     6  3.472200094056, 3.592000203927, 3.695399712520, 3.825799825190, CN   
     7  3.920000101159, 4.000000000000,     64*0.0D+00/                 CN   
      DATA  K_CN/                                                       071215
     1 -1.06579053D-04, 1.59318519D-01, 3.28271241D-01, 8.27655297D-01, CN   
     2  1.96205058D+00, 3.12337722D+00, 4.19687527D+00, 4.85839772D+00, CN   
     3  5.44075535D+00, 6.33954132D+00, 7.22174578D+00, 7.96274106D+00, CN   
     4  8.69426224D+00, 9.26940622D+00, 9.75117932D+00, 1.01638667D+01, CN   
     5  1.04543604D+01, 1.07641665D+01, 1.09591287D+01, 1.11847434D+01, CN   
     6  1.12974527D+01, 1.13948254D+01, 1.14668435D+01, 1.15495592D+01, CN   
     7  1.16099680D+01, 1.16634898D+01,     64*0.0D+00/                 CN   
      DATA TK_CO/                                                       071215
     1  0.699999789529, 0.709700041799, 0.720099799225, 0.752300046791, CO   
     2  0.834000010828, 0.933199906847, 1.043400029734, 1.121699922148, CO   
     3  1.199599810053, 1.339800159627, 1.502600053304, 1.689300055244, CO   
     4  1.801300147671, 1.912899902079, 2.027499911287, 2.132400146174, CO   
     5  2.382500180309, 2.624999988787, 2.808800389852, 3.015700353783, CO   
     6  3.137300287782, 3.263300299672, 3.447299737392, 3.670500071131, CO   
     7  3.808600376226, 3.930600280088, 3.972600084024, 4.000000000000, CO   
     8      62*0.0D+00/                                                 CO   
      DATA  K_CO/                                                       071215
     1 -3.34051951D-04, 1.67284757D-01, 3.43395631D-01, 8.65968750D-01, CO   
     2  2.05313838D+00, 3.27059135D+00, 4.39615280D+00, 5.08033585D+00, CO   
     3  5.68305666D+00, 6.60778540D+00, 7.47700611D+00, 8.26986811D+00, CO   
     4  8.66779633D+00, 9.02199622D+00, 9.35239232D+00, 9.63176251D+00, CO   
     5  1.02289607D+01, 1.07346492D+01, 1.10787480D+01, 1.14238504D+01, CO   
     6  1.16023141D+01, 1.17670379D+01, 1.19734915D+01, 1.21853158D+01, CO   
     7  1.23070927D+01, 1.24127592D+01, 1.24479441D+01, 1.24700742D+01, CO   
     8      62*0.0D+00/                                                 CO   
      DATA TK_CF/                                                       071215
     1  0.699999789529, 0.710200044546, 0.722199849517, 0.757000168598, CF   
     2  0.794000042819, 0.845900006017, 0.951399830436, 1.077399988702, CF   
     3  1.223799881972, 1.372099926670, 1.541099999033, 1.723899879400, CF   
     4  1.905299948986, 2.096000290053, 2.264000322373, 2.431400185652, CF   
     5  2.595500277603, 2.819399743864, 3.005900131657, 3.215300146848, CF   
     6  3.558100365584, 3.684699969909, 3.818599786041, 3.927000256680, CF   
     7  3.971400170248, 4.000000000000,     64*0.0D+00/                 CF   
      DATA  K_CF/                                                       071215
     1 -1.10490697D-04, 1.57395208D-01, 3.38653181D-01, 8.40961168D-01, CF   
     2  1.33944368D+00, 1.98305795D+00, 3.12111023D+00, 4.24218371D+00, CF   
     3  5.30061822D+00, 6.17249768D+00, 6.97890341D+00, 7.67884373D+00, CF   
     4  8.24196157D+00, 8.73535287D+00, 9.11610010D+00, 9.46549020D+00, CF   
     5  9.78658017D+00, 1.01827339D+01, 1.04632375D+01, 1.07204363D+01, CF   
     6  1.10335560D+01, 1.11262424D+01, 1.12134460D+01, 1.12730985D+01, CF   
     7  1.12938126D+01, 1.13058979D+01,     64*0.0D+00/                 CF   
      DATA TK_SiC/                                                      071215
     1  0.699999789529, 0.710700031907, 0.724399902204, 0.762200188320, SiC  
     2  0.802800114178, 0.859400122722, 0.977999956979, 1.106000048634, SiC  
     3  1.219799804366, 1.336500083959, 1.452299858868, 1.558299979405, SiC  
     4  1.740199784555, 1.920000075664, 2.070699700785, 2.309400414643, SiC  
     5  2.535499847909, 2.699299792087, 2.859799657358, 3.093500223498, SiC  
     6  3.372499950555, 3.519000433385, 3.742599822577, 3.830699941976, SiC  
     7  3.915699990677, 3.966800199129, 3.986899707217, 4.000000000000, SiC  
     8      62*0.0D+00/                                                 SiC  
      DATA  K_SiC/                                                      071215
     1  1.57295457D-04, 1.51470561D-01, 3.40603027D-01, 8.37356105D-01, SiC  
     2  1.33287542D+00, 1.96520255D+00, 3.10596587D+00, 4.11812328D+00, SiC  
     3  4.87400965D+00, 5.54275578D+00, 6.12684436D+00, 6.61069496D+00, SiC  
     4  7.35860780D+00, 8.01672700D+00, 8.51560992D+00, 9.21287012D+00, SiC  
     5  9.76845551D+00, 1.01067294D+01, 1.03864687D+01, 1.07115292D+01, SiC  
     6  1.09942941D+01, 1.11052287D+01, 1.12481177D+01, 1.13033960D+01, SiC  
     7  1.13606739D+01, 1.14002076D+01, 1.14178055D+01, 1.14301779D+01, SiC  
     8      62*0.0D+00/                                                 SiC  
      DATA TK_CP/                                                       071215
     1  0.699999789529, 0.710300042018, 0.722399854307, 0.757500181557, CP   
     2  0.794900063134, 0.847199974650, 0.953199877533, 1.080399929986, CP   
     3  1.156299936140, 1.230599999746, 1.382100047842, 1.549300167240, CP   
     4  1.738399818895, 1.918000026767, 2.119199865002, 2.298200150884, CP   
     5  2.471700128406, 2.602700176924, 2.744199863763, 2.893200111822, CP   
     6  3.059200419395, 3.239599764868, 3.359499671848, 3.473799979360, CP   
     7  3.652199673098, 3.732899618064, 3.813100183111, 3.888500320939, CP   
     8  3.963200118450, 4.000000000000,     60*0.0D+00/                 CP   
      DATA  K_CP/                                                       071215
     1 -1.84411538D-05, 1.63069621D-01, 3.50398640D-01, 8.69271185D-01, CP   
     2  1.38464527D+00, 2.04704555D+00, 3.21250696D+00, 4.36243051D+00, CP   
     3  4.94714486D+00, 5.45963764D+00, 6.35358810D+00, 7.15209058D+00, CP   
     4  7.87388457D+00, 8.42772920D+00, 8.94148292D+00, 9.33364184D+00, CP   
     5  9.67349323D+00, 9.90872583D+00, 1.01417182D+01, 1.03608010D+01, CP   
     6  1.05715154D+01, 1.07604140D+01, 1.08601997D+01, 1.09339471D+01, CP   
     7  1.10198560D+01, 1.10579283D+01, 1.11009329D+01, 1.11469392D+01, CP   
     8  1.11968361D+01, 1.12227996D+01,     60*0.0D+00/                 CP   
      DATA TK_CS/                                                       071215
     1  0.699999789529, 0.710100047074, 0.721799839938, 0.756300150457, CS   
     2  0.792900017989, 0.844200047034, 0.948199837999, 1.072900107352, CS   
     3  1.146200013464, 1.218099840754, 1.365799988626, 1.531500194012, CS   
     4  1.708200007251, 1.893599941848, 2.094600254846, 2.263000300530, CS   
     5  2.429900281457, 2.592700217597, 2.804500299288, 3.019200433998, CS   
     6  3.251000024001, 3.582299979332, 3.645199960849, 3.709400021727, CS   
     7  3.781899787085, 3.848400333512, 3.937299797192, 3.975399882833, CS   
     8  3.989199758854, 4.000000000000,     60*0.0D+00/                 CS   
      DATA  K_CS/                                                       071215
     1 -1.27083777D-04, 1.69377666D-01, 3.61456824D-01, 9.02544579D-01, CS   
     2  1.43794734D+00, 2.12809043D+00, 3.34305857D+00, 4.54032922D+00, CS   
     3  5.13989409D+00, 5.66666298D+00, 6.59179770D+00, 7.42932678D+00, CS   
     4  8.14412423D+00, 8.74766709D+00, 9.28496752D+00, 9.67688843D+00, CS   
     5  1.00366417D+01, 1.03684468D+01, 1.07629780D+01, 1.11016933D+01, CS   
     6  1.13927141D+01, 1.17110363D+01, 1.17653564D+01, 1.18194520D+01, CS   
     7  1.18766348D+01, 1.19200842D+01, 1.19484701D+01, 1.19445454D+01, CS   
     8  1.19403115D+01, 1.19359261D+01,     60*0.0D+00/                 CS   
      DATA TK_CCl/                                                      071215
     1  0.699999789529, 0.709800044400, 0.720499808805, 0.753300072707, CCl  
     2  0.836600070431, 0.937099982828, 1.048900151767, 1.132100070281, CCl  
     3  1.215599894267, 1.367799937316, 1.539199996943, 1.722899856205, CCl  
     4  1.839700141311, 1.949299819635, 2.189799603862, 2.379700119417, CCl  
     5  2.514600343902, 2.654599728561, 2.793700061843, 2.983399617622, CCl  
     6  3.140800289963, 3.309900406942, 3.585500054305, 3.715400151868, CCl  
     7  3.852600186198, 4.000000000000,     64*0.0D+00/                 CCl  
      DATA  K_CCl/                                                      071215
     1  1.65825843D-04, 1.54460187D-01, 3.19522657D-01, 8.04369961D-01, CCl  
     2  1.90682215D+00, 3.03165189D+00, 4.07613610D+00, 4.74143451D+00, CCl  
     3  5.33096599D+00, 6.24226708D+00, 7.06859638D+00, 7.77396626D+00, CCl  
     4  8.14711072D+00, 8.45478247D+00, 9.01533243D+00, 9.37495672D+00, CCl  
     5  9.60657096D+00, 9.83895148D+00, 1.00642546D+01, 1.03531332D+01, CCl  
     6  1.05652850D+01, 1.07589111D+01, 1.10068897D+01, 1.11037290D+01, CCl  
     7  1.11966550D+01, 1.12884291D+01,     64*0.0D+00/                 CCl  
      DATA TK_CSe/                                                      071215
     1  0.699999789529, 0.709800044400, 0.720499808805, 0.753200070115, CSe  
     2  0.836300063553, 0.936799976983, 1.048600145111, 1.130400111918, CSe  
     3  1.212099969185, 1.360000137427, 1.531600191453, 1.712799976616, CSe  
     4  1.834900028336, 1.955299936848, 2.088800114537, 2.209800049693, CSe  
     5  2.324599810238, 2.424900186239, 2.688199714822, 2.837800118003, CSe  
     6  3.010900243775, 3.143000125799, 3.322099743478, 3.543000029112, CSe  
     7  3.637800282837, 3.732799615998, 3.889900354170, 3.955699952843, CSe  
     8  3.983099621903, 4.000000000000,     60*0.0D+00/                 CSe  
      DATA  K_CSe/                                                      071215
     1 -2.09326766D-04, 1.63430593D-01, 3.38456797D-01, 8.50833698D-01, CSe  
     2  2.01549861D+00, 3.20407144D+00, 4.30423981D+00, 4.99174506D+00, CSe  
     3  5.59775307D+00, 6.52844808D+00, 7.39650113D+00, 8.12548321D+00, CSe  
     4  8.53372371D+00, 8.88591353D+00, 9.23069398D+00, 9.51101051D+00, CSe  
     5  9.75491178D+00, 9.95343661D+00, 1.04158759D+01, 1.06414233D+01, CSe  
     6  1.08743824D+01, 1.10368765D+01, 1.12394992D+01, 1.14655277D+01, CSe  
     7  1.15573536D+01, 1.16463552D+01, 1.17760127D+01, 1.18176891D+01, CSe  
     8  1.18322691D+01, 1.18404960D+01,     60*0.0D+00/                 CSe  
      DATA TK_CBr/                                                      071215
     1  0.699999789529, 0.710200044546, 0.722199849517, 0.757200173782, CBr  
     2  0.794600056363, 0.846399993953, 0.951999846135, 1.078899949151, CBr  
     3  1.155699951482, 1.229700009114, 1.382200045408, 1.553700091464, CBr  
     4  1.730600006464, 1.958400015078, 2.089600130956, 2.229399589439, CBr  
     5  2.343000231691, 2.454199896753, 2.674600186057, 2.832600000942, CBr  
     6  2.989299750231, 3.218000210292, 3.493499838372, 3.652999691765, CBr  
     7  3.778099879784, 3.875300012858, 3.945899735807, 4.000000000000, CBr  
     8      62*0.0D+00/                                                 CBr  
      DATA  K_CBr/                                                      071215
     1  7.73724810D-05, 1.52555471D-01, 3.28030525D-01, 8.17045017D-01, CBr  
     2  1.30458555D+00, 1.92617568D+00, 3.02945153D+00, 4.12417487D+00, CBr  
     3  4.69053520D+00, 5.18027130D+00, 6.04681579D+00, 6.83735914D+00, CBr  
     4  7.49564173D+00, 8.17240355D+00, 8.49968976D+00, 8.81259183D+00, CBr  
     5  9.04575307D+00, 9.25878436D+00, 9.63865882D+00, 9.87378794D+00, CBr  
     6  1.00760345D+01, 1.03270973D+01, 1.05912933D+01, 1.07368702D+01, CBr  
     7  1.08502619D+01, 1.09383033D+01, 1.10023296D+01, 1.10518464D+01, CBr  
     8      62*0.0D+00/                                                 CBr  
      DATA TK_RhC/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720199801620, 0.752500051974, RhC  
     2  0.834600024583, 0.933799918537, 1.044200047484, 1.123999977450, RhC  
     3  1.203899883696, 1.350499926622, 1.517100042083, 1.695699924962, RhC  
     4  1.825100049954, 1.966399888824, 2.099400375556, 2.238799760734, RhC  
     5  2.427300231944, 2.567199900095, 2.709700048611, 2.965600169678, RhC  
     6  3.064400112077, 3.177800241905, 3.354300064472, 3.590300166363, RhC  
     7  3.687799754001, 3.782899810661, 3.868699863258, 3.933200092696, RhC  
     8  3.974699933131, 4.000000000000,     60*0.0D+00/                 RhC  
      DATA  K_RhC/                                                      071215
     1  3.55537359D-05, 1.62914322D-01, 3.35689219D-01, 8.45068468D-01, RhC  
     2  2.00451082D+00, 3.18869825D+00, 4.28710073D+00, 4.96655683D+00, RhC  
     3  5.56825508D+00, 6.50597601D+00, 7.36514814D+00, 8.10016568D+00, RhC  
     4  8.54074314D+00, 8.95482513D+00, 9.29550607D+00, 9.61406679D+00, RhC  
     5  9.99825978D+00, 1.02553462D+01, 1.04951677D+01, 1.08829825D+01, RhC  
     6  1.10246668D+01, 1.11852846D+01, 1.14307528D+01, 1.17321543D+01, RhC  
     7  1.18398158D+01, 1.19333949D+01, 1.20092276D+01, 1.20628421D+01, RhC  
     8  1.20970330D+01, 1.21182033D+01,     60*0.0D+00/                 RhC  
      DATA TK_IrC/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720299804015, 0.752700057157, IrC  
     2  0.835100036045, 0.934799938019, 1.045600078547, 1.125900023134, IrC  
     3  1.206199932878, 1.353099984316, 1.520899990367, 1.698699847770, IrC  
     4  1.826200019211, 1.953399888900, 2.101500283730, 2.238299748608, IrC  
     5  2.362899715926, 2.470100244507, 2.744499870147, 2.851700258349, IrC  
     6  2.958800013276, 3.152499660399, 3.344500265430, 3.450999817778, IrC  
     7  3.556700336366, 3.657799803766, 3.766700386160, 3.889600347049, IrC  
     8  4.000000000000,     61*0.0D+00/                                 IrC  
      DATA  K_IrC/                                                      071215
     1 -1.70623848D-04, 1.55337714D-01, 3.21879242D-01, 8.09859912D-01, IrC  
     2  1.92209119D+00, 3.06123036D+00, 4.11853671D+00, 4.77569090D+00, IrC  
     3  5.35795973D+00, 6.26541182D+00, 7.10412577D+00, 7.81636965D+00, IrC  
     4  8.24123265D+00, 8.61012533D+00, 8.98634137D+00, 9.29533280D+00, IrC  
     5  9.55225770D+00, 9.75777204D+00, 1.02209864D+01, 1.03765512D+01, IrC  
     6  1.05193742D+01, 1.07580230D+01, 1.09930735D+01, 1.11306432D+01, IrC  
     7  1.12733656D+01, 1.14127718D+01, 1.15588902D+01, 1.17087456D+01, IrC  
     8  1.18297609D+01,     61*0.0D+00/                                 IrC  
      DATA TK_PtC/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720199801620, 0.752500051974, PtC  
     2  0.834500022290, 0.933499912692, 1.043800038609, 1.124599991876, PtC  
     3  1.205599920048, 1.353599995411, 1.522900043561, 1.701699854318, PtC  
     4  1.819800187501, 1.937999996469, 2.158599825612, 2.405199733313, PtC  
     5  2.491299808324, 2.587100091791, 2.700099810193, 2.827599884477, PtC  
     6  2.981299570422, 3.087900098380, 3.183000071501, 3.410699848212, PtC  
     7  3.499399976070, 3.578999901556, 3.675000178242, 3.763900318467, PtC  
     8  3.883700207006, 3.964700152066, 4.000000000000,     59*0.0D+00/ PtC  
      DATA  K_PtC/                                                      071215
     1  6.08658684D-06, 1.65617844D-01, 3.41283338D-01, 8.59128969D-01, PtC  
     2  2.03613430D+00, 3.23686138D+00, 4.35112685D+00, 5.04852709D+00, PtC  
     3  5.66539538D+00, 6.62005765D+00, 7.49689275D+00, 8.23353725D+00, PtC  
     4  8.63831450D+00, 8.99298378D+00, 9.55472795D+00, 1.00858792D+01, PtC  
     5  1.02603596D+01, 1.04524429D+01, 1.06763259D+01, 1.09204061D+01, PtC  
     6  1.11911610D+01, 1.13588649D+01, 1.14937462D+01, 1.17667515D+01, PtC  
     7  1.18572061D+01, 1.19306243D+01, 1.20063587D+01, 1.20599028D+01, PtC  
     8  1.21035983D+01, 1.21187523D+01, 1.21239336D+01,     59*0.0D+00/ PtC  
      DATA TK_CNp/                                                      071215
     1  0.699999789529, 0.710100047074, 0.721699837543, 0.756000142682, CNp  
     2  0.792400006703, 0.844100049447, 0.948799823267, 1.066200095897, CNp  
     3  1.154199989837, 1.244699878182, 1.325500037379, 1.405499960260, CNp  
     4  1.491699818076, 1.572100119558, 1.716399882777, 2.082899993447, CNp  
     5  2.261600269949, 2.425100190048, 2.599300359039, 2.721300187375, CNp  
     6  2.840900187055, 3.000700015620, 3.221400155708, 3.328499882245, CNp  
     7  3.426000213096, 3.568799772786, 3.684000018662, 3.830499937620, CNp  
     8  3.925800230020, 3.971300177434, 4.000000000000,     59*0.0D+00/ CNp  
      DATA  K_CNp/                                                      071215
     1 -1.47682495D-05, 1.76489457D-01, 3.74608945D-01, 9.33082089D-01, CNp  
     2  1.48366302D+00, 2.19788431D+00, 3.43374408D+00, 4.54738268D+00, CNp  
     3  5.23428530D+00, 5.83714925D+00, 6.30485809D+00, 6.71634986D+00, CNp  
     4  7.11467858D+00, 7.45273836D+00, 7.99509325D+00, 9.08864047D+00, CNp  
     5  9.51166008D+00, 9.85582273D+00, 1.01898480D+01, 1.04079672D+01, CNp  
     6  1.06095894D+01, 1.08571885D+01, 1.11499265D+01, 1.12679936D+01, CNp  
     7  1.13584907D+01, 1.14553975D+01, 1.15021726D+01, 1.15359927D+01, CNp  
     8  1.15547881D+01, 1.15648086D+01, 1.15715872D+01,     59*0.0D+00/ CNp  
      DATA TK_COp/                                                      071215
     1  0.699999789529, 0.710000049601, 0.721299827964, 0.755000116765, COp  
     2  0.790599966072, 0.841200119418, 0.944099938662, 1.059799956925, COp  
     3  1.143899962763, 1.228899991874, 1.304600072910, 1.379500088047, COp  
     4  1.460400016380, 1.537600037892, 1.686299985032, 1.865100007730, COp  
     5  2.033599814086, 2.191899627583, 2.379900124190, 2.577499869288, COp  
     6  2.741799812693, 3.076999869821, 3.281499778472, 3.481799574832, COp  
     7  3.578999901556, 3.680600255465, 3.872299944992, 3.949699821606, COp  
     8  3.980799570266, 4.000000000000,     60*0.0D+00/                 COp  
      DATA  K_COp/                                                      071215
     1  3.07819544D-04, 1.71721779D-01, 3.61119494D-01, 9.00130261D-01, COp  
     2  1.42993580D+00, 2.11914850D+00, 3.32048203D+00, 4.41006545D+00, COp  
     3  5.06524190D+00, 5.63478219D+00, 6.07874070D+00, 6.47106133D+00, COp  
     4  6.85312659D+00, 7.18564385D+00, 7.75839518D+00, 8.35753621D+00, COp  
     5  8.85510779D+00, 9.27641569D+00, 9.72815994D+00, 1.01543852D+01, COp  
     6  1.04767047D+01, 1.10443101D+01, 1.13210429D+01, 1.15411175D+01, COp  
     7  1.16319211D+01, 1.17162243D+01, 1.18371298D+01, 1.18657771D+01, COp  
     8  1.18732712D+01, 1.18767098D+01,     60*0.0D+00/                 COp  
      DATA TK_CNm/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718799827163, 0.749299971126, CNm  
     2  0.826400006258, 0.917900038610, 1.022200028857, 1.093600074017, CNm  
     3  1.166099979664, 1.305300058013, 1.452699867562, 1.602799960979, CNm  
     4  1.734399915084, 1.872599948344, 2.011000253161, 2.151899667160, CNm  
     5  2.402199660122, 2.626600025462, 2.793100048256, 2.987999721012, CNm  
     6  3.195399690403, 3.298400154635, 3.396199905011, 3.552100240366, CNm  
     7  3.701299837307, 3.792800036757, 3.870999915584, 3.948799801285, CNm  
     8  4.000000000000,     61*0.0D+00/                                 CNm  
      DATA  K_CNm/                                                      071215
     1 -3.22622800D-04, 1.68716309D-01, 3.38046222D-01, 8.59246930D-01, CNm  
     2  2.03786775D+00, 3.21452385D+00, 4.31711837D+00, 4.95259740D+00, CNm  
     3  5.51595359D+00, 6.41091150D+00, 7.15463294D+00, 7.75768582D+00, CNm  
     4  8.19435383D+00, 8.58589172D+00, 8.92766450D+00, 9.23804225D+00, CNm  
     5  9.72660562D+00, 1.01212186D+01, 1.03948890D+01, 1.06901272D+01, CNm  
     6  1.09634814D+01, 1.10822716D+01, 1.11865836D+01, 1.13455930D+01, CNm  
     7  1.15066138D+01, 1.16166040D+01, 1.17190920D+01, 1.18286779D+01, CNm  
     8  1.19042319D+01,     61*0.0D+00/                                 CNm  
      DATA TK_CSm/                                                      071215
     1  0.699999789529, 0.710200044546, 0.722099847123, 0.756900166007, CSm  
     2  0.793800038305, 0.845600013255, 0.951099822586, 1.076900001885, CSm  
     3  1.223099866888, 1.370699896139, 1.536200073723, 1.722299842288, CSm  
     4  1.911399865406, 2.112899736865, 2.293800041309, 2.466000161546, CSm  
     5  2.616499793608, 2.829899939998, 2.996499918088, 3.186099846684, CSm  
     6  3.367899843566, 3.568499794347, 3.730799574681, 3.873099963090, CSm  
     7  3.950199832747, 4.000000000000,     64*0.0D+00/                 CSm  
      DATA  K_CSm/                                                      071215
     1  6.32762303D-05, 1.54898600D-01, 3.31612936D-01, 8.25483643D-01, CSm  
     2  1.31432033D+00, 1.94614193D+00, 3.06595874D+00, 4.16835989D+00, CSm  
     3  5.21102382D+00, 6.06888589D+00, 6.85294237D+00, 7.56176391D+00, CSm  
     4  8.14298644D+00, 8.65511060D+00, 9.05455077D+00, 9.40409569D+00, CSm  
     5  9.69194102D+00, 1.00644460D+01, 1.03160806D+01, 1.05566644D+01, CSm  
     6  1.07462322D+01, 1.09206128D+01, 1.10455618D+01, 1.11476780D+01, CSm  
     7  1.11999285D+01, 1.12322910D+01,     64*0.0D+00/                 CSm  
      DATA TK_BN/                                                       071215
     1  0.699999789529, 0.709700041799, 0.719999796830, 0.752200044199, BN   
     2  0.833600001659, 0.929499856560, 1.043400029734, 1.190200039249, BN   
     3  1.348899942883, 1.508500192033, 1.670000067731, 1.785000096136, BN   
     4  1.906499918908, 2.030899748561, 2.158399820882, 2.392000194496, BN   
     5  2.504800109415, 2.603500118520, 2.919200056105, 3.074399804708, BN   
     6  3.234999653147, 3.416299995566, 3.566799916525, 3.680100290289, BN   
     7  3.799000166742, 3.928100281119, 3.971800141507, 3.988299738648, BN   
     8  4.000000000000,     61*0.0D+00/                                 BN   
      DATA  K_BN/                                                       071215
     1 -3.90629298D-05, 1.51082146D-01, 3.08393265D-01, 7.80010975D-01, BN   
     2  1.84886435D+00, 2.91606625D+00, 3.96790850D+00, 5.05564770D+00, BN   
     3  5.97475364D+00, 6.70087170D+00, 7.28944853D+00, 7.64272095D+00, BN   
     4  7.97140382D+00, 8.27152263D+00, 8.55014228D+00, 9.00726319D+00, BN   
     5  9.21030151D+00, 9.38025627D+00, 9.86821516D+00, 1.00688701D+01, BN   
     6  1.02469787D+01, 1.04149825D+01, 1.05316918D+01, 1.06092790D+01, BN   
     7  1.06866356D+01, 1.07755255D+01, 1.08095600D+01, 1.08233731D+01, BN   
     8  1.08335909D+01,     61*0.0D+00/                                 BN   
      DATA TK_NO/                                                       071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749599978008, NO   
     2  0.827299984476, 0.919300069952, 1.024200078247, 1.096599997283, NO   
     3  1.170800049414, 1.314800056211, 1.464299920439, 1.612100087610, NO   
     4  1.753700080411, 1.889199884571, 2.005000114812, 2.121299908539, NO   
     5  2.231599586122, 2.349300373396, 2.587600103505, 2.764800349338, NO   
     6  2.947899764310, 3.080599956605, 3.229499567990, 3.382000155122, NO   
     7  3.522400282270, 3.641200243419, 3.753000058282, 3.877200055840, NO   
     8  3.965600172236, 4.000000000000,     60*0.0D+00/                 NO   
      DATA  K_NO/                                                       071215
     1 -4.20134220D-05, 1.66942859D-01, 3.32439762D-01, 8.45284176D-01, NO   
     2  2.00654003D+00, 3.16341704D+00, 4.24810724D+00, 4.87818738D+00, NO   
     3  5.44111582D+00, 6.34136994D+00, 7.07181646D+00, 7.64537292D+00, NO   
     4  8.09234177D+00, 8.45135992D+00, 8.72091309D+00, 8.96837237D+00, NO   
     5  9.18831558D+00, 9.41165192D+00, 9.83748821D+00, 1.01338144D+01, NO   
     6  1.04158529D+01, 1.05994783D+01, 1.07817704D+01, 1.09428429D+01, NO   
     7  1.10708455D+01, 1.11678706D+01, 1.12557110D+01, 1.13584927D+01, NO   
     8  1.14375405D+01, 1.14688065D+01,     60*0.0D+00/                 NO   
      DATA TK_NF/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750399997549, NF   
     2  0.829299936074, 0.923000013288, 1.029400206662, 1.103900000223, NF   
     3  1.179899838233, 1.326300019164, 1.478899988619, 1.632900038944, NF   
     4  1.757600178678, 1.883800029456, 2.149699643185, 2.360799665049, NF   
     5  2.544400056023, 2.711000078595, 2.871799931165, 3.080099946894, NF   
     6  3.275400073608, 3.634900223053, 3.768300424842, 3.887000285335, NF   
     7  3.956899979045, 4.000000000000,     64*0.0D+00/                 NF   
      DATA  K_NF/                                                       071215
     1  1.03447299D-04, 1.61211502D-01, 3.24105676D-01, 8.19440923D-01, NF   
     2  1.94183913D+00, 3.06041216D+00, 4.10254783D+00, 4.71540592D+00, NF   
     3  5.25929223D+00, 6.12213265D+00, 6.82642012D+00, 7.39449947D+00, NF   
     4  7.77921326D+00, 8.11845410D+00, 8.72407769D+00, 9.14630315D+00, NF   
     5  9.49106405D+00, 9.78199987D+00, 1.00335901D+01, 1.03086612D+01, NF   
     6  1.05162242D+01, 1.07910531D+01, 1.08619993D+01, 1.09166255D+01, NF   
     7  1.09480737D+01, 1.09681656D+01,     64*0.0D+00/                 NF   
      DATA TK_AlN/                                                      071215
     1  0.699999789529, 0.710000049601, 0.721299827964, 0.755000116765, AlN  
     2  0.840900126657, 0.940900017229, 1.061899996410, 1.199399814929, AlN  
     3  1.349999915527, 1.516300062953, 1.674899954961, 1.773200048862, AlN  
     4  1.867799939564, 2.064600085741, 2.278599843313, 2.479899533392, AlN  
     5  2.708700023776, 2.909999831178, 3.114699757845, 3.362999713642, AlN  
     6  3.567199887777, 3.671500094934, 3.771700342921, 3.845400266022, AlN  
     7  3.913899944428, 3.965900178959, 3.986599700481, 4.000000000000, AlN  
     8      62*0.0D+00/                                                 AlN  
      DATA  K_AlN/                                                      071215
     1  1.50206406D-05, 1.50415389D-01, 3.16650633D-01, 7.90082877D-01, AlN  
     2  1.85952992D+00, 2.89414664D+00, 3.90650269D+00, 4.81461793D+00, AlN  
     3  5.59332625D+00, 6.27231659D+00, 6.80961302D+00, 7.11078293D+00, AlN  
     4  7.38526502D+00, 7.91527728D+00, 8.42750476D+00, 8.84472131D+00, AlN  
     5  9.23709306D+00, 9.51175498D+00, 9.73412074D+00, 9.94613569D+00, AlN  
     6  1.00864520D+01, 1.01505979D+01, 1.02126057D+01, 1.02642537D+01, AlN  
     7  1.03251715D+01, 1.03873424D+01, 1.04176407D+01, 1.04392965D+01, AlN  
     8      62*0.0D+00/                                                 AlN  
      DATA TK_SiN/                                                      071215
     1  0.699999789529, 0.710200044546, 0.721999844728, 0.756600158232, SiN  
     2  0.793200024761, 0.845300020494, 0.951499833052, 1.069600174561, SiN  
     3  1.162099889084, 1.256800154076, 1.340200159239, 1.421800032670, SiN  
     4  1.498399958609, 1.573700078479, 1.732799953560, 1.858600104371, SiN  
     5  1.972599849961, 2.231199576422, 2.487099708096, 2.647399804201, SiN  
     6  2.789999978054, 2.947699759908, 3.102900170502, 3.341900206949, SiN  
     7  3.580999948875, 3.730299564351, 3.865899799281, 3.950099830563, SiN  
     8  3.980899572511, 4.000000000000,     60*0.0D+00/                 SiN  
      DATA  K_SiN/                                                      071215
     1  1.32443020D-04, 1.60317977D-01, 3.41412156D-01, 8.47661589D-01, SiN  
     2  1.34534431D+00, 1.99277737D+00, 3.12098694D+00, 4.13097310D+00, SiN  
     3  4.78143073D+00, 5.34756092D+00, 5.78053382D+00, 6.15871024D+00, SiN  
     4  6.48296308D+00, 6.78003209D+00, 7.35725737D+00, 7.77733853D+00, SiN  
     5  8.13569802D+00, 8.87556198D+00, 9.50289065D+00, 9.84079222D+00, SiN  
     6  1.01040062D+01, 1.03537346D+01, 1.05601121D+01, 1.08174742D+01, SiN  
     7  1.10310237D+01, 1.11552163D+01, 1.12694279D+01, 1.13456919D+01, SiN  
     8  1.13764931D+01, 1.13970197D+01,     60*0.0D+00/                 SiN  
      DATA TK_PN/                                                       071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749799982595, PN   
     2  0.827799972376, 0.920300078390, 1.025600112820, 1.098099958916, PN   
     3  1.171800026207, 1.313800035750, 1.463899930279, 1.615200010096, PN   
     4  1.734099922298, 1.860200131437, 1.993199845467, 2.133700176054, PN   
     5  2.259300218702, 2.399099672604, 2.514300337633, 2.619199855761, PN   
     6  2.883400188338, 3.026799971532, 3.166899999219, 3.319699719281, PN   
     7  3.464900141423, 3.720300230303, 3.803200257854, 3.883400199885, PN   
     8  3.956399968127, 3.983499630883, 4.000000000000,     59*0.0D+00/ PN   
      DATA  K_PN/                                                       071215
     1  1.08031853D-04, 1.68291927D-01, 3.36688998D-01, 8.54514167D-01, PN   
     2  2.02647256D+00, 3.19386617D+00, 4.28483180D+00, 4.91638035D+00, PN   
     3  5.47606809D+00, 6.36633165D+00, 7.10318695D+00, 7.69472624D+00, PN   
     4  8.08258047D+00, 8.43936924D+00, 8.77011531D+00, 9.08208856D+00, PN   
     5  9.33700618D+00, 9.60106345D+00, 9.80623320D+00, 9.98392810D+00, PN   
     6  1.03859244D+01, 1.05716124D+01, 1.07301578D+01, 1.08799679D+01, PN   
     7  1.10060939D+01, 1.12349224D+01, 1.13269750D+01, 1.14287937D+01, PN   
     8  1.15317826D+01, 1.15721428D+01, 1.15972223D+01,     59*0.0D+00/ PN   
      DATA TK_NS/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750600002732, NS   
     2  0.829799923973, 0.923799993998, 1.030500208120, 1.104700018666, NS   
     3  1.180299842159, 1.325800030548, 1.478399976911, 1.634499999098, NS   
     4  1.762500171492, 1.890199867481, 2.024100162299, 2.190999609554, NS   
     5  2.412699915895, 2.527699884937, 2.627000034631, 2.882800173667, NS   
     6  3.036399896453, 3.185799868440, 3.341500197952, 3.501100015718, NS   
     7  3.618599865645, 3.723799987775, 3.859699685624, 3.952899891703, NS   
     8  4.000000000000,     61*0.0D+00/                                 NS   
      DATA  K_NS/                                                       071215
     1 -2.34257739D-04, 1.65210279D-01, 3.32484849D-01, 8.44275311D-01, NS   
     2  2.00011609D+00, 3.15016842D+00, 4.21993843D+00, 4.84430850D+00, NS   
     3  5.39770061D+00, 6.27475002D+00, 6.99414152D+00, 7.58005844D+00, NS   
     4  7.97962136D+00, 8.32503976D+00, 8.64364711D+00, 8.99252730D+00, NS   
     5  9.39633141D+00, 9.58691439D+00, 9.74262736D+00, 1.01021108D+01, NS   
     6  1.02858513D+01, 1.04417672D+01, 1.05835695D+01, 1.07130604D+01, NS   
     7  1.08034702D+01, 1.08855812D+01, 1.10002706D+01, 1.10857227D+01, NS   
     8  1.11296885D+01,     61*0.0D+00/                                 NS   
      DATA TK_NCl/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, NCl  
     2  0.830699935179, 0.925399955419, 1.032800146667, 1.107800090129, NCl  
     3  1.184099921280, 1.331399967017, 1.487199847367, 1.643999956472, NCl  
     4  1.763600141723, 1.888999889937, 2.011700269636, 2.141700202088, NCl  
     5  2.364299749845, 2.592400211168, 2.963800128293, 3.138600317575, NCl  
     6  3.336500076907, 3.634200208622, 3.778099879784, 3.897599796505, NCl  
     7  3.960900066905, 4.000000000000,     64*0.0D+00/                 NCl  
      DATA  K_NCl/                                                      071215
     1  1.08267583D-04, 1.58499879D-01, 3.21885658D-01, 8.13270175D-01, NCl  
     2  1.92756851D+00, 3.03564833D+00, 4.06538309D+00, 4.66892920D+00, NCl  
     3  5.20309007D+00, 6.05259746D+00, 6.75553502D+00, 7.32025199D+00, NCl  
     4  7.68205734D+00, 8.01442568D+00, 8.30451772D+00, 8.58343899D+00, NCl  
     5  9.01195819D+00, 9.40055273D+00, 9.92608376D+00, 1.01241365D+01, NCl  
     6  1.03109603D+01, 1.05238699D+01, 1.06008167D+01, 1.06572506D+01, NCl  
     7  1.06875249D+01, 1.07074772D+01,     64*0.0D+00/                 NCl  
      DATA TK_TiN/                                                      071215
     1  0.699999789529, 0.709100026195, 0.718199842329, 0.747799936720, TiN  
     2  0.822500100643, 0.911599897571, 1.014100115482, 1.149800092824, TiN  
     3  1.278800193352, 1.419100059636, 1.563700025120, 1.690400061335, TiN  
     4  1.813600032814, 1.941500003892, 2.064700078287, 2.200799811058, TiN  
     5  2.318299832690, 2.460200030731, 2.577899878273, 2.698499774234, TiN  
     6  2.825899843440, 2.961500075411, 3.093100213860, 3.224499930779, TiN  
     7  3.368499859476, 3.493499838372, 3.613899751614, 3.732199603603, TiN  
     8  3.839400131472, 3.937099811607, 3.975399882833, 3.989199758854, TiN  
     9  4.000000000000,     57*0.0D+00/                                 TiN  
      DATA  K_TiN/                                                      071215
     1  2.30446488D-04, 1.58328523D-01, 3.13498718D-01, 7.98095461D-01, TiN  
     2  1.89640455D+00, 3.00474897D+00, 4.05832689D+00, 5.16900301D+00, TiN  
     3  5.99805906D+00, 6.71584884D+00, 7.31017523D+00, 7.74484480D+00, TiN  
     4  8.11550812D+00, 8.46675203D+00, 8.78711314D+00, 9.12978001D+00, TiN  
     5  9.41734304D+00, 9.74889452D+00, 1.00040014D+01, 1.02412415D+01, TiN  
     6  1.04627376D+01, 1.06663342D+01, 1.08358073D+01, 1.09837009D+01, TiN  
     7  1.11337427D+01, 1.12689316D+01, 1.14199831D+01, 1.16058228D+01, TiN  
     8  1.18198469D+01, 1.20587856D+01, 1.21637447D+01, 1.22029913D+01, TiN  
     9  1.22342082D+01,     57*0.0D+00/                                 TiN  
      DATA TK_AsN/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.750099989774, AsN  
     2  0.828499955435, 0.921500049456, 1.027300154802, 1.100299917234, AsN  
     3  1.174599961229, 1.317600113504, 1.468699812198, 1.620899909971, AsN  
     4  1.739699787634, 1.865499997631, 1.997399940914, 2.143400083321, AsN  
     5  2.270100446339, 2.399199665253, 2.627700050676, 2.736999709519, AsN  
     6  2.868599857900, 3.006400142815, 3.152399658121, 3.299500178307, AsN  
     7  3.443499658290, 3.684699969909, 3.773800190954, 3.862699726165, AsN  
     8  3.948799801285, 4.000000000000,     60*0.0D+00/                 AsN  
      DATA  K_AsN/                                                      071215
     1  9.70215206D-06, 1.67608090D-01, 3.37149101D-01, 8.56319132D-01, AsN  
     2  2.02906476D+00, 3.19662221D+00, 4.28613305D+00, 4.91770709D+00, AsN  
     3  5.47766161D+00, 6.36652207D+00, 7.10132181D+00, 7.69062965D+00, AsN  
     4  8.07467244D+00, 8.42790176D+00, 8.75400326D+00, 9.07615582D+00, AsN  
     5  9.33139203D+00, 9.57369322D+00, 9.96361439D+00, 1.01307614D+01, AsN  
     6  1.03131374D+01, 1.04815593D+01, 1.06364504D+01, 1.07715603D+01, AsN  
     7  1.08897525D+01, 1.10946217D+01, 1.11877922D+01, 1.12953620D+01, AsN  
     8  1.14123915D+01, 1.14862868D+01,     60*0.0D+00/                 AsN  
      DATA TK_SeN/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751400023466, SeN  
     2  0.831999964980, 0.927799897551, 1.036200055825, 1.111800094122, SeN  
     3  1.188900021222, 1.337600109181, 1.494599878904, 1.652500031465, SeN  
     4  1.770099971029, 1.894399959346, 2.020000464990, 2.162299904275, SeN  
     5  2.268600422852, 2.373199964304, 2.511700283305, 2.648599717260, SeN  
     6  2.814600086636, 3.004500100417, 3.116199789932, 3.226199807431, SeN  
     7  3.383300188042, 3.549100175014, 3.652999691765, 3.748899962714, SeN  
     8  3.881300150039, 4.000000000000,     60*0.0D+00/                 SeN  
      DATA  K_SeN/                                                      071215
     1  9.44898197D-05, 1.59768041D-01, 3.24350264D-01, 8.18944424D-01, SeN  
     2  1.93990050D+00, 3.05443019D+00, 4.08595602D+00, 4.68904021D+00, SeN  
     3  5.22361171D+00, 6.07190984D+00, 6.77197066D+00, 7.33398414D+00, SeN  
     4  7.68629633D+00, 8.01352340D+00, 8.30859735D+00, 8.61079157D+00, SeN  
     5  8.81934758D+00, 9.01203712D+00, 9.24649306D+00, 9.44992484D+00, SeN  
     6  9.65505269D+00, 9.84599038D+00, 9.94599491D+00, 1.00397086D+01, SeN  
     7  1.01673752D+01, 1.02951371D+01, 1.03717020D+01, 1.04387764D+01, SeN  
     8  1.05212781D+01, 1.05845697D+01,     60*0.0D+00/                 SeN  
      DATA TK_ZrN/                                                      071215
     1  0.699999789529, 0.709100026195, 0.717899849912, 0.747199922958, ZrN  
     2  0.821100134525, 0.908499896328, 1.008800187822, 1.076400015068, ZrN  
     3  1.145099989216, 1.276300133278, 1.413399928587, 1.553600093900, ZrN  
     4  1.678199879013, 1.804900062978, 1.938400005040, 2.140700271950, ZrN  
     5  2.264900342032, 2.384600225445, 2.506300149329, 2.621199901683, ZrN  
     6  2.860899665059, 2.978899621463, 3.095200264460, 3.280699760082, ZrN  
     7  3.478099671114, 3.575099808384, 3.681600185817, 3.774200162008, ZrN  
     8  3.863799751299, 3.947299767417, 3.979799566676, 4.000000000000, ZrN  
     9      58*0.0D+00/                                                 ZrN  
      DATA  K_ZrN/                                                      071215
     1  1.57250880D-04, 1.57532765D-01, 3.06954458D-01, 7.84904818D-01, ZrN  
     2  1.86878813D+00, 2.95608819D+00, 3.99108659D+00, 4.58326652D+00, ZrN  
     3  5.11266199D+00, 5.95807566D+00, 6.66156460D+00, 7.24189107D+00, ZrN  
     4  7.67100794D+00, 8.04550410D+00, 8.38989136D+00, 8.84397612D+00, ZrN  
     5  9.09797571D+00, 9.33533601D+00, 9.57523233D+00, 9.80186924D+00, ZrN  
     6  1.02613938D+01, 1.04704038D+01, 1.06621147D+01, 1.09490576D+01, ZrN  
     7  1.12556432D+01, 1.14151342D+01, 1.16041874D+01, 1.17870128D+01, ZrN  
     8  1.19853941D+01, 1.21908316D+01, 1.22757689D+01, 1.23297867D+01, ZrN  
     9      58*0.0D+00/                                                 ZrN  
      DATA TK_NOp/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749599978008, NOp  
     2  0.827199986897, 0.919300069952, 1.024100075777, 1.096100010072, NOp  
     3  1.169300052128, 1.309999957995, 1.458900002313, 1.610300132619, NOp  
     4  1.733699931917, 1.862200080944, 1.987499831090, 2.115499789747, NOp  
     5  2.352300218142, 2.612199694623, 2.796800132045, 3.048700154719, NOp  
     6  3.261900266055, 3.415099963990, 3.554900298801, 3.676600216326, NOp  
     7  3.793500051433, 3.911499882763, 3.964500147584, 4.000000000000, NOp  
     8      62*0.0D+00/                                                 NOp  
      DATA  K_NOp/                                                      071215
     1  1.27846357D-04, 1.73439556D-01, 3.45171261D-01, 8.77105125D-01, NOp  
     2  2.07885176D+00, 3.27701923D+00, 4.39620724D+00, 5.04267583D+00, NOp  
     3  5.61564643D+00, 6.52483083D+00, 7.27778064D+00, 7.88601930D+00, NOp  
     4  8.29712281D+00, 8.66601200D+00, 8.98216187D+00, 9.27199918D+00, NOp  
     5  9.74712266D+00, 1.02112321D+01, 1.05175474D+01, 1.09010531D+01, NOp  
     6  1.11793859D+01, 1.13489853D+01, 1.14833049D+01, 1.15885021D+01, NOp  
     7  1.16872811D+01, 1.17974142D+01, 1.18532112D+01, 1.18926431D+01, NOp  
     8      62*0.0D+00/                                                 NOp  
      DATA TK_NSp/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749899984889, NSp  
     2  0.827899969956, 0.920500073567, 1.025800117759, 1.098399951243, NSp  
     3  1.172300014604, 1.314400048027, 1.464599913059, 1.616699972589, NSp  
     4  1.742999842448, 1.877700074512, 2.018900439099, 2.162999918139, NSp  
     5  2.397899760811, 2.612299696925, 2.754400095515, 2.904099698501, NSp  
     6  3.071499732082, 3.245499899759, 3.395899926916, 3.530199734336, NSp  
     7  3.652799687098, 3.754200086721, 3.901999664565, 3.962000091557, NSp  
     8  3.985299671295, 4.000000000000,     60*0.0D+00/                 NSp  
      DATA  K_NSp/                                                      071215
     1  1.04480053D-04, 1.77175622D-01, 3.54452706D-01, 9.01172117D-01, NSp  
     2  2.13360461D+00, 3.36108379D+00, 4.50516416D+00, 5.16729987D+00, NSp  
     3  5.75389370D+00, 6.68272158D+00, 7.44859265D+00, 8.06329284D+00, NSp  
     4  8.48532871D+00, 8.87114027D+00, 9.22320246D+00, 9.54243534D+00, NSp  
     5  1.00040819D+01, 1.03818995D+01, 1.06111505D+01, 1.08305309D+01, NSp  
     6  1.10454029D+01, 1.12349656D+01, 1.13743775D+01, 1.14843307D+01, NSp  
     7  1.15804612D+01, 1.16663603D+01, 1.18221616D+01, 1.18999541D+01, NSp  
     8  1.19325843D+01, 1.19538604D+01,     60*0.0D+00/                 NSp  
      DATA TK_LiO/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, LiO  
     2  0.830299926010, 0.925099962653, 1.032500154683, 1.106100050939, LiO  
     3  1.180399844241, 1.322200112519, 1.476999944129, 1.636399951782, LiO  
     4  1.748499956168, 1.871099911235, 1.982699943112, 2.100500355007, LiO  
     5  2.213900130776, 2.322799771262, 2.560700362716, 2.724399963169, LiO  
     6  2.898099747551, 3.077099872326, 3.247699949829, 3.421500122968, LiO  
     7  3.569299736852, 3.685299928120, 3.738699737886, 3.792000019985, LiO  
     8  3.917400034356, 3.967300210335, 3.986999709462, 4.000000000000, LiO  
     9      58*0.0D+00/                                                 LiO  
      DATA  K_LiO/                                                      071215
     1 -8.81020451D-05, 1.46537033D-01, 2.96323850D-01, 7.51628759D-01, LiO  
     2  1.78191503D+00, 2.81391189D+00, 3.77451199D+00, 4.32927388D+00, LiO  
     3  4.81916028D+00, 5.59696759D+00, 6.26821952D+00, 6.82351888D+00, LiO  
     4  7.15656843D+00, 7.48376435D+00, 7.75833639D+00, 8.03126815D+00, LiO  
     5  8.28096477D+00, 8.50930462D+00, 8.96203012D+00, 9.22789880D+00, LiO  
     6  9.46393942D+00, 9.65736145D+00, 9.79870113D+00, 9.90760076D+00, LiO  
     7  9.98060886D+00, 1.00345488D+01, 1.00628417D+01, 1.00972858D+01, LiO  
     8  1.02351869D+01, 1.03279122D+01, 1.03719331D+01, 1.04033054D+01, LiO  
     9      58*0.0D+00/                                                 LiO  
      DATA TK_BeO/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720899818384, 0.754200096032, BeO  
     2  0.838700118571, 0.937499990621, 1.057300012856, 1.189000023304, BeO  
     3  1.330799953260, 1.489299798656, 1.653699999535, 1.852099955265, BeO  
     4  2.037799916013, 2.231199576422, 2.456799954790, 2.619899871875, BeO  
     5  2.796200118458, 3.001200026778, 3.144699998945, 3.215000139798, BeO  
     6  3.280899764679, 3.412399892945, 3.474999893338, 3.541099983667, BeO  
     7  3.598400345345, 3.655299745432, 3.746599911553, 3.831899968113, BeO  
     8  3.934200020622, 3.974299961872, 3.988899752118, 4.000000000000, BeO  
     9      58*0.0D+00/                                                 BeO  
      DATA  K_BeO/                                                      071215
     1  5.92006936D-05, 1.65434419D-01, 3.45199110D-01, 8.64940857D-01, BeO  
     2  2.03392124D+00, 3.16939006D+00, 4.28073362D+00, 5.24545523D+00, BeO  
     3  6.06138350D+00, 6.77697517D+00, 7.36689145D+00, 7.94644293D+00, BeO  
     4  8.41243792D+00, 8.85187797D+00, 9.32036603D+00, 9.63142047D+00, BeO  
     5  9.93789847D+00, 1.02464406D+01, 1.04285726D+01, 1.05067461D+01, BeO  
     6  1.05718902D+01, 1.06679437D+01, 1.06919085D+01, 1.06996271D+01, BeO  
     7  1.06929476D+01, 1.06772756D+01, 1.06449689D+01, 1.06232189D+01, BeO  
     8  1.06294903D+01, 1.06454612D+01, 1.06534908D+01, 1.06604565D+01, BeO  
     9      58*0.0D+00/                                                 BeO  
      DATA TK_BO/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751200018282, BO   
     2  0.830899939764, 0.925099962653, 1.036900037122, 1.176499917136, BO   
     3  1.325300041933, 1.481599977263, 1.644699972982, 1.754000087970, BO   
     4  1.871799928552, 1.978099966028, 2.085400044756, 2.304500300518, BO   
     5  2.572999768202, 2.762000281703, 3.006800151741, 3.160099833649, BO   
     6  3.317499875745, 3.591500192879, 3.693199666995, 3.793900059819, BO   
     7  3.903099689268, 4.000000000000,     64*0.0D+00/                 BO   
      DATA  K_BO/                                                       071215
     1  1.50370581D-04, 1.62489380D-01, 3.31699433D-01, 8.39314107D-01, BO   
     2  1.98684390D+00, 3.13519859D+00, 4.26405740D+00, 5.39611594D+00, BO   
     3  6.34443188D+00, 7.12742205D+00, 7.77647869D+00, 8.14258755D+00, BO   
     4  8.49340596D+00, 8.78168359D+00, 9.05262487D+00, 9.55903898D+00, BO   
     5  1.01113336D+01, 1.04608128D+01, 1.08585349D+01, 1.10697732D+01, BO   
     6  1.12555983D+01, 1.15141960D+01, 1.15931566D+01, 1.16629036D+01, BO   
     7  1.17274097D+01, 1.17757244D+01,     64*0.0D+00/                 BO   
      DATA TK_FO/                                                       071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.750099989774, FO   
     2  0.828399957855, 0.921800042222, 1.027900169619, 1.100899931065, FO   
     3  1.174899954267, 1.316200084858, 1.468599814658, 1.623799974088, FO   
     4  1.737499840538, 1.849999907092, 1.971199820417, 2.087800094013, FO   
     5  2.192599641606, 2.305000312163, 2.497499931873, 2.686299852342, FO   
     6  2.870099895084, 3.068199828732, 3.277899894468, 3.486099670194, FO   
     7  3.696099727005, 3.811800276964, 3.918200054911, 3.968100228263, FO   
     8  4.000000000000,     61*0.0D+00/                                 FO   
      DATA  K_FO/                                                       071215
     1  7.03685605D-05, 1.58557819D-01, 3.18906575D-01, 8.10089268D-01, FO   
     2  1.91929586D+00, 3.03108666D+00, 4.06827936D+00, 4.66879250D+00, FO   
     3  5.20010970D+00, 6.04054631D+00, 6.75305651D+00, 7.33280351D+00, FO   
     4  7.69187711D+00, 8.00905333D+00, 8.32135164D+00, 8.60235517D+00, FO   
     5  8.84375037D+00, 9.09391795D+00, 9.50278685D+00, 9.87035454D+00, FO   
     6  1.01817416D+01, 1.04589796D+01, 1.06916084D+01, 1.08710212D+01, FO   
     7  1.10032266D+01, 1.10534792D+01, 1.10902126D+01, 1.11069484D+01, FO   
     8  1.11182022D+01,     61*0.0D+00/                                 FO   
      DATA TK_NaO/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750299994957, NaO  
     2  0.828999943334, 0.923000013288, 1.029600211601, 1.172100019245, NaO  
     3  1.308000000556, 1.459900024047, 1.617299957586, 1.740399788690, NaO  
     4  1.874700000295, 1.986399856762, 2.098900362982, 2.228599645597, NaO  
     5  2.442899648528, 2.583900016819, 2.745799897810, 2.973899986284, NaO  
     6  3.190699580254, 3.355999936114, 3.461800074354, 3.566399945272, NaO  
     7  3.663999936752, 3.706299951147, 3.749899984958, 3.827899875959, NaO  
     8  3.855799960587, 3.885900259225, 3.921400132264, 3.953699909171, NaO  
     9  3.982999619658, 4.000000000000,     56*0.0D+00/                 NaO  
      DATA  K_NaO/                                                      071215
     1 -4.39058499D-05, 1.44856088D-01, 2.92974227D-01, 7.43570496D-01, NaO  
     2  1.76352189D+00, 2.78808727D+00, 3.74383807D+00, 4.75781649D+00, NaO  
     3  5.51535103D+00, 6.18808155D+00, 6.74812275D+00, 7.11850737D+00, NaO  
     4  7.47708886D+00, 7.75078615D+00, 8.01017697D+00, 8.29121722D+00, NaO  
     5  8.70729347D+00, 8.94201526D+00, 9.17020531D+00, 9.41980339D+00, NaO  
     6  9.59418421D+00, 9.70070872D+00, 9.76211120D+00, 9.82110102D+00, NaO  
     7  9.87999470D+00, 9.91001497D+00, 9.94781440D+00, 1.00525643D+01, NaO  
     8  1.01093873D+01, 1.01860248D+01, 1.02980265D+01, 1.04186359D+01, NaO  
     9  1.05401839D+01, 1.06145339D+01,     56*0.0D+00/                 NaO  
      DATA TK_MgO/                                                      071215
     1  0.699999789529, 0.710100047074, 0.721699837543, 0.755900140090, MgO  
     2  0.791799993159, 0.843100073575, 0.949499806081, 1.068500149111, MgO  
     3  1.147000031100, 1.224599899212, 1.368399921922, 1.528600195166, MgO  
     4  1.701199842554, 1.845100025319, 1.988099817088, 2.166999997362, MgO  
     5  2.326799857876, 2.482099579955, 2.618699844251, 2.732799618700, MgO  
     6  2.861399677581, 2.954799921121, 3.056600351910, 3.197999751337, MgO  
     7  3.273000245583, 3.353900094674, 3.465000143587, 3.597500325458, MgO  
     8  3.667500008198, 3.737799719293, 3.806500330192, 3.868699863258, MgO  
     9  3.947499771933, 3.979799566676, 3.990299783510, 4.000000000000, MgO  
     A      54*0.0D+00/                                                 MgO  
      DATA  K_MgO/                                                      071215
     1 -6.31769194D-05, 1.65787216D-01, 3.51976815D-01, 8.75517016D-01, MgO  
     2  1.38664088D+00, 2.05478384D+00, 3.23893078D+00, 4.30195170D+00, MgO  
     3  4.88369110D+00, 5.38329831D+00, 6.15112918D+00, 6.82563648D+00, MgO  
     4  7.40611671D+00, 7.81412659D+00, 8.17672575D+00, 8.59182190D+00, MgO  
     5  8.93409824D+00, 9.23891029D+00, 9.47895938D+00, 9.65315110D+00, MgO  
     6  9.80660238D+00, 9.87755673D+00, 9.91079271D+00, 9.89781925D+00, MgO  
     7  9.87675455D+00, 9.85179798D+00, 9.82167802D+00, 9.79961446D+00, MgO  
     8  9.79662990D+00, 9.80240002D+00, 9.82031635D+00, 9.85239076D+00, MgO  
     9  9.92975476D+00, 9.97989017D+00, 9.99917506D+00, 1.00183965D+01, MgO  
     A      54*0.0D+00/                                                 MgO  
      DATA TK_AlO/                                                      071215
     1  0.699999789529, 0.710000049601, 0.721399830359, 0.755300124540, AlO  
     2  0.841600109767, 0.942399980401, 1.064000044996, 1.200899819545, AlO  
     3  1.346200010028, 1.509400213195, 1.595200024243, 1.683799926523, AlO  
     4  1.822200131004, 1.957099982272, 2.259800230343, 2.638400302425, AlO  
     5  2.781499762697, 3.006900153972, 3.136600271740, 3.271800331570, AlO  
     6  3.361999687127, 3.450499805687, 3.706999967084, 3.788799949762, AlO  
     7  3.870499904273, 3.949899826122, 3.980799570266, 3.990599790206, AlO  
     8  4.000000000000,     61*0.0D+00/                                 AlO  
      DATA  K_AlO/                                                      071215
     1 -1.11822434D-04, 1.65792986D-01, 3.50722534D-01, 8.75433475D-01, AlO  
     2  2.05672714D+00, 3.19900219D+00, 4.30784789D+00, 5.28842929D+00, AlO  
     3  6.10179867D+00, 6.82103531D+00, 7.14298964D+00, 7.44873255D+00, AlO  
     4  7.89145459D+00, 8.29691334D+00, 9.12756149D+00, 9.97911579D+00, AlO  
     5  1.02363312D+01, 1.05672902D+01, 1.07177181D+01, 1.08388582D+01, AlO  
     6  1.08955971D+01, 1.09316240D+01, 1.09550799D+01, 1.09512429D+01, AlO  
     7  1.09512363D+01, 1.09675654D+01, 1.09827525D+01, 1.09890087D+01, AlO  
     8  1.09957215D+01,     61*0.0D+00/                                 AlO  
      DATA TK_SiO/                                                      071215
     1  0.699999789529, 0.710100047074, 0.721599835148, 0.755700134907, SiO  
     2  0.791599988645, 0.842800080814, 0.947699850275, 1.064600058878, SiO  
     3  1.153300012850, 1.241499792069, 1.397600037167, 1.474799892614, SiO  
     4  1.551900135313, 1.706099957842, 1.822300128210, 1.935999953612, SiO  
     5  2.033099801951, 2.128700063746, 2.272000311545, 2.422000131012, SiO  
     6  2.674200177000, 2.850600339965, 2.988099723260, 3.115499774958, SiO  
     7  3.443099649964, 3.560900340553, 3.669700053106, 3.759100202850, SiO  
     8  3.844200239026, 3.937099811607, 3.975099904389, 3.989099756609, SiO  
     9  4.000000000000,     57*0.0D+00/                                 SiO  
      DATA  K_SiO/                                                      071215
     1 -3.18827157D-04, 1.68714833D-01, 3.56853661D-01, 8.88973192D-01, SiO  
     2  1.40996953D+00, 2.08961609D+00, 3.28082682D+00, 4.34926348D+00, SiO  
     3  5.01682978D+00, 5.58403027D+00, 6.41189371D+00, 6.76198541D+00, SiO  
     4  7.08484768D+00, 7.67494236D+00, 8.08593979D+00, 8.46880365D+00, SiO  
     5  8.78348807D+00, 9.08228140D+00, 9.50650216D+00, 9.91518628D+00, SiO  
     6  1.05113761D+01, 1.08558162D+01, 1.10826557D+01, 1.12629613D+01, SiO  
     7  1.16281105D+01, 1.17382342D+01, 1.18345874D+01, 1.19094500D+01, SiO  
     8  1.19720187D+01, 1.20152112D+01, 1.20203057D+01, 1.20200862D+01, SiO  
     9  1.20191609D+01,     57*0.0D+00/                                 SiO  
      DATA TK_PO/                                                       071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749699980301, PO   
     2  0.827399982056, 0.920000085623, 1.025300105412, 1.096799992167, PO   
     3  1.168900043070, 1.306400034605, 1.456599952324, 1.610800120116, PO   
     4  1.736499864585, 1.860800116289, 1.989999772746, 2.106999891710, PO   
     5  2.351700262750, 2.564900063791, 2.696899738527, 2.833300016701, PO   
     6  2.985799671564, 3.138200308408, 3.279199801316, 3.407899787323, PO   
     7  3.656099764099, 3.721300161009, 3.792600032564, 3.923900187807, PO   
     8  3.969200252915, 4.000000000000,     60*0.0D+00/                 PO   
      DATA  K_PO/                                                       071215
     1 -2.33563751D-05, 1.57868924D-01, 3.15983552D-01, 8.00810454D-01, PO   
     2  1.89880959D+00, 2.99979645D+00, 4.02947165D+00, 4.61902988D+00, PO   
     3  5.13945765D+00, 5.96529405D+00, 6.67708348D+00, 7.26113908D+00, PO   
     4  7.66030317D+00, 8.00846221D+00, 8.33710889D+00, 8.61442276D+00, PO   
     5  9.14858490D+00, 9.56850837D+00, 9.80461189D+00, 1.00255654D+01, PO   
     6  1.02423766D+01, 1.04277547D+01, 1.05744683D+01, 1.06919417D+01, PO   
     7  1.09057423D+01, 1.09677026D+01, 1.10406339D+01, 1.11842962D+01, PO   
     8  1.12332555D+01, 1.12653028D+01,     60*0.0D+00/                 PO   
      DATA TK_SO/                                                       071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749699980301, SO   
     2  0.827299984476, 0.919800081146, 1.025100100473, 1.097499974263, SO   
     3  1.170900047093, 1.310899976411, 1.461399991780, 1.614400030099, SO   
     4  1.728599988419, 1.841500112180, 1.965199920067, 2.083099997551, SO   
     5  2.185899883515, 2.298700163336, 2.394600003380, 2.489299764479, SO   
     6  2.668400046539, 2.835900075231, 2.999599990639, 3.257700170233, SO   
     7  3.462500089499, 3.680200283324, 3.757400162560, 3.837200083553, SO   
     8  3.937499782778, 3.975299890018, 4.000000000000,     59*0.0D+00/ SO   
      DATA  K_SO/                                                       071215
     1  1.58945603D-04, 1.62872394D-01, 3.25802115D-01, 8.25308098D-01, SO   
     2  1.95468194D+00, 3.08670158D+00, 4.14558279D+00, 4.75836417D+00, SO   
     3  5.30063668D+00, 6.15740494D+00, 6.88120694D+00, 7.46865719D+00, SO   
     4  7.83822740D+00, 8.16297505D+00, 8.48667390D+00, 8.77438479D+00, SO   
     5  9.01437874D+00, 9.26992919D+00, 9.48220078D+00, 9.68704050D+00, SO   
     6  1.00550181D+01, 1.03638991D+01, 1.06238646D+01, 1.09480455D+01, SO   
     7  1.11352619D+01, 1.12767628D+01, 1.13173820D+01, 1.13558404D+01, SO   
     8  1.13976915D+01, 1.14104979D+01, 1.14176897D+01,     59*0.0D+00/ SO   
      DATA TK_ClO/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751400023466, ClO  
     2  0.831799960396, 0.927499904784, 1.035800066512, 1.108600108571, ClO  
     3  1.181699871309, 1.322500105688, 1.484499909996, 1.650900074038, ClO  
     4  1.751600027497, 1.849299923982, 1.950999828335, 2.035399857769, ClO  
     5  2.279899751085, 2.620299881053, 2.803000267696, 3.144699998945, ClO  
     6  3.334400024524, 3.520600413121, 3.649499657085, 3.811800276964, ClO  
     7  3.920900121155, 3.968900246192, 4.000000000000,     63*0.0D+00/ ClO  
      DATA  K_ClO/                                                      071215
     1 -4.08848518D-05, 1.66862246D-01, 3.38881747D-01, 8.55729626D-01, ClO  
     2  2.02381947D+00, 3.18609621D+00, 4.26068515D+00, 4.86696324D+00, ClO  
     3  5.39858179D+00, 6.24720409D+00, 7.00862940D+00, 7.62459966D+00, ClO  
     4  7.94044226D+00, 8.21754749D+00, 8.48242062D+00, 8.68713484D+00, ClO  
     5  9.21148525D+00, 9.79247691D+00, 1.00490858D+01, 1.04414652D+01, ClO  
     6  1.06123820D+01, 1.07505660D+01, 1.08296634D+01, 1.09062868D+01, ClO  
     7  1.09388087D+01, 1.09478507D+01, 1.09522407D+01,     63*0.0D+00/ ClO  
      DATA TK_KO/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, KO   
     2  0.830299926010, 0.924899967475, 1.032200162698, 1.104000002529, KO   
     3  1.176199924098, 1.306100040989, 1.392499915309, 1.478999990961, KO   
     4  1.634799991627, 1.801700138261, 1.950999828335, 2.054200295045, KO   
     5  2.170900079261, 2.264900342032, 2.365999791031, 2.624499977326, KO   
     6  2.762800301027, 2.905899738979, 3.074999819734, 3.256200137495, KO   
     7  3.385900253881, 3.474399936349, 3.561700283058, 3.638500297268, KO   
     8  3.714600134611, 3.785799879033, 3.822399742991, 3.859899671523, KO   
     9  3.900099621896, 3.941699640976, 3.980999574756, 4.000000000000, KO   
     A      54*0.0D+00/                                                 KO   
      DATA  K_KO/                                                       071215
     1  1.74263223D-05, 1.59188628D-01, 3.21753089D-01, 8.15642718D-01, KO   
     2  1.93167587D+00, 3.04481070D+00, 4.07938809D+00, 4.66208099D+00, KO   
     3  5.17505449D+00, 5.94832864D+00, 6.37724241D+00, 6.75265929D+00, KO   
     4  7.32344770D+00, 7.82812242D+00, 8.21803207D+00, 8.46073674D+00, KO   
     5  8.70696898D+00, 8.88037924D+00, 9.04074226D+00, 9.34644959D+00, KO   
     6  9.46813645D+00, 9.57637539D+00, 9.68961614D+00, 9.80000031D+00, KO   
     7  9.87504530D+00, 9.92637504D+00, 9.98043076D+00, 1.00364742D+01, KO   
     8  1.01121158D+01, 1.02223161D+01, 1.03022958D+01, 1.04037122D+01, KO   
     9  1.05334823D+01, 1.06862014D+01, 1.08415481D+01, 1.09187141D+01, KO   
     A      54*0.0D+00/                                                 KO   
      DATA TK_CaO/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720799815989, 0.753900088257, CaO  
     2  0.838000102524, 0.934899939967, 1.054700071024, 1.190800024619, CaO  
     3  1.334400035806, 1.484899900717, 1.642799928169, 1.764400120072, CaO  
     4  1.895199976845, 2.013300307295, 2.137500263394, 2.293900043799, CaO  
     5  2.474699910718, 2.646899840426, 2.833800027956, 3.075599834760, CaO  
     6  3.151299633068, 3.226699771152, 3.348200348654, 3.401199651822, CaO  
     7  3.456899960449, 3.530499741163, 3.602300214247, 3.706799962531, CaO  
     8  3.805400306079, 3.861799705601, 3.914999972691, 3.967000203611, CaO  
     9  3.986899707217, 4.000000000000,     56*0.0D+00/                 CaO  
      DATA  K_CaO/                                                      071215
     1  5.18575316D-05, 1.62666323D-01, 3.37849991D-01, 8.46129965D-01, CaO  
     2  1.99179427D+00, 3.09134726D+00, 4.19128339D+00, 5.17481122D+00, CaO  
     3  5.98642118D+00, 6.65789048D+00, 7.22403141D+00, 7.59303520D+00, CaO  
     4  7.94620943D+00, 8.23863650D+00, 8.52696583D+00, 8.86657221D+00, CaO  
     5  9.22383636D+00, 9.52133871D+00, 9.79269460D+00, 1.00711463D+01, CaO  
     6  1.01434066D+01, 1.02072381D+01, 1.02834158D+01, 1.03009557D+01, CaO  
     7  1.03063334D+01, 1.02936595D+01, 1.02675210D+01, 1.02386080D+01, CaO  
     8  1.02714105D+01, 1.03380874D+01, 1.04427146D+01, 1.05849989D+01, CaO  
     9  1.06486019D+01, 1.06927880D+01,     56*0.0D+00/                 CaO  
      DATA TK_ScO/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718499834746, 0.748699957364, ScO  
     2  0.824800044980, 0.915699989358, 1.019799979305, 1.156399933583, ScO  
     3  1.284600095453, 1.429599826155, 1.583700009793, 1.701699854318, ScO  
     4  1.813200022834, 1.917700019433, 2.020300442842, 2.130000091012, ScO  
     5  2.235499680704, 2.521100375053, 2.652899686900, 2.786299884310, ScO  
     6  2.935799885162, 3.094000235546, 3.249399988520, 3.402499678114, ScO  
     7  3.508600190569, 3.668200022487, 3.755100108051, 3.851900235550, ScO  
     8  3.949599819348, 3.980399561285, 4.000000000000,     59*0.0D+00/ ScO  
      DATA  K_ScO/                                                      071215
     1 -1.88867057D-04, 1.59042413D-01, 3.16963504D-01, 8.08955191D-01, ScO  
     2  1.92014847D+00, 3.03966458D+00, 4.09584750D+00, 5.19746358D+00, ScO  
     3  6.00988864D+00, 6.74013198D+00, 7.35932867D+00, 7.75842604D+00, ScO  
     4  8.09777249D+00, 8.39651951D+00, 8.67984838D+00, 8.97511827D+00, ScO  
     5  9.25074651D+00, 9.93366293D+00, 1.02062814D+01, 1.04498453D+01, ScO  
     6  1.06840440D+01, 1.08912262D+01, 1.10605588D+01, 1.12023891D+01, ScO  
     7  1.12938112D+01, 1.14467870D+01, 1.15513109D+01, 1.16886794D+01, ScO  
     8  1.18444807D+01, 1.18961397D+01, 1.19295102D+01,     59*0.0D+00/ ScO  
      DATA TK_TiO/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718299839802, 0.747999941308, TiO  
     2  0.823100086123, 0.912499917719, 1.014700101148, 1.083099996392, TiO  
     3  1.152200040978, 1.284000111983, 1.427199889698, 1.570900150367, TiO  
     4  1.721799830690, 1.842800080814, 1.963999951310, 2.091400174372, TiO  
     5  2.201499829619, 2.309400414643, 2.546200096387, 2.659299843744, TiO  
     6  2.790099980318, 2.932700102711, 3.078999919908, 3.255300117852, TiO  
     7  3.449899791515, 3.622399952994, 3.730199562285, 3.837800096622, TiO  
     8  3.937299797192, 3.975499875647, 4.000000000000,     59*0.0D+00/ TiO  
      DATA  K_TiO/                                                      071215
     1  1.90603056D-04, 1.61496506D-01, 3.18072439D-01, 8.08695766D-01, TiO  
     2  1.92237341D+00, 3.04328075D+00, 4.10208526D+00, 4.70209837D+00, TiO  
     3  5.23472582D+00, 6.08332664D+00, 6.81207918D+00, 7.39315843D+00, TiO  
     4  7.88583088D+00, 8.22237727D+00, 8.52911048D+00, 8.83500317D+00, TiO  
     5  9.09284966D+00, 9.34110591D+00, 9.85748455D+00, 1.00805386D+01, TiO  
     6  1.03131417D+01, 1.05340784D+01, 1.07276904D+01, 1.09267518D+01, TiO  
     7  1.11317760D+01, 1.13364296D+01, 1.14872237D+01, 1.16567819D+01, TiO  
     8  1.18275724D+01, 1.18958011D+01, 1.19401091D+01,     59*0.0D+00/ TiO  
      DATA TK_VO/                                                       071215
     1  0.699999789529, 0.709100026195, 0.717899849912, 0.747099920664, VO   
     2  0.820700144206, 0.908599894023, 1.009800209165, 1.142199925287, VO   
     3  1.265200094643, 1.403000021284, 1.550500169418, 1.662999902651, VO   
     4  1.772800038819, 1.884999997259, 1.998099956822, 2.106499927349, VO   
     5  2.219700243473, 2.323099777758, 2.425000188143, 2.568799786219, VO   
     6  2.707699998941, 2.831799982933, 3.031599783068, 3.269400446144, VO   
     7  3.442699641637, 3.587400098819, 3.726199821470, 3.808800380610, VO   
     8  3.884800233116, 3.955599950659, 3.983199624148, 4.000000000000, VO   
     9      58*0.0D+00/                                                 VO   
      DATA  K_VO/                                                       071215
     1 -9.88797237D-05, 1.52185835D-01, 2.96782465D-01, 7.57828683D-01, VO   
     2  1.80351404D+00, 2.86330613D+00, 3.87491906D+00, 4.93452112D+00, VO   
     3  5.71482109D+00, 6.41656472D+00, 7.02199681D+00, 7.41386246D+00, VO   
     4  7.75928674D+00, 8.09066917D+00, 8.41434541D+00, 8.72091007D+00, VO   
     5  9.03920488D+00, 9.32699937D+00, 9.60447408D+00, 9.97665394D+00, VO   
     6  1.03045350D+01, 1.05674594D+01, 1.09378683D+01, 1.13111596D+01, VO   
     7  1.15425559D+01, 1.17154243D+01, 1.18762024D+01, 1.19771055D+01, VO   
     8  1.20772862D+01, 1.21781407D+01, 1.22193573D+01, 1.22449171D+01, VO   
     9      58*0.0D+00/                                                 VO   
      DATA TK_CrO/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718699829691, 0.749099966539, CrO  
     2  0.825800020779, 0.917500029655, 1.021900021448, 1.091600125173, CrO  
     3  1.161599877761, 1.294800054256, 1.441300012201, 1.593300076862, CrO  
     4  1.725899925791, 1.871499921131, 2.005900135478, 2.144999971540, CrO  
     5  2.264600335479, 2.376300038281, 2.646699854916, 2.778799814723, CrO  
     6  2.913199909414, 3.049500170648, 3.168000026002, 3.295200085770, CrO  
     7  3.430900229510, 3.517500398356, 3.600400351751, 3.711600069900, CrO  
     8  3.846800297517, 3.936599847644, 3.975499875647, 3.989199758854, CrO  
     9  4.000000000000,     57*0.0D+00/                                 CrO  
      DATA  K_CrO/                                                      071215
     1  1.76241156D-05, 1.52044408D-01, 3.02780076D-01, 7.70465412D-01, CrO  
     2  1.82801506D+00, 2.89438051D+00, 3.89524694D+00, 4.46048322D+00, CrO  
     3  4.95942881D+00, 5.75419748D+00, 6.44863676D+00, 7.02673204D+00, CrO  
     4  7.44796646D+00, 7.85066996D+00, 8.18610257D+00, 8.50877070D+00, CrO  
     5  8.77046454D+00, 9.00221713D+00, 9.50318573D+00, 9.71086197D+00, CrO  
     6  9.89562510D+00, 1.00574463D+01, 1.01802023D+01, 1.02993953D+01, CrO  
     7  1.04253328D+01, 1.05132257D+01, 1.06068701D+01, 1.07515218D+01, CrO  
     8  1.09735888D+01, 1.11637556D+01, 1.12579573D+01, 1.12927182D+01, CrO  
     9  1.13206725D+01,     57*0.0D+00/                                 CrO  
      DATA TK_MnO/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749699980301, MnO  
     2  0.827299984476, 0.920100083212, 1.025500110351, 1.096699994725, MnO  
     3  1.168600036276, 1.305600051629, 1.455299924070, 1.608300097823, MnO  
     4  1.728599988419, 1.847799960174, 1.968799826337, 2.087400085803, MnO  
     5  2.212100095801, 2.331399960353, 2.617099807420, 2.754300093053, MnO  
     6  2.889300332598, 3.235699670148, 3.358099777555, 3.481999579267, MnO  
     7  3.602600192536, 3.713500110884, 3.794500072398, 3.876500040005, MnO  
     8  3.951999872051, 3.981699590472, 4.000000000000,     59*0.0D+00/ MnO  
      DATA  K_MnO/                                                      071215
     1 -4.93740226D-05, 1.55338593D-01, 3.10947069D-01, 7.88102484D-01, MnO  
     2  1.86755216D+00, 2.95382233D+00, 3.96858986D+00, 4.54698643D+00, MnO  
     3  5.05874194D+00, 5.87115432D+00, 6.57290509D+00, 7.14751656D+00, MnO  
     4  7.52806758D+00, 7.86227031D+00, 8.17142944D+00, 8.45394740D+00, MnO  
     5  8.73422524D+00, 8.98770214D+00, 9.52737336D+00, 9.74581303D+00, MnO  
     6  9.93338794D+00, 1.03073349D+01, 1.04107547D+01, 1.05048958D+01, MnO  
     7  1.05927718D+01, 1.06857083D+01, 1.07773109D+01, 1.09058917D+01, MnO  
     8  1.10642415D+01, 1.11372862D+01, 1.11851447D+01,     59*0.0D+00/ MnO  
      DATA TK_FeO/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718399837274, 0.748299948189, FeO  
     2  0.823900066761, 0.913399937868, 1.016500058145, 1.086100070177, FeO  
     3  1.156399933583, 1.290499957883, 1.432799879129, 1.580299924277, FeO  
     4  1.698099863209, 1.822200131004, 1.934899930040, 2.059200410192, FeO  
     5  2.185999876344, 2.320699725790, 2.421000111969, 2.515200356440, FeO  
     6  2.763000305858, 2.888600315483, 3.020200438191, 3.206599946114, FeO  
     7  3.400399635643, 3.558900382279, 3.687699760966, 3.879400105609, FeO  
     8  3.952599885152, 3.981999597207, 4.000000000000,     59*0.0D+00/ FeO  
      DATA  K_FeO/                                                      071215
     1 -1.68298602D-04, 1.52062173D-01, 3.01451455D-01, 7.67557253D-01, FeO  
     2  1.82546961D+00, 2.88504071D+00, 3.89434156D+00, 4.47114567D+00, FeO  
     3  4.98282420D+00, 5.79824055D+00, 6.48622026D+00, 7.06052847D+00, FeO  
     4  7.44644150D+00, 7.80395368D+00, 8.09842241D+00, 8.40126227D+00, FeO  
     5  8.69543380D+00, 8.99900961D+00, 9.22164369D+00, 9.42789333D+00, FeO  
     6  9.94067216D+00, 1.01717596D+01, 1.03864423D+01, 1.06401957D+01, FeO  
     7  1.08510939D+01, 1.10047337D+01, 1.11359610D+01, 1.13853170D+01, FeO  
     8  1.15104428D+01, 1.15666007D+01, 1.16027418D+01,     59*0.0D+00/ FeO  
      DATA TK_NiO/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718799827163, 0.749299971126, NiO  
     2  0.826300008678, 0.918300047565, 1.023000048613, 1.092900091921, NiO  
     3  1.163199913993, 1.297300110287, 1.446399893013, 1.598299938393, NiO  
     4  1.727199955945, 1.851899950677, 1.964899927878, 2.072099735881, NiO  
     5  2.227499722815, 2.379500114644, 2.571599736754, 2.816199972378, NiO  
     6  3.023800183650, 3.326899847553, 3.538299918669, 3.631300148838, NiO  
     7  3.726299814540, 3.811900269745, 3.893700083893, 3.958000003064, NiO  
     8  3.984099644354, 4.000000000000,     60*0.0D+00/                 NiO  
      DATA  K_NiO/                                                      071215
     1  8.77248413D-05, 1.75418024D-01, 3.51041092D-01, 8.91548809D-01, NiO  
     2  2.11194108D+00, 3.33726323D+00, 4.48151231D+00, 5.12469629D+00, NiO  
     3  5.69086205D+00, 6.58857337D+00, 7.37181934D+00, 8.00379814D+00, NiO  
     4  8.44852304D+00, 8.82564352D+00, 9.13902843D+00, 9.42091662D+00, NiO  
     5  9.80997105D+00, 1.01654272D+01, 1.05684168D+01, 1.10052474D+01, NiO  
     6  1.13232730D+01, 1.17179463D+01, 1.19439623D+01, 1.20288099D+01, NiO  
     7  1.21055371D+01, 1.21677450D+01, 1.22259259D+01, 1.22763317D+01, NiO  
     8  1.22993533D+01, 1.23143806D+01,     60*0.0D+00/                 NiO  
      DATA TK_CuO/                                                      071215
     1  0.699999789529, 0.709700041799, 0.719999796830, 0.752000039016, CuO  
     2  0.833199992489, 0.929899846916, 1.039199975669, 1.112700070761, CuO  
     3  1.186499971251, 1.328299973624, 1.491899822271, 1.660799850769, CuO  
     4  1.759600229072, 1.863500048124, 1.964099948707, 2.059500417101, CuO  
     5  2.365499778918, 2.539899954924, 2.710500067329, 2.964300139789, CuO  
     6  3.220400228265, 3.306900340919, 3.393800080248, 3.612499717648, CuO  
     7  3.696399733213, 3.769600456272, 3.826499842113, 3.880400128676, CuO  
     8  3.952899891703, 3.982199601697, 4.000000000000,     59*0.0D+00/ CuO  
      DATA  K_CuO/                                                      071215
     1  4.50561297D-05, 1.54278670D-01, 3.14686603D-01, 7.91719897D-01, CuO  
     2  1.87037119D+00, 2.94449943D+00, 3.93737784D+00, 4.49862197D+00, CuO  
     3  4.99153107D+00, 5.77905923D+00, 6.49190543D+00, 7.07598226D+00, CuO  
     4  7.36804516D+00, 7.64689290D+00, 7.89511730D+00, 8.11359241D+00, CuO  
     5  8.71099147D+00, 8.98272479D+00, 9.20421981D+00, 9.46592129D+00, CuO  
     6  9.66761217D+00, 9.72497422D+00, 9.77878552D+00, 9.91497332D+00, CuO  
     7  9.97618440D+00, 1.00357890D+01, 1.00861544D+01, 1.01376804D+01, CuO  
     8  1.02163258D+01, 1.02530561D+01, 1.02773661D+01,     59*0.0D+00/ CuO  
      DATA TK_GaO/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718799827163, 0.749299971126, GaO  
     2  0.826500003838, 0.918400049804, 1.023000048613, 1.094000063786, GaO  
     3  1.165899975135, 1.303100104831, 1.452399861042, 1.604700008252, GaO  
     4  1.714799924483, 1.833399993031, 1.934899930040, 2.036099874757, GaO  
     5  2.281499774394, 2.535099838180, 2.703999907050, 2.817399886685, GaO  
     6  2.939599618490, 3.090000139165, 3.228499640548, 3.438999656211, GaO  
     7  3.616399812269, 3.738999744084, 3.800700203052, 3.859299713825, GaO  
     8  3.944499704197, 3.978499660086, 3.989999776814, 4.000000000000, GaO  
     9      58*0.0D+00/                                                 GaO  
      DATA  K_GaO/                                                      071215
     1  1.00343974D-04, 1.54352069D-01, 3.08897439D-01, 7.84782707D-01, GaO  
     2  1.86344004D+00, 2.94523500D+00, 3.95946878D+00, 4.54076725D+00, GaO  
     3  5.05650446D+00, 5.87595997D+00, 6.58050324D+00, 7.15623792D+00, GaO  
     4  7.50851445D+00, 7.84581580D+00, 8.10972871D+00, 8.35658415D+00, GaO  
     5  8.90764581D+00, 9.42787531D+00, 9.75692139D+00, 9.96713598D+00, GaO  
     6  1.01778654D+01, 1.04082901D+01, 1.05902448D+01, 1.08155256D+01, GaO  
     7  1.09643368D+01, 1.10489805D+01, 1.10880172D+01, 1.11258778D+01, GaO  
     8  1.11937148D+01, 1.12304459D+01, 1.12447572D+01, 1.12580938D+01, GaO  
     9      58*0.0D+00/                                                 GaO  
      DATA TK_GeO/                                                      071215
     1  0.699999789529, 0.709100026195, 0.718199842329, 0.747799936720, GeO  
     2  0.822400103064, 0.910799879661, 1.012100163264, 1.145900006851, GeO  
     3  1.275200106846, 1.421100051204, 1.569000149925, 1.678699867506, GeO  
     4  1.796300102192, 1.906799911388, 2.008100185995, 2.282999804796, GeO  
     5  2.369399873405, 2.464400125459, 2.635600236082, 2.753100063508, GeO  
     6  2.888200305702, 3.071999744604, 3.267800407725, 3.425300199076, GeO  
     7  3.576399839441, 3.697999766323, 3.821799728486, 3.868299854118, GeO  
     8  3.917900047203, 3.963600127415, 4.000000000000,     59*0.0D+00/ GeO  
      DATA  K_GeO/                                                      071215
     1 -4.83004633D-05, 1.52144416D-01, 3.01528983D-01, 7.68130020D-01, GeO  
     2  1.82481118D+00, 2.88585518D+00, 3.89345886D+00, 4.95777762D+00, GeO  
     3  5.76803860D+00, 6.49402706D+00, 7.08296672D+00, 7.45090225D+00, GeO  
     4  7.79873628D+00, 8.09422895D+00, 8.34657425D+00, 8.98768439D+00, GeO  
     5  9.18896952D+00, 9.41510811D+00, 9.83185237D+00, 1.01154219D+01, GeO  
     6  1.04258169D+01, 1.08015357D+01, 1.11310398D+01, 1.13492245D+01, GeO  
     7  1.15296218D+01, 1.16576442D+01, 1.17651681D+01, 1.17954276D+01, GeO  
     8  1.18194642D+01, 1.18342397D+01, 1.18426347D+01,     59*0.0D+00/ GeO  
      DATA TK_AsO/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, AsO  
     2  0.830299926010, 0.925099962653, 1.032500154683, 1.106300055550, AsO  
     3  1.180899854652, 1.323300087473, 1.478399976911, 1.637299929369, AsO  
     4  1.758600203875, 1.878500094304, 2.009300213550, 2.134400192142, AsO  
     5  2.243199857785, 2.352100233011, 2.567699864509, 2.668500048749, AsO  
     6  2.779399769708, 2.907399772710, 3.039299964956, 3.329499903927, AsO  
     7  3.455999938685, 3.584700035562, 3.663599928587, 3.742499820353, AsO  
     8  3.824699798596, 3.903199691514, 3.960900066905, 3.985099666805, AsO  
     9  4.000000000000,     57*0.0D+00/                                 AsO  
      DATA  K_AsO/                                                      071215
     1 -1.11834435D-04, 1.65260024D-01, 3.34146236D-01, 8.47166742D-01, AsO  
     2  2.00593217D+00, 3.16298744D+00, 4.23561164D+00, 4.85406349D+00, AsO  
     3  5.39863202D+00, 6.25766802D+00, 6.98998944D+00, 7.58622959D+00, AsO  
     4  7.96865261D+00, 8.30434625D+00, 8.63774299D+00, 8.93408945D+00, AsO  
     5  9.17761638D+00, 9.40886868D+00, 9.82211874D+00, 9.98840394D+00, AsO  
     6  1.01474963D+01, 1.03002616D+01, 1.04279217D+01, 1.06387824D+01, AsO  
     7  1.07159802D+01, 1.07979858D+01, 1.08549556D+01, 1.09190029D+01, AsO  
     8  1.09919619D+01, 1.10626591D+01, 1.11108615D+01, 1.11292988D+01, AsO  
     9  1.11400167D+01,     57*0.0D+00/                                 AsO  
      DATA TK_SeO/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, SeO  
     2  0.830299926010, 0.924699972297, 1.031700176058, 1.106600062466, SeO  
     3  1.182899896294, 1.330799953260, 1.485399889119, 1.639699869601, SeO  
     4  1.766300068652, 1.898200042463, 2.030299734000, 2.159099837437, SeO  
     5  2.282099786555, 2.395499937225, 2.616399791306, 2.790899998435, SeO  
     6  2.914399938752, 3.047900138790, 3.176700218325, 3.336000064435, SeO  
     7  3.503700076333, 3.655599752432, 3.806100321424, 3.918100052341, SeO  
     8  3.967500214817, 4.000000000000,     60*0.0D+00/                 SeO  
      DATA  K_SeO/                                                      071215
     1  5.51697394D-06, 1.63437982D-01, 3.30345935D-01, 8.37373687D-01, SeO  
     2  1.98272029D+00, 3.12222781D+00, 4.18024709D+00, 4.80149215D+00, SeO  
     3  5.35154353D+00, 6.22798378D+00, 6.94308974D+00, 7.51217865D+00, SeO  
     4  7.90122296D+00, 8.25104925D+00, 8.55859078D+00, 8.82818391D+00, SeO  
     5  9.06536170D+00, 9.27026825D+00, 9.63356991D+00, 9.88433431D+00, SeO  
     6  1.00439134D+01, 1.02036253D+01, 1.03474171D+01, 1.05103137D+01, SeO  
     7  1.06606326D+01, 1.07768975D+01, 1.08714106D+01, 1.09248464D+01, SeO  
     8  1.09432205D+01, 1.09537971D+01,     60*0.0D+00/                 SeO  
      DATA TK_BrO/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720999820779, 0.754400101215, BrO  
     2  0.839400134618, 0.937099982828, 1.058399988246, 1.197599858818, BrO  
     3  1.344100062252, 1.497099931342, 1.658199879798, 1.772000018733, BrO  
     4  1.895399981219, 2.007600174514, 2.120999902247, 2.281399772367, BrO  
     5  2.507800189242, 2.594400254029, 2.690299591235, 2.869599882944, BrO  
     6  3.058000388248, 3.184599955466, 3.313100188672, 3.490899777691, BrO  
     7  3.697599758045, 3.833900011675, 3.953799911355, 4.000000000000, BrO  
     8      62*0.0D+00/                                                 BrO  
      DATA  K_BrO/                                                      071215
     1 -8.51813461D-05, 1.64382431D-01, 3.44772011D-01, 8.63100312D-01, BrO  
     2  2.03177204D+00, 3.14858820D+00, 4.26782544D+00, 5.27436998D+00, BrO  
     3  6.09866696D+00, 6.77576664D+00, 7.34736992D+00, 7.69062777D+00, BrO  
     4  8.02445253D+00, 8.30388938D+00, 8.56946414D+00, 8.92075915D+00, BrO  
     5  9.35802047D+00, 9.50032623D+00, 9.63874632D+00, 9.84523018D+00, BrO  
     6  1.00060979D+01, 1.00960604D+01, 1.01811471D+01, 1.02949127D+01, BrO  
     7  1.04193160D+01, 1.04928123D+01, 1.05523682D+01, 1.05750080D+01, BrO  
     8      62*0.0D+00/                                                 BrO  
      DATA TK_RbO/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, RbO  
     2  0.830299926010, 0.925199960242, 1.032700149339, 1.105400034803, RbO  
     3  1.178499870723, 1.318100123735, 1.472799845782, 1.631100083770, RbO  
     4  1.735599886227, 1.837500089531, 1.937699990040, 2.033899821366, RbO  
     5  2.173600146707, 2.300100198039, 2.399799621149, 2.502800056197, RbO  
     6  2.632000150784, 2.772300302387, 2.876500030920, 2.981799581660, RbO  
     7  3.135100237364, 3.279999743991, 3.412299890313, 3.491699796362, RbO  
     8  3.569399729665, 3.638900305514, 3.706799962531, 3.780699758793, RbO  
     9  3.815000045942, 3.850800313104, 3.895099980728, 3.947999783222, RbO  
     A  4.000000000000,     53*0.0D+00/                                 RbO  
      DATA  K_RbO/                                                      071215
     1 -7.01058611D-05, 1.58033504D-01, 3.19509883D-01, 8.10108982D-01, RbO  
     2  1.91880399D+00, 3.02797807D+00, 4.05725541D+00, 4.64268742D+00, RbO  
     3  5.15732625D+00, 5.97405152D+00, 6.68628051D+00, 7.26760213D+00, RbO  
     4  7.59383548D+00, 7.88073487D+00, 8.14071366D+00, 8.37429398D+00, RbO  
     5  8.68772174D+00, 8.93969077D+00, 9.11049570D+00, 9.25751561D+00, RbO  
     6  9.40123610D+00, 9.51716138D+00, 9.58557551D+00, 9.64572042D+00, RbO  
     7  9.72478316D+00, 9.79560414D+00, 9.86026665D+00, 9.90129856D+00, RbO  
     8  9.94656600D+00, 9.99711408D+00, 1.00667812D+01, 1.01875853D+01, RbO  
     9  1.02669137D+01, 1.03674032D+01, 1.05147018D+01, 1.07148849D+01, RbO  
     A  1.09243834D+01,     53*0.0D+00/                                 RbO  
      DATA TK_SrO/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720799815989, 0.753900088257, SrO  
     2  0.838200107109, 0.935199945812, 1.055300057600, 1.192199990484, SrO  
     3  1.336500083959, 1.487499840408, 1.645199984775, 1.759400224032, SrO  
     4  1.883800029456, 1.993899861375, 2.105100027135, 2.269800449065, SrO  
     5  2.462800089372, 2.622499931482, 2.854600043180, 3.037599924799, SrO  
     6  3.141200260115, 3.236699694435, 3.372099942346, 3.446699724903, SrO  
     7  3.533799816262, 3.618399860792, 3.709100014897, 3.795900101749, SrO  
     8  3.861299694176, 3.916700016370, 3.967700219299, 3.987099711707, SrO  
     9  4.000000000000,     57*0.0D+00/                                 SrO  
      DATA  K_SrO/                                                      071215
     1  1.86739090D-04, 1.62552435D-01, 3.37468727D-01, 8.44976833D-01, SrO  
     2  1.99141927D+00, 3.08995593D+00, 4.19015502D+00, 5.17638009D+00, SrO  
     3  5.98855359D+00, 6.65910025D+00, 7.22210869D+00, 7.56874123D+00, SrO  
     4  7.90653943D+00, 8.18145328D+00, 8.44256128D+00, 8.80368924D+00, SrO  
     5  9.18448592D+00, 9.45784205D+00, 9.78346524D+00, 9.98667944D+00, SrO  
     6  1.00840940D+01, 1.01624284D+01, 1.02442729D+01, 1.02650498D+01, SrO  
     7  1.02622221D+01, 1.02418342D+01, 1.02313557D+01, 1.02682590D+01, SrO  
     8  1.03436421D+01, 1.04443548D+01, 1.05652355D+01, 1.06171981D+01, SrO  
     9  1.06532825D+01,     57*0.0D+00/                                 SrO  
      DATA TK_YO/                                                       071215
     1  0.699999789529, 0.709200028796, 0.718299839802, 0.748199945895, YO   
     2  0.823600074022, 0.913199933390, 1.015700077257, 1.085000043123, YO   
     3  1.154899971938, 1.288399990758, 1.431699854157, 1.577099991186, YO   
     4  1.695799922389, 1.814200047783, 1.946899876329, 2.082099977028, YO   
     5  2.205599938330, 2.319999710632, 2.433900011170, 2.535699852773, YO   
     6  2.666500004549, 2.788099929915, 2.942899654252, 3.090800158441, YO   
     7  3.246099913414, 3.386600271607, 3.480899554872, 3.583300002761, YO   
     8  3.672000106835, 3.730599570549, 3.786399893179, 3.909899841977, YO   
     9  3.964600149825, 4.000000000000,     56*0.0D+00/                 YO   
      DATA  K_YO/                                                       071215
     1 -1.82747836D-04, 1.58405819D-01, 3.12343292D-01, 7.97844304D-01, YO   
     2  1.89635747D+00, 2.99965201D+00, 4.04235730D+00, 4.63884547D+00, YO   
     3  5.16693906D+00, 6.00859669D+00, 6.72450348D+00, 7.30775014D+00, YO   
     4  7.70716181D+00, 8.05686126D+00, 8.40890971D+00, 8.73904785D+00, YO   
     5  9.02495422D+00, 9.28240562D+00, 9.53459320D+00, 9.75556355D+00, YO   
     6  1.00275725D+01, 1.02612738D+01, 1.05244209D+01, 1.07381320D+01, YO   
     7  1.09265468D+01, 1.10723771D+01, 1.11643970D+01, 1.12709111D+01, YO   
     8  1.13804494D+01, 1.14655490D+01, 1.15560783D+01, 1.17786117D+01, YO   
     9  1.18803190D+01, 1.19453105D+01,     56*0.0D+00/                 YO   
      DATA TK_ZrO/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749499975714, ZrO  
     2  0.826899994157, 0.918700056520, 1.023900070838, 1.096000012630, ZrO  
     3  1.169000045334, 1.308099998428, 1.454699911030, 1.607800085383, ZrO  
     4  1.729200002336, 1.855200026377, 1.977299949146, 2.105599991497, ZrO  
     5  2.239499777710, 2.346900319413, 2.440699600663, 2.525600040883, ZrO  
     6  2.622999942943, 2.706399966655, 2.789299960319, 2.870899912063, ZrO  
     7  3.020000452333, 3.170400083274, 3.323199767329, 3.485899665759, ZrO  
     8  3.639400315821, 3.762100274950, 3.859699685624, 3.933700056659, ZrO  
     9  4.000000000000,     57*0.0D+00/                                 ZrO  
      DATA  K_ZrO/                                                      071215
     1  1.71619479D-04, 1.71536553D-01, 3.41345001D-01, 8.65700662D-01, ZrO  
     2  2.05172240D+00, 3.23434940D+00, 4.34716435D+00, 4.98829950D+00, ZrO  
     3  5.55442984D+00, 6.44677447D+00, 7.18558111D+00, 7.79998948D+00, ZrO  
     4  8.20650111D+00, 8.57614794D+00, 8.89935984D+00, 9.21387511D+00, ZrO  
     5  9.52401556D+00, 9.76432992D+00, 9.96913116D+00, 1.01482929D+01, ZrO  
     6  1.03405549D+01, 1.04876199D+01, 1.06134427D+01, 1.07168755D+01, ZrO  
     7  1.08641019D+01, 1.09847015D+01, 1.11085619D+01, 1.12605962D+01, ZrO  
     8  1.14269523D+01, 1.15759153D+01, 1.17032411D+01, 1.18031886D+01, ZrO  
     9  1.18934504D+01,     57*0.0D+00/                                 ZrO  
      DATA TK_NbO/                                                      071215
     1  0.699999789529, 0.708900020993, 0.717399862551, 0.745699888552, NbO  
     2  0.817000097626, 0.901800050770, 1.000100002134, 1.126800044774, NbO  
     3  1.245499899711, 1.379000077143, 1.516700052518, 1.638599896994, NbO  
     4  1.748399954100, 1.851099932325, 1.956799974701, 2.066799921758, NbO  
     5  2.181100227702, 2.412999923170, 2.744799876530, 2.888500313038, NbO  
     6  3.046500110914, 3.208699993175, 3.348400353152, 3.614699771024, NbO  
     7  3.713200104413, 3.814900053161, 4.000000000000,     63*0.0D+00/ NbO  
      DATA  K_NbO/                                                      071215
     1  2.11947035D-04, 1.42480604D-01, 2.75994333D-01, 7.03704551D-01, NbO  
     2  1.67700728D+00, 2.66549722D+00, 3.62202898D+00, 4.61806308D+00, NbO  
     3  5.36511291D+00, 6.04452798D+00, 6.61395576D+00, 7.03950294D+00, NbO  
     4  7.38208097D+00, 7.68310302D+00, 7.98451991D+00, 8.29603307D+00, NbO  
     5  8.62022423D+00, 9.27728538D+00, 1.01714862D+01, 1.05152986D+01, NbO  
     6  1.08508286D+01, 1.11481432D+01, 1.13708041D+01, 1.17557907D+01, NbO  
     7  1.19042782D+01, 1.20672152D+01, 1.23831207D+01,     63*0.0D+00/ NbO  
      DATA TK_InO/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749799982595, InO  
     2  0.827599977216, 0.920300078390, 1.025700115290, 1.096799992167, InO  
     3  1.168400031747, 1.304800068654, 1.453699889296, 1.607300072942, InO  
     4  1.724399890998, 1.852599966735, 1.965199920067, 2.080399942137, InO  
     5  2.265900363875, 2.476099809130, 2.653199694252, 2.776399994784, InO  
     6  2.901399637785, 3.053400268852, 3.184699948214, 3.346200303668, InO  
     7  3.477499714125, 3.605100011609, 3.709500024004, 3.795800099653, InO  
     8  3.881100145292, 3.952899891703, 3.982199601697, 3.990999799133, InO  
     9  4.000000000000,     57*0.0D+00/                                 InO  
      DATA  K_InO/                                                      071215
     1  2.38171124D-05, 1.57091123D-01, 3.14378002D-01, 7.98178753D-01, InO  
     2  1.89144362D+00, 2.98723335D+00, 4.01178442D+00, 4.59480018D+00, InO  
     3  5.10924432D+00, 5.92596333D+00, 6.63090581D+00, 7.21292727D+00, InO  
     4  7.58670071D+00, 7.94802936D+00, 8.23681672D+00, 8.51289363D+00, InO  
     5  8.92545103D+00, 9.34385318D+00, 9.64765870D+00, 9.83301576D+00, InO  
     6  1.00058983D+01, 1.02077713D+01, 1.03804452D+01, 1.05860814D+01, InO  
     7  1.07405442D+01, 1.08757739D+01, 1.09751826D+01, 1.10514692D+01, InO  
     8  1.11282397D+01, 1.12063970D+01, 1.12458960D+01, 1.12589309D+01, InO  
     9  1.12728499D+01,     57*0.0D+00/                                 InO  
      DATA TK_SnO/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718499834746, 0.748499952777, SnO  
     2  0.824300057081, 0.914499962494, 1.017400036643, 1.154399984723, SnO  
     3  1.285600067901, 1.430299822375, 1.577099991186, 1.697299883793, SnO  
     4  1.816600107662, 1.937499985755, 2.051400230563, 2.245799912993, SnO  
     5  2.413899944996, 2.541499990993, 2.666600006759, 2.870699907819, SnO  
     6  2.966000178875, 3.059900437564, 3.301900230882, 3.431800165810, SnO  
     7  3.566899909338, 3.714500132454, 3.836100059594, 3.923900187807, SnO  
     8  4.000000000000,     61*0.0D+00/                                 SnO  
      DATA  K_SnO/                                                      071215
     1 -1.91993806D-04, 1.51821354D-01, 3.02597005D-01, 7.69403682D-01, SnO  
     2  1.82768136D+00, 2.89195115D+00, 3.89527639D+00, 4.95985222D+00, SnO  
     3  5.76059366D+00, 6.46290675D+00, 7.03564894D+00, 7.42969450D+00, SnO  
     4  7.77380045D+00, 8.08935915D+00, 8.36547406D+00, 8.80252907D+00, SnO  
     5  9.14878421D+00, 9.38981127D+00, 9.60834302D+00, 9.94710894D+00, SnO  
     6  1.01092733D+01, 1.02752753D+01, 1.07128836D+01, 1.09338416D+01, SnO  
     7  1.11434999D+01, 1.13459515D+01, 1.14866649D+01, 1.15689916D+01, SnO  
     8  1.16293841D+01,     61*0.0D+00/                                 SnO  
      DATA TK_SbO/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720799815989, 0.753800085665, SbO  
     2  0.837700095647, 0.934499932174, 1.054200082210, 1.189900042043, SbO  
     3  1.332799999119, 1.481699974944, 1.640999885715, 1.777900166867, SbO  
     4  1.921700033058, 2.051100223654, 2.188799675568, 2.367899837063, SbO  
     5  2.598600344038, 2.861899690103, 3.014800333157, 3.182300122266, SbO  
     6  3.347100323912, 3.534799839019, 3.615199783155, 3.696599737352, SbO  
     7  3.774700125825, 3.856999875983, 3.937799761155, 3.975299890018, SbO  
     8  4.000000000000,     61*0.0D+00/                                 SbO  
      DATA  K_SbO/                                                      071215
     1 -7.67819929D-05, 1.70276781D-01, 3.53780293D-01, 8.84535716D-01, SbO  
     2  2.08127822D+00, 3.23050652D+00, 4.37878977D+00, 5.40154238D+00, SbO  
     3  6.24214182D+00, 6.93215501D+00, 7.52270526D+00, 7.94646910D+00, SbO  
     4  8.33693705D+00, 8.65683868D+00, 8.97416950D+00, 9.35624078D+00, SbO  
     5  9.79042648D+00, 1.01853178D+01, 1.03577311D+01, 1.04984769D+01, SbO  
     6  1.05996045D+01, 1.07012655D+01, 1.07505239D+01, 1.08065703D+01, SbO  
     7  1.08646451D+01, 1.09256567D+01, 1.09796533D+01, 1.10019474D+01, SbO  
     8  1.10158904D+01,     61*0.0D+00/                                 SbO  
      DATA TK_TeO/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751100015691, TeO  
     2  0.831099944349, 0.926599926485, 1.034700095902, 1.108400103961, TeO  
     3  1.182699892130, 1.324300064703, 1.479700007352, 1.640599876281, TeO  
     4  1.759700231591, 1.888899892620, 2.002900066591, 2.118399848731, TeO  
     5  2.284299831146, 2.464200120948, 2.560400384067, 2.659699853547, TeO  
     6  2.876400028797, 3.067399888384, 3.238299733295, 3.405899746875, TeO  
     7  3.521600340426, 3.633400192130, 3.733599632526, 3.822499745409, TeO  
     8  3.923700183363, 3.969200252915, 4.000000000000,     59*0.0D+00/ TeO  
      DATA  K_TeO/                                                      071215
     1  1.11994920D-04, 1.72265176D-01, 3.51573527D-01, 8.86909686D-01, TeO  
     2  2.09845790D+00, 3.30699412D+00, 4.42355194D+00, 5.06117889D+00, TeO  
     3  5.62075595D+00, 6.50137877D+00, 7.25603778D+00, 7.87395422D+00, TeO  
     4  8.25705577D+00, 8.62406618D+00, 8.91878803D+00, 9.19727977D+00, TeO  
     5  9.56719645D+00, 9.92196657D+00, 1.00865559D+01, 1.02360739D+01, TeO  
     6  1.04958315D+01, 1.06681519D+01, 1.07978959D+01, 1.09191573D+01, TeO  
     7  1.10053001D+01, 1.10921122D+01, 1.11720904D+01, 1.12423198D+01, TeO  
     8  1.13168813D+01, 1.13475642D+01, 1.13674706D+01,     59*0.0D+00/ TeO  
      DATA TK_IO/                                                       071215
     1  0.699999789529, 0.710000049601, 0.721199825569, 0.754700108990, IO   
     2  0.840100145959, 0.938500010103, 1.060099954764, 1.198999824682, IO   
     3  1.345300032409, 1.500099994521, 1.665099952175, 1.797000116582, IO   
     4  1.938800013612, 2.067599862127, 2.193299655628, 2.495899899990, IO   
     5  2.892700148993, 3.022600268497, 3.151799644456, 3.362699705687, IO   
     6  3.515100342308, 3.659199836433, 3.807600354305, 3.890000356544, IO   
     7  3.959900044552, 4.000000000000,     64*0.0D+00/                 IO   
      DATA  K_IO/                                                       071215
     1  4.88599562D-06, 1.64770030D-01, 3.45249085D-01, 8.60665750D-01, IO   
     2  2.02429579D+00, 3.13803858D+00, 4.24829590D+00, 5.24267907D+00, IO   
     3  6.05861139D+00, 6.73748823D+00, 7.31630687D+00, 7.70623744D+00, IO   
     4  8.07867908D+00, 8.38917662D+00, 8.67334342D+00, 9.28115182D+00, IO   
     5  9.86539819D+00, 9.99518420D+00, 1.00936508D+01, 1.02042295D+01, IO   
     6  1.02648847D+01, 1.03198031D+01, 1.03756725D+01, 1.04047012D+01, IO   
     7  1.04282242D+01, 1.04421743D+01,     64*0.0D+00/                 IO   
      DATA TK_BaO/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718599832218, 0.748899961951, BaO  
     2  0.825400030459, 0.916500007268, 1.020899996753, 1.157799897784, BaO  
     3  1.287900004534, 1.432099863238, 1.585900065128, 1.730200016083, BaO  
     4  1.884200018724, 2.015500359075, 2.149199678116, 2.312500249121, BaO  
     5  2.482699595332, 2.616099784400, 2.748399953136, 2.970600227065, BaO  
     6  3.128500083097, 3.214000116300, 3.297200128811, 3.396099912313, BaO  
     7  3.461400065700, 3.519000433385, 3.577299860942, 3.654499726765, BaO  
     8  3.736599694503, 3.805700312655, 3.903499698251, 3.959500035818, BaO  
     9  3.984599655579, 4.000000000000,     56*0.0D+00/                 BaO  
      DATA  K_BaO/                                                      071215
     1  1.76986566D-04, 1.52864850D-01, 3.02659493D-01, 7.70950092D-01, BaO  
     2  1.83080226D+00, 2.89592352D+00, 3.90258149D+00, 4.95369801D+00, BaO  
     3  5.73966006D+00, 6.43392943D+00, 7.02720990D+00, 7.48762941D+00, BaO  
     4  7.91135146D+00, 8.23721208D+00, 8.54622608D+00, 8.89681718D+00, BaO  
     5  9.22677937D+00, 9.45541531D+00, 9.65445386D+00, 9.93072063D+00, BaO  
     6  1.00901598D+01, 1.01675097D+01, 1.02415416D+01, 1.03402236D+01, BaO  
     7  1.04216802D+01, 1.05092998D+01, 1.06132298D+01, 1.07669803D+01, BaO  
     8  1.09343930D+01, 1.10694748D+01, 1.12568730D+01, 1.13690125D+01, BaO  
     9  1.14211204D+01, 1.14536361D+01,     56*0.0D+00/                 BaO  
      DATA TK_LaO/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718299839802, 0.747999941308, LaO  
     2  0.823000088543, 0.912199911003, 1.014200113093, 1.149400084006, LaO  
     3  1.278800193352, 1.421100051204, 1.566900100474, 1.691500033031, LaO  
     4  1.825900027595, 1.940500027515, 2.062700227363, 2.218900227929, LaO  
     5  2.390700290054, 2.650099618280, 2.849600375805, 3.022600268497, LaO  
     6  3.133200193822, 3.237999726009, 3.381200134863, 3.492999826702, LaO  
     7  3.575599820329, 3.647999763049, 3.805300303887, 3.906499765622, LaO  
     8  4.000000000000,     61*0.0D+00/                                 LaO  
      DATA  K_LaO/                                                      071215
     1 -4.05868144D-04, 1.58087912D-01, 3.11933554D-01, 7.93999993D-01, LaO  
     2  1.88693429D+00, 2.98643454D+00, 4.02600442D+00, 5.12466621D+00, LaO  
     3  5.95081498D+00, 6.67305552D+00, 7.26645517D+00, 7.68937158D+00, LaO  
     4  8.08496054D+00, 8.38809101D+00, 8.68756547D+00, 9.04406485D+00, LaO  
     5  9.40690727D+00, 9.90203414D+00, 1.02517680D+01, 1.05452356D+01, LaO  
     6  1.07342874D+01, 1.09157324D+01, 1.11617911D+01, 1.13451318D+01, LaO  
     7  1.14729180D+01, 1.15799686D+01, 1.18070593D+01, 1.19584820D+01, LaO  
     8  1.21034556D+01,     61*0.0D+00/                                 LaO  
      DATA TK_TbO/                                                      071215
     1  0.699999789529, 0.709000023594, 0.717699854968, 0.746699911489, TbO  
     2  0.819600152678, 0.905999953956, 1.005900125926, 1.072300123172, TbO  
     3  1.139699884138, 1.267200047663, 1.400200089630, 1.538100025096, TbO  
     4  1.670600053922, 1.820500178517, 1.965399914860, 2.108299799051, TbO  
     5  2.222200094866, 2.330099929537, 2.566099978385, 2.666199997919, TbO  
     6  2.779199784713, 2.898799695512, 3.018100408788, 3.268000412527, TbO  
     7  3.587700105848, 3.695699718728, 3.799100168838, 3.878400082987, TbO  
     8  3.954499926640, 4.000000000000,     60*0.0D+00/                 TbO  
      DATA  K_TbO/                                                      071215
     1 -1.62938853D-04, 1.73626850D-01, 3.38586781D-01, 8.66816402D-01, TbO  
     2  2.06092000D+00, 3.26122469D+00, 4.41064611D+00, 5.05818317D+00, TbO  
     3  5.63607990D+00, 6.54986546D+00, 7.30797574D+00, 7.93926722D+00, TbO  
     4  8.43819720D+00, 8.91469333D+00, 9.31926031D+00, 9.69049350D+00, TbO  
     5  9.97940165D+00, 1.02520123D+01, 1.08317023D+01, 1.10588710D+01, TbO  
     6  1.12954307D+01, 1.15224700D+01, 1.17288829D+01, 1.21280379D+01, TbO  
     7  1.26713372D+01, 1.28850883D+01, 1.31103162D+01, 1.32957576D+01, TbO  
     8  1.34804626D+01, 1.35920639D+01,     60*0.0D+00/                 TbO  
      DATA TK_LuO/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718699829691, 0.749099966539, LuO  
     2  0.825800020779, 0.917100020700, 1.021200004161, 1.090800145635, LuO  
     3  1.160899861910, 1.294300043050, 1.440300035571, 1.591100137788, LuO  
     4  1.717899843678, 1.857100069962, 1.982399950113, 2.109199734902, LuO  
     5  2.291599986522, 2.492999842200, 2.654299721209, 2.773100242367, LuO  
     6  2.886800271471, 3.054400294808, 3.422400140993, 3.506800148605, LuO  
     7  3.592500214976, 3.687099802754, 3.769600456272, 3.835500046525, LuO  
     8  3.900299626388, 3.961000069146, 3.984999664560, 4.000000000000, LuO  
     9      58*0.0D+00/                                                 LuO  
      DATA  K_LuO/                                                      071215
     1  1.70670455D-05, 1.60647627D-01, 3.19895647D-01, 8.13866005D-01, LuO  
     2  1.93005177D+00, 3.04965861D+00, 4.10100881D+00, 4.69479744D+00, LuO  
     3  5.21962437D+00, 6.05336072D+00, 6.77554115D+00, 7.37180755D+00, LuO  
     4  7.78998520D+00, 8.18926816D+00, 8.51268113D+00, 8.81630617D+00, LuO  
     5  9.22126580D+00, 9.62571850D+00, 9.91193002D+00, 1.01010516D+01, LuO  
     6  1.02687445D+01, 1.05024455D+01, 1.09813499D+01, 1.10820824D+01, LuO  
     7  1.11809878D+01, 1.12890576D+01, 1.13876460D+01, 1.14744130D+01, LuO  
     8  1.15709404D+01, 1.16734339D+01, 1.17170791D+01, 1.17451557D+01, LuO  
     9      58*0.0D+00/                                                 LuO  
      DATA TK_HfO/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718499834746, 0.748499952777, HfO  
     2  0.824400054661, 0.914599964732, 1.017600031865, 1.086400077556, HfO  
     3  1.155599954039, 1.287700010044, 1.431399847347, 1.577999968079, HfO  
     4  1.703899906080, 1.828899943750, 1.958600020125, 2.076899856213, HfO  
     5  2.235599683129, 2.364499754690, 2.540699973053, 2.668800055379, HfO  
     6  2.776399994784, 2.876200024552, 3.106499910289, 3.214600130399, HfO  
     7  3.321599732637, 3.480199539348, 3.659999855099, 3.735299667646, HfO  
     8  3.810600363598, 3.927200261124, 3.971400170248, 4.000000000000, HfO  
     9      58*0.0D+00/                                                 HfO  
      DATA  K_HfO/                                                      071215
     1  1.80245752D-04, 1.68152516D-01, 3.34726920D-01, 8.50237924D-01, HfO  
     2  2.01900430D+00, 3.19030011D+00, 4.29240061D+00, 4.91429536D+00, HfO  
     3  5.46316975D+00, 6.33708160D+00, 7.08817310D+00, 7.69982351D+00, HfO  
     4  8.13619226D+00, 8.51248600D+00, 8.86173458D+00, 9.15499910D+00, HfO  
     5  9.52063432D+00, 9.79748275D+00, 1.01458923D+01, 1.03744837D+01, HfO  
     6  1.05495425D+01, 1.06996649D+01, 1.10214221D+01, 1.11724132D+01, HfO  
     7  1.13271813D+01, 1.15666120D+01, 1.18480877D+01, 1.19684853D+01, HfO  
     8  1.20894212D+01, 1.22730329D+01, 1.23397599D+01, 1.23817072D+01, HfO  
     9      58*0.0D+00/                                                 HfO  
      DATA TK_TaO/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718399837274, 0.748399950483, TaO  
     2  0.824100061921, 0.914099953539, 1.017000046199, 1.154299987280, TaO  
     3  1.285600067901, 1.428799847336, 1.575500032265, 1.697099888939, TaO  
     4  1.830599927128, 1.948599836171, 2.074099786020, 2.207699994012, TaO  
     5  2.342100211448, 2.504100090789, 2.660799878579, 2.822299756538, TaO  
     6  3.057700380461, 3.218500222041, 3.380900127267, 3.592500214976, TaO  
     7  3.767500405501, 3.859699685624, 3.947199765159, 4.000000000000, TaO  
     8      62*0.0D+00/                                                 TaO  
      DATA  K_TaO/                                                      071215
     1 -1.07340253D-04, 1.53529430D-01, 3.04294261D-01, 7.76186097D-01, TaO  
     2  1.84467197D+00, 2.91846758D+00, 3.93302534D+00, 5.01117421D+00, TaO  
     3  5.82021013D+00, 6.52213059D+00, 7.10023420D+00, 7.50206512D+00, TaO  
     4  7.88705906D+00, 8.19380555D+00, 8.49655256D+00, 8.79912371D+00, TaO  
     5  9.08599359D+00, 9.40820676D+00, 9.69219031D+00, 9.95596537D+00, TaO  
     6  1.03060241D+01, 1.05376345D+01, 1.07744115D+01, 1.10931606D+01, TaO  
     7  1.13653272D+01, 1.15123172D+01, 1.16535919D+01, 1.17386781D+01, TaO  
     8      62*0.0D+00/                                                 TaO  
      DATA TK_WO/                                                       071215
     1  0.699999789529, 0.708900020993, 0.717499860023, 0.745999895433, WO   
     2  0.817900116682, 0.903100020804, 1.001300027746, 1.128400083244, WO   
     3  1.249900018117, 1.383900004034, 1.518799997734, 1.636099959253, WO   
     4  1.754500100568, 1.873499970609, 1.989399786748, 2.114399767374, WO   
     5  2.242699847168, 2.447099739907, 2.597300316178, 2.744399868019, WO   
     6  2.887000276361, 2.956699964895, 3.030099747635, 3.115999785654, WO   
     7  3.220800199242, 3.312900202896, 3.418700058718, 3.512700286261, WO   
     8  3.607999801735, 3.711000056958, 3.813500154234, 3.924800207802, WO   
     9  3.970900206175, 4.000000000000,     56*0.0D+00/                 WO   
      DATA  K_WO/                                                       071215
     1 -2.16325388D-05, 1.48760544D-01, 2.90001924D-01, 7.40152258D-01, WO   
     2  1.76453695D+00, 2.79912603D+00, 3.79308057D+00, 4.83039891D+00, WO   
     3  5.62037749D+00, 6.32119495D+00, 6.89300135D+00, 7.30997811D+00, WO   
     4  7.67578553D+00, 8.00375281D+00, 8.29667795D+00, 8.59194508D+00, WO   
     5  8.87770500D+00, 9.30081401D+00, 9.58683788D+00, 9.85248466D+00, WO   
     6  1.01165799D+01, 1.02569458D+01, 1.04161333D+01, 1.06159879D+01, WO   
     7  1.08707216D+01, 1.10936021D+01, 1.13385249D+01, 1.15435814D+01, WO   
     8  1.17433061D+01, 1.19603371D+01, 1.21904898D+01, 1.24658151D+01, WO   
     9  1.25873168D+01, 1.26657149D+01,     56*0.0D+00/                 WO   
      DATA TK_PtO/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718699829691, 0.749099966539, PtO  
     2  0.825900018358, 0.917200022939, 1.021300006631, 1.091900117499, PtO  
     3  1.163299916258, 1.299700164077, 1.447699862632, 1.598499932854, PtO  
     4  1.711600007895, 1.837000077762, 1.951399838429, 2.058100384860, PtO  
     5  2.379200107485, 2.541099982023, 2.678800281158, 2.806800347729, PtO  
     6  2.925700199083, 3.071999744604, 3.229299582502, 3.540699974099, PtO  
     7  3.665099959206, 3.805100299503, 3.903799704988, 4.000000000000, PtO  
     8      62*0.0D+00/                                                 PtO  
      DATA  K_PtO/                                                      071215
     1  2.38334156D-04, 1.73097152D-01, 3.44445713D-01, 8.75795205D-01, PtO  
     2  2.07688496D+00, 3.27816742D+00, 4.40392396D+00, 5.04690260D+00, PtO  
     3  5.61525532D+00, 6.51560284D+00, 7.28218880D+00, 7.90255793D+00, PtO  
     4  8.29186546D+00, 8.66995938D+00, 8.98024151D+00, 9.24854338D+00, PtO  
     5  9.97536261D+00, 1.03138777D+01, 1.05930297D+01, 1.08409705D+01, PtO  
     6  1.10547131D+01, 1.12890384D+01, 1.15036711D+01, 1.18376894D+01, PtO  
     7  1.19494130D+01, 1.20619342D+01, 1.21297373D+01, 1.21880029D+01, PtO  
     8      62*0.0D+00/                                                 PtO  
      DATA TK_PbO/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750600002732, PbO  
     2  0.829799923973, 0.924299981942, 1.031300186745, 1.104300009445, PbO  
     3  1.177799886967, 1.318000121689, 1.471899824708, 1.631200081280, PbO  
     4  1.751200017419, 1.882500064336, 2.001400032147, 2.123199948390, PbO  
     5  2.303000265582, 2.445699709447, 2.592800219740, 2.806400339305, PbO  
     6  3.061400335769, 3.236099679863, 3.399399671363, 3.485999667977, PbO  
     7  3.567799844655, 3.631600155023, 3.699399795293, 3.777099952149, PbO  
     8  3.865799796996, 3.946499749354, 3.979299602603, 4.000000000000, PbO  
     9      58*0.0D+00/                                                 PbO  
      DATA  K_PbO/                                                      071215
     1 -5.33413254D-05, 1.56233018D-01, 3.14262272D-01, 7.97863689D-01, PbO  
     2  1.89069880D+00, 2.98455579D+00, 4.00037531D+00, 4.58354660D+00, PbO  
     3  5.09669706D+00, 5.90980203D+00, 6.61282301D+00, 7.19377938D+00, PbO  
     4  7.56298082D+00, 7.92145084D+00, 8.21816256D+00, 8.50277870D+00, PbO  
     5  8.89352412D+00, 9.17699506D+00, 9.43951175D+00, 9.76125056D+00, PbO  
     6  1.00595072D+01, 1.02226546D+01, 1.03653564D+01, 1.04485925D+01, PbO  
     7  1.05356789D+01, 1.06066456D+01, 1.06787074D+01, 1.07474019D+01, PbO  
     8  1.07996191D+01, 1.08286328D+01, 1.08396790D+01, 1.08476440D+01, PbO  
     9      58*0.0D+00/                                                 PbO  
      DATA TK_BiO/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720999820779, 0.754300098624, BiO  
     2  0.839100127740, 0.937099982828, 1.057700003907, 1.193599956348, BiO  
     3  1.337500106888, 1.491799820174, 1.655299956962, 1.791900011739, BiO  
     4  1.926899902733, 2.068999757774, 2.191299615563, 2.555700313233, BiO  
     5  2.749699980799, 2.936299850074, 3.082599995447, 3.222600068638, BiO  
     6  3.393700087549, 3.542200009977, 3.739799760611, 3.876700044529, BiO  
     7  3.954099917906, 3.982499608432, 4.000000000000,     63*0.0D+00/ BiO  
      DATA  K_BiO/                                                      071215
     1  9.68561972D-06, 1.62685750D-01, 3.41114020D-01, 8.52347369D-01, BiO  
     2  2.00626414D+00, 3.11514616D+00, 4.21740055D+00, 5.19470031D+00, BiO  
     3  6.00429246D+00, 6.68785848D+00, 7.26705655D+00, 7.67338770D+00, BiO  
     4  8.02972000D+00, 8.37302271D+00, 8.64941176D+00, 9.36756505D+00, BiO  
     5  9.67032374D+00, 9.90765867D+00, 1.00617380D+01, 1.01869824D+01, BiO  
     6  1.03124843D+01, 1.03949175D+01, 1.04688706D+01, 1.05062021D+01, BiO  
     7  1.05275813D+01, 1.05365417D+01, 1.05425947D+01,     63*0.0D+00/ BiO  
      DATA TK_ThO/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718499834746, 0.748699957364, ThO  
     2  0.824700047400, 0.915199978165, 1.018500010363, 1.087400102151, ThO  
     3  1.156599928468, 1.288699982493, 1.432799879129, 1.580499929307, ThO  
     4  1.705099934314, 1.843700059098, 1.971199820417, 2.102000248092, ThO  
     5  2.215200156035, 2.328899903348, 2.531699755486, 2.620899894806, ThO  
     6  2.679200290215, 2.699899805478, 2.711900098875, 2.726799789591, ThO  
     7  2.767900424219, 2.864399752714, 2.947399753304, 2.981599577164, ThO  
     8  2.989999765965, 2.999699992979, 3.008000178518, 3.013500303363, ThO  
     9  3.021900317991, 3.034099842122, 3.090000139165, 3.144200036255, ThO  
     A  3.159999831214, 3.175600194744, 3.181700165779, 3.190299570880, ThO  
     B  3.199999798209, 3.236699694435, 3.265800359701, 3.280999766978, ThO  
     C  3.289799969268, 3.299500178307, 3.309200391537, 3.316299961088, ThO  
     D  3.345700292422, 3.376900040855, 3.388500319721, 3.398099766283, ThO  
     E  3.412799903470, 3.421300118962, 3.434100003021, 3.458700003975, ThO  
     F  3.468800225801, 3.479699556418, 3.492799822034, 3.501100015718, ThO  
     G  3.515400349314, 3.528499838828, 3.544600067381, 3.558500373931, ThO  
     H  3.571199715211, 3.588400122248, 3.600200366225, 3.611399690960, ThO  
     I  3.623299973012, 3.639300313760, 3.651199649764, 3.663199920421, ThO  
     J  3.675200183002, 3.689099663459, 3.699399795293, 3.710800052644, ThO  
     K  3.725599863046, 3.741199791436, 3.750600001402, 3.762400282203, ThO  
     L  3.779499778472, 3.794500072398, 3.814100110917, 3.831899968113, ThO  
     M  3.848500335761, 3.861799705601, 3.875200010596, 3.901699657828, ThO  
     N  3.948799801285, 4.000000000000/                                 ThO  
      DATA  K_ThO/                                                      071215
     1 -4.47770202D-04, 1.66434789D-01, 3.31930433D-01, 8.47423669D-01, ThO  
     2  2.00966428D+00, 3.17645397D+00, 4.27343329D+00, 4.89143659D+00, ThO  
     3  5.43614947D+00, 6.30374372D+00, 7.05167606D+00, 7.66342301D+00, ThO  
     4  8.09257016D+00, 8.50528541D+00, 8.84414384D+00, 9.16421513D+00, ThO  
     5  9.42409567D+00, 9.67096886D+00, 1.00743009D+01, 1.02345709D+01, ThO  
     6  1.03331149D+01, 1.03669875D+01, 1.03868573D+01, 1.04112094D+01, ThO  
     7  1.04766573D+01, 1.06203101D+01, 1.07329419D+01, 1.07765195D+01, ThO  
     8  1.07869497D+01, 1.07988551D+01, 1.08109923D+01, 1.08192104D+01, ThO  
     9  1.08313876D+01, 1.08490162D+01, 1.09270100D+01, 1.09982549D+01, ThO  
     A  1.10181773D+01, 1.10379517D+01, 1.10467793D+01, 1.10603492D+01, ThO  
     B  1.10753650D+01, 1.11307054D+01, 1.11728200D+01, 1.11941621D+01, ThO  
     C  1.12062955D+01, 1.12194149D+01, 1.12351312D+01, 1.12468867D+01, ThO  
     D  1.12941312D+01, 1.13420555D+01, 1.13592527D+01, 1.13735072D+01, ThO  
     E  1.13981835D+01, 1.14122490D+01, 1.14331295D+01, 1.14718868D+01, ThO  
     F  1.14872530D+01, 1.15040996D+01, 1.15260496D+01, 1.15396248D+01, ThO  
     G  1.15625298D+01, 1.15829107D+01, 1.16074638D+01, 1.16299695D+01, ThO  
     H  1.16501758D+01, 1.16766022D+01, 1.16940694D+01, 1.17116941D+01, ThO  
     I  1.17302387D+01, 1.17543546D+01, 1.17717201D+01, 1.17900029D+01, ThO  
     J  1.18081504D+01, 1.18285261D+01, 1.18432786D+01, 1.18604890D+01, ThO  
     K  1.18822593D+01, 1.19045421D+01, 1.19186046D+01, 1.19358879D+01, ThO  
     L  1.19604725D+01, 1.19826532D+01, 1.20110806D+01, 1.20373411D+01, ThO  
     M  1.20616812D+01, 1.20814328D+01, 1.21012458D+01, 1.21407876D+01, ThO  
     N  1.22116486D+01, 1.22889728D+01/                                 ThO  
      DATA TK_BOp/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749599978008, BOp  
     2  0.827099989317, 0.919400072191, 1.024400083186, 1.095700020303, BOp  
     3  1.167800018160, 1.305300058013, 1.454899915377, 1.608300097823, BOp  
     4  1.721599826051, 1.845600013255, 1.956699972178, 2.068099824858, BOp  
     5  2.283899823038, 2.556200324780, 2.738999752766, 2.937299779897, BOp  
     6  3.083100005158, 3.235399662862, 3.405799744853, 3.589900157392, BOp  
     7  3.804200279774, 3.928200283341, 3.971500163063, 4.000000000000, BOp  
     8      62*0.0D+00/                                                 BOp  
      DATA  K_BOp/                                                      071215
     1  7.89712194D-05, 1.59495997D-01, 3.17484379D-01, 8.07018201D-01, BOp  
     2  1.91274996D+00, 3.02085182D+00, 4.05770608D+00, 4.65142533D+00, BOp  
     3  5.17689304D+00, 6.01023917D+00, 6.72532035D+00, 7.31129406D+00, BOp  
     4  7.67628559D+00, 8.02978284D+00, 8.31792476D+00, 8.58775911D+00, BOp  
     5  9.07103053D+00, 9.62082142D+00, 9.95552138D+00, 1.02821734D+01, BOp  
     6  1.04929045D+01, 1.06847416D+01, 1.08677333D+01, 1.10348020D+01, BOp  
     7  1.12036621D+01, 1.12990493D+01, 1.13345485D+01, 1.13591718D+01, BOp  
     8      62*0.0D+00/                                                 BOp  
      DATA TK_SiOp/                                                     071215
     1  0.699999789529, 0.709200028796, 0.718299839802, 0.748199945895, SiOp 
     2  0.823500076442, 0.913099931151, 1.015500082035, 1.085600057880, SiOp 
     3  1.156499931026, 1.291999991501, 1.435699944965, 1.577599978349, SiOp 
     4  1.701799856671, 1.832299967141, 1.941500003892, 2.044700076995, SiOp 
     5  2.144400013458, 2.236299700105, 2.438299704081, 2.620199878761, SiOp 
     6  2.824599812059, 3.010600236899, 3.188999636371, 3.359899641647, SiOp 
     7  3.516800382009, 3.765000345061, 3.915799993246, 3.967300210335, SiOp 
     8  4.000000000000,     61*0.0D+00/                                 SiOp 
      DATA  K_SiOp/                                                     071215
     1 -4.99615999D-06, 1.53996099D-01, 3.03488548D-01, 7.75026553D-01, SiOp 
     2  1.84100061D+00, 2.91379761D+00, 3.92762103D+00, 4.51498568D+00, SiOp 
     3  5.03598028D+00, 5.86591009D+00, 6.56366049D+00, 7.11975161D+00, SiOp 
     4  7.52909895D+00, 7.90626572D+00, 8.19640604D+00, 8.46066801D+00, SiOp 
     5  8.71366953D+00, 8.94720124D+00, 9.45507459D+00, 9.88710334D+00, SiOp 
     6  1.03234493D+01, 1.06639403D+01, 1.09363242D+01, 1.11515712D+01, SiOp 
     7  1.13166904D+01, 1.15340017D+01, 1.16527393D+01, 1.16928497D+01, SiOp 
     8  1.17185080D+01,     61*0.0D+00/                                 SiOp 
      DATA TK_POp/                                                      071215
     1  0.699999789529, 0.709100026195, 0.717999847385, 0.747299925252, POp  
     2  0.821000136946, 0.908999884802, 1.010400203878, 1.142299927492, POp  
     3  1.264400113435, 1.403800001756, 1.553800089028, 1.667099999341, POp  
     4  1.764900106541, 1.858200095195, 1.953599893947, 2.056100338801, POp  
     5  2.164399945867, 2.353700114056, 2.539699950060, 2.719600272383, POp  
     6  2.898099747551, 3.072599759630, 3.207399964042, 3.335100041985, POp  
     7  3.593600239282, 3.803200257854, 4.000000000000,     63*0.0D+00/ POp  
      DATA  K_POp/                                                      071215
     1 -1.08732311D-04, 1.46459779D-01, 2.87204023D-01, 7.32428515D-01, POp  
     2  1.74039965D+00, 2.76237040D+00, 3.73950716D+00, 4.75904206D+00, POp  
     3  5.50965808D+00, 6.19881235D+00, 6.79617383D+00, 7.17991543D+00, POp  
     4  7.48329555D+00, 7.76103454D+00, 8.04147018D+00, 8.34300865D+00, POp  
     5  8.66215749D+00, 9.21345880D+00, 9.72778604D+00, 1.01767517D+01, POp  
     6  1.05595027D+01, 1.08691221D+01, 1.10680780D+01, 1.12294703D+01, POp  
     7  1.15010847D+01, 1.16974362D+01, 1.18802761D+01,     63*0.0D+00/ POp  
      DATA TK_SOp/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, SOp  
     2  0.829599928813, 0.923599998821, 1.030300213463, 1.102599970255, SOp  
     3  1.175499940343, 1.314900058258, 1.469199799898, 1.628900086844, SOp  
     4  1.741899819704, 1.863100058222, 1.970699809865, 2.072199738388, SOp  
     5  2.282299790608, 2.473599990537, 2.637500281101, 2.931900158853, SOp  
     6  3.213800111600, 3.342300215946, 3.467300193348, 3.572899755825, SOp  
     7  3.681500192782, 3.803500264430, 3.934999962962, 3.973500019355, SOp  
     8  4.000000000000,     61*0.0D+00/                                 SOp  
      DATA  K_SOp/                                                      071215
     1  7.40339447D-05, 1.66147965D-01, 3.34057146D-01, 8.46189265D-01, SOp  
     2  2.00520919D+00, 3.15987812D+00, 4.23379194D+00, 4.84551204D+00, SOp  
     3  5.38422607D+00, 6.23837996D+00, 6.97961782D+00, 7.58796049D+00, SOp  
     4  7.95006361D+00, 8.29457783D+00, 8.57248700D+00, 8.81528524D+00, SOp  
     5  9.26378633D+00, 9.61438609D+00, 9.88021062D+00, 1.02896184D+01, SOp  
     6  1.05992632D+01, 1.07168206D+01, 1.08205150D+01, 1.09035241D+01, SOp  
     7  1.09902156D+01, 1.10969903D+01, 1.12247026D+01, 1.12626603D+01, SOp  
     8  1.12883422D+01,     61*0.0D+00/                                 SOp  
      DATA TK_AsOp/                                                     071215
     1  0.699999789529, 0.709100026195, 0.718199842329, 0.747799936720, AsOp 
     2  0.822500100643, 0.910999884138, 1.012400156096, 1.147600044326, AsOp 
     3  1.277500162114, 1.419900078029, 1.564200036894, 1.672200017099, AsOp 
     4  1.778800189464, 1.883600034822, 1.980299999122, 2.141800195101, AsOp 
     5  2.295000071193, 2.424900186239, 2.581899969961, 2.712800119156, AsOp 
     6  2.857399835430, 3.140100342197, 3.254900109121, 3.376700036751, AsOp 
     7  3.671200087793, 3.770700415286, 3.873999983450, 3.950599841481, AsOp 
     8  3.981299581491, 4.000000000000,     60*0.0D+00/                 AsOp 
      DATA  K_AsOp/                                                     071215
     1  1.26023995D-04, 1.53004684D-01, 3.03061283D-01, 7.71752591D-01, AsOp 
     2  1.83442727D+00, 2.90095374D+00, 3.91329995D+00, 4.99113058D+00, AsOp 
     3  5.80537914D+00, 6.51566271D+00, 7.09442289D+00, 7.46014301D+00, AsOp 
     4  7.77995320D+00, 8.06525910D+00, 8.30997470D+00, 8.69048545D+00, AsOp 
     5  9.02674552D+00, 9.29758802D+00, 9.61672187D+00, 9.88513806D+00, AsOp 
     6  1.01879956D+01, 1.07697262D+01, 1.09838184D+01, 1.11889116D+01, AsOp 
     7  1.15967085D+01, 1.17136256D+01, 1.18267504D+01, 1.19044337D+01, AsOp 
     8  1.19336001D+01, 1.19506867D+01,     60*0.0D+00/                 AsOp 
      DATA TK_TaOp/                                                     071215
     1  0.699999789529, 0.709100026195, 0.717899849912, 0.747099920664, TaOp 
     2  0.820700144206, 0.907999907853, 1.008200175016, 1.140199881199, TaOp 
     3  1.266500064106, 1.404699979788, 1.545900097495, 1.672000021702, TaOp 
     4  1.811699985409, 1.933899908611, 2.054700306560, 2.183800034097, TaOp 
     5  2.323099777758, 2.540099959599, 2.654099716308, 2.864499755218, TaOp 
     6  3.048300146754, 3.185699875692, 3.391600240881, 3.595000270217, TaOp 
     7  3.676800221086, 3.760300231432, 3.900099621896, 3.961200073629, TaOp 
     8  4.000000000000,     61*0.0D+00/                                 TaOp 
      DATA  K_TaOp/                                                     071215
     1  1.12102871D-05, 1.62184377D-01, 3.16153279D-01, 8.06961178D-01, TaOp 
     2  1.91936138D+00, 3.03832094D+00, 4.10293998D+00, 5.22437691D+00, TaOp 
     3  6.06989105D+00, 6.80700566D+00, 7.41164490D+00, 7.85992927D+00, TaOp 
     4  8.28576517D+00, 8.61651027D+00, 8.91749534D+00, 9.21831533D+00, TaOp 
     5  9.52347541D+00, 9.96719173D+00, 1.01896447D+01, 1.05883004D+01, TaOp 
     6  1.09294251D+01, 1.11868419D+01, 1.15816473D+01, 1.19748067D+01, TaOp 
     7  1.21328156D+01, 1.22932740D+01, 1.25550141D+01, 1.26644897D+01, TaOp 
     8  1.27319890D+01,     61*0.0D+00/                                 TaOp 
      DATA TK_FeOm/                                                     071215
     1  0.699999789529, 0.709200028796, 0.718399837274, 0.748399950483, FeOm 
     2  0.824200059501, 0.914199955777, 1.017100043810, 1.086900089853, FeOm 
     3  1.157599902898, 1.293000013914, 1.436699967666, 1.581899964520, FeOm 
     4  1.697999865782, 1.825100049954, 1.957299987319, 2.100900326496, FeOm 
     5  2.214900150206, 2.336300076506, 2.451499836485, 2.570599714290, FeOm 
     6  2.687999729298, 2.820099703431, 2.967300208765, 3.122499933536, FeOm 
     7  3.265400350096, 3.416800008723, 3.585800061333, 3.730899576746, FeOm 
     8  3.889700349423, 3.956899979045, 3.983599633128, 4.000000000000, FeOm 
     9      58*0.0D+00/                                                 FeOm 
      DATA  K_FeOm/                                                     071215
     1 -1.67030310D-04, 1.69487524D-01, 3.35937356D-01, 8.56690762D-01, FeOm 
     2  2.03577744D+00, 3.21646923D+00, 4.32881818D+00, 4.96551124D+00, FeOm 
     3  5.52994574D+00, 6.42752422D+00, 7.17786319D+00, 7.78341551D+00, FeOm 
     4  8.18756329D+00, 8.57127169D+00, 8.92439291D+00, 9.27306894D+00, FeOm 
     5  9.53372575D+00, 9.80163517D+00, 1.00486575D+01, 1.02953021D+01, FeOm 
     6  1.05259575D+01, 1.07653782D+01, 1.10023668D+01, 1.12174731D+01, FeOm 
     7  1.13885295D+01, 1.15524831D+01, 1.17343753D+01, 1.19102530D+01, FeOm 
     8  1.21562782D+01, 1.22869227D+01, 1.23441168D+01, 1.23808081D+01, FeOm 
     9      58*0.0D+00/                                                 FeOm 
      DATA TK_LiF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750299994957, LiF  
     2  0.828899945754, 0.922400027755, 1.028600186906, 1.102699972560, LiF  
     3  1.178299875364, 1.323600080642, 1.474999897297, 1.627900064735, LiF  
     4  1.754800108127, 1.883300042871, 2.020400435459, 2.154999740474, LiF  
     5  2.266200370428, 2.373099961918, 2.573299774941, 2.733099625187, LiF  
     6  2.888700317928, 3.288199932488, 3.451199822614, 3.610599671551, LiF  
     7  3.720700202585, 3.784699853099, 3.841500178285, 3.934699984585, LiF  
     8  3.974399954687, 3.988899752118, 4.000000000000,     59*0.0D+00/ LiF  
      DATA  K_LiF/                                                      071215
     1  2.92465602D-05, 1.57801784D-01, 3.19048111D-01, 8.09379160D-01, LiF  
     2  1.91668054D+00, 3.02284826D+00, 4.05437505D+00, 4.65938375D+00, LiF  
     3  5.19688428D+00, 6.04903253D+00, 6.74589975D+00, 7.30946098D+00, LiF  
     4  7.70045712D+00, 8.04469521D+00, 8.37005635D+00, 8.66010109D+00, LiF  
     5  8.88505046D+00, 9.09238116D+00, 9.45780931D+00, 9.72072450D+00, LiF  
     6  9.94518403D+00, 1.03769265D+01, 1.05006875D+01, 1.05975553D+01, LiF  
     7  1.06561084D+01, 1.06928441D+01, 1.07340420D+01, 1.08429527D+01, LiF  
     8  1.09133168D+01, 1.09431474D+01, 1.09674939D+01,     59*0.0D+00/ LiF  
      DATA TK_BeF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.750099989774, BeF  
     2  0.828399957855, 0.921400051867, 1.027200152332, 1.100999933371, BeF  
     3  1.176199924098, 1.320900142120, 1.471899824708, 1.624399987353, BeF  
     4  1.748599958236, 1.874299990400, 2.140600278937, 2.357199853840, BeF  
     5  2.549700174872, 2.695199700588, 2.839900165278, 2.997799948512, BeF  
     6  3.150299610293, 3.515400349314, 3.637500276652, 3.761000248356, BeF  
     7  3.905699747657, 3.964300143102, 3.985999687011, 4.000000000000, BeF  
     8      62*0.0D+00/                                                 BeF  
      DATA  K_BeF/                                                      071215
     1 -1.11726526D-04, 1.48310070D-01, 2.98486870D-01, 7.58596125D-01, BeF  
     2  1.79824045D+00, 2.83747017D+00, 3.81075645D+00, 4.38292144D+00, BeF  
     3  4.89182675D+00, 5.70257017D+00, 6.37000690D+00, 6.91276279D+00, BeF  
     4  7.28469936D+00, 7.61422684D+00, 8.20911577D+00, 8.63626839D+00, BeF  
     5  8.99500182D+00, 9.25083593D+00, 9.48473800D+00, 9.71027792D+00, BeF  
     6  9.89697936D+00, 1.02350348D+01, 1.03217489D+01, 1.04063846D+01, BeF  
     7  1.05231921D+01, 1.05813012D+01, 1.06048899D+01, 1.06208074D+01, BeF  
     8      62*0.0D+00/                                                 BeF  
      DATA TK_BF/                                                       071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751500026057, BF   
     2  0.831799960396, 0.925499953008, 1.038499994372, 1.184799935855, BF   
     3  1.341900116962, 1.498199954414, 1.656799917049, 1.778700186953, BF   
     4  1.902300024183, 2.029999726719, 2.156799783043, 2.267700403193, BF   
     5  2.377500066917, 2.566099978385, 2.772300302387, 2.959000017884, BF   
     6  3.141100267577, 3.331699957174, 3.472100101224, 3.605000018846, BF   
     7  3.708700005790, 3.803900273198, 3.918200054911, 3.967000203611, BF   
     8  4.000000000000,     61*0.0D+00/                                 BF   
      DATA  K_BF/                                                       071215
     1  1.47858073D-04, 1.68028529D-01, 3.41157217D-01, 8.63658989D-01, BF   
     2  2.04490350D+00, 3.21100999D+00, 4.37404302D+00, 5.57390791D+00, BF   
     3  6.57292228D+00, 7.34845918D+00, 7.97445529D+00, 8.37490182D+00, BF   
     4  8.72773743D+00, 9.05010914D+00, 9.33971336D+00, 9.57592721D+00, BF   
     5  9.79921885D+00, 1.01647481D+01, 1.05340971D+01, 1.08276543D+01, BF   
     6  1.10688611D+01, 1.12753306D+01, 1.14012136D+01, 1.15025890D+01, BF   
     7  1.15701657D+01, 1.16218675D+01, 1.16656830D+01, 1.16769814D+01, BF   
     8  1.16824399D+01,     61*0.0D+00/                                 BF   
      DATA TK_NaF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750299994957, NaF  
     2  0.828999943334, 0.922500025344, 1.028700189375, 1.101699949507, NaF  
     3  1.175899931060, 1.318500131920, 1.469499792518, 1.622699949768, NaF  
     4  1.752000037576, 1.892899926537, 2.022800258274, 2.155899761759, NaF  
     5  2.261300263396, 2.365799786186, 2.520100449313, 2.684899953673, NaF  
     6  2.899699628605, 3.126500033243, 3.353300139976, 3.479699556418, NaF  
     7  3.601100301092, 3.647599791306, 3.692799658717, 3.736599694503, NaF  
     8  3.778399858074, 3.843400221029, 3.871999938206, 3.901199646599, NaF  
     9  3.934799977377, 3.962000091557, 3.985599678030, 4.000000000000, NaF  
     A      54*0.0D+00/                                                 NaF  
      DATA  K_NaF/                                                      071215
     1 -7.28506848D-05, 1.57237320D-01, 3.18012875D-01, 8.06923852D-01, NaF  
     2  1.91238023D+00, 3.01529591D+00, 4.04391408D+00, 4.63892723D+00, NaF  
     3  5.16680627D+00, 6.00649585D+00, 6.70582512D+00, 7.27362937D+00, NaF  
     4  7.67332169D+00, 8.04909847D+00, 8.35520024D+00, 8.63968017D+00, NaF  
     5  8.84855730D+00, 9.04283203D+00, 9.30587741D+00, 9.55304246D+00, NaF  
     6  9.82148110D+00, 1.00450703D+01, 1.02177262D+01, 1.02948542D+01, NaF  
     7  1.03581726D+01, 1.03810155D+01, 1.04042515D+01, 1.04305028D+01, NaF  
     8  1.04635577D+01, 1.05472545D+01, 1.06033453D+01, 1.06755804D+01, NaF  
     9  1.07774839D+01, 1.08729776D+01, 1.09632361D+01, 1.10208835D+01, NaF  
     A      54*0.0D+00/                                                 NaF  
      DATA TK_MgF/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, MgF  
     2  0.830699935179, 0.925499953008, 1.032900143995, 1.108200099350, MgF  
     3  1.184999940019, 1.333100005998, 1.488399819532, 1.643799951755, MgF  
     4  1.759100216473, 1.882000077751, 2.013400309648, 2.155499752299, MgF  
     5  2.260200239368, 2.362999718349, 2.579999925446, 2.736299694382, MgF  
     6  2.888700317928, 3.306200325514, 3.478299656777, 3.625200015273, MgF  
     7  3.743599844821, 3.798100147873, 3.856299925335, 3.901299648845, MgF  
     8  3.943699686134, 3.977999696013, 3.989899774569, 4.000000000000, MgF  
     9      58*0.0D+00/                                                 MgF  
      DATA  K_MgF/                                                      071215
     1  8.57643181D-05, 1.48005757D-01, 3.00611597D-01, 7.59720152D-01, MgF  
     2  1.80176563D+00, 2.84066681D+00, 3.80711798D+00, 4.37700798D+00, MgF  
     3  4.88340319D+00, 5.68989849D+00, 6.35514473D+00, 6.89072108D+00, MgF  
     4  7.22760231D+00, 7.54490064D+00, 7.84865601D+00, 8.14807829D+00, MgF  
     5  8.35517060D+00, 8.54957407D+00, 8.92856885D+00, 9.16831107D+00, MgF  
     6  9.37130091D+00, 9.78486480D+00, 9.90589432D+00, 9.99008739D+00, MgF  
     7  1.00500931D+01, 1.00793860D+01, 1.01165992D+01, 1.01539928D+01, MgF  
     8  1.02012832D+01, 1.02521890D+01, 1.02731389D+01, 1.02924042D+01, MgF  
     9      58*0.0D+00/                                                 MgF  
      DATA TK_AlF/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720899818384, 0.754200096032, AlF  
     2  0.838800120863, 0.937099982828, 1.057000019567, 1.193599956348, AlF  
     3  1.338200122939, 1.491399811784, 1.654199986231, 1.784000120828, AlF  
     4  1.907799886322, 2.070599698278, 2.234799663728, 2.412299906195, AlF  
     5  2.652999689350, 2.773200234865, 2.885900249465, 3.315799996648, AlF  
     6  3.488899732291, 3.649599650021, 3.724999904622, 3.799200170935, AlF  
     7  3.862599723880, 3.925800230020, 3.970200256473, 3.987899729667, AlF  
     8  4.000000000000,     61*0.0D+00/                                 AlF  
      DATA  K_AlF/                                                      071215
     1  1.04372306D-04, 1.69030721D-01, 3.52648893D-01, 8.83489086D-01, AlF  
     2  2.07843621D+00, 3.23179828D+00, 4.36702614D+00, 5.38284458D+00, AlF  
     3  6.22153087D+00, 6.92299270D+00, 7.53156657D+00, 7.95616000D+00, AlF  
     4  8.32899611D+00, 8.78069979D+00, 9.19628127D+00, 9.60492530D+00, AlF  
     5  1.00880432D+01, 1.02942419D+01, 1.04654992D+01, 1.09485810D+01, AlF  
     6  1.10842211D+01, 1.11864428D+01, 1.12261219D+01, 1.12593003D+01, AlF  
     7  1.12837611D+01, 1.13082388D+01, 1.13296330D+01, 1.13400961D+01, AlF  
     8  1.13480976D+01,     61*0.0D+00/                                 AlF  
      DATA TK_SiF/                                                      071215
     1  0.699999789529, 0.710200044546, 0.722299851912, 0.757300176373, SiF  
     2  0.794500054105, 0.847199974650, 0.954599914164, 1.073900080985, SiF  
     3  1.168800040805, 1.266900054710, 1.352599973221, 1.436099954045, SiF  
     4  1.517200039474, 1.594700038090, 1.756500150962, 1.969799800301, SiF  
     5  2.208800023178, 2.391400238599, 2.578499891751, 2.707399991490, SiF  
     6  2.840900187055, 3.088700113917, 3.333399999580, 3.446299716576, SiF  
     7  3.563000189628, 3.675400187763, 3.779199800182, 3.860599678182, SiF  
     8  3.916700016370, 3.960300053459, 3.985199669050, 4.000000000000, SiF  
     9      58*0.0D+00/                                                 SiF  
      DATA  K_SiF/                                                      071215
     1  3.24461244D-05, 1.53780635D-01, 3.31982579D-01, 8.23188238D-01, SiF  
     2  1.30809914D+00, 1.93548483D+00, 3.02784138D+00, 4.00425636D+00, SiF  
     3  4.64231467D+00, 5.20188535D+00, 5.62623635D+00, 5.99622637D+00, SiF  
     4  6.32534127D+00, 6.61938589D+00, 7.18734741D+00, 7.86699788D+00, SiF  
     5  8.55450270D+00, 9.02848441D+00, 9.46151259D+00, 9.72340187D+00, SiF  
     6  9.96089233D+00, 1.03156551D+01, 1.05792378D+01, 1.06818465D+01, SiF  
     7  1.07797578D+01, 1.08671001D+01, 1.09402466D+01, 1.09906642D+01, SiF  
     8  1.10215593D+01, 1.10443419D+01, 1.10576191D+01, 1.10659002D+01, SiF  
     9      58*0.0D+00/                                                 SiF  
      DATA TK_PF/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750600002732, PF   
     2  0.829699926393, 0.923699996409, 1.030400210791, 1.105100027887, PF   
     3  1.181199860898, 1.327899982732, 1.480899993501, 1.635199981666, PF   
     4  1.759600229072, 1.886399959696, 2.020300442842, 2.150999645876, PF   
     5  2.259700228015, 2.363199723195, 2.566699935681, 2.690099586771, PF   
     6  2.819599729581, 2.964500144387, 3.137200285491, 3.284699852032, PF   
     7  3.413299916627, 3.525700042375, 3.632200167392, 3.768300424842, PF   
     8  3.897299818612, 3.958100005248, 4.000000000000,     59*0.0D+00/ PF   
      DATA  K_PF/                                                       071215
     1 -9.88690845D-05, 1.59268244D-01, 3.20406217D-01, 8.13477743D-01, PF   
     2  1.92613940D+00, 3.03550278D+00, 4.06845854D+00, 4.67580236D+00, PF   
     3  5.21412855D+00, 6.06906992D+00, 6.76769607D+00, 7.33130939D+00, PF   
     4  7.71185394D+00, 8.05021001D+00, 8.36763923D+00, 8.64969217D+00, PF   
     5  8.87005244D+00, 9.07106885D+00, 9.44173459D+00, 9.64617424D+00, PF   
     6  9.84016775D+00, 1.00308276D+01, 1.02235438D+01, 1.03616128D+01, PF   
     7  1.04644228D+01, 1.05440549D+01, 1.06157705D+01, 1.07099987D+01, PF   
     8  1.08010987D+01, 1.08409030D+01, 1.08662511D+01,     59*0.0D+00/ PF   
      DATA TK_SF/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, SF   
     2  0.830399928302, 0.924699972297, 1.031900170714, 1.108400103961, SF   
     3  1.186399969169, 1.334500038099, 1.486699858965, 1.642699925811, SF   
     4  1.753700080411, 1.860400126388, 2.030699743707, 2.217000191011, SF   
     5  2.305000312163, 2.397999753460, 2.597400318321, 2.709800051095, SF   
     6  2.837400108999, 2.977999687131, 3.127900068141, 3.313300174448, SF   
     7  3.556800338453, 3.703799894227, 3.852300207349, 3.929600314445, SF   
     8  4.000000000000,     61*0.0D+00/                                 SF   
      DATA  K_SF/                                                       071215
     1 -1.15598151D-04, 1.58592304D-01, 3.20686545D-01, 8.13166371D-01, SF   
     2  1.92740983D+00, 3.03421025D+00, 4.06577480D+00, 4.68286678D+00, SF   
     3  5.22863316D+00, 6.08053491D+00, 6.76699075D+00, 7.33035106D+00, SF   
     4  7.66838326D+00, 7.95559608D+00, 8.35349617D+00, 8.72622578D+00, SF   
     5  8.88783198D+00, 9.05222905D+00, 9.38451338D+00, 9.55640454D+00, SF   
     6  9.73451637D+00, 9.90849340D+00, 1.00695745D+01, 1.02394122D+01, SF   
     7  1.04293071D+01, 1.05354109D+01, 1.06390773D+01, 1.06910573D+01, SF   
     8  1.07372934D+01,     61*0.0D+00/                                 SF   
      DATA TK_KF/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750399997549, KF   
     2  0.829199938494, 0.922800018110, 1.029100199253, 1.101899954118, KF   
     3  1.175699935701, 1.317700115550, 1.469299797438, 1.623599969666, KF   
     4  1.755200118206, 1.893499939661, 2.021200376397, 2.151599660065, KF   
     5  2.351600270185, 2.548500147962, 2.733099625187, 3.079599934934, KF   
     6  3.340400173210, 3.453899887904, 3.566799916525, 3.657699801432, KF   
     7  3.735799677975, 3.805700312655, 3.839800140184, 3.874299990236, KF   
     8  3.952899891703, 3.982699612923, 4.000000000000,     59*0.0D+00/ KF   
      DATA  K_KF/                                                       071215
     1  2.14931542D-05, 1.57723792D-01, 3.17181875D-01, 8.02117826D-01, KF   
     2  1.89998263D+00, 2.99509167D+00, 4.01636838D+00, 4.60517763D+00, KF   
     3  5.12657833D+00, 5.95796816D+00, 6.65648207D+00, 7.22513620D+00, KF   
     4  7.62914574D+00, 7.99570588D+00, 8.29491793D+00, 8.57049999D+00, KF   
     5  8.94491818D+00, 9.26121110D+00, 9.51159080D+00, 9.87149113D+00, KF   
     6  1.00644013D+01, 1.01307084D+01, 1.01888371D+01, 1.02375677D+01, KF   
     7  1.02980717D+01, 1.03924644D+01, 1.04600265D+01, 1.05445115D+01, KF   
     8  1.07886987D+01, 1.08938155D+01, 1.09564812D+01,     59*0.0D+00/ KF   
      DATA TK_CaF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.749999987182, CaF  
     2  0.828099965115, 0.921100059100, 1.026700139985, 1.099099933338, CaF  
     3  1.172500009963, 1.313600031657, 1.463099949959, 1.614300032600, CaF  
     4  1.742399830043, 1.880300123363, 2.011900274343, 2.149499657157, CaF  
     5  2.264100324557, 2.379100105099, 2.531199743325, 2.698299769770, CaF  
     6  2.946099724689, 3.290099976017, 3.385900253881, 3.486099670194, CaF  
     7  3.594300254749, 3.704399907888, 3.873699976663, 3.911399880194, CaF  
     8  3.947199765159, 3.980099554550, 3.990399785742, 4.000000000000, CaF  
     9      58*0.0D+00/                                                 CaF  
      DATA  K_CaF/                                                      071215
     1 -1.23257460D-04, 1.44927092D-01, 2.91700770D-01, 7.40033560D-01, CaF  
     2  1.75439945D+00, 2.77175346D+00, 3.72364252D+00, 4.27487670D+00, CaF  
     3  4.76439337D+00, 5.54763237D+00, 6.20551409D+00, 6.74256339D+00, CaF  
     4  7.12500264D+00, 7.48281196D+00, 7.78623907D+00, 8.07519227D+00, CaF  
     5  8.29933112D+00, 8.51103136D+00, 8.76946980D+00, 9.02119709D+00, CaF  
     6  9.32887705D+00, 9.64336639D+00, 9.71272234D+00, 9.77998157D+00, CaF  
     7  9.85446078D+00, 9.95261741D+00, 1.02311914D+01, 1.03284573D+01, CaF  
     8  1.04344433D+01, 1.05420229D+01, 1.05773315D+01, 1.06108032D+01, CaF  
     9      58*0.0D+00/                                                 CaF  
      DATA TK_ScF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749499975714, ScF  
     2  0.826899994157, 0.919100065475, 1.024600088125, 1.095100035650, ScF  
     3  1.165799972870, 1.299700164077, 1.447999855620, 1.605000015717, ScF  
     4  1.725599918833, 1.842500088052, 1.960600039833, 2.075899831144, ScF  
     5  2.192599641606, 2.305900333125, 2.536599874662, 2.680200293854, ScF  
     6  2.823599787919, 2.982799604136, 3.194199662280, 3.335900061941, ScF  
     7  3.457699979794, 3.555600313410, 3.673600144919, 3.763400306379, ScF  
     8  3.865699794711, 3.951199854582, 3.981199579246, 4.000000000000, ScF  
     9      58*0.0D+00/                                                 ScF  
      DATA  K_ScF/                                                      071215
     1 -1.06483189D-04, 1.66214232D-01, 3.31033190D-01, 8.40042649D-01, ScF  
     2  1.99174771D+00, 3.14541648D+00, 4.22952795D+00, 4.83981237D+00, ScF  
     3  5.37591674D+00, 6.22216342D+00, 6.96170333D+00, 7.58401023D+00, ScF  
     4  7.98465990D+00, 8.33105437D+00, 8.65451414D+00, 8.95382543D+00, ScF  
     5  9.24418454D+00, 9.51429168D+00, 1.00190669D+01, 1.02933322D+01, ScF  
     6  1.05319592D+01, 1.07564959D+01, 1.09975576D+01, 1.11292794D+01, ScF  
     7  1.12270743D+01, 1.12992724D+01, 1.13897272D+01, 1.14704682D+01, ScF  
     8  1.15793203D+01, 1.16826399D+01, 1.17213051D+01, 1.17461690D+01, ScF  
     9      58*0.0D+00/                                                 ScF  
      DATA TK_MnF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.750099989774, MnF  
     2  0.828499955435, 0.922000037400, 1.028100174558, 1.101199937981, MnF  
     3  1.175599938022, 1.318700136012, 1.469799785138, 1.621599925448, MnF  
     4  1.738099826109, 1.866899962286, 2.004900112516, 2.137900272587, MnF  
     5  2.266700381350, 2.384500223295, 2.652699681998, 2.783099803235, MnF  
     6  2.905299725487, 3.255600124399, 3.371799936189, 3.489899754468, MnF  
     7  3.605999946476, 3.712800095785, 3.794100064012, 3.875200010596, MnF  
     8  3.951399858949, 3.981399583736, 4.000000000000,     59*0.0D+00/ MnF  
      DATA  K_MnF/                                                      071215
     1  1.79853405D-05, 1.52818122D-01, 3.07416977D-01, 7.81018957D-01, MnF  
     2  1.85204439D+00, 2.92572165D+00, 3.92708915D+00, 4.50821135D+00, MnF  
     3  5.02474261D+00, 5.84770967D+00, 6.53237206D+00, 7.08475305D+00, MnF  
     4  7.44208748D+00, 7.78765406D+00, 8.11527964D+00, 8.40128061D+00, MnF  
     5  8.65825259D+00, 8.87884322D+00, 9.32493076D+00, 9.50939745D+00, MnF  
     6  9.66249731D+00, 1.00101113D+01, 1.01018536D+01, 1.01867525D+01, MnF  
     7  1.02678924D+01, 1.03548044D+01, 1.04448319D+01, 1.05699449D+01, MnF  
     8  1.07268839D+01, 1.07990360D+01, 1.08464357D+01,     59*0.0D+00/ MnF  
      DATA TK_NiF/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718599832218, 0.748899961951, NiF  
     2  0.825300032879, 0.916200000552, 1.019799979305, 1.090200160982, NiF  
     3  1.161499875497, 1.298500137182, 1.444699932742, 1.591100137788, NiF  
     4  1.716699874957, 1.843600061511, 1.967299865391, 2.090900161798, NiF  
     5  2.204999922421, 2.313800155783, 2.570399709797, 2.727999702802, NiF  
     6  2.946899742299, 3.239199755153, 3.404099710472, 3.572399743880, NiF  
     7  3.708399998959, 3.837900098800, 3.936099883681, 3.974999911574, NiF  
     8  3.989099756609, 4.000000000000,     60*0.0D+00/                 NiF  
      DATA  K_NiF/                                                      071215
     1 -2.37047283D-04, 1.72816742D-01, 3.42552679D-01, 8.72913881D-01, NiF  
     2  2.07005585D+00, 3.26931471D+00, 4.39361213D+00, 5.03740435D+00, NiF  
     3  5.60738812D+00, 6.51500275D+00, 7.27541633D+00, 7.88207025D+00, NiF  
     4  8.31384921D+00, 8.69111306D+00, 9.01971279D+00, 9.32313464D+00, NiF  
     5  9.58900785D+00, 9.83309736D+00, 1.03709688D+01, 1.06695821D+01, NiF  
     6  1.10426699D+01, 1.14689629D+01, 1.16730190D+01, 1.18556196D+01, NiF  
     7  1.19884796D+01, 1.21117910D+01, 1.22147296D+01, 1.22611169D+01, NiF  
     8  1.22790460D+01, 1.22933721D+01,     60*0.0D+00/                 NiF  
      DATA TK_CuF/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, CuF  
     2  0.830399928302, 0.924999965064, 1.032300160026, 1.107200076297, CuF  
     3  1.183399906705, 1.330299941795, 1.484799903037, 1.639799867110, CuF  
     4  1.761700193143, 1.893599941848, 2.022300295187, 2.159099837437, CuF  
     5  2.272300290261, 2.382200173861, 2.608799731587, 2.844100256480, CuF  
     6  3.124899993360, 3.242699836033, 3.354000087123, 3.565699995581, CuF  
     7  3.667200002074, 3.789999978054, 3.871099917846, 3.941999647750, CuF  
     8  3.977899703198, 3.989799772324, 4.000000000000,     59*0.0D+00/ CuF  
      DATA  K_CuF/                                                      071215
     1 -9.45903239D-05, 1.62362178D-01, 3.28275650D-01, 8.32295574D-01, CuF  
     2  1.97225086D+00, 3.10718497D+00, 4.16137284D+00, 4.77859278D+00, CuF  
     3  5.32452988D+00, 6.19074812D+00, 6.90284968D+00, 7.47279765D+00, CuF  
     4  7.84819452D+00, 8.20176861D+00, 8.50789350D+00, 8.80324307D+00, CuF  
     5  9.03037982D+00, 9.23810689D+00, 9.62387722D+00, 9.95444175D+00, CuF  
     6  1.02571154D+01, 1.03590176D+01, 1.04446512D+01, 1.05946516D+01, CuF  
     7  1.06700609D+01, 1.07644875D+01, 1.08250025D+01, 1.08783741D+01, CuF  
     8  1.09081042D+01, 1.09187664D+01, 1.09283272D+01,     59*0.0D+00/ CuF  
      DATA TK_ZnF/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751300020874, ZnF  
     2  0.831499953518, 0.926899919251, 1.034900090559, 1.110700122675, ZnF  
     3  1.187900000401, 1.336900093131, 1.493899864221, 1.651200066056, ZnF  
     4  1.767900025351, 1.893899948410, 2.017600408502, 2.149199678116, ZnF  
     5  2.254000095305, 2.354800032274, 2.586300073048, 2.721000209072, ZnF  
     6  2.856599894787, 3.163099906695, 3.399899634855, 3.617099829252, ZnF  
     7  3.767100395831, 3.828199883212, 3.885700254478, 3.958000003064, ZnF  
     8  3.983899639864, 4.000000000000,     60*0.0D+00/                 ZnF  
      DATA  K_ZnF/                                                      071215
     1 -4.26109000D-06, 1.49546583D-01, 3.05357128D-01, 7.73687971D-01, ZnF  
     2  1.83098816D+00, 2.88519365D+00, 3.86375633D+00, 4.44065365D+00, ZnF  
     3  4.95203252D+00, 5.76608549D+00, 6.43931390D+00, 6.98068219D+00, ZnF  
     4  7.32062278D+00, 7.64442700D+00, 7.92976899D+00, 8.20784445D+00, ZnF  
     5  8.41540006D+00, 8.60520716D+00, 9.00219038D+00, 9.20388629D+00, ZnF  
     6  9.38335948D+00, 9.70941301D+00, 9.90176363D+00, 1.00445413D+01, ZnF  
     7  1.01267653D+01, 1.01577312D+01, 1.01876879D+01, 1.02322971D+01, ZnF  
     8  1.02525589D+01, 1.02670474D+01,     60*0.0D+00/                 ZnF  
      DATA TK_GaF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749899984889, GaF  
     2  0.827899969956, 0.920500073567, 1.025800117759, 1.097099984494, GaF  
     3  1.169300052128, 1.309499968635, 1.463799932739, 1.619599900075, GaF  
     4  1.737799833323, 1.860400126388, 1.980699989787, 2.096700307657, GaF  
     5  2.337200097840, 2.445099696393, 2.551200209311, 2.701699849929, GaF  
     6  2.820599715501, 2.941599625637, 3.120899893654, 3.301300217677, GaF  
     7  3.488599725638, 3.656199766432, 3.727899703670, 3.799500177225, GaF  
     8  3.860499675897, 3.919800096021, 3.968300232745, 3.987299716197, GaF  
     9  4.000000000000,     57*0.0D+00/                                 GaF  
      DATA  K_GaF/                                                      071215
     1  2.09658143D-04, 1.60752369D-01, 3.21512727D-01, 8.17505524D-01, GaF  
     2  1.93692690D+00, 3.05420013D+00, 4.09851102D+00, 4.69455567D+00, GaF  
     3  5.22271379D+00, 6.07347470D+00, 6.80755567D+00, 7.39641774D+00, GaF  
     4  7.76981357D+00, 8.10794220D+00, 8.40311314D+00, 8.66235694D+00, GaF  
     5  9.14866302D+00, 9.35375038D+00, 9.55199094D+00, 9.82892184D+00, GaF  
     6  1.00400731D+01, 1.02412830D+01, 1.05038456D+01, 1.07210995D+01, GaF  
     7  1.09004673D+01, 1.10247657D+01, 1.10674092D+01, 1.11035613D+01, GaF  
     8  1.11303623D+01, 1.11562677D+01, 1.11820440D+01, 1.11945011D+01, GaF  
     9  1.12038552D+01,     57*0.0D+00/                                 GaF  
      DATA TK_GeF/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718699829691, 0.748999964245, GeF  
     2  0.825700023199, 0.916500007268, 1.020099976997, 1.091800120057, GeF  
     3  1.164799950225, 1.304100083550, 1.448399846272, 1.602399951027, GeF  
     4  1.705699948431, 1.803400098266, 1.979900004013, 2.134000182949, GeF  
     5  2.218500220157, 2.305300319151, 2.390900275352, 2.471600135663, GeF  
     6  2.650699632984, 2.819399743864, 3.016100362951, 3.207199959560, GeF  
     7  3.356799875711, 3.500800008724, 3.649699642957, 3.765300352314, GeF  
     8  3.871999938206, 3.923900187807, 3.963700129656, 3.986199691501, GeF  
     9  4.000000000000,     57*0.0D+00/                                 GeF  
      DATA  K_GeF/                                                      071215
     1 -2.32012964D-04, 1.46134720D-01, 2.91268622D-01, 7.40208834D-01, GeF  
     2  1.75928893D+00, 2.77849642D+00, 3.73981388D+00, 4.30263963D+00, GeF  
     3  4.80455841D+00, 5.60178808D+00, 6.25801680D+00, 6.82214247D+00, GeF  
     4  7.14277309D+00, 7.41363402D+00, 7.84361976D+00, 8.17814996D+00, GeF  
     5  8.35434299D+00, 8.53531769D+00, 8.71657045D+00, 8.89037400D+00, GeF  
     6  9.27641029D+00, 9.61774936D+00, 9.96175960D+00, 1.02320763D+01, GeF  
     7  1.04058299D+01, 1.05501514D+01, 1.06809174D+01, 1.07680278D+01, GeF  
     8  1.08329424D+01, 1.08590928D+01, 1.08783253D+01, 1.08897867D+01, GeF  
     9  1.08973612D+01,     57*0.0D+00/                                 GeF  
      DATA TK_AsF/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751100015691, AsF  
     2  0.831299948934, 0.926099938541, 1.033700122621, 1.110600125271, AsF  
     3  1.189800039961, 1.265100096992, 1.346300007541, 1.424199969127, AsF  
     4  1.505400119142, 1.653500004857, 1.866699967335, 2.029099793163, AsF  
     5  2.206599964845, 2.309200409984, 2.407599791865, 2.615599772890, AsF  
     6  2.741399804181, 2.868099845377, 3.171600108998, 3.375000001862, AsF  
     7  3.575599820329, 3.648099755985, 3.725799849187, 3.797700139487, AsF  
     8  3.866499812990, 3.946399747096, 3.979399595417, 4.000000000000, AsF  
     9      58*0.0D+00/                                                 AsF  
      DATA  K_AsF/                                                      071215
     1 -1.23998712D-04, 1.58665174D-01, 3.24080048D-01, 8.18107268D-01, AsF  
     2  1.93984671D+00, 3.05033009D+00, 4.08280144D+00, 4.70092420D+00, AsF  
     3  5.25238888D+00, 5.70961226D+00, 6.14147103D+00, 6.50581556D+00, AsF  
     4  6.84198406D+00, 7.36201988D+00, 7.94772479D+00, 8.29982002D+00, AsF  
     5  8.62745696D+00, 8.80107737D+00, 8.96019486D+00, 9.27213714D+00, AsF  
     6  9.44121121D+00, 9.59451271D+00, 9.89488753D+00, 1.00515686D+01, AsF  
     7  1.01864396D+01, 1.02355437D+01, 1.02893535D+01, 1.03380213D+01, AsF  
     8  1.03793425D+01, 1.04145825D+01, 1.04241031D+01, 1.04285802D+01, AsF  
     9      58*0.0D+00/                                                 AsF  
      DATA TK_SeF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.749999987182, SeF  
     2  0.828199962695, 0.921100059100, 1.026700139985, 1.099699917991, SeF  
     3  1.173999975153, 1.317000101227, 1.467199849098, 1.618399930081, SeF  
     4  1.739399794848, 1.866299977434, 2.014700340246, 2.157799806693, SeF  
     5  2.274400141279, 2.377700071690, 2.613399722247, 2.793800064108, SeF  
     6  3.020700402838, 3.169500062525, 3.336400074413, 3.496599910722, SeF  
     7  3.653199696431, 3.842700205281, 3.963600127415, 4.000000000000, SeF  
     8      62*0.0D+00/                                                 SeF  
      DATA  K_SeF/                                                      071215
     1  2.87672005D-05, 1.56653572D-01, 3.15115094D-01, 7.98985799D-01, SeF  
     2  1.89398949D+00, 2.98772256D+00, 4.00984877D+00, 4.60484400D+00, SeF  
     3  5.13338241D+00, 5.97502115D+00, 6.67054665D+00, 7.23194552D+00, SeF  
     4  7.60882759D+00, 7.95333107D+00, 8.30803876D+00, 8.61604668D+00, SeF  
     5  8.85044784D+00, 9.04856499D+00, 9.46490979D+00, 9.74278041D+00, SeF  
     6  1.00462123D+01, 1.02257982D+01, 1.04129078D+01, 1.05795050D+01, SeF  
     7  1.07332066D+01, 1.09114598D+01, 1.10216629D+01, 1.10547787D+01, SeF  
     8      62*0.0D+00/                                                 SeF  
      DATA TK_BrF/                                                      071215
     1  0.699999789529, 0.709900047001, 0.721099823174, 0.754500103807, BrF  
     2  0.789599963486, 0.839700141495, 0.942999965670, 1.059199970348, BrF  
     3  1.135499987006, 1.211299986309, 1.354200008725, 1.506700149709, BrF  
     4  1.661699871994, 1.774100071459, 1.893899948410, 2.016500382612, BrF  
     5  2.148699713048, 2.252100051068, 2.354000091751, 2.562900206136, BrF  
     6  2.728499666640, 2.891200260505, 3.213100095152, 3.383900203235, BrF  
     7  3.457799982212, 3.528199860636, 3.649799635892, 3.702899873736, BrF  
     8  3.757800172040, 3.827199859036, 3.906799772360, 3.962200096039, BrF  
     9  3.985499675785, 4.000000000000,     56*0.0D+00/                 BrF  
      DATA  K_BrF/                                                      071215
     1 -9.77214354D-05, 1.76012825D-01, 3.70861947D-01, 9.25498518D-01, BrF  
     2  1.46815449D+00, 2.17755046D+00, 3.43049219D+00, 4.56469068D+00, BrF  
     3  5.18328584D+00, 5.71736483D+00, 6.54981153D+00, 7.24830122D+00, BrF  
     4  7.81593944D+00, 8.16272314D+00, 8.48793271D+00, 8.78502103D+00, BrF  
     5  9.07607200D+00, 9.28826161D+00, 9.48656819D+00, 9.85904697D+00, BrF  
     6  1.01162292D+01, 1.03331153D+01, 1.06741462D+01, 1.08237803D+01, BrF  
     7  1.08827470D+01, 1.09335082D+01, 1.09962859D+01, 1.10075037D+01, BrF  
     8  1.10058989D+01, 1.09853693D+01, 1.09437141D+01, 1.09100825D+01, BrF  
     9  1.08961206D+01, 1.08877071D+01,     56*0.0D+00/                 BrF  
      DATA TK_RbF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750299994957, RbF  
     2  0.828999943334, 0.922500025344, 1.028700189375, 1.101299940286, RbF  
     3  1.174899954267, 1.316400088950, 1.467399844178, 1.622299940924, RbF  
     4  1.754400098048, 1.891399893728, 2.020800405928, 2.155099742839, RbF  
     5  2.397699775512, 2.582099974647, 2.761200262379, 3.090700156032, RbF  
     6  3.339300146752, 3.449799789433, 3.557300348888, 3.639500317883, RbF  
     7  3.683700039557, 3.721700133292, 3.790799994826, 3.827299861454, RbF  
     8  3.861899707886, 3.963400122932, 4.000000000000,     59*0.0D+00/ RbF  
      DATA  K_RbF/                                                      071215
     1 -6.87031206D-06, 1.56134416D-01, 3.15717755D-01, 8.01019535D-01, RbF  
     2  1.89841988D+00, 2.99348399D+00, 4.01501411D+00, 4.60304518D+00, RbF  
     3  5.12394110D+00, 5.95425234D+00, 6.65200307D+00, 7.22433895D+00, RbF  
     4  7.63049535D+00, 7.99397504D+00, 8.29650720D+00, 8.57718395D+00, RbF  
     5  9.01293973D+00, 9.28969854D+00, 9.51687435D+00, 9.84056781D+00, RbF  
     6  1.00184643D+01, 1.00821979D+01, 1.01384903D+01, 1.01848882D+01, RbF  
     7  1.02167314D+01, 1.02524405D+01, 1.03508533D+01, 1.04263360D+01, RbF  
     8  1.05143036D+01, 1.08440251D+01, 1.09769250D+01,     59*0.0D+00/ RbF  
      DATA TK_SrF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.749999987182, SrF  
     2  0.828299960275, 0.921700044633, 1.027800167149, 1.100199914928, SrF  
     3  1.173699982115, 1.314900058258, 1.465399893379, 1.617899942583, SrF  
     4  1.744899881733, 1.884500010674, 2.013300307295, 2.141900188115, SrF  
     5  2.332999998281, 2.536099862502, 2.725799861915, 2.892600156427, SrF  
     6  3.073899792186, 3.266800383713, 3.364799761369, 3.464300128442, SrF  
     7  3.573099760603, 3.673100133017, 3.769000441766, 3.855299995839, SrF  
     8  3.899899627020, 3.940599616140, 3.977399739125, 4.000000000000, SrF  
     9      58*0.0D+00/                                                 SrF  
      DATA  K_SrF/                                                      071215
     1 -1.22560561D-04, 1.44810709D-01, 2.91466189D-01, 7.39439766D-01, SrF  
     2  1.75539938D+00, 2.77554633D+00, 3.72969364D+00, 4.27941471D+00, SrF  
     3  4.76825755D+00, 5.54990425D+00, 6.21002272D+00, 6.74932146D+00, SrF  
     4  7.12697763D+00, 7.48788108D+00, 7.78384570D+00, 8.05320218D+00, SrF  
     5  8.41589319D+00, 8.75504114D+00, 9.02548530D+00, 9.22605860D+00, SrF  
     6  9.40843455D+00, 9.56829619D+00, 9.63837930D+00, 9.70456433D+00, SrF  
     7  9.77990059D+00, 9.87031693D+00, 9.99902538D+00, 1.01654276D+01, SrF  
     8  1.02726536D+01, 1.03823761D+01, 1.04895792D+01, 1.05581874D+01, SrF  
     9      58*0.0D+00/                                                 SrF  
      DATA TK_YF/                                                       071215
     1  0.699999789529, 0.709300031396, 0.718599832218, 0.748899961951, YF   
     2  0.825500028039, 0.916300002791, 1.019999974527, 1.091700122615, YF   
     3  1.164699947961, 1.304700070782, 1.449099829913, 1.596399991011, YF   
     4  1.723799877081, 1.844800032558, 2.141300230033, 2.245699910869, YF   
     5  2.359199705145, 2.461400057796, 2.555700313233, 2.771500362407, YF   
     6  2.981499574917, 3.160199836084, 3.322499752151, 3.461600070027, YF   
     7  3.604100083980, 3.673000130637, 3.740899784763, 3.888500320939, YF   
     8  3.955699952843, 3.983199624148, 4.000000000000,     59*0.0D+00/ YF   
      DATA  K_YF/                                                       071215
     1  2.20958065D-04, 1.64145348D-01, 3.24942248D-01, 8.27479288D-01, YF   
     2  1.96524807D+00, 3.10191620D+00, 4.17133854D+00, 4.79457721D+00, YF   
     3  5.34856077D+00, 6.22771711D+00, 6.94191562D+00, 7.52572304D+00, YF   
     4  7.94528523D+00, 8.28978944D+00, 8.99053895D+00, 9.21037821D+00, YF   
     5  9.44386345D+00, 9.65145806D+00, 9.83971831D+00, 1.02427511D+01, YF   
     6  1.05745317D+01, 1.08040762D+01, 1.09759858D+01, 1.11057504D+01, YF   
     7  1.12449009D+01, 1.13243445D+01, 1.14132030D+01, 1.16273451D+01, YF   
     8  1.17214323D+01, 1.17576071D+01, 1.17789372D+01,     59*0.0D+00/ YF   
      DATA TK_AgF/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719799801886, 0.751600028649, AgF  
     2  0.832399974150, 0.928399883083, 1.037100031778, 1.113200057782, AgF  
     3  1.190800024619, 1.340600149291, 1.499299977487, 1.659199853190, AgF  
     4  1.777600159335, 1.905699938960, 2.023400213978, 2.145299950581, AgF  
     5  2.252700065038, 2.355799957926, 2.606699884900, 2.742799833972, AgF  
     6  2.879600096715, 3.108799744043, 3.354900019169, 3.502000036700, AgF  
     7  3.671800102074, 3.738199727557, 3.804100277582, 3.859599692674, AgF  
     8  3.909899841977, 3.965600172236, 3.986499698236, 4.000000000000, AgF  
     9      58*0.0D+00/                                                 AgF  
      DATA  K_AgF/                                                      071215
     1 -1.01563759D-04, 1.62089778D-01, 3.30896766D-01, 8.34651912D-01, AgF  
     2  1.97497687D+00, 3.10740200D+00, 4.15516188D+00, 4.76931865D+00, AgF  
     3  5.31291011D+00, 6.17461364D+00, 6.88604091D+00, 7.45640721D+00, AgF  
     4  7.81119218D+00, 8.14782030D+00, 8.42459226D+00, 8.68651278D+00, AgF  
     5  8.90054181D+00, 9.09277881D+00, 9.50565554D+00, 9.69527959D+00, AgF  
     6  9.86183742D+00, 1.00925198D+01, 1.02847863D+01, 1.03765429D+01, AgF  
     7  1.04600624D+01, 1.04855250D+01, 1.05074169D+01, 1.05248486D+01, AgF  
     8  1.05424302D+01, 1.05689278D+01, 1.05823744D+01, 1.05925085D+01, AgF  
     9      58*0.0D+00/                                                 AgF  
      DATA TK_CdF/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751300020874, CdF  
     2  0.831599955811, 0.926999916840, 1.035100085215, 1.110500127866, CdF  
     3  1.187299987908, 1.335400058736, 1.492199828564, 1.650300090003, CdF  
     4  1.771800013712, 1.905499943973, 2.029399771015, 2.156099766488, CdF  
     5  2.263900320189, 2.363899740154, 2.619899871875, 2.776200009789, CdF  
     6  2.917800021877, 3.341000186705, 3.583300002761, 3.685699900261, CdF  
     7  3.787099909682, 3.857699826631, 3.922100147816, 3.969400257397, CdF  
     8  3.987599722932, 4.000000000000,     60*0.0D+00/                 CdF  
      DATA  K_CdF/                                                      071215
     1 -3.81281666D-05, 1.47580653D-01, 3.01382621D-01, 7.63705089D-01, CdF  
     2  1.80882523D+00, 2.84981641D+00, 3.81734334D+00, 4.38461343D+00, CdF  
     3  4.88807761D+00, 5.69020181D+00, 6.35791781D+00, 6.89858262D+00, CdF  
     4  7.24951427D+00, 7.58863165D+00, 7.87029735D+00, 8.13450053D+00, CdF  
     5  8.34433808D+00, 8.52788047D+00, 8.94621742D+00, 9.16198953D+00, CdF  
     6  9.33188451D+00, 9.72405654D+00, 9.89657594D+00, 9.96289808D+00, CdF  
     7  1.00274376D+01, 1.00746888D+01, 1.01237651D+01, 1.01675875D+01, CdF  
     8  1.01872842D+01, 1.02019062D+01,     60*0.0D+00/                 CdF  
      DATA TK_InF/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, InF  
     2  0.829599928813, 0.923599998821, 1.030200216135, 1.103900000223, InF  
     3  1.178699866081, 1.322900096580, 1.475599911347, 1.630100108674, InF  
     4  1.754400098048, 1.887399932866, 2.007600174514, 2.130000091012, InF  
     5  2.320999732286, 2.491999822273, 2.655599753068, 2.788499940050, InF  
     6  2.914399938752, 3.008100180750, 3.101800250011, 3.335700056952, InF  
     7  3.561500297432, 3.668500028611, 3.762400282203, 3.854400059292, InF  
     8  3.929800318889, 3.972800069653, 3.988499743138, 4.000000000000, InF  
     9      58*0.0D+00/                                                 InF  
      DATA  K_InF/                                                      071215
     1  2.24697975D-04, 1.60741390D-01, 3.23039508D-01, 8.18110323D-01, InF  
     2  1.93888218D+00, 3.05614996D+00, 4.09535527D+00, 4.69922332D+00, InF  
     3  5.23343617D+00, 6.08410795D+00, 6.79068539D+00, 7.36168867D+00, InF  
     4  7.74576330D+00, 8.10264441D+00, 8.38943507D+00, 8.65563295D+00, InF  
     5  9.03165896D+00, 9.33177721D+00, 9.58513698D+00, 9.76849075D+00, InF  
     6  9.93032017D+00, 1.00481122D+01, 1.01658932D+01, 1.04533662D+01, InF  
     7  1.06923778D+01, 1.07831475D+01, 1.08478927D+01, 1.08984157D+01, InF  
     8  1.09369269D+01, 1.09630105D+01, 1.09741105D+01, 1.09829318D+01, InF  
     9      58*0.0D+00/                                                 InF  
      DATA TK_SnF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749799982595, SnF  
     2  0.827699974796, 0.920100083212, 1.025300105412, 1.097299979378, SnF  
     3  1.170400058697, 1.311099980503, 1.460500013920, 1.611600100113, SnF  
     4  1.725899925791, 1.847199974650, 1.968899823733, 2.094700257361, SnF  
     5  2.269200435959, 2.441799624595, 2.551900225477, 2.660199865319, SnF  
     6  2.853700109956, 2.946199726891, 3.033299823225, 3.263600306875, SnF  
     7  3.401599659912, 3.530699745715, 3.642800130391, 3.796300110135, SnF  
     8  3.878700089773, 3.956799976862, 3.983299626393, 4.000000000000, SnF  
     9      58*0.0D+00/                                                 SnF  
      DATA  K_SnF/                                                      071215
     1 -1.34264558D-04, 1.46249468D-01, 2.92859192D-01, 7.43963384D-01, SnF  
     2  1.76548708D+00, 2.78660196D+00, 3.74518517D+00, 4.29953089D+00, SnF  
     3  4.79264578D+00, 5.58257432D+00, 6.24705926D+00, 6.78886949D+00, SnF  
     4  7.13556925D+00, 7.45940064D+00, 7.74981163D+00, 8.02319602D+00, SnF  
     5  8.36995451D+00, 8.68221436D+00, 8.86534310D+00, 9.03448373D+00, SnF  
     6  9.32378994D+00, 9.46374169D+00, 9.59838483D+00, 9.95250449D+00, SnF  
     7  1.01464726D+01, 1.03088712D+01, 1.04335999D+01, 1.05741174D+01, SnF  
     8  1.06321330D+01, 1.06788528D+01, 1.06946862D+01, 1.07051896D+01, SnF  
     9      58*0.0D+00/                                                 SnF  
      DATA TK_SbF/                                                      071215
     1  0.699999789529, 0.710000049601, 0.721199825569, 0.754800111582, SbF  
     2  0.840400138721, 0.943599950938, 1.059699959162, 1.135999974760, SbF  
     3  1.211699977747, 1.354700019820, 1.507800175574, 1.663599916801, SbF  
     4  1.773400053884, 1.892499917788, 2.006400146959, 2.116299806018, SbF  
     5  2.221100172084, 2.320799727955, 2.541499990993, 2.640100333090, SbF  
     6  2.739799770065, 3.006400142815, 3.222800054127, 3.367799840915, SbF  
     7  3.488399721202, 3.601100301092, 3.726799779893, 3.843500223279, SbF  
     8  3.912099898179, 3.972700076838, 4.000000000000,     59*0.0D+00/ SbF  
      DATA  K_SbF/                                                      071215
     1 -1.56792362D-04, 1.65472968D-01, 3.46896106D-01, 8.66484136D-01, SbF  
     2  2.03831643D+00, 3.20606340D+00, 4.26606261D+00, 4.84622354D+00, SbF  
     3  5.34763551D+00, 6.13358008D+00, 6.79858287D+00, 7.34277070D+00, SbF  
     4  7.66820592D+00, 7.98096137D+00, 8.25025683D+00, 8.48937912D+00, SbF  
     5  8.70270631D+00, 8.89378759D+00, 9.26852681D+00, 9.40901113D+00, SbF  
     6  9.53248569D+00, 9.78384163D+00, 9.93128347D+00, 1.00149004D+01, SbF  
     7  1.00823897D+01, 1.01499617D+01, 1.02338514D+01, 1.03151494D+01, SbF  
     8  1.03608124D+01, 1.03992686D+01, 1.04164969D+01,     59*0.0D+00/ SbF  
      DATA TK_IF/                                                       071215
     1  0.699999789529, 0.710000049601, 0.721399830359, 0.755200121949, IF   
     2  0.790799970587, 0.841500112180, 0.946099889558, 1.063600035742, IF   
     3  1.141199903243, 1.218399834333, 1.364200029675, 1.520499979728, IF   
     4  1.680299844609, 1.798400145363, 1.926599910252, 2.044300067839, IF   
     5  2.169600048856, 2.278499850407, 2.381600160965, 2.647599789711, IF   
     6  2.796800132045, 2.953299886563, 3.203399874403, 3.373299966973, IF   
     7  3.442799643719, 3.513100295602, 3.642100179841, 3.690899619400, IF   
     8  3.743199835924, 3.815600002625, 3.905399740920, 3.962400100522, IF   
     9  3.985499675785, 4.000000000000,     56*0.0D+00/                 IF   
      DATA  K_IF/                                                       071215
     1 -1.18501366D-04, 1.76472979D-01, 3.73288015D-01, 9.29968954D-01, IF   
     2  1.47536905D+00, 2.18602075D+00, 3.43983916D+00, 4.57087348D+00, IF   
     3  5.19011337D+00, 5.72449300D+00, 6.55678788D+00, 7.25610295D+00, IF   
     4  7.82615956D+00, 8.18016280D+00, 8.51727230D+00, 8.79454524D+00, IF   
     5  9.06498341D+00, 9.28416942D+00, 9.47977504D+00, 9.92718257D+00, IF   
     6  1.01376784D+01, 1.03270661D+01, 1.05716105D+01, 1.07055976D+01, IF   
     7  1.07538887D+01, 1.07971278D+01, 1.08452313D+01, 1.08447928D+01, IF   
     8  1.08290256D+01, 1.07814930D+01, 1.06937011D+01, 1.06325851D+01, IF   
     9  1.06086383D+01, 1.05941888D+01,     56*0.0D+00/                 IF   
      DATA TK_CsF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749899984889, CsF  
     2  0.827999967535, 0.920800066334, 1.026200127637, 1.097799966589, CsF  
     3  1.170300061017, 1.309599966507, 1.458900002313, 1.612100087610, CsF  
     4  1.738899806871, 1.869499896645, 1.993699856830, 2.125099988240, CsF  
     5  2.402299662561, 2.608399760790, 2.802100248741, 3.080699958547, CsF  
     6  3.322299747815, 3.417700032405, 3.514900337638, 3.607599830683, CsF  
     7  3.686199865437, 3.754600096201, 3.788399940332, 3.820399694639, CsF  
     8  3.920600114490, 3.967300210335, 4.000000000000,     59*0.0D+00/ CsF  
      DATA  K_CsF/                                                      071215
     1  2.27607029D-04, 1.56064128D-01, 3.12120952D-01, 7.93663949D-01, CsF  
     2  1.88215754D+00, 2.96990736D+00, 3.98601292D+00, 4.56818504D+00, CsF  
     3  5.08419530D+00, 5.90822902D+00, 6.60533722D+00, 7.17793498D+00, CsF  
     4  7.57327918D+00, 7.92594638D+00, 8.22193847D+00, 8.50142981D+00, CsF  
     5  8.99834448D+00, 9.29813967D+00, 9.53173988D+00, 9.79676769D+00, CsF  
     6  9.97116823D+00, 1.00291138D+01, 1.00863219D+01, 1.01492644D+01, CsF  
     7  1.02279429D+01, 1.03385635D+01, 1.04135554D+01, 1.04975592D+01, CsF  
     8  1.08227063D+01, 1.09900506D+01, 1.11077269D+01,     59*0.0D+00/ CsF  
      DATA TK_BaF/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718799827163, 0.749299971126, BaF  
     2  0.826300008678, 0.918100043088, 1.022500036265, 1.163299916258, BaF  
     3  1.298800143905, 1.444399939753, 1.595800007627, 1.749599978912, BaF  
     4  1.918300034102, 2.068099824858, 2.215900169637, 2.379300109872, BaF  
     5  2.539699950060, 2.680600264902, 2.815800000943, 3.001100024546, BaF  
     6  3.164099931043, 3.236099679863, 3.307000343120, 3.404399716539, BaF  
     7  3.493299833704, 3.586000066019, 3.653099694098, 3.724399946199, BaF  
     8  3.790199982247, 3.888400318566, 3.949699821606, 4.000000000000, BaF  
     9      58*0.0D+00/                                                 BaF  
      DATA  K_BaF/                                                      071215
     1  4.60585425D-05, 1.42997526D-01, 2.86243145D-01, 7.27484568D-01, BaF  
     2  1.72619685D+00, 2.73136884D+00, 3.67546539D+00, 4.68908786D+00, BaF  
     3  5.45396252D+00, 6.10925983D+00, 6.65920362D+00, 7.12071835D+00, BaF  
     4  7.54862552D+00, 7.88072783D+00, 8.17627540D+00, 8.47182703D+00, BaF  
     5  8.73057455D+00, 8.93124133D+00, 9.10049277D+00, 9.29821958D+00, BaF  
     6  9.44450340D+00, 9.50374691D+00, 9.56258983D+00, 9.65493848D+00, BaF  
     7  9.76631277D+00, 9.92025998D+00, 1.00518092D+01, 1.02016038D+01, BaF  
     8  1.03441355D+01, 1.05647461D+01, 1.07093530D+01, 1.08311756D+01, BaF  
     9      58*0.0D+00/                                                 BaF  
      DATA TK_LaF/                                                      071215
     1  0.699999789529, 0.709100026195, 0.718099844857, 0.747599932133, LaF  
     2  0.821900115164, 0.909999861751, 1.010900191932, 1.146400017873, LaF  
     3  1.276700142890, 1.415099967672, 1.557699994021, 1.694499955839, LaF  
     4  1.839900146018, 1.975299906939, 2.114799775509, 2.250500013816, LaF  
     5  2.413599937721, 2.722100129515, 2.879200088225, 3.054400294808, LaF  
     6  3.178900265485, 3.364799761369, 3.494299857043, 3.627600068656, LaF  
     7  3.731499589142, 3.830199931085, 3.929000301115, 4.000000000000, LaF  
     8      62*0.0D+00/                                                 LaF  
      DATA  K_LaF/                                                      071215
     1 -7.79595630D-05, 1.59754292D-01, 3.14917625D-01, 8.03296909D-01, LaF  
     2  1.90841116D+00, 3.01833418D+00, 4.07102137D+00, 5.19793741D+00, LaF  
     3  6.04724921D+00, 6.76480160D+00, 7.35899356D+00, 7.82856348D+00, LaF  
     4  8.24984267D+00, 8.59008904D+00, 8.90335548D+00, 9.18191090D+00, LaF  
     5  9.49016805D+00, 1.00157408D+01, 1.02683933D+01, 1.05505100D+01, LaF  
     6  1.07577163D+01, 1.10791892D+01, 1.13101478D+01, 1.15583607D+01, LaF  
     7  1.17661178D+01, 1.19803088D+01, 1.22114265D+01, 1.23843531D+01, LaF  
     8      62*0.0D+00/                                                 LaF  
      DATA TK_HoF/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718499834746, 0.748699957364, HoF  
     2  0.824800044980, 0.915499984881, 1.018900000807, 1.088600131665, HoF  
     3  1.158999867100, 1.293600027361, 1.436299958586, 1.584300024885, HoF  
     4  1.733699931917, 1.893099930912, 2.047000129646, 2.198499759798, HoF  
     5  2.380700141621, 2.565600013971, 2.740799791413, 2.958700010972, HoF  
     6  3.068999769081, 3.174700175452, 3.264900338091, 3.352100230582, HoF  
     7  3.536399875431, 3.617399836531, 3.698099768392, 3.772700270556, HoF  
     8  3.879300103347, 3.951899869867, 3.981799592717, 4.000000000000, HoF  
     9      58*0.0D+00/                                                 HoF  
      DATA  K_HoF/                                                      071215
     1  3.78917885D-05, 1.75086162D-01, 3.48665580D-01, 8.89239982D-01, HoF  
     2  2.10887586D+00, 3.33287240D+00, 4.48041774D+00, 5.13237336D+00, HoF  
     3  5.70852213D+00, 6.62301099D+00, 7.38632015D+00, 8.01607957D+00, HoF  
     4  8.53247400D+00, 8.99037177D+00, 9.36922823D+00, 9.70105608D+00, HoF  
     5  1.00615584D+01, 1.03872441D+01, 1.06543028D+01, 1.09287831D+01, HoF  
     6  1.10459854D+01, 1.11478372D+01, 1.12302508D+01, 1.13113459D+01, HoF  
     7  1.15282128D+01, 1.16649423D+01, 1.18344749D+01, 1.20159654D+01, HoF  
     8  1.22913321D+01, 1.24706999D+01, 1.25399140D+01, 1.25804733D+01, HoF  
     9      58*0.0D+00/                                                 HoF  
      DATA TK_YbF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750199992366, YbF  
     2  0.828799948174, 0.922200032577, 1.028200177028, 1.100799928760, YbF  
     3  1.174499963549, 1.316200084858, 1.466899856478, 1.619699897575, YbF  
     4  1.750099989702, 1.892999928725, 2.024400140151, 2.155199745204, YbF  
     5  2.343500242938, 2.538799928170, 2.719300265623, 2.905799736730, YbF  
     6  3.115999785654, 3.328199875740, 3.445099691596, 3.556100323845, YbF  
     7  3.648099755985, 3.727399738317, 3.799400175128, 3.834400022566, YbF  
     8  3.869699886107, 3.950499839297, 3.981999597207, 4.000000000000, YbF  
     9      58*0.0D+00/                                                 YbF  
      DATA  K_YbF/                                                      071215
     1  7.57268190D-05, 1.45943538D-01, 2.95047798D-01, 7.47210638D-01, YbF  
     2  1.77278303D+00, 2.79819666D+00, 3.75628653D+00, 4.31017647D+00, YbF  
     3  4.80244811D+00, 5.58957591D+00, 6.25227547D+00, 6.79367751D+00, YbF  
     4  7.18130314D+00, 7.54944238D+00, 7.84981490D+00, 8.12214642D+00, YbF  
     5  8.47746357D+00, 8.80238340D+00, 9.06034542D+00, 9.28405631D+00, YbF  
     6  9.48981717D+00, 9.65679307D+00, 9.73422159D+00, 9.80043264D+00, YbF  
     7  9.85683504D+00, 9.92315324D+00, 1.00217375D+01, 1.00897296D+01, YbF  
     8  1.01730192D+01, 1.04107167D+01, 1.05150554D+01, 1.05761259D+01, YbF  
     9      58*0.0D+00/                                                 YbF  
      DATA TK_LuF/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749899984889, LuF  
     2  0.827899969956, 0.920600071156, 1.026000122698, 1.098299953800, LuF  
     3  1.171700028528, 1.312800015288, 1.462199972099, 1.614000040101, LuF  
     4  1.739199799657, 1.872899955765, 2.003200073480, 2.136200233514, LuF  
     5  2.231799590973, 2.325099821065, 2.509100223834, 2.658999836392, LuF  
     6  2.884900225014, 3.023600197791, 3.144600006407, 3.270900396060, LuF  
     7  3.521600340426, 3.613199734631, 3.704599912441, 3.787699923828, LuF  
     8  3.869199874682, 3.946799756128, 3.979699573861, 4.000000000000, LuF  
     9      58*0.0D+00/                                                 LuF  
      DATA  K_LuF/                                                      071215
     1  2.09607965D-04, 1.66385915D-01, 3.32776752D-01, 8.46067430D-01, LuF  
     2  2.00405327D+00, 3.16010151D+00, 4.23922504D+00, 4.86189837D+00, LuF  
     3  5.41344216D+00, 6.29012724D+00, 7.01847339D+00, 7.60843369D+00, LuF  
     4  8.01307074D+00, 8.38564075D+00, 8.70520858D+00, 8.99961546D+00, LuF  
     5  9.19651121D+00, 9.37882082D+00, 9.71043382D+00, 9.95107786D+00, LuF  
     6  1.02697296D+01, 1.04512036D+01, 1.06063312D+01, 1.07647630D+01, LuF  
     7  1.10573831D+01, 1.11529273D+01, 1.12406436D+01, 1.13153070D+01, LuF  
     8  1.13891137D+01, 1.14674071D+01, 1.15044501D+01, 1.15285671D+01, LuF  
     9      58*0.0D+00/                                                 LuF  
      DATA TK_HgF/                                                      071215
     1  0.699999789529, 0.710000049601, 0.721199825569, 0.754800111582, HgF  
     2  0.840500136308, 0.938600012051, 1.060099954764, 1.199499812491, HgF  
     3  1.347799970238, 1.501900036845, 1.665399959250, 1.842200095290, HgF  
     4  2.036399882037, 2.205399933027, 2.360999669894, 2.557600357111, HgF  
     5  2.741999816949, 3.131000143404, 3.286199886513, 3.537399898188, HgF  
     6  3.669800055148, 3.825699822772, 3.883100192764, 3.936399862059, HgF  
     7  3.975199897204, 3.989199758854, 4.000000000000,     63*0.0D+00/ HgF  
      DATA  K_HgF/                                                      071215
     1 -4.85334699D-05, 1.56321555D-01, 3.27624777D-01, 8.18370456D-01, HgF  
     2  1.92724524D+00, 2.98300080D+00, 4.04018419D+00, 4.99343020D+00, HgF  
     3  5.78464469D+00, 6.43320309D+00, 6.98589429D+00, 7.47603493D+00, HgF  
     4  7.92900158D+00, 8.27445534D+00, 8.56152097D+00, 8.88287522D+00, HgF  
     5  9.14004150D+00, 9.55078537D+00, 9.67362030D+00, 9.83034501D+00, HgF  
     6  9.89339993D+00, 9.96121004D+00, 9.98851054D+00, 1.00173103D+01, HgF  
     7  1.00419108D+01, 1.00519245D+01, 1.00601968D+01,     63*0.0D+00/ HgF  
      DATA TK_TlF/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751100015691, TlF  
     2  0.831099944349, 0.926199936130, 1.034000114605, 1.108900115487, TlF  
     3  1.185099942101, 1.331999980775, 1.487699835769, 1.645199984775, TlF  
     4  1.772300026265, 1.910099833623, 2.034699840781, 2.161399886450, TlF  
     5  2.374599997713, 2.651299647689, 2.867699835360, 3.097700324698, TlF  
     6  3.251800041461, 3.397799788187, 3.484599636928, 3.566299952459, TlF  
     7  3.750900008512, 3.826099832443, 3.901199646599, 3.960700062423, TlF  
     8  3.984999664560, 4.000000000000,     60*0.0D+00/                 TlF  
      DATA  K_TlF/                                                      071215
     1  5.75531567D-05, 1.59965348D-01, 3.26543197D-01, 8.24028446D-01, TlF  
     2  1.95092803D+00, 3.07257126D+00, 4.11332288D+00, 4.71992496D+00, TlF  
     3  5.25678400D+00, 6.10941899D+00, 6.81625494D+00, 7.38635369D+00, TlF  
     4  7.77076287D+00, 8.13244233D+00, 8.42291242D+00, 8.69139233D+00, TlF  
     5  9.09407351D+00, 9.52981843D+00, 9.80152443D+00, 1.00310946D+01, TlF  
     6  1.01580230D+01, 1.02671568D+01, 1.03313470D+01, 1.03931220D+01, TlF  
     7  1.05345870D+01, 1.05905847D+01, 1.06493630D+01, 1.07059783D+01, TlF  
     8  1.07340366D+01, 1.07532548D+01,     60*0.0D+00/                 TlF  
      DATA TK_PbF/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, PbF  
     2  0.830599932887, 0.925299957830, 1.032700149339, 1.106900069382, PbF  
     3  1.182199881719, 1.327399994117, 1.482099965666, 1.638399901975, PbF  
     4  1.771300001158, 1.916599992539, 2.044400070128, 2.178100259117, PbF  
     5  2.346700314915, 2.518800431664, 2.677800258515, 2.870699907819, PbF  
     6  3.113899740732, 3.211300052855, 3.300700204473, 3.419000066612, PbF  
     7  3.505900127623, 3.577099856164, 3.643100109198, 3.777499923203, PbF  
     8  3.867399833554, 3.941299631945, 3.977499731940, 4.000000000000, PbF  
     9      58*0.0D+00/                                                 PbF  
      DATA  K_PbF/                                                      071215
     1 -8.68244515D-05, 1.40529246D-01, 2.85616475D-01, 7.22215953D-01, PbF  
     2  1.71273918D+00, 2.70192673D+00, 3.62477250D+00, 4.16248810D+00, PbF  
     3  4.63976431D+00, 5.40444214D+00, 6.04928258D+00, 6.57529233D+00, PbF  
     4  6.95207357D+00, 7.31102852D+00, 7.59344700D+00, 7.86427038D+00, PbF  
     5  8.17657696D+00, 8.46245205D+00, 8.69455093D+00, 8.93383635D+00, PbF  
     6  9.17599873D+00, 9.25808562D+00, 9.32909389D+00, 9.42398476D+00, PbF  
     7  9.50130577D+00, 9.57265976D+00, 9.64521077D+00, 9.80300454D+00, PbF  
     8  9.90767372D+00, 9.99287311D+00, 1.00363087D+01, 1.00646725D+01, PbF  
     9      58*0.0D+00/                                                 PbF  
      DATA TK_LiNa/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751200018282, LiNa 
     2  0.831399951226, 0.926199936130, 1.033700122621, 1.108900115487, LiNa 
     3  1.185599952512, 1.333800022049, 1.488599814893, 1.645299987133, LiNa 
     4  1.754700105607, 1.863800040550, 2.067099899396, 2.322199758270, LiNa 
     5  2.617899825835, 2.901599642283, 3.060400410334, 3.200599811655, LiNa 
     6  3.292400025514, 3.392300189771, 3.488099714549, 3.586000066019, LiNa 
     7  3.676000202044, 3.715200147554, 3.754500093831, 3.833600005141, LiNa 
     8  3.861999710171, 3.892900142845, 3.928400287785, 3.957499992147, LiNa 
     9  3.983899639864, 3.991499810292, 4.000000000000,     55*0.0D+00/ LiNa 
      DATA  K_LiNa/                                                     071215
     1 -5.10605503D-06, 1.42111886D-01, 2.90194056D-01, 7.34046366D-01, LiNa 
     2  1.74006820D+00, 2.73866545D+00, 3.66970372D+00, 4.21827930D+00, LiNa 
     3  4.70650262D+00, 5.48734356D+00, 6.13132248D+00, 6.65751677D+00, LiNa 
     4  6.97028190D+00, 7.24740755D+00, 7.68713696D+00, 8.11684349D+00, LiNa 
     5  8.47696580D+00, 8.72389073D+00, 8.83087742D+00, 8.90552978D+00, LiNa 
     6  8.94152271D+00, 8.96890745D+00, 8.98854528D+00, 9.01244339D+00, LiNa 
     7  9.05179188D+00, 9.07973385D+00, 9.11916611D+00, 9.25926217D+00, LiNa 
     8  9.33991732D+00, 9.45157449D+00, 9.61281613D+00, 9.77047316D+00, LiNa 
     9  9.93101323D+00, 9.97995977D+00, 1.00358900D+01,     55*0.0D+00/ LiNa 
      DATA TK_AsP/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750600002732, AsP  
     2  0.829799923973, 0.923899991587, 1.030700202776, 1.104000002529, AsP  
     3  1.178299875364, 1.321400130735, 1.474599887931, 1.630200106183, AsP  
     4  1.767600033470, 1.916599992539, 2.048600166273, 2.180700256385, AsP  
     5  2.327299868703, 2.487199710659, 2.704999931886, 3.229699553478, AsP  
     6  3.318899776177, 3.411399866632, 3.604100083980, 3.687099802754, AsP  
     7  3.775700053460, 3.909799839731, 3.959200029267, 4.000000000000, AsP  
     8      62*0.0D+00/                                                 AsP  
      DATA  K_AsP/                                                      071215
     1  1.19410171D-04, 1.68427858D-01, 3.38588711D-01, 8.59149013D-01, AsP  
     2  2.03442907D+00, 3.20436169D+00, 4.29160406D+00, 4.91797106D+00, AsP  
     3  5.47109239D+00, 6.35022442D+00, 7.08611925D+00, 7.67980371D+00, AsP  
     4  8.11267656D+00, 8.51312337D+00, 8.82424379D+00, 9.10475484D+00, AsP  
     5  9.38592078D+00, 9.65788386D+00, 9.96924835D+00, 1.04866644D+01, AsP  
     6  1.05523140D+01, 1.06180281D+01, 1.07706857D+01, 1.08560581D+01, AsP  
     7  1.09661787D+01, 1.11632919D+01, 1.12400019D+01, 1.13033095D+01, AsP  
     8      62*0.0D+00/                                                 AsP  
      DATA TK_SbP/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, SbP  
     2  0.830199923717, 0.924699972297, 1.031800173386, 1.105200030192, SbP  
     3  1.179499847516, 1.322800098857, 1.477099946471, 1.634200006569, SbP  
     4  1.770399978561, 1.918600041436, 2.048300159406, 2.173100134217, SbP  
     5  2.465100141247, 2.603000155023, 2.734899664109, 3.062300268662, SbP  
     6  3.214100118650, 3.299600180459, 3.384500218429, 3.556300328018, SbP  
     7  3.653399701098, 3.749599978285, 3.874900003810, 3.925000212246, SbP  
     8  3.969000248433, 4.000000000000,     60*0.0D+00/                 SbP  
      DATA  K_SbP/                                                      071215
     1 -9.10510657D-05, 1.66954935D-01, 3.37546161D-01, 8.55715166D-01, SbP  
     2  2.02456593D+00, 3.18960811D+00, 4.27028552D+00, 4.89191584D+00, SbP  
     3  5.44020928D+00, 6.31315403D+00, 7.04792974D+00, 7.64170400D+00, SbP  
     4  8.06719829D+00, 8.46301882D+00, 8.76710509D+00, 9.03047700D+00, SbP  
     5  9.54879161D+00, 9.74746865D+00, 9.91204680D+00, 1.02335980D+01, SbP  
     6  1.03515933D+01, 1.04126666D+01, 1.04721398D+01, 1.06058910D+01, SbP  
     7  1.07026515D+01, 1.08201367D+01, 1.10000913D+01, 1.10764932D+01, SbP  
     8  1.11444386D+01, 1.11929694D+01,     60*0.0D+00/                 SbP  
      DATA TK_BeS/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720899818384, 0.754100093440, BeS  
     2  0.838700118571, 0.936099963346, 1.056900021805, 1.195199917336, BeS  
     3  1.341900116962, 1.492099826466, 1.651900047430, 1.842800080814, BeS  
     4  2.051700237472, 2.236999717081, 2.402099657682, 2.663499938249, BeS  
     5  2.780399734828, 2.934899948322, 3.089400127512, 3.273400216920, BeS  
     6  3.447399739474, 3.514300323626, 3.586100068362, 3.669100040858, BeS  
     7  3.749699980509, 3.884900235489, 3.956399968127, 3.983399628638, BeS  
     8  4.000000000000,     61*0.0D+00/                                 BeS  
      DATA  K_BeS/                                                      071215
     1  5.04320546D-05, 1.61131189D-01, 3.36239451D-01, 8.41122107D-01, BeS  
     2  1.98190807D+00, 3.07510412D+00, 4.17101990D+00, 5.15636030D+00, BeS  
     3  5.97104292D+00, 6.62898078D+00, 7.19100968D+00, 7.73516721D+00, BeS  
     4  8.22837811D+00, 8.61484588D+00, 8.93978053D+00, 9.42300268D+00, BeS  
     5  9.61812456D+00, 9.84840066D+00, 1.00455667D+01, 1.02381014D+01, BeS  
     6  1.03701382D+01, 1.04031635D+01, 1.04254948D+01, 1.04349739D+01, BeS  
     7  1.04320091D+01, 1.04260932D+01, 1.04383759D+01, 1.04483395D+01, BeS  
     8  1.04563483D+01,     61*0.0D+00/                                 BeS  
      DATA TK_BS/                                                       071215
     1  0.699999789529, 0.709600039198, 0.719799801886, 0.751600028649, BS   
     2  0.832199969565, 0.926599926485, 1.039799959638, 1.186899979579, BS   
     3  1.345000039870, 1.502300046250, 1.661599869635, 1.783200140581, BS   
     4  1.904899959013, 2.160499868625, 2.261200261212, 2.364599757113, BS   
     5  2.545000069478, 2.752900058583, 2.955299932640, 3.172400126147, BS   
     6  3.365699785233, 3.473100029539, 3.581199953561, 3.774400147535, BS   
     7  3.832199974647, 3.894300039680, 4.000000000000,     63*0.0D+00/ BS   
      DATA  K_BS/                                                       071215
     1  2.04835206D-04, 1.59769559D-01, 3.25957161D-01, 8.22702269D-01, BS   
     2  1.95062568D+00, 3.06908688D+00, 4.17998528D+00, 5.33308894D+00, BS   
     3  6.29653684D+00, 7.04672422D+00, 7.65332915D+00, 8.04044030D+00, BS   
     4  8.37895072D+00, 8.98298920D+00, 9.19678443D+00, 9.40880782D+00, BS   
     5  9.76532343D+00, 1.01465291D+01, 1.04679537D+01, 1.07489256D+01, BS   
     6  1.09496185D+01, 1.10452532D+01, 1.11313640D+01, 1.12493140D+01, BS   
     7  1.12711710D+01, 1.12863441D+01, 1.12970493D+01,     63*0.0D+00/ BS   
      DATA TK_MgS/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, MgS  
     2  0.830399928302, 0.924999965064, 1.032300160026, 1.107100073992, MgS  
     3  1.183199902541, 1.329999934916, 1.484699905356, 1.639699869601, MgS  
     4  1.756900161040, 1.882900053603, 2.005400123997, 2.138000274886, MgS  
     5  2.252000048740, 2.362599708658, 2.611499678510, 2.756500147219, MgS  
     6  2.904199700750, 3.291900014754, 3.466000165222, 3.605299997135, MgS  
     7  3.715200147554, 3.770800408050, 3.824499793761, 3.880600133424, MgS  
     8  3.934200020622, 3.973799997799, 3.988799749873, 4.000000000000, MgS  
     9      58*0.0D+00/                                                 MgS  
      DATA  K_MgS/                                                      071215
     1 -4.92854724D-05, 1.53937422D-01, 3.11217563D-01, 7.89124588D-01, MgS  
     2  1.87075612D+00, 2.94894042D+00, 3.95208252D+00, 4.53971864D+00, MgS  
     3  5.06044439D+00, 5.88953603D+00, 6.57528891D+00, 7.12602615D+00, MgS  
     4  7.47731643D+00, 7.80906423D+00, 8.09751084D+00, 8.38222691D+00, MgS  
     5  8.60995860D+00, 8.81827918D+00, 9.24014500D+00, 9.45153618D+00, MgS  
     6  9.63907572D+00, 1.00128532D+01, 1.01374438D+01, 1.02244658D+01, MgS  
     7  1.02897571D+01, 1.03245827D+01, 1.03623267D+01, 1.04106466D+01, MgS  
     8  1.04728932D+01, 1.05355720D+01, 1.05642491D+01, 1.05877153D+01, MgS  
     9      58*0.0D+00/                                                 MgS  
      DATA TK_AlS/                                                      071215
     1  0.699999789529, 0.709900047001, 0.720999820779, 0.754400101215, AlS  
     2  0.839300132325, 0.938000000362, 1.058299990483, 1.195299914898, AlS  
     3  1.341500126910, 1.498099952317, 1.659499845207, 1.782200165272, AlS  
     4  1.900400071808, 2.192899647615, 2.371999935668, 2.609699665882, AlS  
     5  2.728699652175, 2.846900317228, 3.072499757126, 3.294500070706, AlS  
     6  3.448399760290, 3.594400256959, 3.683500053486, 3.767800412754, AlS  
     7  3.838700116225, 3.907999799308, 3.963700129656, 3.985899684766, AlS  
     8  4.000000000000,     61*0.0D+00/                                 AlS  
      DATA  K_AlS/                                                      071215
     1  3.80702804D-05, 1.61736903D-01, 3.39095032D-01, 8.48747063D-01, AlS  
     2  1.99684737D+00, 3.10600065D+00, 4.19777161D+00, 5.17577687D+00, AlS  
     3  5.99093933D+00, 6.68072837D+00, 7.26429386D+00, 7.65642815D+00, AlS  
     4  8.00708676D+00, 8.78415320D+00, 9.20191731D+00, 9.68550479D+00, AlS  
     5  9.89332711D+00, 1.00761872D+01, 1.03644911D+01, 1.05833895D+01, AlS  
     6  1.07066095D+01, 1.08068098D+01, 1.08593727D+01, 1.09016286D+01, AlS  
     7  1.09324191D+01, 1.09636642D+01, 1.09975519D+01, 1.10154257D+01, AlS  
     8  1.10285493D+01,     61*0.0D+00/                                 AlS  
      DATA TK_SiS/                                                      071215
     1  0.699999789529, 0.710100047074, 0.721799839938, 0.756200147865, SiS  
     2  0.792500008960, 0.844200047034, 0.949799798715, 1.067400123660, SiS  
     3  1.157899895227, 1.249099996588, 1.330099937209, 1.409599860182, SiS  
     4  1.486599861284, 1.563000008636, 1.640999885715, 1.720499800536, SiS  
     5  1.854600012613, 1.983099933777, 2.121699916929, 2.256300148854, SiS  
     6  2.396399871070, 2.549000159175, 2.680000308330, 2.827399879649, SiS  
     7  2.977699709020, 3.129000095560, 3.440699600004, 3.551300223671, SiS  
     8  3.662699910215, 3.755200110421, 3.852700179147, 3.895399958621, SiS  
     9  3.936599847644, 3.974599940316, 4.000000000000,     55*0.0D+00/ SiS  
      DATA  K_SiS/                                                      071215
     1  2.36527912D-05, 1.66095385D-01, 3.54120550D-01, 8.81204000D-01, SiS  
     2  1.39817457D+00, 2.07117673D+00, 3.24625331D+00, 4.29890275D+00, SiS  
     3  4.96507935D+00, 5.53704939D+00, 5.97872910D+00, 6.36513419D+00, SiS  
     4  6.70547451D+00, 7.01831885D+00, 7.31816587D+00, 7.60754792D+00, SiS  
     5  8.06561685D+00, 8.47536717D+00, 8.88996637D+00, 9.26733750D+00, SiS  
     6  9.63336177D+00, 9.99703693D+00, 1.02754712D+01, 1.05483573D+01, SiS  
     7  1.07835973D+01, 1.09815995D+01, 1.13041354D+01, 1.14049615D+01, SiS  
     8  1.15044938D+01, 1.15842720D+01, 1.16574544D+01, 1.16818689D+01, SiS  
     9  1.16989897D+01, 1.17088417D+01, 1.17127549D+01,     55*0.0D+00/ SiS  
      DATA TK_PS/                                                       071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, PS   
     2  0.830399928302, 0.924899967475, 1.032000168042, 1.104500014055, PS   
     3  1.177799886967, 1.320100160335, 1.478199972228, 1.638999887033, PS   
     4  1.757500176158, 1.873399968135, 2.002600059702, 2.151899667160, PS   
     5  2.287699900059, 2.412899920745, 2.667000015599, 2.873699971491, PS   
     6  3.139300333617, 3.266900386114, 3.388300314656, 3.617699843809, PS   
     7  3.732699613933, 3.822899755079, 3.899799634389, 3.953699909171, PS   
     8  4.000000000000,     61*0.0D+00/                                 PS   
      DATA  K_PS/                                                       071215
     1  2.35093310D-04, 1.63534334D-01, 3.31972972D-01, 8.38490627D-01, PS   
     2  1.98271665D+00, 3.12225536D+00, 4.18002144D+00, 4.78195086D+00, PS   
     3  5.31348636D+00, 6.16758068D+00, 6.90892292D+00, 7.50539655D+00, PS   
     4  7.87235933D+00, 8.18690818D+00, 8.49690198D+00, 8.81235415D+00, PS   
     5  9.06761157D+00, 9.28186911D+00, 9.66406868D+00, 9.92327427D+00, PS   
     6  1.01919520D+01, 1.02997518D+01, 1.03936774D+01, 1.05733262D+01, PS   
     7  1.06809385D+01, 1.07764104D+01, 1.08623357D+01, 1.09228810D+01, PS   
     8  1.09744595D+01,     61*0.0D+00/                                 PS   
      DATA TK_CaS/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749599978008, CaS  
     2  0.827099989317, 0.919200067714, 1.024100075777, 1.095400027976, CaS  
     3  1.167700015896, 1.306600030349, 1.454299902336, 1.604199995812, CaS  
     4  1.734499912679, 1.874499995348, 2.006100140071, 2.148999692089, CaS  
     5  2.286899883844, 2.413399932870, 2.582699988704, 2.730199562478, CaS  
     6  2.962200091505, 3.295700096530, 3.389400342512, 3.485699661323, CaS  
     7  3.592100206137, 3.698399774600, 3.872299944992, 3.910999869916, CaS  
     8  3.947099762901, 3.979999552305, 3.990399785742, 4.000000000000, CaS  
     9      58*0.0D+00/                                                 CaS  
      DATA  K_CaS/                                                      071215
     1  6.48807422D-05, 1.52109991D-01, 3.02808229D-01, 7.69860146D-01, CaS  
     2  1.82548640D+00, 2.88246491D+00, 3.87443829D+00, 4.44391333D+00, CaS  
     3  4.94991782D+00, 5.75903407D+00, 6.44011762D+00, 6.99608886D+00, CaS  
     4  7.39974076D+00, 7.77375129D+00, 8.08405269D+00, 8.38801941D+00, CaS  
     5  8.65637621D+00, 8.88381639D+00, 9.16000191D+00, 9.37246490D+00, CaS  
     6  9.65266628D+00, 9.95569999D+00, 1.00249075D+01, 1.00928407D+01, CaS  
     7  1.01732264D+01, 1.02803064D+01, 1.06092669D+01, 1.07249469D+01, CaS  
     8  1.08488617D+01, 1.09736730D+01, 1.10150940D+01, 1.10540273D+01, CaS  
     9      58*0.0D+00/                                                 CaS  
      DATA TK_ScS/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718599832218, 0.748799959658, ScS  
     2  0.825000040140, 0.915899993836, 1.020199979466, 1.157899895227, ScS  
     3  1.288099999024, 1.432899881399, 1.585200047521, 1.708900023720, ScS  
     4  1.828199963314, 1.942599977907, 2.055000313469, 2.166399985478, ScS  
     5  2.279599772368, 2.427700239561, 2.598300337609, 2.722600093353, ScS  
     6  2.848300347601, 3.152799667231, 3.276899966124, 3.398999700569, ScS  
     7  3.512100272249, 3.634600216868, 3.728499662094, 3.828999902553, ScS  
     8  3.939799617007, 3.976399810979, 4.000000000000,     59*0.0D+00/ ScS  
      DATA  K_ScS/                                                      071215
     1  2.08427957D-04, 1.57743726D-01, 3.12284770D-01, 7.93804845D-01, ScS  
     2  1.88294609D+00, 2.97919876D+00, 4.01588329D+00, 5.10352691D+00, ScS  
     3  5.91055578D+00, 6.62397202D+00, 7.22438469D+00, 7.63325604D+00, ScS  
     4  7.98356847D+00, 8.29405636D+00, 8.58389499D+00, 8.86000808D+00, ScS  
     5  9.12931807D+00, 9.46084715D+00, 9.80530263D+00, 1.00266467D+01, ScS  
     6  1.02240520D+01, 1.06029975D+01, 1.07250346D+01, 1.08332531D+01, ScS  
     7  1.09313495D+01, 1.10541707D+01, 1.11772210D+01, 1.13470929D+01, ScS  
     8  1.15766752D+01, 1.16602515D+01, 1.17158338D+01,     59*0.0D+00/ ScS  
      DATA TK_TiS/                                                      071215
     1  0.699999789529, 0.709200028796, 0.718299839802, 0.748199945895, TiS  
     2  0.823600074022, 0.913099931151, 1.015500082035, 1.085300050501, TiS  
     3  1.156399933583, 1.292299998225, 1.436199956315, 1.576899996321, TiS  
     4  1.875200012665, 1.996099911371, 2.125800002922, 2.230899569146, TiS  
     5  2.332399984058, 2.572499756971, 2.687799743773, 2.817999843839, TiS  
     6  2.959500029403, 3.111599691531, 3.222100104917, 3.333199994591, TiS  
     7  3.481299563743, 3.631000142653, 3.741599800333, 3.845500268272, TiS  
     8  3.942799665813, 3.977799710384, 4.000000000000,     59*0.0D+00/ TiS  
      DATA  K_TiS/                                                      071215
     1  1.91574573D-04, 1.59579574D-01, 3.14292496D-01, 8.02234941D-01, TiS  
     2  1.90623876D+00, 3.01384404D+00, 4.06085716D+00, 4.66441061D+00, TiS  
     3  5.20310406D+00, 6.05912096D+00, 6.77289041D+00, 7.32736912D+00, TiS  
     4  8.18827480D+00, 8.46606018D+00, 8.74608940D+00, 8.96779616D+00, TiS  
     5  9.17959316D+00, 9.65781781D+00, 9.86528825D+00, 1.00760038D+01, TiS  
     6  1.02754279D+01, 1.04581705D+01, 1.05747320D+01, 1.06845321D+01, TiS  
     7  1.08339988D+01, 1.10097098D+01, 1.11669931D+01, 1.13390463D+01, TiS  
     8  1.15194689D+01, 1.15883579D+01, 1.16330811D+01,     59*0.0D+00/ TiS  
      DATA TK_CrS/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718699829691, 0.748999964245, CrS  
     2  0.825500028039, 0.916700011745, 1.020599989344, 1.091300132846, CrS  
     3  1.162999909464, 1.300300164416, 1.444899928068, 1.592500099017, CrS  
     4  1.730300013678, 1.880400120680, 2.035799867476, 2.190999609554, CrS  
     5  2.302600256266, 2.406899774788, 2.682900098431, 2.812000272304, CrS  
     6  2.941299619034, 3.070899717056, 3.190499575567, 3.311700288239, CrS  
     7  3.441699620821, 3.607499837920, 3.709000012620, 3.845200261523, CrS  
     8  3.936599847644, 3.975399882833, 3.989199758854, 4.000000000000, CrS  
     9      58*0.0D+00/                                                 CrS  
      DATA  K_CrS/                                                      071215
     1  9.26517373D-05, 1.71171964D-01, 3.40759649D-01, 8.64994425D-01, CrS  
     2  2.04988759D+00, 3.23903895D+00, 4.35311676D+00, 4.99178708D+00, CrS  
     3  5.55782553D+00, 6.45607868D+00, 7.20025574D+00, 7.80619681D+00, CrS  
     4  8.27121020D+00, 8.69863757D+00, 9.08015585D+00, 9.42005960D+00, CrS  
     5  9.64715905D+00, 9.84853476D+00, 1.03247699D+01, 1.05141014D+01, CrS  
     6  1.06812597D+01, 1.08278105D+01, 1.09473861D+01, 1.10594152D+01, CrS  
     7  1.11827264D+01, 1.13718909D+01, 1.15171594D+01, 1.17695538D+01, CrS  
     8  1.19907507D+01, 1.20978589D+01, 1.21376588D+01, 1.21693919D+01, CrS  
     9      58*0.0D+00/                                                 CrS  
      DATA TK_CuS/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719799801886, 0.751700031241, CuS  
     2  0.832499976442, 0.928399883083, 1.036900037122, 1.110800120079, CuS  
     3  1.185599952512, 1.331199962432, 1.492399832759, 1.654699972927, CuS  
     4  1.757700181198, 1.867199954712, 1.961500016401, 2.056000336498, CuS  
     5  2.248399968201, 2.395799915173, 2.559000389443, 2.766200383155, CuS  
     6  2.982699601888, 3.234899650718, 3.322199745646, 3.405599740808, CuS  
     7  3.596900312200, 3.687999740071, 3.792100022081, 3.868799865543, CuS  
     8  3.949699821606, 3.980699568021, 3.990499787974, 4.000000000000, CuS  
     9      58*0.0D+00/                                                 CuS  
      DATA  K_CuS/                                                      071215
     1  1.49789107D-05, 1.58487346D-01, 3.23431768D-01, 8.17210527D-01, CuS  
     2  1.93167520D+00, 3.03798376D+00, 4.06181673D+00, 4.64740812D+00, CuS  
     3  5.16453547D+00, 5.99660643D+00, 6.71615478D+00, 7.29107075D+00, CuS  
     4  7.60001667D+00, 7.89276194D+00, 8.12122120D+00, 8.33173625D+00, CuS  
     5  8.70913419D+00, 8.95275351D+00, 9.17867363D+00, 9.41015939D+00, CuS  
     6  9.60339144D+00, 9.78606738D+00, 9.84217406D+00, 9.89462925D+00, CuS  
     7  1.00284860D+01, 1.01095527D+01, 1.02202469D+01, 1.03137476D+01, CuS  
     8  1.04261729D+01, 1.04750956D+01, 1.04915080D+01, 1.05079101D+01, CuS  
     9      58*0.0D+00/                                                 CuS  
      DATA TK_GeS/                                                      071215
     1  0.699999789529, 0.709100026195, 0.718199842329, 0.747899939014, GeS  
     2  0.822800093383, 0.910999884138, 1.012800146540, 1.083099996392, GeS  
     3  1.154799974495, 1.291899989260, 1.429899818212, 1.573700078479, GeS  
     4  1.688500036521, 1.801800135908, 1.989099793750, 2.162499908236, GeS  
     5  2.242999853538, 2.328999905514, 2.415099974097, 2.497099923902, GeS  
     6  2.641400238905, 2.891100267939, 3.046800116888, 3.206499943873, GeS  
     7  3.507400162593, 3.670900080652, 3.800900207436, 3.886800280588, GeS  
     8  3.958200007432, 3.983999642109, 4.000000000000,     59*0.0D+00/ GeS  
      DATA  K_GeS/                                                      071215
     1  2.00106619D-04, 1.50919807D-01, 2.98861281D-01, 7.62487654D-01, GeS  
     2  1.81280571D+00, 2.86096740D+00, 3.86354927D+00, 4.44983070D+00, GeS  
     3  4.97384314D+00, 5.80755305D+00, 6.47513934D+00, 7.03792818D+00, GeS  
     4  7.41552480D+00, 7.74138179D+00, 8.20649550D+00, 8.58649763D+00, GeS  
     5  8.75651986D+00, 8.93924027D+00, 9.12618064D+00, 9.30871793D+00, GeS  
     6  9.63670926D+00, 1.01879704D+01, 1.04928318D+01, 1.07613654D+01, GeS  
     7  1.11627592D+01, 1.13468292D+01, 1.14819540D+01, 1.15641624D+01, GeS  
     8  1.16292246D+01, 1.16533689D+01, 1.16690302D+01,     59*0.0D+00/ GeS  
      DATA TK_AsS/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, AsS  
     2  0.829499931233, 0.923500001232, 1.030100218807, 1.104100004834, AsS  
     3  1.179399849837, 1.324600057872, 1.477599958179, 1.631200081280, AsS  
     4  1.752500050174, 1.883900026773, 2.015600361429, 2.156199768853, AsS  
     5  2.272300290261, 2.386000255535, 2.611999690019, 2.849200367127, AsS  
     6  3.127200050692, 3.248299963484, 3.369199878036, 3.572899755825, AsS  
     7  3.650199626431, 3.729599585871, 3.891700231272, 3.953899913538, AsS  
     8  4.000000000000,     61*0.0D+00/                                 AsS  
      DATA  K_AsS/                                                      071215
     1 -3.95606826D-05, 1.57293782D-01, 3.16379776D-01, 8.01696412D-01, AsS  
     2  1.89938983D+00, 2.99565134D+00, 4.01594533D+00, 4.61145961D+00, AsS  
     3  5.13959278D+00, 5.98061082D+00, 6.67621453D+00, 7.23538378D+00, AsS  
     4  7.60616035D+00, 7.95637279D+00, 8.26794837D+00, 8.56992006D+00, AsS  
     5  8.80217940D+00, 9.01720118D+00, 9.40549136D+00, 9.74597207D+00, AsS  
     6  1.00559690D+01, 1.01658228D+01, 1.02646859D+01, 1.04258725D+01, AsS  
     7  1.04930243D+01, 1.05679608D+01, 1.07319662D+01, 1.07932977D+01, AsS  
     8  1.08369945D+01,     61*0.0D+00/                                 AsS  
      DATA TK_SeS/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.751000013099, SeS  
     2  0.830799937472, 0.925299957830, 1.032400157355, 1.107000071687, SeS  
     3  1.183099900459, 1.330599948674, 1.486299868243, 1.644199961189, SeS  
     4  1.749999987182, 1.850399916268, 2.006800156144, 2.216400179352, SeS  
     5  2.321499743113, 2.421700125299, 2.678400272101, 3.035799882280, SeS  
     6  3.236699694435, 3.365099769324, 3.505900127623, 3.620699915181, SeS  
     7  3.724999904622, 3.826299837278, 3.926300241128, 3.970600227731, SeS  
     8  4.000000000000,     61*0.0D+00/                                 SeS  
      DATA  K_SeS/                                                      071215
     1  5.99124719D-05, 1.62130019D-01, 3.29302808D-01, 8.33576287D-01, SeS  
     2  1.97303201D+00, 3.10338255D+00, 4.15286940D+00, 4.76636394D+00, SeS  
     3  5.31071243D+00, 6.17881491D+00, 6.89421292D+00, 7.47119625D+00, SeS  
     4  7.79615233D+00, 8.06813589D+00, 8.43332491D+00, 8.83779214D+00, SeS  
     5  9.01551491D+00, 9.17385655D+00, 9.53646507D+00, 9.95264050D+00, SeS  
     6  1.01599388D+01, 1.02848569D+01, 1.04161545D+01, 1.05208741D+01, SeS  
     7  1.06146278D+01, 1.07016297D+01, 1.07779866D+01, 1.08076004D+01, SeS  
     8  1.08259741D+01,     61*0.0D+00/                                 SeS  
      DATA TK_SrS/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749499975714, SrS  
     2  0.826999991737, 0.919100065475, 1.023900070838, 1.095000038208, SrS  
     3  1.167000000044, 1.305200060142, 1.452599865388, 1.603499978396, SrS  
     4  1.735899879013, 1.873899980504, 2.004900112516, 2.150199626956, SrS  
     5  2.298900168316, 2.432800087942, 2.768600441128, 3.038099936610, SrS  
     6  3.290499984625, 3.381900152589, 3.474099957854, 3.569699708104, SrS  
     7  3.670500071131, 3.762100274950, 3.855499981738, 3.901299648845, SrS  
     8  3.941599638719, 3.977699717569, 3.989799772324, 4.000000000000, SrS  
     9      58*0.0D+00/                                                 SrS  
      DATA  K_SrS/                                                      071215
     1  2.06019431D-05, 1.52379390D-01, 3.03387909D-01, 7.69920541D-01, SrS  
     2  1.82790012D+00, 2.88718063D+00, 3.88039008D+00, 4.44963373D+00, SrS  
     3  4.95493829D+00, 5.76276192D+00, 6.44513996D+00, 7.00652528D+00, SrS  
     4  7.41700579D+00, 7.78568950D+00, 8.09408947D+00, 8.39998342D+00, SrS  
     5  8.68140016D+00, 8.91123012D+00, 9.39461184D+00, 9.68889593D+00, SrS  
     6  9.90177271D+00, 9.96788638D+00, 1.00328061D+01, 1.01073893D+01, SrS  
     7  1.02149629D+01, 1.03634019D+01, 1.05854415D+01, 1.07235170D+01, SrS  
     8  1.08598851D+01, 1.09918559D+01, 1.10377383D+01, 1.10769414D+01, SrS  
     9      58*0.0D+00/                                                 SrS  
      DATA TK_YS/                                                       071215
     1  0.699999789529, 0.709200028796, 0.718499834746, 0.748499952777, YS   
     2  0.824300057081, 0.914299958016, 1.017200041421, 1.087900114448, YS   
     3  1.159699849201, 1.297300110287, 1.439400028961, 1.584900039976, YS   
     4  1.716799872351, 1.845500015668, 1.989099793750, 2.147499796883, YS   
     5  2.255500130229, 2.370499899873, 2.476099809130, 2.585900063676, YS   
     6  2.695399705051, 2.820799720329, 2.973400022766, 3.119099851967, YS   
     7  3.266000364504, 3.404499718561, 3.499599980738, 3.644200031491, YS   
     8  3.729899565083, 3.823399767167, 3.893700083893, 3.946699753870, YS   
     9  4.000000000000,     57*0.0D+00/                                 YS   
      DATA  K_YS/                                                       071215
     1  2.10416202D-04, 1.55747920D-01, 3.10012404D-01, 7.87574757D-01, YS   
     2  1.86994005D+00, 2.95570098D+00, 3.98108207D+00, 4.57657455D+00, YS   
     3  5.10613547D+00, 5.94975173D+00, 6.63991601D+00, 7.20899415D+00, YS   
     4  7.63774291D+00, 7.99772342D+00, 8.34994178D+00, 8.69754985D+00, YS   
     5  8.92001562D+00, 9.15079164D+00, 9.35959264D+00, 9.57220796D+00, YS   
     6  9.77559856D+00, 9.99215608D+00, 1.02263770D+01, 1.04189342D+01, YS   
     7  1.05857727D+01, 1.07255120D+01, 1.08199964D+01, 1.09871342D+01, YS   
     8  1.11140625D+01, 1.12805947D+01, 1.14212817D+01, 1.15330918D+01, YS   
     9  1.16486317D+01,     57*0.0D+00/                                 YS   
      DATA TK_SnS/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718599832218, 0.748899961951, SnS  
     2  0.825300032879, 0.915999996074, 1.019399988862, 1.089300148881, SnS  
     3  1.159999841529, 1.296200085633, 1.442199991168, 1.589000143099, SnS  
     4  1.703299891963, 1.824600063928, 1.945299914126, 2.071499720840, SnS  
     5  2.267300394456, 2.450199807466, 2.573999790666, 2.695399705051, SnS  
     6  2.891300253071, 2.987899718765, 3.076399854795, 3.206499943873, SnS  
     7  3.315899989536, 3.428800269176, 3.597200318829, 3.706899964808, SnS  
     8  3.839300129294, 3.920300107825, 4.000000000000,     59*0.0D+00/ SnS  
      DATA  K_SnS/                                                      071215
     1  1.65179785D-04, 1.53070727D-01, 3.03078547D-01, 7.72032891D-01, SnS  
     2  1.83208802D+00, 2.89461094D+00, 3.89485709D+00, 4.46709998D+00, SnS  
     3  4.97524687D+00, 5.79181079D+00, 6.48535082D+00, 7.04616292D+00, SnS  
     4  7.41358981D+00, 7.75419257D+00, 8.05489263D+00, 8.33889835D+00, SnS  
     5  8.73534162D+00, 9.06728574D+00, 9.27246228D+00, 9.46160841D+00, SnS  
     6  9.76313552D+00, 9.92056771D+00, 1.00719730D+01, 1.03017877D+01, SnS  
     7  1.04934005D+01, 1.06826678D+01, 1.09432530D+01, 1.10964284D+01, SnS  
     8  1.12533538D+01, 1.13272037D+01, 1.13860187D+01,     59*0.0D+00/ SnS  
      DATA TK_TeS/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750700005324, TeS  
     2  0.830099921425, 0.924499977120, 1.031500181401, 1.105900046329, TeS  
     3  1.181499867144, 1.327399994117, 1.481399981903, 1.636299954272, TeS  
     4  1.760300231031, 1.893999950597, 2.020200450224, 2.153599707364, TeS  
     5  2.275400070334, 2.393600076886, 2.630200108135, 2.838800140515, TeS  
     6  3.045800096976, 3.389600347576, 3.499499978404, 3.618399860792, TeS  
     7  3.711800074215, 3.829899924311, 3.923300174476, 3.968800243951, TeS  
     8  4.000000000000,     61*0.0D+00/                                 TeS  
      DATA  K_TeS/                                                      071215
     1  8.35750184D-05, 1.71257379D-01, 3.46056006D-01, 8.75315114D-01, TeS  
     2  2.07277349D+00, 3.26455645D+00, 4.36959546D+00, 5.01347441D+00, TeS  
     3  5.58197169D+00, 6.48375692D+00, 7.22540481D+00, 7.81768492D+00, TeS  
     4  8.21228463D+00, 8.58017397D+00, 8.88680488D+00, 9.17929460D+00, TeS  
     5  9.42434534D+00, 9.64487539D+00, 1.00359397D+01, 1.03219216D+01, TeS  
     6  1.05538983D+01, 1.08675272D+01, 1.09639916D+01, 1.10714117D+01, TeS  
     7  1.11555037D+01, 1.12492154D+01, 1.13020049D+01, 1.13203210D+01, TeS  
     8  1.13311112D+01,     61*0.0D+00/                                 TeS  
      DATA TK_BaS/                                                      071215
     1  0.699999789529, 0.709300031396, 0.718799827163, 0.749299971126, BaS  
     2  0.826300008678, 0.918000040849, 1.022400033796, 1.092800094479, BaS  
     3  1.163999932109, 1.300000170800, 1.444599935079, 1.597199968856, BaS  
     4  1.750399997261, 1.908099878803, 2.060600383893, 2.223200024668, BaS  
     5  2.464300123203, 2.625199993371, 2.791100002964, 2.973100044655, BaS  
     6  3.146799842244, 3.220300235521, 3.293300044882, 3.393600094851, BaS  
     7  3.483099603662, 3.569099751225, 3.634300210684, 3.697999766323, BaS  
     8  3.766800388578, 3.843900232277, 3.901099644354, 3.959600038001, BaS  
     9  3.984699657824, 4.000000000000,     56*0.0D+00/                 BaS  
      DATA  K_BaS/                                                      071215
     1 -1.32738725D-04, 1.50598088D-01, 3.01622767D-01, 7.66714280D-01, BaS  
     2  1.81868931D+00, 2.87512593D+00, 3.86690734D+00, 4.43242760D+00, BaS  
     3  4.93441032D+00, 5.73481807D+00, 6.41101701D+00, 6.98434426D+00, BaS  
     4  7.45728785D+00, 7.86858055D+00, 8.21399783D+00, 8.53960195D+00, BaS  
     5  8.95820193D+00, 9.19928543D+00, 9.41611355D+00, 9.61840413D+00, BaS  
     6  9.78128786D+00, 9.84398626D+00, 9.90617578D+00, 1.00022962D+01, BaS  
     7  1.01138768D+01, 1.02523187D+01, 1.03731247D+01, 1.04968471D+01, BaS  
     8  1.06323208D+01, 1.07901193D+01, 1.09175092D+01, 1.10604888D+01, BaS  
     9  1.11257136D+01, 1.11664926D+01,     56*0.0D+00/                 BaS  
      DATA TK_LaS/                                                      071215
     1  0.699999789529, 0.709100026195, 0.718099844857, 0.747699934427, LaS  
     2  0.822400103064, 0.910699877422, 1.011900168042, 1.147800048735, LaS  
     3  1.278700190949, 1.418200038944, 1.561699978023, 1.696999891512, LaS  
     4  1.840200143546, 1.970299801424, 2.099200370527, 2.232299603099, LaS  
     5  2.383900210399, 2.590000159734, 2.721300187375, 2.878500073368, LaS  
     6  3.046300106932, 3.165599967566, 3.352800177728, 3.468100210656, LaS  
     7  3.586800084762, 3.687499774895, 3.781099768224, 3.864099758153, LaS  
     8  3.928700294450, 3.972400098394, 4.000000000000,     59*0.0D+00/ LaS  
      DATA  K_LaS/                                                      071215
     1 -1.35798865D-04, 1.54763302D-01, 3.05145982D-01, 7.80093382D-01, LaS  
     2  1.85679064D+00, 2.93486415D+00, 3.95847967D+00, 5.05501557D+00, LaS  
     3  5.88361729D+00, 6.58665709D+00, 7.16870004D+00, 7.62239742D+00, LaS  
     4  8.02990099D+00, 8.35288954D+00, 8.64041535D+00, 8.91160017D+00, LaS  
     5  9.19524579D+00, 9.54841902D+00, 9.76147044D+00, 1.00109907D+01, LaS  
     6  1.02782643D+01, 1.04743235D+01, 1.07939839D+01, 1.09970577D+01, LaS  
     7  1.12132506D+01, 1.14059167D+01, 1.15952535D+01, 1.17714476D+01, LaS  
     8  1.19121306D+01, 1.20076642D+01, 1.20676852D+01,     59*0.0D+00/ LaS  
      DATA TK_PbS/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, PbS  
     2  0.829599928813, 0.923500001232, 1.030100218807, 1.103599993308, PbS  
     3  1.178099880005, 1.321600126181, 1.474399883248, 1.629400097898, PbS  
     4  1.761400201262, 1.902600016664, 2.035799867476, 2.175300189173, PbS  
     5  2.302000242291, 2.422500140534, 2.674400181528, 2.886000251910, PbS  
     6  3.114399751427, 3.199999798209, 3.282199794563, 3.456699955612, PbS  
     7  3.538799930048, 3.619499887480, 3.806000319232, 3.883400199885, PbS  
     8  3.928900298893, 3.962800109486, 3.985999687011, 4.000000000000, PbS  
     9      58*0.0D+00/                                                 PbS  
      DATA  K_PbS/                                                      071215
     1  2.37073613D-04, 1.56216985D-01, 3.13937032D-01, 7.95103959D-01, PbS  
     2  1.88479413D+00, 2.97073192D+00, 3.98286349D+00, 4.57005219D+00, PbS  
     3  5.08950665D+00, 5.91778679D+00, 6.61135102D+00, 7.17438882D+00, PbS  
     4  7.57437078D+00, 7.94399878D+00, 8.25214925D+00, 8.54334123D+00, PbS  
     5  8.78538292D+00, 8.99799562D+00, 9.38751251D+00, 9.65546872D+00, PbS  
     6  9.88883394D+00, 9.96440133D+00, 1.00336884D+01, 1.01893307D+01, PbS  
     7  1.02779366D+01, 1.03778486D+01, 1.06208587D+01, 1.07035523D+01, PbS  
     8  1.07450331D+01, 1.07739742D+01, 1.07937096D+01, 1.08059080D+01, PbS  
     9      58*0.0D+00/                                                 PbS  
      DATA TK_BiS/                                                      071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751500026057, BiS  
     2  0.832199969565, 0.928099890317, 1.036600045137, 1.112500075952, BiS  
     3  1.189800039961, 1.339000141283, 1.496799925049, 1.656999911728, BiS  
     4  1.780100217124, 1.909799836191, 2.036099874757, 2.168000017167, BiS  
     5  2.395399944576, 2.624499977326, 2.810400386561, 3.003900087028, BiS  
     6  3.272000317239, 3.483099603662, 3.604300069506, 3.733299626328, BiS  
     7  3.838400109691, 3.931700200807, 3.973600012170, 4.000000000000, BiS  
     8      62*0.0D+00/                                                 BiS  
      DATA  K_BiS/                                                      071215
     1 -9.97128436D-05, 1.57233467D-01, 3.19406759D-01, 8.08272215D-01, BiS  
     2  1.91406901D+00, 3.01341322D+00, 4.03100906D+00, 4.62783205D+00, BiS  
     3  5.15614543D+00, 5.99537490D+00, 6.68943840D+00, 7.25175309D+00, BiS  
     4  7.61454828D+00, 7.94935219D+00, 8.24002057D+00, 8.51476977D+00, BiS  
     5  8.93371985D+00, 9.29349738D+00, 9.53915316D+00, 9.75234407D+00, BiS  
     6  9.98704542D+00, 1.01356746D+01, 1.02126978D+01, 1.02897051D+01, BiS  
     7  1.03474038D+01, 1.03962596D+01, 1.04194144D+01, 1.04352269D+01, BiS  
     8      62*0.0D+00/                                                 BiS  
      DATA TK_LiCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, LiCl 
     2  0.830499930595, 0.925099962653, 1.032300160026, 1.106600062466, LiCl 
     3  1.182099879637, 1.327899982732, 1.482799949428, 1.639199882052, LiCl 
     4  1.763700139017, 1.895199976845, 2.012200281404, 2.134300189844, LiCl 
     5  2.286599877763, 2.475899823643, 2.647099825936, 2.809900413020, LiCl 
     6  2.973000051951, 3.329899912600, 3.471500144235, 3.612099707943, LiCl 
     7  3.717900205793, 3.769100444183, 3.816799915992, 3.926600247793, LiCl 
     8  3.971100191804, 3.988099734158, 4.000000000000,     59*0.0D+00/ LiCl 
      DATA  K_LiCl/                                                     071215
     1 -1.95924446D-04, 1.55067997D-01, 3.15233930D-01, 7.96977567D-01, LiCl 
     2  1.88716956D+00, 2.97367428D+00, 3.98344429D+00, 4.57171675D+00, LiCl 
     3  5.09274733D+00, 5.92389697D+00, 6.61677436D+00, 7.17630899D+00, LiCl 
     4  7.54984915D+00, 7.89404377D+00, 8.16757963D+00, 8.42796238D+00, LiCl 
     5  8.72396549D+00, 9.05249706D+00, 9.31429389D+00, 9.53458699D+00, LiCl 
     6  9.72877675D+00, 1.00638107D+01, 1.01643586D+01, 1.02496438D+01, LiCl 
     7  1.03101122D+01, 1.03419870D+01, 1.03768361D+01, 1.05018374D+01, LiCl 
     8  1.05817640D+01, 1.06177483D+01, 1.06447642D+01,     59*0.0D+00/ LiCl 
      DATA TK_BeCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, BeCl 
     2  0.830499930595, 0.925199960242, 1.032500154683, 1.105000025581, BeCl 
     3  1.178099880005, 1.320900142120, 1.481399981903, 1.642699925811, BeCl 
     4  1.757600178678, 1.875900029982, 1.991399804561, 2.118399848731, BeCl 
     5  2.233899641901, 2.352100233011, 2.595100269031, 2.890400319978, BeCl 
     6  3.040699995429, 3.203299872162, 3.391100277389, 3.537099891361, BeCl 
     7  3.658999831766, 3.789699970981, 3.920400110046, 3.969300255156, BeCl 
     8  3.987499720687, 4.000000000000,     60*0.0D+00/                 BeCl 
      DATA  K_BeCl/                                                     071215
     1 -2.05551369D-05, 1.45520680D-01, 2.95677584D-01, 7.47452963D-01, BeCl 
     2  1.77070897D+00, 2.79306455D+00, 3.74491521D+00, 4.28752883D+00, BeCl 
     3  4.76712334D+00, 5.54682545D+00, 6.23634925D+00, 6.78967262D+00, BeCl 
     4  7.12293595D+00, 7.42705048D+00, 7.69474033D+00, 7.96376537D+00, BeCl 
     5  8.19074528D+00, 8.40912435D+00, 8.81981195D+00, 9.24741033D+00, BeCl 
     6  9.43158436D+00, 9.60454592D+00, 9.77210796D+00, 9.88067790D+00, BeCl 
     7  9.95966327D+00, 1.00414315D+01, 1.01365677D+01, 1.01795118D+01, BeCl 
     8  1.01970130D+01, 1.02096253D+01,     60*0.0D+00/                 BeCl 
      DATA TK_BCl/                                                      071215
     1  0.699999789529, 0.709700041799, 0.720199801620, 0.752500051974, BCl  
     2  0.834500022290, 0.930599856194, 1.045500076328, 1.194399936842, BCl  
     3  1.355000026477, 1.516900047300, 1.680699853970, 1.795900093969, BCl  
     4  1.916799997429, 2.034799843208, 2.160899876547, 2.272000311545, BCl  
     5  2.384100214698, 2.601200286434, 2.840200171868, 3.042700035252, BCl  
     6  3.315000053544, 3.428400261164, 3.545600091299, 3.667300004115, BCl  
     7  3.810500370817, 3.909699837485, 4.000000000000,     63*0.0D+00/ BCl  
      DATA  K_BCl/                                                      071215
     1 -8.04191104D-05, 1.64548357D-01, 3.39143696D-01, 8.53629425D-01, BCl  
     2  2.02120280D+00, 3.17572856D+00, 4.31484495D+00, 5.48942700D+00, BCl  
     3  6.46987572D+00, 7.23877734D+00, 7.85613889D+00, 8.21975524D+00, BCl  
     4  8.55459360D+00, 8.84582099D+00, 9.12736688D+00, 9.35566294D+00, BCl  
     5  9.57073582D+00, 9.94933230D+00, 1.03100915D+01, 1.05678977D+01, BCl  
     6  1.08459529D+01, 1.09405119D+01, 1.10259053D+01, 1.10990952D+01, BCl  
     7  1.11538261D+01, 1.11630759D+01, 1.11536953D+01,     63*0.0D+00/ BCl  
      DATA TK_NaCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, NaCl 
     2  0.830499930595, 0.924999965064, 1.032000168042, 1.104800020971, NaCl 
     3  1.178399873043, 1.321300133012, 1.478799986278, 1.637499924388, NaCl 
     4  1.743399850719, 1.856800063080, 1.962599987761, 2.059300412495, NaCl 
     5  2.345800294671, 2.566799928564, 2.752500048735, 2.920400084325, NaCl 
     6  3.103800105449, 3.386700274140, 3.502200041363, 3.614999778302, NaCl 
     7  3.659999855099, 3.703599889674, 3.744399862616, 3.782199794158, NaCl 
     8  3.846400288519, 3.874699999285, 3.903299693759, 3.936899826022, NaCl 
     9  3.963700129656, 3.985999687011, 4.000000000000,     55*0.0D+00/ NaCl 
      DATA  K_NaCl/                                                     071215
     1 -1.76362179D-04, 1.56411253D-01, 3.17940068D-01, 8.03765590D-01, NaCl 
     2  1.90308138D+00, 2.99742040D+00, 4.01367168D+00, 4.59558828D+00, NaCl 
     3  5.10977828D+00, 5.93717448D+00, 6.65199162D+00, 7.22511246D+00, NaCl 
     4  7.54827570D+00, 7.85483743D+00, 8.11138998D+00, 8.32487919D+00, NaCl 
     5  8.85427562D+00, 9.17154015D+00, 9.39421136D+00, 9.57042076D+00, NaCl 
     6  9.73858730D+00, 9.94757960D+00, 1.00145352D+01, 1.00704655D+01, NaCl 
     7  1.00915735D+01, 1.01133596D+01, 1.01377635D+01, 1.01678843D+01, NaCl 
     8  1.02516575D+01, 1.03081491D+01, 1.03800100D+01, 1.04833307D+01, NaCl 
     9  1.05786314D+01, 1.06648244D+01, 1.07214218D+01,     55*0.0D+00/ NaCl 
      DATA TK_MgCl/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719899799358, 0.751800033832, MgCl 
     2  0.832899985612, 0.929299861383, 1.038200002388, 1.112700070761, MgCl 
     3  1.188200006647, 1.335100051857, 1.498099952317, 1.662099881427, MgCl 
     4  1.768999995582, 1.883600034822, 1.986499854428, 2.085100038599, MgCl 
     5  2.371699928509, 2.675200199643, 2.835800072980, 3.015500349200, MgCl 
     6  3.384600220961, 3.519300440391, 3.645699925527, 3.748699958265, MgCl 
     7  3.797900143680, 3.845200261523, 3.894800002835, 3.940899622913, MgCl 
     8  3.976699789423, 3.989499765589, 4.000000000000,     59*0.0D+00/ MgCl 
      DATA  K_MgCl/                                                     071215
     1  3.98888105D-05, 1.46683166D-01, 3.00820640D-01, 7.57890724D-01, MgCl 
     2  1.79410890D+00, 2.82541803D+00, 3.78004335D+00, 4.32946530D+00, MgCl 
     3  4.81592291D+00, 5.60035514D+00, 6.28340848D+00, 6.83222638D+00, MgCl 
     4  7.13656311D+00, 7.42803421D+00, 7.66553943D+00, 7.87544550D+00, MgCl 
     5  8.40089245D+00, 8.83571090D+00, 9.02720905D+00, 9.21592047D+00, MgCl 
     6  9.52357888D+00, 9.60951599D+00, 9.67766629D+00, 9.72824209D+00, MgCl 
     7  9.75397748D+00, 9.78261726D+00, 9.82091245D+00, 9.87002240D+00, MgCl 
     8  9.92208577D+00, 9.94446577D+00, 9.96450656D+00,     59*0.0D+00/ MgCl 
      DATA TK_AlCl/                                                     071215
     1  0.699999789529, 0.710000049601, 0.721499832753, 0.755400127132, AlCl 
     2  0.791299981873, 0.842000100116, 0.942299982856, 1.064400054251, AlCl 
     3  1.205399915771, 1.358300099704, 1.520699985047, 1.602399951027, AlCl 
     4  1.685699970990, 1.822000136594, 1.956699972178, 2.282899802770, AlCl 
     5  2.483199608146, 2.725699869148, 2.945699715885, 3.207099957319, AlCl 
     6  3.426100215099, 3.525700042375, 3.609899664231, 3.782399798873, AlCl 
     7  3.864699771863, 3.928000278898, 3.972600084024, 3.988499743138, AlCl 
     8  4.000000000000,     61*0.0D+00/                                 AlCl 
      DATA  K_AlCl/                                                     071215
     1  2.22427612D-04, 1.67478983D-01, 3.55527622D-01, 8.84375190D-01, AlCl 
     2  1.40533532D+00, 2.07863852D+00, 3.22373091D+00, 4.34518960D+00, AlCl 
     3  5.35798234D+00, 6.20796174D+00, 6.91578573D+00, 7.22018171D+00, AlCl 
     4  7.50624492D+00, 7.93715270D+00, 8.32920178D+00, 9.14614847D+00, AlCl 
     5  9.54557433D+00, 9.93514943D+00, 1.02201650D+01, 1.04908831D+01, AlCl 
     6  1.06683799D+01, 1.07354426D+01, 1.07852267D+01, 1.08627307D+01, AlCl 
     7  1.08870167D+01, 1.09050920D+01, 1.09228685D+01, 1.09312287D+01, AlCl 
     8  1.09381646D+01,     61*0.0D+00/                                 AlCl 
      DATA TK_SiCl/                                                     071215
     1  0.699999789529, 0.710400039490, 0.723199873466, 0.759300228206, SiCl 
     2  0.797800128595, 0.852099955265, 0.962399992968, 1.084100020987, SiCl 
     3  1.193699953910, 1.305800047373, 1.409999850418, 1.506800152061, SiCl 
     4  1.583500004763, 1.660599846053, 1.813500030319, 1.954799924230, SiCl 
     5  2.108399791923, 2.259300218702, 2.455599928004, 2.678200267572, SiCl 
     6  2.867999842873, 3.091600177717, 3.289799969268, 3.417300021879, SiCl 
     7  3.545600091299, 3.666999997991, 3.783099815377, 3.866799819845, SiCl 
     8  3.920200105603, 3.961400078111, 3.985499675785, 4.000000000000, SiCl 
     9      58*0.0D+00/                                                 SiCl 
      DATA  K_SiCl/                                                     071215
     1 -2.29853875D-04, 1.60258738D-01, 3.53047028D-01, 8.70166735D-01, SiCl 
     2  1.38103391D+00, 2.03698965D+00, 3.16997030D+00, 4.17043225D+00, SiCl 
     3  4.89971456D+00, 5.51908828D+00, 6.01081567D+00, 6.41690276D+00, SiCl 
     4  6.71411974D+00, 6.99631037D+00, 7.51597325D+00, 7.95322296D+00, SiCl 
     5  8.38110569D+00, 8.75070537D+00, 9.15623855D+00, 9.52722407D+00, SiCl 
     6  9.78686221D+00, 1.00404361D+01, 1.02264223D+01, 1.03319561D+01, SiCl 
     7  1.04302783D+01, 1.05161999D+01, 1.05886514D+01, 1.06317163D+01, SiCl 
     8  1.06547012D+01, 1.06710790D+01, 1.06809097D+01, 1.06872349D+01, SiCl 
     9      58*0.0D+00/                                                 SiCl 
      DATA TK_PCl/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, PCl  
     2  0.830599932887, 0.925299957830, 1.032700149339, 1.107100073992, PCl  
     3  1.182699892130, 1.328699964517, 1.483999921594, 1.640399871564, PCl  
     4  1.765400093009, 1.898900057774, 2.015400356721, 2.131000113996, PCl  
     5  2.261000256843, 2.430100276382, 2.656899784927, 2.847600332414, PCl  
     6  3.029999745273, 3.230899553570, 3.410699848212, 3.649899628828, PCl  
     7  3.760900245938, 3.861099689607, 3.950499839297, 3.981099577001, PCl  
     8  4.000000000000,     61*0.0D+00/                                 PCl  
      DATA  K_PCl/                                                      071215
     1 -8.92326134D-05, 1.60436099D-01, 3.26018610D-01, 8.23981468D-01, PCl  
     2  1.95171007D+00, 3.07433345D+00, 4.11715337D+00, 4.72354493D+00, PCl  
     3  5.25999086D+00, 6.11419735D+00, 6.82506280D+00, 7.39586379D+00, PCl  
     4  7.77733032D+00, 8.13147345D+00, 8.40670284D+00, 8.65569725D+00, PCl  
     5  8.91185349D+00, 9.21083471D+00, 9.55628612D+00, 9.80474451D+00, PCl  
     6  1.00109523D+01, 1.02049924D+01, 1.03532604D+01, 1.05346393D+01, PCl  
     7  1.06242202D+01, 1.07113684D+01, 1.07944816D+01, 1.08243748D+01, PCl  
     8  1.08433680D+01,     61*0.0D+00/                                 PCl  
      DATA TK_KCl/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750600002732, KCl  
     2  0.829899921553, 0.923699996409, 1.030200216135, 1.104300009445, KCl  
     3  1.179799840554, 1.325300041933, 1.477999967545, 1.632300053886, KCl  
     4  1.734499912679, 1.839700141311, 1.944499933024, 2.040099971693, KCl  
     5  2.286999885871, 2.566599942798, 2.677400249458, 2.799100184130, KCl  
     6  2.953599893474, 3.111999700088, 3.359999634096, 3.456199943522, KCl  
     7  3.555300307149, 3.650299628764, 3.732499609801, 3.799400175128, KCl  
     8  3.832299976826, 3.866799819845, 3.906999776851, 3.948499794511, KCl  
     9  3.981499585982, 4.000000000000,     56*0.0D+00/                 KCl  
      DATA  K_KCl/                                                      071215
     1 -1.35830021D-04, 1.55551512D-01, 3.12976331D-01, 7.94736593D-01, KCl  
     2  1.88473228D+00, 2.96701994D+00, 3.97614947D+00, 4.56670911D+00, KCl  
     3  5.09116180D+00, 5.92594497D+00, 6.61423117D+00, 7.17156378D+00, KCl  
     4  7.48483076D+00, 7.77146432D+00, 8.02696818D+00, 8.23707505D+00, KCl  
     5  8.68952573D+00, 9.07377054D+00, 9.20029333D+00, 9.32795690D+00, KCl  
     6  9.47572155D+00, 9.61115692D+00, 9.78790830D+00, 9.84409489D+00, KCl  
     7  9.89703572D+00, 9.95216664D+00, 1.00234483D+01, 1.01227606D+01, KCl  
     8  1.01924767D+01, 1.02823032D+01, 1.04072823D+01, 1.05545527D+01, KCl  
     9  1.06804353D+01, 1.07530213D+01,     56*0.0D+00/                 KCl  
      DATA TK_CaCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750500000141, CaCl 
     2  0.829599928813, 0.923699996409, 1.030400210791, 1.102599970255, CaCl 
     3  1.175599938022, 1.317200105320, 1.473299857490, 1.630400101202, CaCl 
     4  1.733399939131, 1.843000075988, 1.939400026469, 2.034899845635, CaCl 
     5  2.194699683674, 2.344500265430, 2.586100068362, 2.784899848840, CaCl 
     6  3.158699801606, 3.324699799852, 3.412999908733, 3.505800125291, CaCl 
     7  3.609899664231, 3.711300063429, 3.768700434513, 3.827599868707, CaCl 
     8  3.876700044529, 3.914399957275, 3.948799801285, 3.980799570266, CaCl 
     9  4.000000000000,     57*0.0D+00/                                 CaCl 
      DATA  K_CaCl/                                                     071215
     1 -1.49079900D-04, 1.44247376D-01, 2.90278478D-01, 7.35944982D-01, CaCl 
     2  1.74627305D+00, 2.75692576D+00, 3.69984394D+00, 4.23868714D+00, CaCl 
     3  4.71649851D+00, 5.48874266D+00, 6.16139132D+00, 6.70477207D+00, CaCl 
     4  7.00833236D+00, 7.29676148D+00, 7.52669830D+00, 7.73596852D+00, CaCl 
     5  8.04927342D+00, 8.30362598D+00, 8.64436726D+00, 8.87776983D+00, CaCl 
     6  9.23726565D+00, 9.36654938D+00, 9.42840474D+00, 9.49031275D+00, CaCl 
     7  9.56361447D+00, 9.65677094D+00, 9.72964243D+00, 9.82842819D+00, CaCl 
     8  9.93549364D+00, 1.00354411D+01, 1.01400387D+01, 1.02473313D+01, CaCl 
     9  1.03153795D+01,     57*0.0D+00/                                 CaCl 
      DATA TK_ScCl/                                                     071215
     1  0.699999789529, 0.709400033997, 0.719199817052, 0.750299994957, ScCl 
     2  0.828899945754, 0.922800018110, 1.029300204192, 1.099599920549, ScCl 
     3  1.169700061186, 1.302500117599, 1.456599952324, 1.617999940083, ScCl 
     4  1.735099898251, 1.848299948110, 1.961800008590, 2.064200115556, ScCl 
     5  2.234199649177, 2.349500377895, 2.500800002979, 2.653299696702, ScCl 
     6  2.904699711994, 3.082199987679, 3.259300205155, 3.347900341906, ScCl 
     7  3.440699600004, 3.526599976949, 3.643600073877, 3.738399731689, ScCl 
     8  3.841900187284, 3.944799710970, 3.978199681642, 4.000000000000, ScCl 
     9      58*0.0D+00/                                                 ScCl 
      DATA  K_ScCl/                                                     071215
     1 -2.77618265D-05, 1.64097704D-01, 3.31825024D-01, 8.41783717D-01, ScCl 
     2  1.99288572D+00, 3.14632023D+00, 4.21841001D+00, 4.81447902D+00, ScCl 
     3  5.33583299D+00, 6.16100361D+00, 6.91502345D+00, 7.53980702D+00, ScCl 
     4  7.92104616D+00, 8.25208976D+00, 8.55993478D+00, 8.82267485D+00, ScCl 
     5  9.22816086D+00, 9.47723296D+00, 9.76829299D+00, 1.00225106D+01, ScCl 
     6  1.03706173D+01, 1.05720594D+01, 1.07411739D+01, 1.08152846D+01, ScCl 
     7  1.08873294D+01, 1.09521460D+01, 1.10472885D+01, 1.11390077D+01, ScCl 
     8  1.12567154D+01, 1.13881498D+01, 1.14336424D+01, 1.14641578D+01, ScCl 
     9      58*0.0D+00/                                                 ScCl 
      DATA TK_MnCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750600002732, MnCl 
     2  0.829799923973, 0.924099986765, 1.030900197432, 1.103199984087, MnCl 
     3  1.176199924098, 1.318100123735, 1.474999897297, 1.632900038944, MnCl 
     4  1.737599838133, 1.849399921569, 1.951399838429, 2.046100109044, MnCl 
     5  2.348500355402, 2.572899765956, 2.750499999493, 2.976399803873, MnCl 
     6  3.311400309575, 3.408699803502, 3.506500141611, 3.612099707943, MnCl 
     7  3.705099923825, 3.791400007405, 3.873699976663, 3.950699843664, MnCl 
     8  3.981199579246, 4.000000000000,     60*0.0D+00/                 MnCl 
      DATA  K_MnCl/                                                     071215
     1  2.90508997D-05, 1.52643059D-01, 3.06966549D-01, 7.79278627D-01, MnCl 
     2  1.84693187D+00, 2.91411827D+00, 3.90653839D+00, 4.47277186D+00, MnCl 
     3  4.97329781D+00, 5.78158394D+00, 6.48413409D+00, 7.04824360D+00, MnCl 
     4  7.36516503D+00, 7.66581295D+00, 7.91297029D+00, 8.12308484D+00, MnCl 
     5  8.68654559D+00, 9.01281398D+00, 9.23123472D+00, 9.47350067D+00, MnCl 
     6  9.77230164D+00, 9.84752977D+00, 9.91983336D+00, 9.99930380D+00, MnCl 
     7  1.00814463D+01, 1.01838312D+01, 1.03205614D+01, 1.04910291D+01, MnCl 
     8  1.05697908D+01, 1.06212182D+01,     60*0.0D+00/                 MnCl 
      DATA TK_FeCl/                                                     071215
     1  0.699999789529, 0.709400033997, 0.719099819580, 0.749999987182, FeCl 
     2  0.828299960275, 0.921200056689, 1.026800142454, 1.097399976820, FeCl 
     3  1.168600036276, 1.307400013324, 1.463699935199, 1.620999912182, FeCl 
     4  1.739199799657, 1.865100007730, 1.990299779563, 2.115999799917, FeCl 
     5  2.382200173861, 2.489999782419, 2.603700103918, 2.730299564640, FeCl 
     6  2.891700223334, 3.015500349200, 3.150599617126, 3.279699765488, FeCl 
     7  3.445099691596, 3.588100115220, 3.715000143240, 3.885500249731, FeCl 
     8  3.955299944108, 3.982999619658, 4.000000000000,     59*0.0D+00/ FeCl 
      DATA  K_FeCl/                                                     071215
     1 -1.04124388D-04, 1.49750967D-01, 3.01376856D-01, 7.64465734D-01, FeCl 
     2  1.81426539D+00, 2.86245539D+00, 3.84336422D+00, 4.39759471D+00, FeCl 
     3  4.88844666D+00, 5.68639506D+00, 6.39420211D+00, 6.96219557D+00, FeCl 
     4  7.32075657D+00, 7.65542701D+00, 7.95173802D+00, 8.22153396D+00, FeCl 
     5  8.72933464D+00, 8.91853494D+00, 9.11070925D+00, 9.31448782D+00, FeCl 
     6  9.55111368D+00, 9.70788803D+00, 9.85014598D+00, 9.96040073D+00, FeCl 
     7  1.00809186D+01, 1.01875242D+01, 1.03019979D+01, 1.05196896D+01, FeCl 
     8  1.06408334D+01, 1.06951986D+01, 1.07304380D+01,     59*0.0D+00/ FeCl 
      DATA TK_CuCl/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719799801886, 0.751700031241, CuCl 
     2  0.832699981027, 0.928899871028, 1.037600018419, 1.112000088931, CuCl 
     3  1.187499992072, 1.334300033514, 1.496699922952, 1.660399841336, CuCl 
     4  1.768800000994, 1.885499983844, 1.993599854557, 2.094200244787, CuCl 
     5  2.360099648090, 2.675500206436, 2.801900244528, 2.929600283528, CuCl 
     6  3.173400147584, 3.319499733505, 3.452399851632, 3.533899818538, CuCl 
     7  3.614899775876, 3.698499776669, 3.772800263319, 3.867799842694, CuCl 
     8  3.958100005248, 3.983599633128, 4.000000000000,     59*0.0D+00/ CuCl 
      DATA  K_CuCl/                                                     071215
     1  4.69096212D-05, 1.63227249D-01, 3.33061506D-01, 8.41410064D-01, CuCl 
     2  1.99092909D+00, 3.13143349D+00, 4.18409091D+00, 4.78816986D+00, CuCl 
     3  5.32197642D+00, 6.17746176D+00, 6.91387122D+00, 7.50069824D+00, CuCl 
     4  7.82793048D+00, 8.13963091D+00, 8.39868876D+00, 8.61846374D+00, CuCl 
     5  9.11312035D+00, 9.56246183D+00, 9.71224185D+00, 9.85000645D+00, CuCl 
     6  1.00779617D+01, 1.01937838D+01, 1.02892030D+01, 1.03463127D+01, CuCl 
     7  1.04040769D+01, 1.04637949D+01, 1.05127761D+01, 1.05620454D+01, CuCl 
     8  1.05964668D+01, 1.06066903D+01, 1.06141253D+01,     59*0.0D+00/ CuCl 
      DATA TK_ZnCl/                                                     071215
     1  0.699999789529, 0.709700041799, 0.719999796830, 0.752200044199, ZnCl 
     2  0.833900008536, 0.930799860090, 1.040399963170, 1.117499946166, ZnCl 
     3  1.195499910021, 1.344100062252, 1.507100159115, 1.671900024004, ZnCl 
     4  1.772600033798, 1.874800002769, 2.045600097598, 2.349100368897, ZnCl 
     5  2.679000285687, 2.828499906203, 2.975999833059, 3.456299945940, ZnCl 
     6  3.712600091471, 3.829199907388, 3.918400060049, 3.969300255156, ZnCl 
     7  3.987499720687, 4.000000000000,     64*0.0D+00/                 ZnCl 
      DATA  K_ZnCl/                                                     071215
     1  2.37376571D-05, 1.51750085D-01, 3.09554855D-01, 7.81713236D-01, ZnCl 
     2  1.84870057D+00, 2.90655392D+00, 3.88514008D+00, 4.46214481D+00, ZnCl 
     3  4.96957050D+00, 5.76750800D+00, 6.45306439D+00, 7.00604836D+00, ZnCl 
     4  7.29418935D+00, 7.55745117D+00, 7.94606005D+00, 8.50838221D+00, ZnCl 
     5  8.96628733D+00, 9.13768388D+00, 9.29044200D+00, 9.68672140D+00, ZnCl 
     6  9.83860477D+00, 9.89353996D+00, 9.93447990D+00, 9.96276795D+00, ZnCl 
     7  9.97505865D+00, 9.98451057D+00,     64*0.0D+00/                 ZnCl 
      DATA TK_GaCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750700005324, GaCl 
     2  0.830099921425, 0.924199984353, 1.031000194760, 1.106400057855, GaCl 
     3  1.183299904623, 1.330899955553, 1.482899947109, 1.639799867110, GaCl 
     4  1.772100021244, 1.903899984078, 2.036599886891, 2.151999669525, GaCl 
     5  2.258300195419, 2.367799834641, 2.466600175078, 2.614499747568, GaCl 
     6  2.767500414557, 2.898599710380, 3.030399754721, 3.216400172695, GaCl 
     7  3.394600021836, 3.558900382279, 3.709500024004, 3.767400403084, GaCl 
     8  3.826699846948, 3.879000096560, 3.928400287785, 3.971800141507, GaCl 
     9  3.988299738648, 4.000000000000,     56*0.0D+00/                 GaCl 
      DATA  K_GaCl/                                                     071215
     1 -4.59571784D-05, 1.60382848D-01, 3.24229306D-01, 8.20463266D-01, GaCl 
     2  1.94408271D+00, 3.06072099D+00, 4.09993843D+00, 4.71554916D+00, GaCl 
     3  5.26099986D+00, 6.12213133D+00, 6.81686538D+00, 7.38990380D+00, GaCl 
     4  7.79205877D+00, 8.13891947D+00, 8.44617934D+00, 8.68447237D+00, GaCl 
     5  8.88275006D+00, 9.06920121D+00, 9.22684233D+00, 9.45557810D+00, GaCl 
     6  9.69167445D+00, 9.89013199D+00, 1.00779143D+01, 1.03116421D+01, GaCl 
     7  1.04964780D+01, 1.06331265D+01, 1.07284779D+01, 1.07562934D+01, GaCl 
     8  1.07793750D+01, 1.07963019D+01, 1.08122221D+01, 1.08301656D+01, GaCl 
     9  1.08389614D+01, 1.08461164D+01,     56*0.0D+00/                 GaCl 
      DATA TK_GeCl/                                                     071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749499975714, GeCl 
     2  0.826999991737, 0.918900060997, 1.023500060960, 1.094500050997, GeCl 
     3  1.166499988722, 1.305600051629, 1.456899958845, 1.612300082609, GeCl 
     4  1.741299807299, 1.858700106665, 2.050000198322, 2.166099979536, GeCl 
     5  2.266400374797, 2.361399679585, 2.453899890057, 2.542300008932, GeCl 
     6  2.735899685733, 2.896399873931, 3.055500323359, 3.209600013344, GeCl 
     7  3.367299827657, 3.499999990074, 3.618999875349, 3.730599570549, GeCl 
     8  3.849000347010, 3.948399792254, 3.979899559490, 4.000000000000, GeCl 
     9      58*0.0D+00/                                                 GeCl 
      DATA  K_GeCl/                                                     071215
     1 -9.69168195D-05, 1.45985290D-01, 2.90785438D-01, 7.38223664D-01, GeCl 
     2  1.75345848D+00, 2.76892107D+00, 3.72291648D+00, 4.27085782D+00, GeCl 
     3  4.75857633D+00, 5.54434612D+00, 6.22099097D+00, 6.77865623D+00, GeCl 
     4  7.16602040D+00, 7.47479334D+00, 7.91287491D+00, 8.15024254D+00, GeCl 
     5  8.34467765D+00, 8.52467176D+00, 8.69913547D+00, 8.86584286D+00, GeCl 
     6  9.22441231D+00, 9.50165160D+00, 9.74628805D+00, 9.95067597D+00, GeCl 
     7  1.01306712D+01, 1.02653384D+01, 1.03765466D+01, 1.04730912D+01, GeCl 
     8  1.05659975D+01, 1.06396183D+01, 1.06650165D+01, 1.06828247D+01, GeCl 
     9      58*0.0D+00/                                                 GeCl 
      DATA TK_AsCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.751000013099, AsCl 
     2  0.830899939764, 0.925699948186, 1.033100138652, 1.107300078603, AsCl 
     3  1.182699892130, 1.328199975901, 1.483999921594, 1.641899906942, AsCl 
     4  1.767100047002, 1.902500019170, 2.022600273039, 2.132200141577, AsCl 
     5  2.289099928435, 2.408699818702, 2.692199633637, 2.896999829326, AsCl 
     6  3.238199730866, 3.407099771144, 3.637500276652, 3.748499953817, AsCl 
     7  3.852800172097, 3.947499771933, 3.979699573861, 4.000000000000, AsCl 
     8      62*0.0D+00/                                                 AsCl 
      DATA  K_AsCl/                                                     071215
     1  2.17318711D-04, 1.67518334D-01, 3.40076092D-01, 8.60522376D-01, AsCl 
     2  2.03742010D+00, 3.20587705D+00, 4.28869258D+00, 4.91595374D+00, AsCl 
     3  5.47042168D+00, 6.35138847D+00, 7.08705981D+00, 7.67848350D+00, AsCl 
     4  8.06880672D+00, 8.43410883D+00, 8.72055463D+00, 8.95698254D+00, AsCl 
     5  9.25879493D+00, 9.46185080D+00, 9.86501407D+00, 1.01052134D+01, AsCl 
     6  1.04321078D+01, 1.05657021D+01, 1.07429784D+01, 1.08391065D+01, AsCl 
     7  1.09388687D+01, 1.10367113D+01, 1.10716286D+01, 1.10941829D+01, AsCl 
     8      62*0.0D+00/                                                 AsCl 
      DATA TK_SeCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750399997549, SeCl 
     2  0.829299936074, 0.923000013288, 1.029300204192, 1.102399965644, SeCl 
     3  1.176499917136, 1.319400150335, 1.472199831733, 1.626800040415, SeCl 
     4  1.752400047655, 1.883600034822, 1.995299893190, 2.114499769408, SeCl 
     5  2.265800361691, 2.449699796475, 2.675200199643, 2.865699785271, SeCl 
     6  3.041600013349, 3.286599895708, 3.508800195232, 3.658399817766, SeCl 
     7  3.823199762332, 3.898699715447, 3.962700107245, 4.000000000000, SeCl 
     8      62*0.0D+00/                                                 SeCl 
      DATA  K_SeCl/                                                     071215
     1  1.09783196D-04, 1.56714552D-01, 3.15065280D-01, 7.96649563D-01, SeCl 
     2  1.88830372D+00, 2.97696509D+00, 3.99127962D+00, 4.57846931D+00, SeCl 
     3  5.09824630D+00, 5.92844547D+00, 6.62645332D+00, 7.19134807D+00, SeCl 
     4  7.57533451D+00, 7.92419077D+00, 8.18944835D+00, 8.44774668D+00, SeCl 
     5  8.74583510D+00, 9.06850093D+00, 9.41060111D+00, 9.66420477D+00, SeCl 
     6  9.88083169D+00, 1.01637146D+01, 1.03990257D+01, 1.05475230D+01, SeCl 
     7  1.07048616D+01, 1.07749975D+01, 1.08337876D+01, 1.08681893D+01, SeCl 
     8      62*0.0D+00/                                                 SeCl 
      DATA TK_BrCl/                                                     071215
     1  0.699999789529, 0.710000049601, 0.721399830359, 0.755200121949, BrCl 
     2  0.790999975101, 0.841500112180, 0.940900017229, 1.063400031114, BrCl 
     3  1.203699879419, 1.354600017601, 1.515200091649, 1.679299853698, BrCl 
     4  1.800300171197, 1.927599885190, 2.038599935428, 2.140600278937, BrCl 
     5  2.395899907823, 2.606499899501, 2.777799889749, 2.991099791708, BrCl 
     6  3.241299804170, 3.324199789011, 3.405299734741, 3.481799574832, BrCl 
     7  3.554700294627, 3.648299741856, 3.739999764743, 3.814000118136, BrCl 
     8  3.897599796505, 3.960100048977, 3.984799660070, 4.000000000000, BrCl 
     9      58*0.0D+00/                                                 BrCl 
      DATA  K_BrCl/                                                     071215
     1  2.59424041D-05, 1.76496728D-01, 3.73177559D-01, 9.29480810D-01, BrCl 
     2  1.47746368D+00, 2.18469249D+00, 3.38141383D+00, 4.56630361D+00, BrCl 
     3  5.62488246D+00, 6.50404109D+00, 7.23064594D+00, 7.81836063D+00, BrCl 
     4  8.18072387D+00, 8.51446291D+00, 8.77456638D+00, 8.99246966D+00, BrCl 
     5  9.45964198D+00, 9.77131663D+00, 9.98746622D+00, 1.02210503D+01, BrCl 
     6  1.04568124D+01, 1.05281763D+01, 1.05945499D+01, 1.06521325D+01, BrCl 
     7  1.06983782D+01, 1.07349906D+01, 1.07348948D+01, 1.07088473D+01, BrCl 
     8  1.06622281D+01, 1.06245286D+01, 1.06105346D+01, 1.06024297D+01, BrCl 
     9      58*0.0D+00/                                                 BrCl 
      DATA TK_RbCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719499809469, 0.750900010507, RbCl 
     2  0.830499930595, 0.924499977120, 1.031300186745, 1.106600062466, RbCl 
     3  1.183599910869, 1.332499992240, 1.485599884480, 1.639299879562, RbCl 
     4  1.748399954100, 1.857000067668, 2.071599723347, 2.304700305176, RbCl 
     5  2.541499990993, 2.670600095485, 2.805400318243, 2.965100158182, RbCl 
     6  3.126800040721, 3.356699883261, 3.461100059210, 3.564700067450, RbCl 
     7  3.646199890206, 3.726199821470, 3.794100064012, 3.830599939798, RbCl 
     8  3.864699771863, 3.958900022716, 3.984699657824, 4.000000000000, RbCl 
     9      58*0.0D+00/                                                 RbCl 
      DATA  K_RbCl/                                                     071215
     1 -1.02164907D-04, 1.56188092D-01, 3.17410809D-01, 8.02319937D-01, RbCl 
     2  1.89958946D+00, 2.98668531D+00, 4.00042514D+00, 4.60072763D+00, RbCl 
     3  5.13445331D+00, 5.98357434D+00, 6.66778284D+00, 7.21863303D+00, RbCl 
     4  7.54888819D+00, 7.83820537D+00, 8.31440418D+00, 8.71120406D+00, RbCl 
     5  9.01765805D+00, 9.15790135D+00, 9.29119156D+00, 9.43481740D+00, RbCl 
     6  9.56428141D+00, 9.71590383D+00, 9.77026522D+00, 9.81819622D+00, RbCl 
     7  9.86080987D+00, 9.92721306D+00, 1.00283444D+01, 1.01075831D+01, RbCl 
     8  1.01984391D+01, 1.05156110D+01, 1.06125328D+01, 1.06708580D+01, RbCl 
     9      58*0.0D+00/                                                 RbCl 
      DATA TK_SrCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750600002732, SrCl 
     2  0.829699926393, 0.923699996409, 1.030400210791, 1.104000002529, SrCl 
     3  1.178599868402, 1.323000094303, 1.476999944129, 1.632200056376, SrCl 
     4  1.821500150568, 1.995399895463, 2.155199745204, 2.298900168316, SrCl 
     5  2.573899788419, 2.686899808915, 2.808600385640, 3.093800230727, SrCl 
     6  3.307900362927, 3.395699941519, 3.484999645799, 3.578199882444, SrCl 
     7  3.677300232987, 3.767100395831, 3.857999805480, 3.902199669056, SrCl 
     8  3.941899645492, 3.977899703198, 4.000000000000,     59*0.0D+00/ SrCl 
      DATA  K_SrCl/                                                     071215
     1  1.83368546D-04, 1.43607172D-01, 2.88656768D-01, 7.32716124D-01, SrCl 
     2  1.73617668D+00, 2.73916178D+00, 3.67622990D+00, 4.22156030D+00, SrCl 
     3  4.70516293D+00, 5.48279900D+00, 6.13843485D+00, 6.67101206D+00, SrCl 
     4  7.20091886D+00, 7.60579395D+00, 7.92248915D+00, 8.16636040D+00, SrCl 
     5  8.54351903D+00, 8.67392829D+00, 8.80358806D+00, 9.07170370D+00, SrCl 
     6  9.24196358D+00, 9.30493757D+00, 9.36734301D+00, 9.43764859D+00, SrCl 
     7  9.53485936D+00, 9.66264971D+00, 9.84698752D+00, 9.95916856D+00, SrCl 
     8  1.00716938D+01, 1.01816024D+01, 1.02518221D+01,     59*0.0D+00/ SrCl 
      DATA TK_YCl/                                                      071215
     1  0.699999789529, 0.709400033997, 0.718899824635, 0.749599978008, YCl  
     2  0.827199986897, 0.919300069952, 1.024000073308, 1.094000063786, YCl  
     3  1.164799950225, 1.301900130367, 1.454199900163, 1.609700132656, YCl  
     4  1.736799857371, 1.863700043074, 1.988499807752, 2.125900005019, YCl  
     5  2.300000195710, 2.419100071099, 2.532299770079, 2.636600259776, YCl  
     6  2.888500313038, 3.023400211932, 3.161599870172, 3.288999950878, YCl  
     7  3.408499799457, 3.492799822034, 3.590100161944, 3.659299838766, YCl  
     8  3.746999920451, 3.869199874682, 3.960300053459, 4.000000000000, YCl  
     9      58*0.0D+00/                                                 YCl  
      DATA  K_YCl/                                                      071215
     1 -1.15898409D-04, 1.63765049D-01, 3.26170738D-01, 8.29345096D-01, YCl  
     2  1.96692274D+00, 3.10250735D+00, 4.16405417D+00, 4.76293804D+00, YCl  
     3  5.29392154D+00, 6.14901114D+00, 6.89495688D+00, 7.49897568D+00, YCl  
     4  7.90807780D+00, 8.26099924D+00, 8.56681635D+00, 8.86646731D+00, YCl  
     5  9.20419214D+00, 9.41709691D+00, 9.61067794D+00, 9.78241706D+00, YCl  
     6  1.01633230D+01, 1.03414524D+01, 1.05026696D+01, 1.06327635D+01, YCl  
     7  1.07418816D+01, 1.08160538D+01, 1.09077035D+01, 1.09830314D+01, YCl  
     8  1.10955119D+01, 1.12793375D+01, 1.14268344D+01, 1.14925559D+01, YCl  
     9      58*0.0D+00/                                                 YCl  
      DATA TK_AgCl/                                                     071215
     1  0.699999789529, 0.709700041799, 0.719999796830, 0.752200044199, AgCl 
     2  0.834000010828, 0.930899862038, 1.040599967607, 1.118899909826, AgCl 
     3  1.198299841750, 1.349399930448, 1.510600211651, 1.673599984879, AgCl 
     4  1.771500006180, 1.868199929466, 2.031499763122, 2.207799996663, AgCl 
     5  2.340300170960, 2.681500199761, 2.855799954144, 3.061200350682, AgCl 
     6  3.300000189067, 3.479199592261, 3.601300286617, 3.712700093628, AgCl 
     7  3.823399767167, 3.906399763377, 3.965100161031, 3.986299693746, AgCl 
     8  4.000000000000,     61*0.0D+00/                                 AgCl 
      DATA  K_AgCl/                                                     071215
     1  1.21205840D-05, 1.63718957D-01, 3.33956533D-01, 8.43133162D-01, AgCl 
     2  1.99389619D+00, 3.13121438D+00, 4.18156127D+00, 4.80754609D+00, AgCl 
     3  5.35672087D+00, 6.21333741D+00, 6.92455404D+00, 7.49586091D+00, AgCl 
     4  7.78742258D+00, 8.04631707D+00, 8.43086700D+00, 8.78348201D+00, AgCl 
     5  9.01006112D+00, 9.47303637D+00, 9.66406381D+00, 9.86090297D+00, AgCl 
     6  1.00529496D+01, 1.01707610D+01, 1.02366523D+01, 1.02846099D+01, AgCl 
     7  1.03214486D+01, 1.03487110D+01, 1.03766294D+01, 1.03906258D+01, AgCl 
     8  1.04012662D+01,     61*0.0D+00/                                 AgCl 
      DATA TK_CdCl/                                                     071215
     1  0.699999789529, 0.709700041799, 0.720099799225, 0.752300046791, CdCl 
     2  0.834200015413, 0.930899862038, 1.040999976483, 1.119299899443, CdCl 
     3  1.199099822244, 1.350999937717, 1.502500050953, 1.675599938851, CdCl 
     4  1.851199934619, 2.016000370843, 2.183900026926, 2.326099842718, CdCl 
     5  2.530299721435, 2.699699801014, 2.941899632241, 3.142100192957, CdCl 
     6  3.348300350903, 3.659699848099, 3.762400282203, 3.825899827607, CdCl 
     7  3.886100263973, 3.955299944108, 3.983099621903, 4.000000000000, CdCl 
     8      62*0.0D+00/                                                 CdCl 
      DATA  K_CdCl/                                                     071215
     1 -9.32393546D-05, 1.50487973D-01, 3.08607914D-01, 7.77131920D-01, CdCl 
     2  1.83840275D+00, 2.88604800D+00, 3.86163788D+00, 4.44244910D+00, CdCl 
     3  4.95582024D+00, 5.75969349D+00, 6.39206289D+00, 6.97072759D+00, CdCl 
     4  7.45026518D+00, 7.83076428D+00, 8.16256976D+00, 8.40356188D+00, CdCl 
     5  8.69441026D+00, 8.89891173D+00, 9.15271908D+00, 9.33434929D+00, CdCl 
     6  9.49636306D+00, 9.69610943D+00, 9.75058863D+00, 9.78285007D+00, CdCl 
     7  9.81490915D+00, 9.85947564D+00, 9.88229609D+00, 9.89834870D+00, CdCl 
     8      62*0.0D+00/                                                 CdCl 
      DATA TK_InCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751300020874, InCl 
     2  0.831599955811, 0.926399931307, 1.034600098574, 1.110700122675, InCl 
     3  1.188200006647, 1.335800067908, 1.484299914635, 1.649800093269, InCl 
     4  1.812900015349, 1.964799930481, 2.206699967497, 2.362299701390, InCl 
     5  2.561400312895, 2.725599876380, 2.927500238058, 3.016400369826, InCl 
     6  3.102700184958, 3.344200258683, 3.583400005104, 3.774600133062, InCl 
     7  3.862399719310, 3.934899970170, 3.974799925945, 3.989099756609, InCl 
     8  4.000000000000,     61*0.0D+00/                                 InCl 
      DATA  K_InCl/                                                     071215
     1  1.02888148D-04, 1.59543315D-01, 3.25635294D-01, 8.24722847D-01, InCl 
     2  1.95180528D+00, 3.06598399D+00, 4.10697477D+00, 4.72032396D+00, InCl 
     3  5.26257296D+00, 6.11225070D+00, 6.78452711D+00, 7.38215454D+00, InCl 
     4  7.86221493D+00, 8.23992775D+00, 8.73476387D+00, 8.99288406D+00, InCl 
     5  9.26766617D+00, 9.46237160D+00, 9.68780500D+00, 9.78954088D+00, InCl 
     6  9.89132432D+00, 1.01784554D+01, 1.04235411D+01, 1.05605314D+01, InCl 
     7  1.06017822D+01, 1.06329090D+01, 1.06545719D+01, 1.06639767D+01, InCl 
     8  1.06718985D+01,     61*0.0D+00/                                 InCl 
      DATA TK_SnCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750399997549, SnCl 
     2  0.829299936074, 0.922700020521, 1.028800191845, 1.172900000680, SnCl 
     3  1.313900037796, 1.470399789584, 1.628200071368, 1.813300025329, SnCl 
     4  1.974399887947, 2.217300196840, 2.351800255316, 2.456799954790, SnCl 
     5  2.556000320161, 2.731999601401, 2.883900200563, 3.029199801837, SnCl 
     6  3.152299655844, 3.262500280462, 3.392800153263, 3.528799817019, SnCl 
     7  3.763600311214, 3.857199861882, 3.936199876474, 3.975499875647, SnCl 
     8  3.989199758854, 4.000000000000,     60*0.0D+00/                 SnCl 
      DATA  K_SnCl/                                                     071215
     1 -1.64243202D-05, 1.40226690D-01, 2.82066700D-01, 7.13661035D-01, SnCl 
     2  1.69344527D+00, 2.67027252D+00, 3.58529234D+00, 4.57271207D+00, SnCl 
     3  5.32762366D+00, 5.99103612D+00, 6.52875634D+00, 7.04463642D+00, SnCl 
     4  7.42339454D+00, 7.90130362D+00, 8.12267407D+00, 8.27650718D+00, SnCl 
     5  8.40966512D+00, 8.63126411D+00, 8.82820388D+00, 9.03432552D+00, SnCl 
     6  9.21879059D+00, 9.38267683D+00, 9.56538242D+00, 9.73790420D+00, SnCl 
     7  9.98410453D+00, 1.00610245D+01, 1.01199004D+01, 1.01502870D+01, SnCl 
     8  1.01615889D+01, 1.01708875D+01,     60*0.0D+00/                 SnCl 
      DATA TK_SbCl/                                                     071215
     1  0.699999789529, 0.709700041799, 0.719999796830, 0.752200044199, SbCl 
     2  0.834000010828, 0.930999863987, 1.040599967607, 1.115599995485, SbCl 
     3  1.191500007551, 1.339000141283, 1.502700055655, 1.668400029998, SbCl 
     4  1.777300151803, 1.890799880604, 1.996199913643, 2.100000390645, SbCl 
     5  2.454799910147, 2.616199786702, 2.786799896978, 2.999699992979, SbCl 
     6  3.266900386114, 3.380800124734, 3.522700260461, 3.619999899611, SbCl 
     7  3.728299675953, 3.837500090088, 3.902899684777, 3.954899935374, SbCl 
     8  3.983399628638, 4.000000000000,     60*0.0D+00/                 SbCl 
      DATA  K_SbCl/                                                     071215
     1  1.14884800D-04, 1.59468445D-01, 3.25188400D-01, 8.20914054D-01, SbCl 
     2  1.94167308D+00, 3.05107914D+00, 4.07454068D+00, 4.66167834D+00, SbCl 
     3  5.17927191D+00, 6.00934183D+00, 6.72763719D+00, 7.30372297D+00, SbCl 
     4  7.62340065D+00, 7.91946632D+00, 8.16711965D+00, 8.38906153D+00, SbCl 
     5  8.99709864D+00, 9.19796624D+00, 9.36579195D+00, 9.53061441D+00, SbCl 
     6  9.70063512D+00, 9.76925372D+00, 9.86175760D+00, 9.93553200D+00, SbCl 
     7  1.00304820D+01, 1.01375507D+01, 1.02047004D+01, 1.02590020D+01, SbCl 
     8  1.02892533D+01, 1.03071520D+01,     60*0.0D+00/                 SbCl 
      DATA TK_ICl/                                                      071215
     1  0.699999789529, 0.710100047074, 0.721699837543, 0.756000142682, ICl  
     2  0.792099999931, 0.843600061511, 0.949599803626, 1.068400146797, ICl  
     3  1.146900028895, 1.224899905677, 1.373099948478, 1.535000104435, ICl  
     4  1.700899835496, 1.806200032394, 1.915799972980, 2.021100383780, ICl  
     5  2.116599812120, 2.367499827372, 2.576099837839, 2.757000159529, ICl  
     6  2.927600240223, 3.188699658127, 3.263500304474, 3.340200168711, ICl  
     7  3.411199861369, 3.481199561525, 3.554000280018, 3.630800138530, ICl  
     8  3.720500216444, 3.787299914398, 3.859999664473, 3.945399724518, ICl  
     9  3.979099616974, 4.000000000000,     56*0.0D+00/                 ICl  
      DATA  K_ICl/                                                      071215
     1  9.10604240D-05, 1.77231338D-01, 3.76066239D-01, 9.36562944D-01, ICl  
     2  1.48475964D+00, 2.19939169D+00, 3.45479242D+00, 4.58219127D+00, ICl  
     3  5.19881938D+00, 5.72978113D+00, 6.56040210D+00, 7.26870826D+00, ICl  
     4  7.84465084D+00, 8.15390218D+00, 8.44086012D+00, 8.68928245D+00, ICl  
     5  8.89461960D+00, 9.35315172D+00, 9.65730917D+00, 9.87990028D+00, ICl  
     6  1.00640007D+01, 1.03036055D+01, 1.03632968D+01, 1.04201355D+01, ICl  
     7  1.04673364D+01, 1.05050396D+01, 1.05280099D+01, 1.05257975D+01, ICl  
     8  1.04860642D+01, 1.04383103D+01, 1.03810770D+01, 1.03210000D+01, ICl  
     9  1.03021636D+01, 1.02924075D+01,     56*0.0D+00/                 ICl  
      DATA TK_CsCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719299814524, 0.750399997549, CsCl 
     2  0.829299936074, 0.922500025344, 1.028500184436, 1.102799974865, CsCl 
     3  1.178699866081, 1.325200044210, 1.475199901980, 1.628500078000, CsCl 
     4  1.754600103088, 1.877200062143, 2.195399697697, 2.335700062283, CsCl 
     5  2.553700267046, 2.691299613552, 2.833400018952, 3.011800264401, CsCl 
     6  3.196499716183, 3.338700131785, 3.437099790688, 3.537099891361, CsCl 
     7  3.615499790433, 3.693899681480, 3.757600167300, 3.792800036757, CsCl 
     8  3.826199834860, 3.921300130042, 3.967200208093, 4.000000000000, CsCl 
     9      58*0.0D+00/                                                 CsCl 
      DATA  K_CsCl/                                                     071215
     1 -2.95760310D-06, 1.56010194D-01, 3.13763775D-01, 7.93539172D-01, CsCl 
     2  1.88114027D+00, 2.96057931D+00, 3.96972610D+00, 4.56476277D+00, CsCl 
     3  5.09418675D+00, 5.93701655D+00, 6.61555600D+00, 7.17215598D+00, CsCl 
     4  7.55424389D+00, 7.87460277D+00, 8.51850402D+00, 8.73184978D+00, CsCl 
     5  9.00254533D+00, 9.14785131D+00, 9.28443208D+00, 9.43960704D+00, CsCl 
     6  9.58006583D+00, 9.67237690D+00, 9.72850414D+00, 9.78364175D+00, CsCl 
     7  9.83523932D+00, 9.91423754D+00, 1.00198135D+01, 1.00998674D+01, CsCl 
     8  1.01905528D+01, 1.05083166D+01, 1.06778563D+01, 1.08001089D+01, CsCl 
     9      58*0.0D+00/                                                 CsCl 
      DATA TK_BaCl/                                                     071215
     1  0.699999789529, 0.709400033997, 0.718999822108, 0.749699980301, BaCl 
     2  0.827499979636, 0.919300069952, 1.024600088125, 1.096300004956, BaCl 
     3  1.169000045334, 1.307300015452, 1.449999808880, 1.610300132619, BaCl 
     4  1.770599983583, 1.928699857621, 2.152199674255, 2.267500398825, BaCl 
     5  2.384300218997, 2.623099945235, 2.745999902066, 2.873399965124, BaCl 
     6  3.038799953145, 3.188599665379, 3.292000016906, 3.417000013985, BaCl 
     7  3.462000078681, 3.510500234884, 3.571899731935, 3.642100179841, BaCl 
     8  3.702099855522, 3.768500429678, 3.890200341806, 3.954099917906, BaCl 
     9  3.982699612923, 4.000000000000,     56*0.0D+00/                 BaCl 
      DATA  K_BaCl/                                                     071215
     1 -6.06971025D-06, 1.41593536D-01, 2.83421416D-01, 7.18505599D-01, BaCl 
     2  1.70631937D+00, 2.69001087D+00, 3.62180996D+00, 4.15857621D+00, BaCl 
     3  4.63618457D+00, 5.39522919D+00, 6.02122419D+00, 6.58851050D+00, BaCl 
     4  7.05446135D+00, 7.44261224D+00, 7.89639675D+00, 8.09198968D+00, BaCl 
     5  8.26656001D+00, 8.56607676D+00, 8.69942389D+00, 8.82728559D+00, BaCl 
     6  8.97944517D+00, 9.10501340D+00, 9.18958063D+00, 9.30913107D+00, BaCl 
     7  9.36358722D+00, 9.43192356D+00, 9.53331049D+00, 9.66631433D+00, BaCl 
     8  9.78961961D+00, 9.93288922D+00, 1.02185800D+01, 1.03852703D+01, BaCl 
     9  1.04632875D+01, 1.05112400D+01,     56*0.0D+00/                 BaCl 
      DATA TK_YbCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719399811997, 0.750800007916, YbCl 
     2  0.830399928302, 0.924599974709, 1.031400184073, 1.106000048634, YbCl 
     3  1.181899875473, 1.328499969070, 1.481899970305, 1.636599946801, YbCl 
     4  1.734999900656, 1.835000030689, 2.019800460282, 2.158399820882, YbCl 
     5  2.297000121000, 2.575099815375, 2.685599903007, 2.803800284545, YbCl 
     6  3.146199887016, 3.365799787884, 3.460700050556, 3.556300328018, YbCl 
     7  3.649199678278, 3.731299585010, 3.797400133197, 3.865799796996, YbCl 
     8  3.948599796769, 3.981299581491, 4.000000000000,     59*0.0D+00/ YbCl 
      DATA  K_YbCl/                                                     071215
     1  6.55705035D-05, 1.43768274D-01, 2.90566175D-01, 7.36768003D-01, YbCl 
     2  1.74758552D+00, 2.75298357D+00, 3.69090986D+00, 4.24309324D+00, YbCl 
     3  4.73365036D+00, 5.51868915D+00, 6.16762221D+00, 6.69591997D+00, YbCl 
     4  6.98372583D+00, 7.24650186D+00, 7.66860093D+00, 7.93843220D+00, YbCl 
     5  8.17190778D+00, 8.54993284D+00, 8.67610996D+00, 8.80091082D+00, YbCl 
     6  9.11392731D+00, 9.27754436D+00, 9.33897937D+00, 9.39613355D+00, YbCl 
     7  9.45400546D+00, 9.52602909D+00, 9.62113567D+00, 9.77320483D+00, YbCl 
     8  1.00275250D+01, 1.01417787D+01, 1.02089331D+01,     59*0.0D+00/ YbCl 
      DATA TK_AuCl/                                                     071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751300020874, AuCl 
     2  0.831699958103, 0.927099914429, 1.035000087887, 1.108500106266, AuCl 
     3  1.182899896294, 1.327599989563, 1.487899831130, 1.649700090910, AuCl 
     4  1.761300203968, 1.878900104199, 1.986399856762, 2.093700232213, AuCl 
     5  2.228399659637, 2.362599708658, 2.569199757750, 2.711900098875, AuCl 
     6  2.894700000311, 3.156299746945, 3.316199968200, 3.512400279255, AuCl 
     7  3.582199976990, 3.649199678278, 3.815300024283, 3.892200194427, AuCl 
     8  3.963800131897, 4.000000000000,     60*0.0D+00/                 AuCl 
      DATA  K_AuCl/                                                     071215
     1  3.60367640D-05, 1.57059146D-01, 3.20638053D-01, 8.12207974D-01, AuCl 
     2  1.92380393D+00, 3.02789032D+00, 4.05030772D+00, 4.63562127D+00, AuCl 
     3  5.15282175D+00, 5.98493823D+00, 6.70549912D+00, 7.28268822D+00, AuCl 
     4  7.61795879D+00, 7.93001116D+00, 8.18588007D+00, 8.41776207D+00, AuCl 
     5  8.67883265D+00, 8.90780649D+00, 9.20751041D+00, 9.38563291D+00, AuCl 
     6  9.58863338D+00, 9.83793079D+00, 9.97019997D+00, 1.01287515D+01, AuCl 
     7  1.01898905D+01, 1.02520321D+01, 1.04125416D+01, 1.04838206D+01, AuCl 
     8  1.05476980D+01, 1.05803149D+01,     60*0.0D+00/                 AuCl 
      DATA TK_HgCl/                                                     071215
     1  0.699999789529, 0.710100047074, 0.721599835148, 0.755800137498, HgCl 
     2  0.843100073575, 0.948099840454, 1.066000091269, 1.145199991420, HgCl 
     3  1.224399894902, 1.372699939755, 1.525000099416, 1.690600056189, HgCl 
     4  1.839500136604, 1.977799959697, 2.168100019148, 2.318999782431, HgCl 
     5  2.565600013971, 2.677900260779, 2.798800177337, 3.124599985882, HgCl 
     6  3.379200088058, 3.492299810365, 3.601500272143, 3.728399669023, HgCl 
     7  3.819199742725, 3.882100169028, 3.940499613882, 3.976599796608, HgCl 
     8  3.989499765589, 4.000000000000,     60*0.0D+00/                 HgCl 
      DATA  K_HgCl/                                                     071215
     1 -4.02890578D-06, 1.53051633D-01, 3.23447001D-01, 8.06993856D-01, HgCl 
     2  1.89861732D+00, 2.98279042D+00, 3.96516496D+00, 4.51449863D+00, HgCl 
     3  4.99206122D+00, 5.73328785D+00, 6.33874767D+00, 6.87161680D+00, HgCl 
     4  7.27193612D+00, 7.59362682D+00, 7.97030391D+00, 8.22000925D+00, HgCl 
     5  8.55203283D+00, 8.68067207D+00, 8.80828114D+00, 9.10604492D+00, HgCl 
     6  9.29069267D+00, 9.35703875D+00, 9.41203177D+00, 9.46888653D+00, HgCl 
     7  9.50882160D+00, 9.53773826D+00, 9.56683111D+00, 9.58708398D+00, HgCl 
     8  9.59504738D+00, 9.60192462D+00,     60*0.0D+00/                 HgCl 
      DATA TK_TlCl/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719899799358, 0.751800033832, TlCl 
     2  0.832899985612, 0.928799873439, 1.037500021091, 1.113900039612, TlCl 
     3  1.191899997798, 1.342700097068, 1.501700032142, 1.662099881427, TlCl 
     4  1.763100155255, 1.866199979958, 2.063200190094, 2.178400266611, TlCl 
     5  2.310900363998, 2.535099838180, 2.728299681105, 2.883200183448, TlCl 
     6  3.141600230267, 3.277899894468, 3.412699900839, 3.506400139279, TlCl 
     7  3.588900133963, 3.678800268691, 3.755000105681, 3.825299813102, TlCl 
     8  3.893300113369, 3.958200007432, 3.984099644354, 4.000000000000, TlCl 
     9      58*0.0D+00/                                                 TlCl 
      DATA  K_TlCl/                                                     071215
     1  2.92614823D-05, 1.61014637D-01, 3.30193195D-01, 8.31644457D-01, TlCl 
     2  1.96705130D+00, 3.08922325D+00, 4.12895870D+00, 4.74081693D+00, TlCl 
     3  5.28283207D+00, 6.14275986D+00, 6.84924001D+00, 7.41663011D+00, TlCl 
     4  7.71958630D+00, 7.99538046D+00, 8.44602452D+00, 8.66907527D+00, TlCl 
     5  8.89262647D+00, 9.20422394D+00, 9.42576740D+00, 9.58280472D+00, TlCl 
     6  9.81107013D+00, 9.91621689D+00, 1.00148712D+01, 1.00839558D+01, TlCl 
     7  1.01463074D+01, 1.02149369D+01, 1.02721862D+01, 1.03246244D+01, TlCl 
     8  1.03797499D+01, 1.04468635D+01, 1.04808679D+01, 1.05044554D+01, TlCl 
     9      58*0.0D+00/                                                 TlCl 
      DATA TK_PbCl/                                                     071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751400023466, PbCl 
     2  0.831799960396, 0.926899919251, 1.034700095902, 1.111000114888, PbCl 
     3  1.188400010811, 1.335200054150, 1.491699818076, 1.649900095627, PbCl 
     4  1.837400087177, 1.997199936369, 2.165499967653, 2.310200414257, PbCl 
     5  2.560000412536, 2.662899924989, 2.777299927261, 2.993999859579, PbCl 
     6  3.154999717337, 3.249099981692, 3.330099917263, 3.430500257821, PbCl 
     7  3.518000410032, 3.626900053086, 3.771200379104, 3.861899707886, PbCl 
     8  3.944699708712, 3.978599652901, 4.000000000000,     59*0.0D+00/ PbCl 
      DATA  K_PbCl/                                                     071215
     1  1.10907338D-04, 1.41359503D-01, 2.86987859D-01, 7.24871012D-01, PbCl 
     2  1.71664465D+00, 2.70171303D+00, 3.61983669D+00, 4.16686506D+00, PbCl 
     3  4.65062567D+00, 5.41136165D+00, 6.05278311D+00, 6.57623841D+00, PbCl 
     4  7.08505737D+00, 7.45022102D+00, 7.77933944D+00, 8.02083209D+00, PbCl 
     5  8.36232837D+00, 8.48252263D+00, 8.60642903D+00, 8.81853052D+00, PbCl 
     6  8.95882025D+00, 9.03544680D+00, 9.10060328D+00, 9.18564592D+00, PbCl 
     7  9.26952885D+00, 9.39187775D+00, 9.57771052D+00, 9.69949652D+00, PbCl 
     8  9.81401106D+00, 9.86385718D+00, 9.89695416D+00,     59*0.0D+00/ PbCl 
      DATA TK_AlSe/                                                     071215
     1  0.699999789529, 0.710000049601, 0.721199825569, 0.754800111582, AlSe 
     2  0.840300141134, 0.940200034416, 1.060699968646, 1.196499885639, AlSe 
     3  1.345500027436, 1.512700156867, 1.669900065373, 1.820500178517, AlSe 
     4  1.964999925274, 2.228399659637, 2.364699759536, 2.508600210530, AlSe 
     5  2.700999832545, 2.912799899634, 3.067799858558, 3.242199824653, AlSe 
     6  3.489699750033, 3.615599792859, 3.693899681480, 3.769300449019, AlSe 
     7  3.829999926729, 3.886200266346, 3.956699974678, 3.983499630883, AlSe 
     8  4.000000000000,     61*0.0D+00/                                 AlSe 
      DATA  K_AlSe/                                                     071215
     1 -1.45426307D-05, 1.56919456D-01, 3.28839010D-01, 8.21341460D-01, AlSe 
     2  1.93175007D+00, 3.00914788D+00, 4.05884542D+00, 4.99225394D+00, AlSe 
     3  5.79359384D+00, 6.50058820D+00, 7.04854887D+00, 7.51341923D+00, AlSe 
     4  7.92345503D+00, 8.58284030D+00, 8.87250275D+00, 9.13732998D+00, AlSe 
     5  9.43021641D+00, 9.69102241D+00, 9.85810084D+00, 1.00330798D+01, AlSe 
     6  1.02608295D+01, 1.03673250D+01, 1.04301801D+01, 1.04883878D+01, AlSe 
     7  1.05349187D+01, 1.05808070D+01, 1.06511858D+01, 1.06849411D+01, AlSe 
     8  1.07084427D+01,     61*0.0D+00/                                 AlSe 
      DATA TK_SiSe/                                                     071215
     1  0.699999789529, 0.710200044546, 0.722099847123, 0.756900166007, SiSe 
     2  0.793700036047, 0.846100001191, 0.952899869683, 1.071600141629, SiSe 
     3  1.164199936638, 1.258800193273, 1.342700097068, 1.424799953241, SiSe 
     4  1.504800105034, 1.582599982126, 1.741299807299, 1.886199965062, SiSe 
     5  2.025100088472, 2.358099786927, 2.490999802346, 2.634900219497, SiSe 
     6  2.762800301027, 2.914699946087, 3.046600112905, 3.193499645875, SiSe 
     7  3.397499810092, 3.534899841295, 3.657799803766, 3.781299772939, SiSe 
     8  3.837900098800, 3.895599943884, 3.941499636461, 3.982799615168, SiSe 
     9  4.000000000000,     57*0.0D+00/                                 SiSe 
      DATA  K_SiSe/                                                     071215
     1 -2.86676368D-05, 1.67972913D-01, 3.59471379D-01, 8.93081974D-01, SiSe 
     2  1.41707609D+00, 2.09823622D+00, 3.28284227D+00, 4.33972707D+00, SiSe 
     3  5.01599872D+00, 5.60225382D+00, 6.05273187D+00, 6.44513614D+00, SiSe 
     4  6.79287465D+00, 7.10655283D+00, 7.69263882D+00, 8.18071861D+00, SiSe 
     5  8.61477405D+00, 9.51798328D+00, 9.81406286D+00, 1.00900270D+01, SiSe 
     6  1.02996873D+01, 1.05142065D+01, 1.06793196D+01, 1.08482198D+01, SiSe 
     7  1.10657873D+01, 1.12061550D+01, 1.13288415D+01, 1.14430186D+01, SiSe 
     8  1.14877398D+01, 1.15255603D+01, 1.15497268D+01, 1.15683790D+01, SiSe 
     9  1.15759824D+01,     57*0.0D+00/                                 SiSe 
      DATA TK_GeSe/                                                     071215
     1  0.699999789529, 0.709200028796, 0.718499834746, 0.748599955070, GeSe 
     2  0.824600049820, 0.914699966971, 1.017500034254, 1.155999943811, GeSe 
     3  1.290399955641, 1.437299981287, 1.588500130523, 1.721099814453, GeSe 
     4  1.846199998778, 2.147299810856, 2.261400265580, 2.372599949986, GeSe 
     5  2.476899751080, 2.572599759217, 2.965400165080, 3.143300103413, GeSe 
     6  3.312800210008, 3.556300328018, 3.664899955124, 3.775300082406, GeSe 
     7  3.871899935943, 3.958500013982, 3.983799637619, 4.000000000000, GeSe 
     8      62*0.0D+00/                                                 GeSe 
      DATA  K_GeSe/                                                     071215
     1 -1.56020518D-04, 1.53455949D-01, 3.05814219D-01, 7.79023722D-01, GeSe 
     2  1.85061539D+00, 2.92389863D+00, 3.93576406D+00, 5.02040987D+00, GeSe 
     3  5.84336272D+00, 6.55534773D+00, 7.14126469D+00, 7.56704133D+00, GeSe 
     4  7.91362637D+00, 8.60096724D+00, 8.82737273D+00, 9.04129982D+00, GeSe 
     5  9.24102954D+00, 9.42497621D+00, 1.01612410D+01, 1.04613284D+01, GeSe 
     6  1.07179232D+01, 1.10408897D+01, 1.11708529D+01, 1.12924711D+01, GeSe 
     7  1.13868695D+01, 1.14635268D+01, 1.14865137D+01, 1.15020669D+01, GeSe 
     8      62*0.0D+00/                                                 GeSe 
      DATA TK_KBr/                                                      071215
     1  0.699999789529, 0.709500036598, 0.719599806941, 0.751300020874, KBr  
     2  0.831599955811, 0.926399931307, 1.034000114605, 1.110700122675, KBr  
     3  1.189400031633, 1.263400136925, 1.341400129397, 1.494399874709, KBr  
     4  1.649100076759, 1.771600008690, 1.890299869668, 2.163599930022, KBr  
     5  2.523600189403, 2.683100083955, 2.849200367127, 3.142800140723, KBr  
     6  3.386600271607, 3.472800051045, 3.565100038702, 3.654499726765, KBr  
     7  3.730199562285, 3.802000231549, 3.835000035635, 3.867099826700, KBr  
     8  3.912699913596, 3.955299944108, 3.983499630883, 4.000000000000, KBr  
     9      58*0.0D+00/                                                 KBr  
      DATA  K_KBr/                                                      071215
     1 -6.64667582D-05, 1.56166069D-01, 3.18923060D-01, 8.08034102D-01, KBr  
     2  1.91287510D+00, 3.00557730D+00, 4.02205280D+00, 4.62946876D+00, KBr  
     3  5.16998626D+00, 5.61412740D+00, 6.02593686D+00, 6.70033242D+00, KBr  
     4  7.24761158D+00, 7.61057532D+00, 7.91539255D+00, 8.47362172D+00, KBr  
     5  8.97169445D+00, 9.13171652D+00, 9.27124418D+00, 9.47201135D+00, KBr  
     6  9.61206545D+00, 9.65609105D+00, 9.70138320D+00, 9.75089611D+00, KBr  
     7  9.81344681D+00, 9.91637626D+00, 9.98550540D+00, 1.00677177D+01, KBr  
     8  1.02082165D+01, 1.03592189D+01, 1.04661202D+01, 1.05302942D+01, KBr  
     9      58*0.0D+00/                                                 KBr  
      DATA TK_SiTe/                                                     071215
     1  0.699999789529, 0.710200044546, 0.722299851912, 0.757400178965, SiTe 
     2  0.794700058620, 0.847599964999, 0.955399935096, 1.075000051982, SiTe 
     3  1.169700061186, 1.267500040616, 1.353199986535, 1.436899972207, SiTe 
     4  1.518200013387, 1.596599985472, 1.760700220206, 1.904999956506, SiTe 
     5  2.036899894172, 2.213100115231, 2.341600200201, 2.471300157431, SiTe 
     6  2.597900329036, 2.850500347385, 3.042900039234, 3.384500218429, SiTe 
     7  3.489099736726, 3.595000270217, 3.704699914718, 3.812700211989, SiTe 
     8  3.889000332808, 3.956599972495, 3.983599633128, 4.000000000000, SiTe 
     9      58*0.0D+00/                                                 SiTe 
      DATA  K_SiTe/                                                     071215
     1 -8.76790059D-05, 1.60827624D-01, 3.47312078D-01, 8.62619580D-01, SiTe 
     2  1.37090514D+00, 2.02867532D+00, 3.17190998D+00, 4.19005453D+00, SiTe 
     3  4.85094443D+00, 5.42904158D+00, 5.86785794D+00, 6.25033993D+00, SiTe 
     4  6.58955383D+00, 6.89436866D+00, 7.48160472D+00, 7.95507292D+00, SiTe 
     5  8.35852754D+00, 8.85140722D+00, 9.17120660D+00, 9.45559071D+00, SiTe 
     6  9.69556737D+00, 1.00755310D+01, 1.02961237D+01, 1.06166651D+01, SiTe 
     7  1.07154836D+01, 1.08212703D+01, 1.09364541D+01, 1.10513357D+01, SiTe 
     8  1.11299237D+01, 1.11978294D+01, 1.12257514D+01, 1.12434709D+01, SiTe 
     9      58*0.0D+00/                                                 SiTe 
      DATA TK_GeTe/                                                     071215
     1  0.699999789529, 0.709300031396, 0.718699829691, 0.748999964245, GeTe 
     2  0.825700023199, 0.916500007268, 1.020099976997, 1.090700148193, GeTe 
     3  1.162399895877, 1.299900168559, 1.445799907035, 1.599199913468, GeTe 
     4  1.722899856205, 1.845000027732, 2.039899966977, 2.127600040674, GeTe 
     5  2.209900052345, 2.394200032783, 2.490799798360, 2.584400028533, GeTe 
     6  2.729299608780, 2.917500014543, 3.049900178612, 3.179700282634, GeTe 
     7  3.358099777555, 3.508200181244, 3.590900179621, 3.671500094934, GeTe 
     8  3.802900251277, 3.887700301950, 3.961100071388, 3.984799660070, GeTe 
     9  4.000000000000,     57*0.0D+00/                                 GeTe 
      DATA  K_GeTe/                                                     071215
     1  1.61633431D-04, 1.54636111D-01, 3.07792386D-01, 7.81430955D-01, GeTe 
     2  1.85582291D+00, 2.92902244D+00, 3.93958978D+00, 4.52175104D+00, GeTe 
     3  5.03990255D+00, 5.86691285D+00, 6.56131296D+00, 7.14583120D+00, GeTe 
     4  7.53866851D+00, 7.87552088D+00, 8.33463651D+00, 8.51646331D+00, GeTe 
     5  8.67681345D+00, 9.01666839D+00, 9.19237076D+00, 9.36340199D+00, GeTe 
     6  9.62771855D+00, 9.95732207D+00, 1.01690296D+01, 1.03565700D+01, GeTe 
     7  1.05885637D+01, 1.07723326D+01, 1.08720848D+01, 1.09680484D+01, GeTe 
     8  1.11163077D+01, 1.12012012D+01, 1.12688301D+01, 1.12913312D+01, GeTe 
     9  1.13065609D+01,     57*0.0D+00/                                 GeTe 
      DATA TK_KI/                                                       071215
     1  0.699999789529, 0.709600039198, 0.719699804414, 0.751500026057, KI   
     2  0.832099967273, 0.927199912018, 1.035200082543, 1.112800068165, KI   
     3  1.192599980731, 1.268200024173, 1.347599975212, 1.501700032142, KI   
     4  1.655699946319, 1.785200091198, 1.908999856244, 2.254400104618, KI   
     5  2.538399918442, 2.842000210920, 3.117299813463, 3.375500012124, KI   
     6  3.467600199839, 3.561800275871, 3.651799663764, 3.729099620518, KI   
     7  3.801700224972, 3.834900033457, 3.867199828985, 3.913599936720, KI   
     8  3.955499948475, 3.983399628638, 4.000000000000,     59*0.0D+00/ KI   
      DATA  K_KI/                                                       071215
     1  4.51837132D-05, 1.57019370D-01, 3.18823306D-01, 8.06579914D-01, KI   
     2  1.90864073D+00, 2.99748212D+00, 4.01046828D+00, 4.62006000D+00, KI   
     3  5.16290209D+00, 5.61145134D+00, 6.02504932D+00, 6.69460801D+00, KI   
     4  7.23182524D+00, 7.60789258D+00, 7.91539332D+00, 8.55828208D+00, KI   
     5  8.91087455D+00, 9.17896922D+00, 9.36117882D+00, 9.49298710D+00, KI   
     6  9.53203431D+00, 9.57043858D+00, 9.61363242D+00, 9.67248957D+00, KI   
     7  9.77286747D+00, 9.84111723D+00, 9.92280253D+00, 1.00647913D+01, KI   
     8  1.02127308D+01, 1.03182931D+01, 1.03829221D+01,     59*0.0D+00/ KI   
C
C Length of idividual temperature grids
C
      DATA MTQ/ 23, 19, 20, 27, 21, 21, 20, 22, 21, 22, 21, 18, 22, 23,
     *  20, 19, 21, 23, 20, 18, 24, 23, 26, 26, 20, 22, 25, 17, 21, 23,
     *  21, 23, 25, 23, 22, 21, 20, 22, 20, 25, 20, 23, 19, 20, 21, 23,
     *  21, 23, 22, 23, 23, 22, 20, 22, 22, 22, 22, 24, 22, 19, 19, 20,
     *  22, 22, 22, 23, 19, 23, 19, 22, 25, 24, 23, 20, 23, 20, 26, 23,
     *  24, 23, 22, 22, 20, 18, 18, 22, 22, 22, 22, 20, 22, 21, 23, 19,
     *  23, 23, 21, 23, 24, 21, 19, 19, 19, 19, 20, 21, 22, 21, 19, 21,
     *  20, 20, 20, 19, 19, 18, 22, 22, 19, 20, 19, 19, 19, 21, 18, 20,
     *  24, 18, 18, 20, 20, 23, 20, 21, 17, 23, 19, 21, 21, 19, 20, 20,
     *  22, 22, 23, 21, 20, 22, 19, 20, 20, 19, 21, 23, 20, 19, 19, 21,
     *  20, 24, 19, 20, 20, 20, 21, 19, 19, 20, 16, 21, 22, 20, 20, 20,
     *  20, 20, 20, 19, 15, 17, 21, 19, 20, 18, 21, 20, 22, 19, 20, 22,
     *  20, 21, 17, 19, 18, 19, 19, 18, 19, 18, 19, 22, 22, 17, 20, 21,
     *  21, 19, 17, 17, 19, 21, 21, 23, 20, 21, 16, 18, 19, 20, 19, 17,
     *  21, 21, 19, 19, 22, 20, 20, 19, 19, 20, 19, 18, 22, 17, 22, 19,
     *  20, 19, 19, 17, 19, 20, 21, 19, 21, 22, 19, 19, 21, 17, 19, 18,
     *  20, 17, 20, 20, 21, 18, 19, 21, 19, 19, 18, 16, 18, 23, 20, 20,
     *  20, 19, 16, 20, 21, 19, 23, 19, 19, 20, 19, 18, 18, 17, 18, 19,
     *  19, 19, 19, 20, 19/
      DATA MTK/ 29, 32, 30, 33, 32, 28, 33, 36, 34, 35, 32, 31, 32, 33,
     *  35, 30, 27, 31, 28, 32, 35, 38, 32, 28, 28, 32, 29, 23, 30, 29,
     *  27, 30, 30, 28, 26, 24, 25, 28, 28, 33, 33, 29, 27, 26, 27, 27,
     *  35, 33, 30, 31, 31, 31, 30, 29, 31, 28, 33, 33, 26, 27, 27, 35,
     *  34, 31, 31, 34, 30, 29, 28, 29, 31, 36, 34, 31, 29, 29, 32, 30,
     *  27, 26, 27, 29, 30, 25, 24, 27, 28, 28, 29, 28, 28, 27, 27, 26,
     *  29, 28, 27, 31, 27, 27, 26, 28, 26, 28, 30, 30, 26, 30, 28, 30,
     *  29, 31, 31, 30, 29, 26, 29, 30, 26, 28, 30, 31, 29, 26, 33, 30,
     *  30, 32, 28, 30, 32, 32, 26, 29, 34, 36, 29, 33, 30, 31, 27, 36,
     *  34, 31, 31, 32, 33, 31, 31, 30, 31, 32, 31, 33, 30, 28, 37, 33,
     *  34, 33, 27, 33, 29, 29, 31, 26, 34, 29, 30, 32, 32, 28, 34, 28,
     *  32, 27, 90, 28, 29, 27, 29, 30, 29, 32, 31, 28, 29, 36, 32, 29,
     *  32, 31, 29, 31, 32, 32, 31, 30, 31, 30, 33, 33, 32, 28, 34, 31,
     *  32, 31, 32, 30, 32, 32, 31, 34, 31, 32, 28, 32, 32, 32, 27, 30,
     *  32, 35, 28, 30, 29, 27, 32, 29, 35, 29, 32, 31, 31, 32, 32, 31,
     *  29, 29, 32, 33, 31, 29, 34, 31, 32, 28, 31, 30, 27, 35, 31, 29,
     *  32, 29, 34, 33, 32, 30, 31, 31, 26, 34, 32, 28, 28, 32, 32, 31,
     *  32, 29, 28, 29, 30, 30, 34, 32, 34, 31, 30, 30, 32, 31, 29, 33,
     *  28, 32, 32, 33, 31/
C
      DATA FIRST/.TRUE./
C
C Compute 2nd derivatives for spline interpolation
C
      IF(FIRST) THEN
        DO 1 I=1,MSPEC
          CALL SPL_INIT(TQ(1,I),Q(1,I),Q2(1,I),U,MTQ(I))
          CALL SPL_INIT(TK(1,I),K(1,I),K2(1,I),U,MTK(I))
   1    CONTINUE
        FIRST=.FALSE.
      ENDIF
C
C Fits are made in log10 of temperatures
C
      TLOG=LOG10(TEMP)
C
C Find species name
C
      DO 4 II=1,MSPEC
        ISPEC=II
        IF(SPLIST(II).EQ.SPNAME) THEN
C
C  The species is in Barklem's list.
C  Find the braketing temperatures for the partition functions.
C
          KHI=MTQ(ISPEC)
          KLO=1
   2      CONTINUE
          I=(KLO+KHI)/2
          A=TQ(I,ISPEC)
          IF(A.GT.TLOG) THEN
            KHI=I
          ELSE IF(A.LE.TLOG) THEN
            KLO=I
          END IF
          IF(KHI-KLO.GT.1) GO TO 2
C
C Do the interpolation of the partition functions
C
          Q_spln=SPL_INTERP(KLO,KHI,TQ(1,ISPEC),Q(1,ISPEC),Q2(1,ISPEC),
     *                      MTQ(ISPEC),TLOG)
C  Find the braketing temperatures for the partition functions.
C
          KHI=MTK(ISPEC)
          KLO=1
   3      CONTINUE
          I=(KLO+KHI)/2
          A=TK(I,ISPEC)
          IF(A.GT.TLOG) THEN
            KHI=I
          ELSE IF(A.LE.TLOG) THEN
            KLO=I
          END IF
          IF(KHI-KLO.GT.1) GO TO 3
C
C Do the interpolation of the equilibrium constants
C
          K_spln=SPL_INTERP(KLO,KHI,TK(1,ISPEC),K(1,ISPEC),K2(1,ISPEC),
     *                      MTK(ISPEC),TLOG)
C
C The "+1" converts from pascals (N/m^2 as in Barklem tables) to
C dynes/cm^2 as required by the EOS.
C
          K_spln=K_spln+1.D0
          BARKLEM=.TRUE.
          RETURN
        ENDIF
   4  CONTINUE
C
C Species was not found
C
      BARKLEM=.FALSE.
      RETURN
C
C End of computer-generated subroutine KP_Q_SPLN
      END


C=========================================================================
C=========================================================================
C
C NEGION: Returns partition function and ionization equilibrium for
C         a given negative ion and temperature.
C
C Inputs:
C   ANUM  [integer] atomic number.
C   TEMP  [real]    temperature (in K)
C   PARTN [real]    partition function of neutral atom
C
C                                     (3/2)              Eaffin
C   1     P(A)*P(e)        (2*Pi*m*kT)       2*U(A)    - ----
C  -- =   --------- = kT * -----------    *  ------ * e   kT
C  IT       P(A-)              h^3           U(A-)
C
C  U(A) is passed in as PARTN
C
C                           (3/2)
C  Const = k*(2*Pi*m_e*k/h^2)
C
C History:
C  10-dec-2007: First version written by N. Piskunov including 7 ions.
C               Partition functions tabulated by P. Barklem, resampled
C               for optimal spline interpolation and converted to Fortran
C               DATA statements by J. Valenti
C
C  15-dec-2007: Second version includes the same 7 negative ions tabulated
C               vs alog10(T) on adaptive grid similar to molecular species.
C
C  26-apr-2019: Subroutine data modified and the subroutine text generated
C               by IDL program qk_spl_nodes_f77.pro with errthr=0.000100
C
C Outputs:
C   Q_spln [real*8] partition functions at temperature T,
C          interpolated from Paul Barklem's tables;
C   IT     [real*8] computed according to the formula above.
C
C To obtain partition functions,Q:
C
C   D2 = SPL_INIT(TQ_<species>,Q_<species>)
C   Q(T) = SPL_INTERP(TQ_<species>,Q_<species>,D2,TLOG)
C
C Note that NEGION returns log10(Q)
C
C Reference:
C   Paul Barklem, Remo Collet, 2016, A&A 588, 96.
C
      SUBROUTINE NEGION(ANUM,TEMP,PARTN,IT,Q_atom,POTION,BARKLEM)
C
      IMPLICIT NONE
      INTEGER ANUM
      REAL TEMP,POTION
      REAL*8 PARTN,IT,Q_atom
      LOGICAL BARKLEM
C
C  Local variables
C
      LOGICAL FIRST
      INTEGER MSPEC,NTQ,KLO,KHI,I,II,ISPEC
      PARAMETER(MSPEC=7, NTQ=14)
      INTEGER MTQ(MSPEC)
      REAL*8 TLOG,A,U(14),SPL_INTERP,Const,TkeV,kBoleV
      PARAMETER(Const=0.3333984D0,kBoleV=8.6173175D-5)
C
      REAL*8 TQ(NTQ,MSPEC),Q(NTQ+1,MSPEC),Q2(NTQ,MSPEC)
      REAL*8           TQ_Hm   (NTQ  ),TQ_Cm   (NTQ  ),TQ_Om   (NTQ  ),
     * TQ_Fm   (NTQ  ),TQ_Sim  (NTQ  ),TQ_Sm   (NTQ  ),TQ_Clm  (NTQ  )
      REAL*8            Q_Hm   (NTQ+1), Q_Cm   (NTQ+1), Q_Om   (NTQ+1),
     *  Q_Fm   (NTQ+1), Q_Sim  (NTQ+1), Q_Sm   (NTQ+1), Q_Clm  (NTQ+1)
      EQUIVALENCE (TQ(1,  1),TQ_Hm   ),(TQ(1,  2),TQ_Cm   )
      EQUIVALENCE (TQ(1,  3),TQ_Om   ),(TQ(1,  4),TQ_Fm   )
      EQUIVALENCE (TQ(1,  5),TQ_Sim  ),(TQ(1,  6),TQ_Sm   )
      EQUIVALENCE (TQ(1,  7),TQ_Clm  )
      EQUIVALENCE ( Q(1,  1), Q_Hm   ),( Q(1,  2), Q_Cm   )
      EQUIVALENCE ( Q(1,  3), Q_Om   ),( Q(1,  4), Q_Fm   )
      EQUIVALENCE ( Q(1,  5), Q_Sim  ),( Q(1,  6), Q_Sm   )
      EQUIVALENCE ( Q(1,  7), Q_Clm  )
C
      INTEGER ATLIST(MSPEC)
      SAVE ATLIST,TQ,Q,Q2,FIRST,KHI,KLO
C
C                   H-  C-  O-  F- Si-  S- Cl-
      DATA ATLIST/  1,  6,  8,  9, 14, 16, 17/
C
C Tables of log10(T) and log10(Q)
C
      DATA TQ_Hm/                                                       071215
     1  0.699999789529, 4.000000000000,     12*0.0D+00/                 Hm   
      DATA  Q_Hm/                                                       071215
     1  0.00000000D+00, 0.00000000D+00, 7.54199982D-01,     12*0.0D+00/ Hm   
      DATA TQ_Cm/                                                       071215
     1  0.699999789529, 2.262800296161, 2.956799967199, 3.208699993175, Cm   
     2  3.317899847297, 3.424600185056, 3.571599724768, 3.704599912441, Cm   
     3  3.811100327501, 3.903499698251, 4.000000000000,      3*0.0D+00/ Cm   
      DATA  Q_Cm/                                                       071215
     1  6.02059991D-01, 6.02059991D-01, 6.02060146D-01, 6.02219450D-01, Cm   
     2  6.03193224D-01, 6.07091424D-01, 6.25081501D-01, 6.62530204D-01, Cm   
     3  7.07809633D-01, 7.54553991D-01, 8.06182861D-01, 1.26199996D+00, Cm   
     4       3*0.0D+00/                                                 Cm   
      DATA TQ_Om/                                                       071215
     1  0.699999789529, 1.320900142120, 1.584200022369, 1.862700068321, Om   
     2  1.987999819421, 2.119399869070, 2.306300342441, 2.516800389873, Om   
     3  2.718500247596, 2.916299985204, 3.416600003460, 3.793500051433, Om   
     4  4.000000000000,      1*0.0D+00/                                 Om   
      DATA  Q_Om/                                                       071215
     1  6.02059991D-01, 6.02061119D-01, 6.02344829D-01, 6.08602488D-01, Om   
     2  6.17605419D-01, 6.32336670D-01, 6.59747806D-01, 6.92087919D-01, Om   
     3  7.18403708D-01, 7.37869385D-01, 7.64474877D-01, 7.72298923D-01, Om   
     4  7.74494613D-01, 1.46000004D+00,      1*0.0D+00/                 Om   
      DATA TQ_Fm/                                                       071215
     1  0.699999789529, 4.000000000000,     12*0.0D+00/                 Fm   
      DATA  Q_Fm/                                                       071215
     1  0.00000000D+00, 0.00000000D+00, 3.40109992D+00,     12*0.0D+00/ Fm   
      DATA TQ_Sim/                                                      071215
     1  0.699999789529, 2.199899787843, 2.598800348324, 2.863699735183, Sim  
     2  3.103100156045, 3.213700109251, 3.322099743478, 3.452599856468, Sim  
     3  3.582399981675, 3.691199625608, 3.831499959401, 3.933100099903, Sim  
     4  3.973899990614, 4.000000000000/                                 Sim  
      DATA  Q_Sim/                                                      071215
     1  6.02059991D-01, 6.02059991D-01, 6.02059991D-01, 6.02061225D-01, Sim  
     2  6.02469956D-01, 6.04496833D-01, 6.11581979D-01, 6.35197979D-01, Sim  
     3  6.83836825D-01, 7.44193693D-01, 8.37573130D-01, 9.06940668D-01, Sim  
     4  9.33786269D-01, 9.50480774D-01, 1.38900006D+00/                 Sim  
      DATA TQ_Sm/                                                       071215
     1  0.699999789529, 1.580599931822, 1.959900052931, 2.292500008935, Sm   
     2  2.418900066249, 2.546500103114, 2.909499819934, 3.226499785663, Sm   
     3  3.734799657316, 3.895799929146, 4.000000000000,      3*0.0D+00/ Sm   
      DATA  Q_Sm/                                                       071215
     1  6.02059991D-01, 6.02059994D-01, 6.02165531D-01, 6.08270893D-01, Sm   
     2  6.17116269D-01, 6.31152923D-01, 6.85656158D-01, 7.26190657D-01, Sm   
     3  7.60386090D-01, 7.65723946D-01, 7.68312725D-01, 2.07699990D+00, Sm   
     4       3*0.0D+00/                                                 Sm   
      DATA TQ_Clm/                                                      071215
     1  0.699999789529, 4.000000000000,     12*0.0D+00/                 Clm  
      DATA  Q_Clm/                                                      071215
     1  0.00000000D+00, 0.00000000D+00, 3.61700010D+00,     12*0.0D+00/ Clm  
C
C Length of idividual temperature grids
C
      DATA MTQ/  2, 11, 13,  2, 14, 11,  2/
C
      DATA FIRST/.TRUE./
C
C Compute 2nd derivatives for spline interpolation
C
      IF(FIRST) THEN
        DO 1 I=1,MSPEC
          CALL SPL_INIT(TQ(1,I),Q(1,I),Q2(1,I),U,MTQ(I))
   1    CONTINUE
        FIRST=.FALSE.
      ENDIF
C
C Fits are made in log10 of temperatures
C
      TLOG=LOG10(TEMP)
C
C Find species name
C
      DO 3 II=1,MSPEC
        ISPEC=II
        IF(ANUM.EQ.ATLIST(II)) THEN
C
C  The species is in Barklem's list.
C  Find the braketing temperatures for the partition functions.
C
          KHI=MTQ(ISPEC)
          KLO=1
   2      CONTINUE
          I=(KLO+KHI)/2
          A=TQ(I,ISPEC)
          IF(A.GT.TLOG) THEN
            KHI=I
          ELSE IF(A.LE.TLOG) THEN
            KLO=I
          END IF
          IF(KHI-KLO.GT.1) GO TO 2
C
C Do the interpolation of the partition functions
C
          Q_atom=SPL_INTERP(KLO,KHI,TQ(1,ISPEC),Q(1,ISPEC),Q2(1,ISPEC),
     *                      MTQ(ISPEC),TLOG)
          TkeV=kBoleV*TEMP
          Q_atom=10.d0**Q_atom
          POTION=-Q(MTQ(ISPEC)+1,ISPEC)
          IT=Const*(2.d0*PARTN)/Q_atom*EXP(-POTION/TkeV)*SQRT(TEMP)*
     *       TEMP*TEMP
          IT=1.D0/IT
          BARKLEM=.TRUE.
          RETURN
        ENDIF
   3  CONTINUE
C
C Species was not found
C
      Q_atom=1.D0
      IT=1.D-50
      BARKLEM=.FALSE.
      RETURN
C
C End of computer-generated subroutine NEGION
      END
