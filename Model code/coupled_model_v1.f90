!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Cimatoribus coupled to CC model v1.0
!   Start: 26-10-2022
!   Author: Amber Boot
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION CO2, Ct, Cts, Cn, Cs, Cd, At, Ats, An, As, Ad, Pt, Pts, Pn, Ps, Pd, TC, TA, TP, dCtot
      DOUBLE PRECISION Ht, Hts, Hn, Hs, Hd, CO3t, CO3ts, CO3n, CO3s, CO3d, pCO2t, pCO2ts, pCO2n, pCO2s
      DOUBLE PRECISION K0t, K1t,K2t, K0ts, K1ts, K2ts, K0n, K1n,K2n, K0s, K1s, K2s, K0d, K1d, K2d, BTt, BTts, BTn, BTs, BTd
      DOUBLE PRECISION Kt0, Kts0, Kn0, Ks0, Kd0, K11t, K11ts, K11n, K11s, K11d, Kt1, Kts1, Kn1, Ks1, Kd1
      DOUBLE PRECISION dVCt, dVCts, dVCn, dVCs, dVCd, dKCt, dKCts, dKCn, dKCs, dKCd, Kspt, Kspts, Kspn
      DOUBLE PRECISION Ksps, Kspd, Kspds, Ksprest, Ksprests, Kspresn, Kspress, Kspresd
      DOUBLE PRECISION T0, Tt, Tts, Tn, Ts, Td, S0, St, Sts, Sn, Ss, Sd, rho0, rhot, rhots, rhon, rhos, rhod, alpha, beta
      DOUBLE PRECISION V0, Vt, Vts, Vn, Vs, Vd, VAt, SAt, SAts, SAn, SAs, dt, dts, dn, ds, D
      DOUBLE PRECISION d0, LxA, LxS , Ly, A, Ar, t_dim, S_dim, D_dim,Zt2,Zs2,Zts2,Zn2
      DOUBLE PRECISION Cphys_t, Cphys_ts, Cphys_n, Cphys_s, Aphys_t, Aphys_ts, Aphys_n, Aphys_s, Pphys_t, Pphys_ts, Pphys_n, Pphys_s
      DOUBLE PRECISION Ccarb_t, Ccarb_ts, Ccarb_n, Ccarb_s, Ccarb_d, Acarb_t, Acarb_ts, Acarb_n, Acarb_s
      DOUBLE PRECISION Ccarb_t1, Ccarb_ts1, Ccarb_n1, Ccarb_s1, Ccarb_d1, Ccarb_t2, Ccarb_ts2, Ccarb_n2, Ccarb_s2, Ccarb_d2
      DOUBLE PRECISION Cbio_t, Cbio_ts, Cbio_n, Cbio_s, Pbio_t, Pbio_ts, Pbio_n, Pbio_s
      DOUBLE PRECISION Cair_t, Cair_ts, Cair_n, Cair_s, Criver_t, Ariver_t, Priver_t, EC
      DOUBLE PRECISION qe, qs, qu, qn, qEk, QN1, QS1, QS2, rn, rs, eta, kappa, Agm, tau, Ea, Es, fs
      DOUBLE PRECISION Zt, Zts, Zn, Zs, b, rcp, FCAt, FCAts, FCAn, FCAs, Cat, Cats, Can, Cas, Cad
      DOUBLE PRECISION OmegaC_t, OmegaC_ts, OmegaC_n, OmegaC_s, OmegaC_d, OmegaC_ds,epstol
      DOUBLE PRECISION PVt, PVts, PVn, PVs, A1_C, A2_C, A3_C, B1_C, B2_C, B3_C, n, kCa, PerC, DC, Priver, Csi, Vsi, Vca
      DOUBLE PRECISION RAt, RAts, RAn, RAs, RAd, DMt, DMts, DMn, DMs, DMd, Kd, dtot, Pbart, Pbarts,Pbarn, Pbars
      DOUBLE PRECISION TC1, TA1, TP1, TS1, b1, b2, b3, b4, b5, b6, b7, Kt2, Kts2, Kn2, Ks2, Kd2, Pbard, SAd
      DOUBLE PRECISION CO2_0, Fca_t, Fca_ts, Fca_n, Fca_s, eta_c, Tadjust, PV_t, PV_ts, PV_n, PV_s
      DOUBLE PRECISION epst_base, epst, epsts_base, epsts, epsn_base, epsn, epss_base, epss
      DOUBLE PRECISION Z_t, Z_ts, Z_n, Z_s, PV, Cdil_t, Cdil_n, Cdil_s, Adil_t,Adil_n,Adil_s
           
        ! Switches

        ! Start with defining the state variables
        St = U(1)         ! Salinity thermocline box
        Sts = U(2)         ! Salinity ts box
        Sn = U(3)         ! Salinity northern box
        Ss = U(4)         ! Salinity southern box
        D = U(5)         ! Thermocline depth

        Ct = U(6)         ! DIC thermocline box
        Cts = U(7)         ! DIC ts box
        Cn = U(8)         ! DIC northern box
        Cs = U(9)         ! DIC southern box

        At = U(10)         ! Alkalinity thermocline box
        Ats = U(11)         ! Alkalinity ts box
        An = U(12)         ! Alkalinity northern box
        As = U(13)         ! Alkalinity southern box

        Pt = U(14)         ! Phosphate thermocline box
        Pts = U(15)         ! Phosphate ts box
        Pn = U(16)         ! Phosphate northern box
        Ps = U(17)         ! Phosphate southern box

        CO2 = U(18)         ! CO2 concentration atmosphere
        dCtot = U(19)         ! CO2 concentration atmosphere
        
        t_dim = 100*365*86400
        D_dim = 1000
        S_dim = 35

        ! Parameters physical model (not from Cimatoribus: dn, ds, Tt, Ts)
        V0 = 3d17       ! Total volume of the basin [m3]
        Vn = 3d15       ! Volume northern box [m3]
        Vs = 9d15       ! Volume southern box [m3]
        A = 1d14        ! Horizontal area of Atlantic pycnocline [m2]
        LxA = 1d7       ! Zonal exten of the Atlantic at its southern end [m]
        LxS = 3d7       ! Zonal extent of the Southern Ocean [m]
        Ly = 1d6        ! Meridional extent of frontal region Southern Ocean [m]
        Ar = 4d14       ! Area of thermocline + ts [m2]
        Vts = LxA*Ly*D*D_dim/2          ! Volume ts box [m3]
        Vt = A*D*D_dim           ! Volume thermocline box

        dt = D*D_dim  ! Depth thermocline box
        dts = D*D_dim/2   ! Depth ts box
        dn = 300   ! Depth northern box [m] (first guess for now)
        ds = 300   ! Depth southern box [m] (first guess for now)
        SAn = Vn/dn      ! Surface area northern box [m2]
        SAs = Vs/ds      ! Surfacer area southern box [m2]
        SAt = Vt/dt      ! Surdace area thermocline box [m2]
        SAts = LxA*Ly   ! Surface area ts box [m2]
        SAd = SAn + SAs + SAt + SAts

        Vd = V0-Vn-Vs-Vts-Vt            ! Volume deep box [m3]
        dtot = 4000     ! Total depth [m]
        EC = par(12)*(1d15/12)	! Extra carbon [PgC]
   		
   		CO2_0 = 300d-6 ! 275d-6
   		eta_c = 1
   		
   		S0 = 35         ! Reference salinity [psu]
        tau = 0.1       ! Average zonal wind stress amplitude [Nm-2]
        Agm = 1700      ! Eddy diffusivity [m2s-1]
        kappa = 1d-5    ! Vertical diffusivity [m2s-1]
        eta = 3d4       ! Hydraulic constant [ms-1]
        fs = -1d-4      ! Coriolis paramter [s-1]
        !Es = (1-par(5))*0.56d6 + par(5)*(0.00012*CO2*1d6+0.562)*1d6   ! Symmetric freshwater flux [m3s-1]
        Es = (1-par(5))*0.56d6 + par(5)*(0.099 +0.079 *log(CO2*1d6))*1d6   ! Symmetric freshwater flux [m3s-1]
        !Es = (1-par(5))*0.56d6 + par(5)*((log((CO2+250*1d-6)/CO2_0)*0.75d6+0.1d6))
        rn = 5d6        ! Transport by northern subtropical gyre [m3s-1]
        rs = 1d7        ! Transport by southern subtropical gyre
        Ea = par(1) + (par(15)*0.1*(1+log(CO2/CO2_0))) !* (43*CO2*1d6-64800))/1d6     ! Asymmetric freshater flux [m3s-1]

		Tadjust = par(7)*0.53*5.45*log(CO2/CO2_0)
        rho0 = 1027.5   ! Reference density [kgm-3]
        T0 = 5          ! Reference temperature [degree C]
        Tn = 5 + Tadjust         ! Temperature northern box [degree C]
        Tts = 10 + Tadjust        ! Temperature ts box [degree C]
        alpha = 2d-4    ! Thermal expansion coefficient [K-1]
        beta = 8d-4     ! Haline contraction coefficient [psu-1]

        Tt = 23.44 + Tadjust      ! Temperature thermocline box [degree C] (SCP-M, box 1)
        Ts = 0.93 + Tadjust   ! Temperature southern box [degree C] (SCP-M, box 5)
        Td =  1.8  ! Temperature deep box [degree C] (SCP-M, box 6)

        ! Parameters for the carbon cycle model

        Z_t =  2.1/(365*86400) ! Biological production thermocline box [mol C m-2 s-1]
        Z_ts = 2.1/(365*86400)  ! Biological production ts box [mol C m-2 s-1]
        Z_n = 1.9/(365*86400)  ! Biological production northern box [mol C m-2 s-1]
        Z_s = 1.1/(365*86400)  ! Biological production southern box [mol C m-2 s-1]

		epst_base=0.50
		epsts_base=0.30
		epsn_base=0.10
		epss_base=0.10
		
		epst=par(10)*Tadjust*(-0.1)+epst_base
		epsts=par(10)*Tadjust*(-0.1)+epsts_base
        epsn=par(10)*Tadjust*(-0.1)+epsn_base
        epss=par(10)*Tadjust*(-0.1)+epss_base

        d0 = 100    ! Base depth for export production [m]
        b = 0.75 ! Power law coefficient for export production [-]
        rcp=130 ! Redfield P:C ratio [mol P per mol C]

        Fca_t = 0.15  ! Rain ratio thermocline box [-]
        Fca_ts = 0.15  ! Rain ratio ts box [-]
        Fca_n = 0.15  ! Rain ratio norhtern box [-]
        Fca_s = 0.15  ! Rain ratio southern box [-]
        
        PV_t = 3*SAt/86400     ! Piston velocity thermocline box [-]
        PV_ts = 3*SAts/86400    ! Piston velocity ts box [-]
        PV_n = 3*SAn/86400     ! Piston velocity norhtern box [-]
        PV_s = 3*SAs/86400     ! Piston velocity southern box [-]

		PVt=(PV_t*par(8)*(((2116.8-136.25*Tt+4.7353*Tt**2-0.092307*Tt**3+0.0007555*Tt**4)/660)**-0.5)+PV_t*(1-par(8)))
		PVts=(PV_ts*par(8)*(((2116.8-136.25*Tts+4.7353*Tts**2-0.092307*Tts**3+0.0007555*Tts**4)/660)**-0.5)+PV_ts*(1-par(8)))
		PVn=(PV_n*par(8)*(((2116.8-136.25*Tn+4.7353*Tn**2-0.092307*Tn**3+0.0007555*Tn**4)/660)**-0.5)+PV_n*(1-par(8)))
		PVs=(PV_s*par(8)*(((2116.8-136.25*Ts+4.7353*Ts**2-0.092307*Ts**3+0.0007555*Ts**4)/660)**-0.5)+PV_s*(1-par(8)))
		
        kCa = 0.38/86400    ! Constant CaCO3 dissolution rate [s-1]
        PerC = 0.12 ! Percentage of C in CaCO3 [-]
        n = 1   ! Order of CaCO3 dissolution kinetics [-]
        DC = 2.75d-13   ! Constant CaCO3 dissolution rate [mol m-3 s-1]

        VAt = 1.76d20          ! Volume of the atmosphere [m3]

        Priver = 15356.44572697869*0.2       ! River influx of PO4

        Pbart = dt/2/10
        Pbarts = dts/2/10
        Pbarn = dn/2/10
        Pbars = ds/2/10
        Pbard = (dtot-(Vd/SAd))/2/10
        
        Csi =2.378234398782344*(1d-12) 
        Vsi = 1.585489599188229*(1d-8)  
        Vca = 6.341958396752917*(1d-8)

        ! For K0x
        A1_C = -60.3409
        A2_C = 93.4517
        A3_C = 23.3585
        B1_C = 0.023517
        B2_C = -0.023656
        B3_C = 0.0047036

        ! Determine ocean circulation
        qEk = tau*LxS/(rho0*ABS(fs))
        qe = Agm*D*D_dim*LxA/Ly
        qs = qEk - qe
        qu = kappa*A/(D*D_dim)
        qn = eta*(D*D_dim*D*D_dim)*(alpha*(Tts-Tn)+beta*(Sn*S_dim-Sts*S_dim))

        ! Conservation laws for S, C, A, and P
        Sd= ((S0*V0-(St*Vt+Sts*Vts+Sn*Vn+Ss*Vs)*S_dim)/Vd)/S_dim !-(b1 + b2*Sn + b5*D*St + b6*D*Sts + b7*Ss)/(b3+b4*D)

        TC = CO2*VAt+Ct*Vt+Cts*Vts+Cn*Vn+Cs*Vs   ! Sum of carbon without eliminated variable
        TC1 = dCtot*(1d22)*par(2)+EC                      ! Total carbon in system
        Cd = 1/Vd*(TC1-TC)                               ! Eliminated variable

        TA = At*Vt+Ats*Vts+An*Vn+As*Vs! Sum of alkalinity without eliminated variable
        TA1 = ((697249.347221776*(1d12)-(685539.5393413941*(1d12)-dCtot*(1d22))*2)*par(2))           ! Total alkalinity in system
        Ad = (1/Vd)*(TA1-TA)! Eliminated variable

        TP = Pt*Vt+Pts*Vts+Pn*Vn+Ps*Vs      ! Sum of phosphate without eliminated variable
        TP1 = 2.8733503182616036*(1d15)*0.2*par(2)        ! Total phosphate in system
        Pd = (1/Vd)*(TP1-TP)                          ! Eliminated variable

        ! Carbonate chemistry
        ! Thermocline box
        BTt = 1.179e-5*St*S_dim ! Total Boron
        Kt0 = A1_C+A2_C*(100.0/((Tt+273.16)))+A3_C*log((Tt+273.16)/100.0)
        K0t = exp((Kt0+(St*S_dim)*(B1_C+B2_C*((Tt+273.16)/100.0)+B3_C*(((Tt+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K1t = 10**-((3633.86/(Tt+273.16))-61.2172+(9.67770*log(Tt+273.16))-(0.011555*(St*S_dim))+(0.0001152*((St*S_dim)**2))) !Lueker et al (2000) in mol kg-1
        K2t = 10**-((471.78/((Tt+273.16)))+25.9290-(3.16967*log((Tt+273.16)))-(0.01781*(St*S_dim))+(0.0001122*((St*S_dim)**2)))! Lueker et al (2000) in mol kg-1
        
        K11t = -395.8293+(6537.773/(Tt+273.16))+71.595*log((Tt+273.16))-0.17959*(Tt+273.16)
        Kt2 = (-1.78938+410.64/(Tt+273.16)+0.0065453*(Tt+273.16))*(St*S_dim)**0.5-0.17755*(St*S_dim)+0.0094979*(St*S_dim)**(3.0/2.0) ! Mucci (1983)
        Kt1 = exp(K11t+Kt2)
        dVCt = -65.28+0.397*Tt-0.005155*Tt**2+(19.816-0.0441*Tt-0.00017*Tt**2)*(St)**0.5  ! St no S_dim
        dKCt = 0.01847+0.0001956*Tt-0.000002212*Tt**2+(-0.03217-0.0000711*Tt-0.000002212*Tt**2)*(St)**0.5 ! St no S_dim
        Ksprest = (-dVCt+0.5*dKCt*Pbart)*Pbart/(83.144621*(Tt+273.16))
        Kspt = Kt1*exp(Ksprest)

        ! ts box
        BTts = 1.179e-5*Sts*S_dim ! Total Boron
        Kts0 = A1_C+A2_C*(100.0/((Tts+273.16)))+A3_C*log((Tts+273.16)/100.0)
        K0ts = exp((Kts0+(Sts*S_dim)*(B1_C+B2_C*((Tts+273.16)/100.0)+B3_C*(((Tts+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K1ts = 10**-((3633.86/(Tts+273.16))-61.2172+(9.67770*log(Tts+273.16))-(0.011555*(Sts*S_dim))+(0.0001152*((Sts*S_dim)**2))) !Lueker et al (2000) in mol kg-1
        K2ts = 10**-((471.78/((Tts+273.16)))+25.9290-(3.16967*log((Tts+273.16)))-(0.01781*(Sts*S_dim))+(0.0001122*((Sts*S_dim)**2)))! Lueker et al (2000) in mol kg-1
        
        K11ts = -395.8293+(6537.773/(Tts+273.16))+71.595*log((Tts+273.16))-0.17959*(Tts+273.16)
    Kts2=(-1.78938+410.64/(Tts+273.16)+0.0065453*(Tts+273.16))*(Sts*S_dim)**0.5-0.17755*(Sts*S_dim)+0.0094979*(Sts*S_dim)**(3.0/2.0) ! Mucci (1983)
        Kts1 = exp(K11ts+Kts2)
        dVCts = -65.28+0.397*Tts-0.005155*Tts**2+(19.816-0.0441*Tts-0.00017*Tts**2)*(Sts)**0.5 ! Sts no S_dim
        dKCts = 0.01847+0.0001956*Tts-0.000002212*Tts**2+(-0.03217-0.0000711*Tts-0.000002212*Tts**2)*(Sts)**0.5 ! Sts no S_dim
        Ksprests = (-dVCts+0.5*dKCts*Pbarts)*Pbarts/(83.144621*(Tts+273.16))
        Kspts = Kts1*exp(Ksprests)
		
        ! Northern box
        BTn = 1.179e-5*Sn*S_dim ! Total Boron
        Kn0 = A1_C+A2_C*(100.0/((Tn+273.16)))+A3_C*log((Tn+273.16)/100.0)
        K0n = exp((Kn0+(Sn*S_dim)*(B1_C+B2_C*((Tn+273.16)/100.0)+B3_C*(((Tn+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K1n = 10**-((3633.86/(Tn+273.16))-61.2172+(9.67770*log(Tn+273.16))-(0.011555*(Sn*S_dim))+(0.0001152*((Sn*S_dim)**2))) !Lueker et al (2000) in mol kg-1
        K2n = 10**-((471.78/((Tn+273.16)))+25.9290-(3.16967*log((Tn+273.16)))-(0.01781*(Sn*S_dim))+(0.0001122*((Sn*S_dim)**2)))! Lueker et al (2000) in mol kg-1
        
        K11n = -395.8293+(6537.773/(Tn+273.16))+71.595*log((Tn+273.16))-0.17959*(Tn+273.16)
        Kn2 = (-1.78938+410.64/(Tn+273.16)+0.0065453*(Tn+273.16))*(Sn*S_dim)**0.5-0.17755*(Sn*S_dim)+0.0094979*(Sn*S_dim)**(3.0/2.0) ! Mucci (1983)
        Kn1 = exp(K11n+Kn2)
        dVCn = -65.28+0.397*Tn-0.005155*Tn**2+(19.816-0.0441*Tn-0.00017*Tn**2)*(Sn)**0.5 ! Sn no S_dim
        dKCn = 0.01847+0.0001956*Tn-0.000002212*Tn**2+(-0.03217-0.0000711*Tn-0.000002212*Tn**2)*(Sn)**0.5 ! Sn no S_di,
        Kspresn = (-dVCn+0.5*dKCn*Pbarn)*Pbarn/(83.144621*(Tn+273.16))
        Kspn = Kn1*exp(Kspresn)

        ! Southern box
        BTs = 1.179e-5*Ss*S_dim ! Total Boron
        Ks0 = A1_C+A2_C*(100.0/((Ts+273.16)))+A3_C*log((Ts+273.16)/100.0)
        K0s = exp((Ks0+(Ss*S_dim)*(B1_C+B2_C*((Ts+273.16)/100.0)+B3_C*(((Ts+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K1s = 10**-((3633.86/(Ts+273.16))-61.2172+(9.67770*log(Ts+273.16))-(0.011555*(Ss*S_dim))+(0.0001152*((Ss*S_dim)**2))) !Lueker et al (2000) in mol kg-1
        K2s = 10**-((471.78/((Ts+273.16)))+25.9290-(3.16967*log((Ts+273.16)))-(0.01781*(Ss*S_dim))+(0.0001122*((Ss*S_dim)**2)))! Lueker et al (2000) in mol kg-1
        
        K11s=-395.8293+(6537.773/(Ts+273.16))+71.595*log((Ts+273.16))-0.17959*(Ts+273.16)
        Ks2=(-1.78938+410.64/(Ts+273.16)+0.0065453*(Ts+273.16))*(Ss*S_dim)**0.5-0.17755*(Ss*S_dim)+0.0094979*(Ss*S_dim)**(3.0/2.0) ! Mucci (1983)
        Ks1 = exp(K11s+Ks2)
        dVCs=-65.28+0.397*Ts-0.005155*Ts**2+(19.816-0.0441*Ts-0.00017*Ts**2)*(Ss)**0.5 ! Ss no S_dim
        dKCs=0.01847+0.0001956*Ts-0.000002212*Ts**2+(-0.03217-0.0000711*Ts-0.000002212*Ts**2)*(Ss)**0.5 ! Ss no S_dim
        Kspress=(-dVCs+0.5*dKCs*Pbars)*Pbars/(83.144621*(Ts+273.16))
        Ksps=Ks1*exp(Kspress)

        ! Deep box
        BTd = 1.179e-5*Sd*S_dim ! Total Boron
        Kd0 = A1_C+A2_C*(100.0/((Td+273.16)))+A3_C*log((Td+273.16)/100.0)
        K0d = exp((Kd0+(Sd*S_dim)*(B1_C+B2_C*((Td+273.16)/100.0)+B3_C*(((Td+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K1d = 10**-((3633.86/(Td+273.16))-61.2172+(9.67770*log(Td+273.16))-(0.011555*(Sd*S_dim))+(0.0001152*((Sd*S_dim)**2))) !Lueker et al (2000) in mol kg-1
        K2d = 10**-((471.78/((Td+273.16)))+25.9290-(3.16967*log((Td+273.16)))-(0.01781*(Sd*S_dim))+(0.0001122*((Sd*S_dim)**2)))! Lueker et al (2000) in mol kg-1

        K11d=-395.8293+(6537.773/(Td+273.16))+71.595*log((Td+273.16))-0.17959*(Td+273.16)
        Kd2=(-1.78938+410.64/(Td+273.16)+0.0065453*(Td+273.16))*(Sd*S_dim)**0.5-0.17755*(Sd*S_dim)+0.0094979*(Sd*S_dim)**(3.0/2.0) ! Mucci (1983)
        Kd1=exp(K11d+Kd2)
        dVCd=-65.28+0.397*Td-0.005155*Td**2+(19.816-0.0441*Td-0.00017*Td**2)*(Sd)**0.5 ! Sd no S_dim
        dKCd=0.01847+0.0001956*Td-0.000002212*Td**2+(-0.03217-0.0000711*Td-0.000002212*Td**2)*(Sd)**0.5 ! Sd no S_dim
        Kspresd=(-dVCd+0.5*dKCd*Pbard)*Pbard/(83.144621*(Td+273.16))
        Kspd=Kd1*exp(Kspresd)
        Kspds=Kd1*exp((-dVCd+0.5*dKCd*dtot/10)*(dtot/10)/(83.144621*(Td+273.16)))
		
		epstol = 1.0d-06
        IF (abs(At).gt.epstol) then
           RAt = (Ct/rho0)/(At/rho0)
        ELSE 
           RAt = 0.0
        ENDIF 
        IF (abs(Ats).gt.epstol) then
           RAts = (Cts/rho0)/(Ats/rho0)
        ELSE 
           RAts = 0.0
        ENDIF 
        IF (abs(An).gt.epstol) then
           RAn = (Cn/rho0)/(An/rho0)
        ELSE 
           RAn = 0.0
        ENDIF 
        IF (abs(As).gt.epstol) then
           RAs = (Cs/rho0)/(As/rho0)
        ELSE 
           RAs = 0.0
        ENDIF 
        IF (abs(Ad).gt.epstol) then
           RAd = (Cd/rho0)/(Ad/rho0)
        ELSE 
           RAd = 0.0
        ENDIF 
 
        !RAt = (Ct/rho0)/(At/rho0)       ! C to A ratio thermocline box
        !RAts = (Cts/rho0)/(Ats/rho0)    ! C to A ratio ts box
        !RAn = (Cn/rho0)/(An/rho0)       ! C to A ratio northern box
        !RAs = (Cs/rho0)/(As/rho0)       ! C to A ratio southern box
        !RAd = (Cd/rho0)/(Ad/rho0)       ! C to A ratio deep box

        DMt = (1-RAt)**2*K1t**2-4*K1t*K2t*(1-2*RAt)
        DMts = (1-RAts)**2*K1ts**2-4*K1ts*K2ts*(1-2*RAts)
        DMn = (1-RAn)**2*K1n**2-4*K1n*K2n*(1-2*RAn)
        DMs = (1-RAs)**2*K1s**2-4*K1s*K2s*(1-2*RAs)
        DMd =(1-RAd)**2*K1d**2-4*K1d*K2d*(1-2*RAd)

        Ht = 0.5*((RAt-1)*K1t+DMt**0.5)         ! [H+] concentration thermocline box
        Hts = 0.5*((RAts-1)*K1ts+DMts**0.5)     ! [H+] concentration ts box
        Hn = 0.5*((RAn-1)*K1n+DMn**0.5)         ! [H+] concentration northern box
        Hs = 0.5*((RAs-1)*K1s+DMs**0.5)         ! [H+] concentration southern box
        Hd = 0.5*((RAd-1)*K1d+DMd**0.5)         ! [H+] concentration deep box

        CO3t = (Ct/rho0/K0t)*(Ht**2/(Ht**2+Ht*K1t+K1t*K2t))*K0t*K1t/Ht*K2t/Ht                   ! [CO32-] concentration thermocline box
        CO3ts = (Cts/rho0/K0ts)*(Hts**2/(Hts**2+Hts*K1ts+K1ts*K2ts))*K0ts*K1ts/Hts*K2ts/Hts     ! [CO32-] concentration ts box
        CO3n = (Cn/rho0/K0n)*(Hn**2/(Hn**2+Hn*K1n+K1n*K2n))*K0n*K1n/Hn*K2n/Hn                   ! [CO32-] concentration northern box
        CO3s = (Cs/rho0/K0s)*(Hs**2/(Hs**2+Hs*K1s+K1s*K2s))*K0s*K1s/Hs*K2s/Hs                   ! [CO32-] concentration southern box
        CO3d = (Cd/rho0/K0d)*(Hd**2/(Hd**2+Hd*K1d+K1d*K2d))*K0d*K1d/Hd*K2d/Hd                   ! [CO32-] concentration deep box

        pCO2t = par(3)*(Ct/rho0/K0t)*((Ht**2)/((Ht**2)+Ht*K1t+K1t*K2t))+(1-par(3))*Ct            ! pCO2 thermocline box
        pCO2ts = par(3)*(Cts/rho0/K0ts)*((Hts**2)/((Hts**2)+Hts*K1ts+K1ts*K2ts))+(1-par(3))*Cts    ! pCO2 ts box
        pCO2n = par(3)*(Cn/rho0/K0n)*((Hn**2)/((Hn**2)+Hn*K1n+K1n*K2n))+(1-par(3))*Cn            ! pCO2 northern box
        pCO2s = par(3)*(Cs/rho0/K0s)*((Hs**2)/((Hs**2)+Hs*K1s+K1s*K2s))+(1-par(3))*Cs            ! pCO2 southern box

        Cat = 0.01028*St         ! [Ca2+] in thermocline box [mol m-3] (no S_dim)
        Cats = 0.01028*Sts      ! [Ca2+] in ts box [mol m-3]
        Can = 0.01028*Sn        ! [Ca2+] in northern box [mol m-3]
        Cas = 0.01028*Ss        ! [Ca2+] in southern box [mol m-3]
        Cad = 0.01028*Sd        ! [Ca2+] in deep box [mol m-3]

        OmegaC_t = min((CO3t*Cat/Kspt),n)   ! Saturation state in thermocline box [-]
        OmegaC_ts = min((CO3ts*Cats/Kspts),n) ! Saturation state in ts box [-]
        OmegaC_n = min((CO3n*Can/Kspn),n) ! Saturation state in northern box [-]
        OmegaC_s = min((CO3s*Cas/Ksps),n) ! Saturation state in southern box [-]
        OmegaC_d = min((CO3d*Cad/Kspd),n) ! Saturation state in deep box [-]
        OmegaC_ds = min((CO3d*Cad/Kspds),n) ! Saturation state in sediment deep box [-]

		IF (Cat*CO3t/Kspt-1>0) then
           Fcat=par(9)*0.022*(Cat*CO3ts/Kspt-1)**eta_c+(1-par(9))*Fca_t
        ELSE 
           Fcat = Fca_t
        ENDIF 
        
        IF (Cats*CO3ts/Kspts-1>0) then
           Fcats=par(9)*0.022*(Cats*CO3t/Kspts-1)**eta_c+(1-par(9))*Fca_ts
        ELSE 
           Fcats = Fca_ts
        ENDIF 
        
        IF (Can*CO3n/Kspn-1>0) then
           Fcan=par(9)*0.022*(Can*CO3n/Kspn-1)**eta_c+(1-par(9))*Fca_n
        ELSE 
           Fcan = Fca_n
        ENDIF 
        
        IF (Cas*CO3s/Ksps-1>0) then
           Fcas=par(9)*0.022*(Cas*CO3s/Ksps-1)**eta_c+(1-par(9))*Fca_s
        ELSE 
           Fcas = Fca_s
        ENDIF 
        
		! Heaviside functions
        IF (qn.lt.0) THEN
            QN1 = 0
        ELSE
            QN1 = qn
        ENDIF

        IF (qs.gt.0) THEN
            QS2 = 0
            QS1 = qs
        ELSE
            QS2 = qs
            QS1 = 0
        ENDIF
        
        Zt=((1-par(6))*Z_t+par(6)*(QS1*Pts+rs*Pts+qu*Pd+Priver+rn*Pn)*epst*rcp/SAt)
        Zts=((1-par(6))*Z_ts+par(6)*((-QS2*Pt+qEk*Ps+rs*Pt)*epsts*rcp/SAts))
        Zn=((1-par(6))*Z_n+par(6)*((QN1*Pt+rn*Pt)*epsn*rcp/SAn))
        Zs=((1-par(6))*Z_s+par(6)*(QS1*Pd+qe*Pts)*epss*rcp/SAs)
            
        ! Determine fluxes carbon cycle model
        Cphys_t = (-(-QS2+QN1+rs+rn)*Ct+(QS1+rs)*Cts+qu*Cd+rn*Cn)/Vt  ! Advective DIC flux thermocline box
        Cphys_ts = (-(QS1+qe+rs)*Cts+(-QS2+rs)*Ct+qEk*Cs)/Vts  ! Advective DIC flux ts box
        Cphys_n = (-(QN1+rn)*Cn+(QN1+rn)*Ct)/Vn  ! Advective DIC flux northern box
        Cphys_s = (-(qEk-QS2)*Cs+QS1*Cd+qe*Cts)/Vs  ! Advective DIC flux southern box

        Cbio_t = -(Zt*SAt/(Vt))*((D*D_dim)/d0)**(-b) ! Biological DIC flux thermocline box
        Cbio_ts = -(Zts*SAts/(Vts))*((D*D_dim*0.5)/d0)**(-b)!par(4)*((-Zts*SAts/Vts)*(0.5*D*D_dim/d0)**(-b))  ! Biological DIC flux ts box
        Cbio_n = -(Zn*SAn/(Vn))*((dn)/d0)**(-b)!par(4)*((-Zn*SAn/Vn)*(dn/d0)**(-b)) ! Biological DIC flux northern box
        Cbio_s = -(Zs*SAs/(Vs))*((ds)/d0)**(-b)! 0*par(4)*((-Zs*SAs/Vs)*(ds/d0)**(-b))   ! Biological DIC flux southern box

        Ccarb_t1 = -((Zt*SAt)*(FCAt))/Vt   ! Carbonate DIC flux related to biology thermocline box
        Ccarb_ts1 = -((Zts*SAts)*(FCAts))/Vts ! Carbonate DIC flux related to biology ts box
        Ccarb_n1 = -((Zn*SAn)*(FCAn))/Vn      ! Carbonate DIC flux related to biology northern box
        Ccarb_s1 = -((Zs*SAs)*(FCAs))/Vs      ! Carbonate DIC flux related to biology southern box
      
        Ccarb_t2 = (CO3t*Cat)*(rho0)*kCa*(1-OmegaC_t)**n*PerC+DC      ! Carbonate DIC flux related to dissolution thermocline box
        Ccarb_ts2 = (CO3ts*Cats)*(rho0)*kCa*(1-OmegaC_ts)**n*PerC+DC  ! Carbonate DIC flux related to dissolution ts box
        Ccarb_n2 = (CO3n*Can)*(rho0)*kCa*(1-OmegaC_n)**n*PerC+DC      ! Carbonate DIC flux related to dissolution northern box
        Ccarb_s2 = (CO3s*Cas)*(rho0)*kCa*(1-OmegaC_s)**n*PerC+DC      ! Carbonate DIC flux related to dissolution southern box
        
        Ccarb_d1 = (CO3d*Cad)*(rho0)*kCa*(1-OmegaC_d)**n*PerC+DC      ! Carbonate DIC flux related to water dissolution deep box
        Ccarb_d2 = (CO3d*Cad)*(rho0)*kCa*(1-OmegaC_ds)**n*PerC+DC      ! Carbonate DIC flux related to sediment dissolution deep box
           !1.2845340708965197e-12!     
        Ccarb_t = (Ccarb_t1+Ccarb_t2)  ! Carbonate DIC flux thermocline box
        Ccarb_ts =(Ccarb_ts1+Ccarb_ts2) ! Carbonate DIC flux ts box
        Ccarb_n = (Ccarb_n1+Ccarb_n2) ! Carbonate DIC flux northern box
        Ccarb_s = (Ccarb_s1+Ccarb_s2) ! Carbonate DIC flux southern box
        Ccarb_d = (Ccarb_d1+Ccarb_d2) ! Carbonate DIC flux deep box

        Cair_t = K0t*PVt*rho0*(CO2-pCO2t)/Vt ! Air-sea CO2 flux thermocline box
        Cair_ts = K0ts*PVts*rho0*(CO2-pCO2ts)/Vts ! Air-sea CO2 flux ts box
        Cair_n = K0n*PVn*rho0*(CO2-pCO2n)/Vn ! Air-sea CO2 flux northern box
        Cair_s = K0s*PVs*rho0*(CO2-pCO2s)/Vs ! Air-sea CO2 flux southern box

        Aphys_t = (-(-QS2+QN1+rs+rn)*At+(QS1+rs)*Ats+qu*Ad+rn*An)/Vt  ! Advective alkalinity flux thermocline box
        Aphys_ts = (-(QS1+qe+rs)*Ats+(-QS2+rs)*At+qEk*As)/Vts  ! Advective alkalinity flux ts box
        Aphys_n = (-(QN1+rn)*An+(QN1+rn)*At)/Vn  ! Advective alkalinity flux northern box
        Aphys_s = (-(qEk-QS2)*As+QS1*Ad+qe*Ats)/Vs  ! Advective alkalinity flux southern box

        Acarb_t = 2*Ccarb_t                    ! Carbonate alkalinity flux thermocline box
        Acarb_ts = 2*Ccarb_ts                  ! Carbonate alkalinity flux ts box
        Acarb_n = 2*Ccarb_n                    ! Carbonate alkalinity flux northern box
        Acarb_s = 2*Ccarb_s                    ! Carbonate alkalinity flux southern box

        Pphys_t = (-(-QS2+QN1+rs+rn)*Pt+(QS1+rs)*Pts+qu*Pd+rn*Pn)/Vt   ! Advective phosphate flux thermocline box
        Pphys_ts = (-(QS1+qe+rs)*Pts+(-QS2+rs)*Pt+qEk*Ps)/Vts! Advective phosphate flux ts box
        Pphys_n = (-(QN1+rn)*Pn+(QN1+rn)*Pt)/Vn  ! Advective phosphate flux northern box
        Pphys_s = (-(qEk-QS2)*Ps+QS1*Pd+qe*Pts)/Vs  ! Advective phosphate flux southern box

        Pbio_t = Cbio_t/rcp                     ! Biological phosphate flux thermocline box
        Pbio_ts = Cbio_ts/rcp                   ! Biological phosphate flux ts box
        Pbio_n = Cbio_n/rcp                      ! Biological phosphate flux northern box
        Pbio_s = Cbio_s/rcp                    ! Biological phosphate flux southern box

        Criver_t = ((Csi + (Vsi+Vca)*CO2)*0.2)          ! DIC river flux
        Ariver_t = (2*Criver_t)                   ! Alkalinity river flux
        Priver_t = (Priver/Vt)                    ! Phosphate river flux

		Cdil_t = par(13)*(2*Es*Ct/Vt)
		Adil_t = par(13)*(2*Es*At/Vt)
		
		Cdil_n = par(13)*((-Es-Ea)*Cn/Vn)
		Adil_n = par(13)*((-Es-Ea)*An/Vn)
		
		Cdil_s = par(13)*((-Es+Ea)*Cs/Vs)
		Adil_s = par(13)*((-Es+Ea)*As/Vs)
		
        ! State the ODEs
        F(1) = ((QS1*Sts+QS2*St)+qu*Sd-QN1*St+rs*(Sts-St)+rn*(Sn-St)+2*Es)*t_dim/Vt-(St*(qu+qs-QN1)*t_dim/(Ar*D_dim))/D                                     ! Salinity thermocline box
        F(2) = (qEk*Ss-qe*Sts-(QS1*Sts+QS2*St)+rs*(St-Sts))*t_dim/Vts-(Sts*(qu+qs-QN1)*t_dim/(Ar*D_dim))/D                                            ! Salinity ts box
        F(3) = ((QN1+rn)*(St-Sn)-(Es+Ea*1e6))*t_dim/Vn            ! Salinity northern box
        F(4) = ((QS1*Sd+QS2*Ss)+qe*Sts-qEk*Ss-(Es-Ea*1e6))*t_dim/Vs   ! Salinity southern box
        F(5) = (qu+qs-QN1)*t_dim/(Ar*D_dim)                  ! Thermocline depth

        F(6) = (Cphys_t+Cair_t+par(4)*(Cbio_t+Ccarb_t+Criver_t+Cdil_t))*t_dim       ! DIC thermocline box
        F(7) = (Cphys_ts+Cair_ts+par(4)*(Cbio_ts+Ccarb_ts))*t_dim            ! DIC ts box
        F(8) = (Cphys_n+Cair_n+par(4)*(Cbio_n+Ccarb_n+Cdil_n))*t_dim                ! DIC northern box
        F(9) = (Cphys_s+Cair_s+par(4)*(Cbio_s+Ccarb_s+Cdil_s))*t_dim                ! DIC southern box

        F(10) = (Aphys_t+par(4)*(Acarb_t+Ariver_t+Adil_t))*t_dim                    ! Alkalinity thermocline box
        F(11) = (Aphys_ts+par(4)*(Acarb_ts))*t_dim                           ! Alkalinity ts box
        F(12) = (Aphys_n+par(4)*(Acarb_n+Adil_n))*t_dim                             ! Alkalinity northern box
        F(13) = (Aphys_s+par(4)*(Acarb_s+Adil_s))*t_dim                      ! Alkalinity southern box

        F(14) = (Pphys_t+par(4)*(Pbio_t+Priver_t))*t_dim                     ! Phosphate thermocline box
        F(15) = (Pphys_ts+par(4)*(Pbio_ts))*t_dim                            ! Phosphate ts box
        F(16) = (Pphys_n+par(4)*(Pbio_n))*t_dim                             ! Phosphate northern box
        F(17) = (Pphys_s+par(4)*(Pbio_s))*t_dim                              ! Phosphate southern box

        F(18) = (-(Cair_t*Vt+Cair_ts*Vts+Cair_n*Vn+Cair_s*Vs)/VAt)*t_dim         ! Atmospheric CO2 concentration

        F(19)=((Ccarb_t)*Vt+(Ccarb_ts)*Vts+(Ccarb_n)*Vn+(Ccarb_s)*Vs+(Ccarb_d)*Vd+Criver_t*Vt)/(1d22)*t_dim                              ! Total carbon
        
      END SUBROUTINE FUNC
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      INTEGER I 
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      

! Parameters
        PAR(1) = -0.75 ! EA [Sv]
        PAR(2) = 1.0 ! Homotopy parameter
        PAR(3) = 1.0 ! Homotopy parameter
        PAR(4) = 1.0 ! Homotopy parameter
        PAR(5) = 1.0 ! ES
        PAR(6) = 1.0 ! BIO
        PAR(7) = 1.0 ! TEMP
        PAR(8) = 0.0 ! PV
        PAR(9) = 0.0 ! FCA
        PAR(10) = 0.0 ! ET  
        PAR(12) = 000.0 ! TC  
        PAR(13) = 0.0 ! Virtual fluxes
        PAR(15) = 0.0 ! EA - feedback  
 
! Special solution from transient code 
        U(1) = 1.0058376418610575    ! Salinity thermocline
        U(2) = 0.994840665484994    ! Salinity ts
        U(3) = 0.9984467757381668     ! Salinity northern
        U(4) = 0.9910742010761921     ! Salinity southern
        U(5) = 0.7382511281434748     ! Thermocline depth

        U(6) = 1.9105267705      ! C thermocline
        U(7) = 1.9879092113    ! C ts
        U(8) = 1.9496357353   ! C northern
        U(9) = 2.0426263730     ! C southern

        U(10) = 2.1653666478    ! A thermocline
        U(11) = 2.1576922880    ! A ts
        U(12) = 2.1567531722   ! A northern
        U(13) =  2.1630908657      ! A southern

        U(14) =  9.6051503852E-004   ! P thermocline
        U(15) =  1.0414013836E-003   ! P ts
        U(16) =  7.5116742147E-004   ! P northern
    	U(17) =  1.5651325565E-003   ! P southern

        U(18) = 2.5804155211E-004 ! CO2
        U(19) = 6.6331354527E-005    ! Total C
        
        !U(1) = 9.96811E-01    ! Salinity thermocline
        !U(2) = 1.00423E+00    ! Salinity ts
        !U(3) = 9.12811E-01     ! Salinity northern
        !U(4) = 1.00692E+00     ! Salinity southern
        !U(5) = 1.75107E+00     ! Thermocline depth

        !U(6) = 5.05238E+00      ! C thermocline
        !U(7) = 5.02118E+00    ! C ts
        !U(8) = 5.06684E+00   ! C northern
        !U(9) = 5.02233E+00     ! C southern

        !U(10) = 8.67955E+00    ! A thermocline
        !U(11) = 8.65307E+00    ! A ts
        !U(12) = 8.67988E+00   ! A northern
        !U(13) =  8.65059E+00      ! A southern

        !U(14) =  8.08467E-04   ! P thermocline
        !U(15) =  4.15473E-04   ! P ts
        !U(16) =  8.08467E-04   ! P northern
        !U(17) =  3.51681E-04   ! P southern

        !U(18) = 1.05661E-05 ! CO2
        !U(19) = 1.64754E-04    ! Total C
        
!     Trivial Solution
      !DO I = 6, 19
       ! U(I) = 0.0
      !ENDDO
      END SUBROUTINE STPNT
!----------------------------------------------------------------------

      SUBROUTINE BCND
      END SUBROUTINE BCND

    SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
      !---------- ----

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
    DOUBLE PRECISION, INTENT(IN) :: PAR(*)
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
    DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
    DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
    

    END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
