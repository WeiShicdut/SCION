%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%              110111010                                                                        %
%           111010-1-----101                                                                    %
%        1011111---------101111                                                                 %
%      11011------------------101         SCION: Spatial Continuous Integration                 %
%     111-----------------10011011        Earth Evolution Model                                 %
%    1--10---------------1111011111                                                             %
%    1---1011011---------1010110111       Coded by Benjamin J. W. Mills                         %
%    1---1011000111----------010011       email: b.mills@leeds.ac.uk                            %
%    1----1111011101----------10101                                                             %
%     1----1001111------------0111        Model initialiser                                     %
%      1----1101-------------1101         call this script to perform single runs               %
%        1--111----------------1                                                                %
%           1---------------1                                                                   %
%               111011011                                                                       %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run = SCION_initialise(runcontrol)

    %%%%%%% remove structures from pervious runs 
    clear stepnumber
    clear pars
    clear forcings
    clear workingstate
    clear switches
    clear state
    clear rawoutput
    clear options
    clear geoldata
    clear rawoutput
    clear resample
    %%%%%%% set up global structures
    global stepnumber
    global pars
    global forcings
    global workingstate
    global state 
    global gridstate 
    global INTERPSTACK
    global sensanal
    global plotrun
    global sensparams
    %%%% global tuning variables
    global Gtune
    global Ctune
    global PYRtune
    global GYPtune
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Check for sensitivity analysis   %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if runcontrol >= 1
        sensanal = 1 ;
        plotrun = 0 ;
        pars.telltime = 0 ;
    else
        sensanal = 0 ;
        plotrun = 1 ;
        pars.telltime = 1 ;
    end
    pars.runcontrol = runcontrol ;
    
    %%%%%%% starting to load params
    if sensanal == 0 
        fprintf('setting parameters... \t')
        tic
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Flux values at present   %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% reductant input
    pars.k_reductant_input = 0.4e12 ;  %%%% schopf and klein 1992
    
    %%%%%% water reservoir sizes in m3 (m=margins, s=surface, h= hi-lat, d=deep)
    pars.vol_p  = 2.6e15 ;  %%%% approx volume of all shelves and slope to depth 100m, area pecentage 5%
    pars.vol_di = 5.4e15 ;  %%%% approx volume of all shelves and slope in depth 100-1000m, area pecentage 5%
    pars.vol_s  = 2.75e16 ; %%%% approx volume of suface water to depth 100m, area pecentage 76.5%
    pars.vol_h  = 1.22e16 ; %%%% approx volume of hi-lat to depth 250m, area pecentage 13.5%
    pars.vol_d  = 1.35e18 ;
    pars.vol_ocean = pars.vol_p + pars.vol_di + pars.vol_s + pars.vol_h + pars.vol_d;

    %%%% mixing coefficient (Sv)
    pars.mixcoeff_dip = 30.28;
    pars.mixcoeff_ds  = 46.33;
    pars.mixcoeff_dh  = 54.9;

    %%%%%% initial inorganic carbon reservoirs in moles C
    pars.CO2_a_0  = 5e16 ;
    pars.DIC_p_0  = 5.2e15 ;
    pars.DIC_di_0 = 1.08e16 ;
    pars.DIC_s_0  = 5.37e16 ;
    pars.DIC_h_0  = 2.71e16 ;
    pars.DIC_d_0  = 3e18 ;
    pars.ALK_p_0  = 5.2e15 ;
    pars.ALK_di_0 = 1.08e16 ;
    pars.ALK_s_0  = 5.37e16 ;
    pars.ALK_h_0  = 2.71e16 ;
    pars.ALK_d_0  = 3e18 ;
    %%%%%% initial C isotope composition
    pars.d13c_atm_0    = -7 ;
    pars.d13c_DIC_p_0  = 0.1 ;
    pars.d13c_DIC_di_0 = 0.1 ;
    pars.d13c_DIC_s_0  = 0.1 ;
    pars.d13c_DIC_h_0  = 0.1 ;
    pars.d13c_DIC_d_0  = 0.1 ; 
   %%%%%% initial POC reservoirs in moles C
    pars.POC_p_0  = 630e12 ;
    pars.POC_di_0 = 250e12 ;
    pars.POC_s_0  = 2329e12 ;
    pars.POC_h_0  = 1084e12 ;
    pars.POC_d_0  = 56000e12 ;
    %%%%%% initial dissolved phosphate in moles
    pars.DP_p_0  = 1.82e12 ;
    pars.DP_di_0 = 7.56e12 ;
    pars.DP_s_0  = 0.55e12 ;
    pars.DP_h_0  = 16.27e12 ;
    pars.DP_d_0  = 2970e12 ;
    %%%%%% initial d13C of POC
    pars.d13c_POC_p_0  = -26;
    pars.d13c_POC_di_0 = -26;
    pars.d13c_POC_s_0  = -26 ;
    pars.d13c_POC_h_0  = -26 ;
    pars.d13c_POC_d_0  = -26 ;
    %%%%%% initial amount of O2 in moles C
    pars.O2_a_0  = 3.7e19 ;
    pars.O2_p_0  = 6.705e14 ;
    pars.O2_di_0 = 8.964e14 ;
    pars.O2_s_0  = 9.139e15 ;
    pars.O2_h_0  = 4.02e15 ;
    pars.O2_d_0  = 1.823e17 ;
    %%%%%% initial amount of FeIII in moles
    pars.FeIII_p_0  = 1.56e9 ;
    pars.FeIII_di_0 = 3.24e9 ;
    pars.FeIII_s_0  = 9.625e9 ;
    pars.FeIII_h_0  = 3.66e9 ;
    pars.FeIII_d_0  = 810e9 ;
    %%%%%% initial amount of sulfate in moles
    pars.SO4_p_0  = 7.28e16;
    pars.SO4_di_0 = 1.512e17;
    pars.SO4_s_0  = 7.7e17;
    pars.SO4_h_0  = 3.416e17;
    pars.SO4_d_0  = 3.78e19;
    %%%%%% initial amount of FeII in moles
    pars.FeII_p_0  = 0;
    pars.FeII_di_0 = 0;
    pars.FeII_s_0  = 0;
    pars.FeII_h_0  = 0;
    pars.FeII_d_0  = 0;
    %%%%%% initial amount of H2S in moles
    pars.H2S_p_0  = 0;
    pars.H2S_di_0 = 0;
    pars.H2S_s_0  = 0;
    pars.H2S_h_0  = 0;
    pars.H2S_d_0  = 0; 
    
    %%%%%% initial size of land S reservoirs in moles
    pars.GYP_0 = 9.375e19;
    pars.PYR_0 = 1.875e20;
   
    %%%% org C cycle
    pars.k_locb = 2.5e12 ;
    pars.k_mocb = 2.5e12 ; % 7e12 in MBOX frontend
    pars.k_ocdeg = 1.25e12 ;
    %%%% fluxes calculated for steady state
    pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg - pars.k_reductant_input ;

    %%%% carb C cycle
    pars.k_ccdeg = 12e12 ;
    pars.k_carbw = 8e12 ; % 12e12 in MBOX frontend
%     pars.k_sfw = 1.75e12 ; %%% delete?
    pars.k_mccb = pars.k_carbw + pars.k_ccdeg - pars.k_sfw ; % mass balance minus sfw
    pars.k_silw = pars.k_mccb - pars.k_carbw ; % used in line 703, P weathering
    pars.basfrac = 0.3 ;
    pars.k_granw = pars.k_silw * (1-pars.basfrac) ;
    pars.k_basw = pars.k_silw * pars.basfrac ;

    %%%% S cycle
%     pars.k_mpsb = 0.7e12 ; % stop
    pars.k_mgsb = 1.25e12 ;
    pars.k_pyrw = 1.85e12 ;
    pars.k_gypw = 1.25e12 ;
    pars.k_pyrdeg = 0 ; % need a new value
    pars.k_gypdeg = 0 ; % need a new value
    
    %%%% P cycle
    pars.k_capb = 2e10 ;
    pars.k_fepb = 1e10 ;
    pars.k_mopb = 1e10 ;
    pars.k_phosw = 4.25e10 ; % pars.k_phosw = 0.0967e12;
    pars.k_landfrac = 0.0588 ;
    
%     %%%% N cycle
%     pars.k_nfix = 8.67e12 ;
%     pars.k_denit = 4.3e12 ;
%     %%%% Sr cycle
%     pars.k_Sr_sedw = 17e9 ;
%     pars.k_Sr_mantle = 7.3e9 ;
%     pars.k_Sr_silw = 13e9 ;
%     pars.k_Sr_granw = pars.k_Sr_silw * (1 - pars.basfrac) ;
%     pars.k_Sr_basw = pars.k_Sr_silw * pars.basfrac ;
%     pars.total_Sr_removal = pars.k_Sr_granw + pars.k_Sr_basw + pars.k_Sr_sedw + pars.k_Sr_mantle ;
%     pars.k_Sr_sfw = pars.total_Sr_removal * ( pars.k_sfw / (pars.k_sfw + pars.k_mccb) ) ;
%     pars.k_Sr_sedb = pars.total_Sr_removal * ( pars.k_mccb / (pars.k_sfw + pars.k_mccb) ) ;
%     pars.k_Sr_metam = 13e9 ;

%     %%%% others
%     pars.k_oxfrac = 0.9975 ;
%     Pconc0 = 2.2 ;
%     Nconc0 = 30.9 ;
%     pars.newp0 = 117 * min(Nconc0/16,Pconc0) ;
    %COPSE constant for calculating pO2 from normalised O2
    pars.copsek16 = 3.762 ;
%     % oxidative weathering dependency on O2 concentration
%     pars.a = 0.5 ;
%     % marine organic carbon burial dependency on new production
%     pars.b = 2 ; 
    %%fire feedback
    pars.kfire = 3 ;

    %%%%%% Redfeild ratio
    pars.Red_C_P = 106;
    pars.Red_C_N = 106 / 16;
    pars.Red_C_O = 106 / 138;
    pars.Red_C_Fe = 106 * 2000;
    pars.Red_Fe_P = pars.Red_C_P / pars.Red_C_Fe;
    %%%%%% Monod constant; mol/m3
    pars.KP = 0.1e-3;
    pars.KFe = 0.1e-6;
    pars.KmO2 = 13e-3;
    pars.KmFeIII = 10;
    pars.KmSO4 = 0.5;
    %%%%%% reaction rate constant; m3 mol-1 yr-1
    pars.kpy = 0.3708/1e3 * 24 * 365.25;
    pars.kironO = 1.4e4;
    pars.ksulfO = 1e2;
    pars.kSironR = 8;
    pars.kSide = 4e3;
    %%%%%% Ksp
    pars.Kspside = 10^(-8.4) * 1e6;%%%mol2 m-6
    pars.KspFeSaq = 10^(-5.08) * 1e3;%%%mol m-3
    pars.STFeSaq = 10^(-5.7) * 1e3;%%%mol m-3
    %%%%%% dust Fe fluxes,mol yr-1
    pars.FeIIIa_s = 1.443e9; %%%
    pars.FeIIIa_h = 0;  %%%% 
    %%%%%% riverine solid FeIII fluxes, mol yr-1
    pars.sFeIII_p = 190e9;
    pars.sFeIII_di = 70e9;
    pars.sFeIII_d =  50e9;
    %%%%% reservoir present day sizes (mol)
    pars.G0 = 1.25 * 10^21 ;
    pars.C0 = 5 * 10^21 ;
    pars.PYR0 = 1.8 * 10^20 ;
    pars.GYP0 = 2 * 10^20 ;
    pars.CAL0 = 1.397e19 ;


    %%%% finished loading params
    if sensanal == 0 
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Load Forcings   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%% starting to load forcings
    if sensanal == 0 
        fprintf('loading forcings... \t')
        tic
    end

    %%%% load INTERPSTACK
    load( 'forcings/INTERPSTACK_oct_2021.mat' ) ;

    %%%% relative contribution from latitude bands
    lat_areas = (cosd(INTERPSTACK.lat))' ;
    for n=1:48
        pars.rel_contrib(:,n) = lat_areas / mean(lat_areas) ;
    end

    %%%% load COPSE reloaded forcing set
    load( 'forcings/COPSE_forcings.mat' ) 
    %%%% new BA 
    forcings.GR_BA = xlsread('forcings/GR_BA.xlsx','','','basic') ;
    forcings.GR_BA(:,1) = forcings.GR_BA(:,1)*1e6 ; %%% correct Myr
    %%%% new GA
    forcings.newGA = xlsread('forcings/GA_revised.xlsx','','','basic') ;
    forcings.newGA(:,1) = forcings.newGA(:,1)*1e6 ; %%% correct Myr
    %%%% degassing rate
    load('forcings/combined_D_force_oct_2021.mat') ;
    forcings.D_force_x = D_force_x ;
    forcings.D_force_mid = D_force_mid ;
    forcings.D_force_min = D_force_min ;
    forcings.D_force_max = D_force_max ;
    
    %%%% load shoreline forcing
    load('forcings/shoreline.mat') ;
    forcings.shoreline_time = shoreline_time ;
    forcings.shoreline_relative = shoreline_relative ;
    
    %%%%% finished loading forcings
    if sensanal == 0 
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Generate sensitivity randoms   %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if sensanal == 1
        %%%% generate random number in [-1 +1]
        sensparams.randminusplus1 = 2*(0.5 - rand) ;
        sensparams.randminusplus2 = 2*(0.5 - rand) ;
        sensparams.randminusplus3 = 2*(0.5 - rand) ;
        sensparams.randminusplus4 = 2*(0.5 - rand) ;
        sensparams.randminusplus5 = 2*(0.5 - rand) ;
        sensparams.randminusplus6 = 2*(0.5 - rand) ;
        sensparams.randminusplus7 = 2*(0.5 - rand) ;    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Initialise solver   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% run beginning
    if sensanal == 0 
        fprintf('Beginning run: \n')
    end

    %%%% if no plot or sensitivity command set to single run
    if isempty(sensanal) == 1
        sensanal = 0 ;
    end
    if isempty(plotrun) == 1
        plotrun = 1 ;
    end

    %%%%%%% model timeframe in years (0 = present day)
    pars.whenstart = - 600e6 ;
    pars.whenend = 0 ;

    %%%% setp up grid stamp times
    if runcontrol == -2
        pars.runstamps = 0 ;
    else
        pars.runstamps = INTERPSTACK.time( INTERPSTACK.time > ( pars.whenstart * 1e-6 ) ) ;
    end
    pars.next_gridstamp = pars.runstamps(1) ;
    pars.gridstamp_number = 1 ;
    pars.finishgrid = 0 ;

    %%%%%%% Show current timestep in command window? (1 = yes, 0 = no)
    pars.telltime = 1;

    %%%%%%% set number of model steps to take before beiling out
    pars.bailnumber = 1e5;

    %%%%%%% display every n model steps whilst running
    pars.display_resolution = 200 ;

    %%%%%%% set maximum step size for solver
    options = odeset('maxstep',1e6) ;

    %%%% set stepnumber to 1
    stepnumber = 1 ;

    %%%%%%% set starting reservoir sizes 
    pars.pstart = pars.P0;
    pars.tempstart = 288;
    pars.CAL_start = pars.CAL0;
    pars.N_start = pars.N0;
    pars.OSr_start = pars.OSr0;
    pars.SSr_start = pars.SSr0;
    pars.delta_A_start = 0 ;
    pars.delta_S_start = 35 ;
    pars.delta_G_start = -27 ;
    pars.delta_C_start = -2 ;
    pars.delta_PYR_start = -5 ;
    pars.delta_GYP_start = 20 ;
    pars.delta_OSr_start = 0.708 ;
    pars.delta_SSr_start = 0.708 ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Initial parameter tuning option  %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(Gtune) == 0
        pars.ostart = pars.O0 * abs( Otune )  ;
        pars.astart = pars.A0 * abs( Atune ) ;
        pars.sstart = pars.S0 * abs( Stune ) ;
        pars.gstart = pars.G0 * abs( Gtune ) ;
        pars.cstart = pars.C0 * abs( Ctune ) ;
        pars.pyrstart = pars.PYR0 * abs( PYRtune ) ;
        pars.gypstart = pars.GYP0 * abs( GYPtune ) ; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% if no tuning use previously tuned values
    if isempty(Gtune) == 1

        outputs = [ 0.55 1 1.2 1 0.1 0.05 3 ] ;
        
        pars.gstart = pars.G0 * outputs(1) ;
        pars.cstart = pars.C0 * outputs(2) ;
        pars.pyrstart = pars.PYR0 * outputs(3) ;
        pars.gypstart = pars.GYP0 * outputs(4) ; 
        pars.ostart = pars.O0 * outputs(5)  ;
        pars.sstart = pars.S0 * outputs(6) ;
        pars.astart = pars.A0 * outputs(7) ;

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% model start state
    pars.startstate(1) = pars.CO2_a_start ;
    pars.startstate(2) = pars.DIC_p_start ;
    pars.startstate(3) = pars.DIC_di_start ;
    pars.startstate(4) = pars.DIC_s_start ;
    pars.startstate(5) = pars.DIC_h_start ;
    pars.startstate(6) = pars.DIC_d_start ;
    pars.startstate(7) = pars.ALK_p_start ;
    pars.startstate(8) = pars.ALK_di_start ;
    pars.startstate(9) = pars.ALK_s_start ;
    pars.startstate(10) = pars.ALK_h_start ;
    pars.startstate(11) = pars.ALK_d_start ;
    pars.startstate(12) = pars.d13c_atm_start * pars.CO2_a_start;
    pars.startstate(13) = pars.d13c_DIC_p_start * pars.DIC_p_start;
    pars.startstate(14) = pars.d13c_DIC_di_start * pars.DIC_di_start;
    pars.startstate(15) = pars.d13c_DIC_s_start * pars.DIC_s_start;
    pars.startstate(16) = pars.d13c_DIC_h_start * pars.DIC_h_start ;
    pars.startstate(17) = pars.d13c_DIC_d_start * pars.DIC_d_start;         
    pars.startstate(18) = pars.POC_p_start ;
    pars.startstate(19) = pars.POC_di_start ;
    pars.startstate(20) = pars.POC_s_start ;
    pars.startstate(21) = pars.POC_h_start ;
    pars.startstate(22) = pars.POC_d_start ;
    pars.startstate(23) = pars.DP_p_start ;
    pars.startstate(24) = pars.DP_di_start ;
    pars.startstate(25) = pars.DP_s_start ;
    pars.startstate(26) = pars.DP_h_start ;
    pars.startstate(27) = pars.DP_d_start ;
    pars.startstate(28) = pars.d13c_POC_p_start * pars.POC_p_start;
    pars.startstate(29) = pars.d13c_POC_di_start * pars.POC_di_start;
    pars.startstate(30) = pars.d13c_POC_s_start * pars.POC_s_start;
    pars.startstate(31) = pars.d13c_POC_h_start * pars.POC_h_start;
    pars.startstate(32) = pars.d13c_POC_d_start * pars.POC_d_start;
    pars.startstate(33) = pars.O2_a_start ;
    pars.startstate(34) = pars.O2_p_start ;
    pars.startstate(35) = pars.O2_di_start ;
    pars.startstate(36) = pars.O2_s_start ;
    pars.startstate(37) = pars.O2_h_start ;
    pars.startstate(38) = pars.O2_d_start ;
    pars.startstate(39) = pars.FeIII_p_start ;
    pars.startstate(40) = pars.FeIII_di_start ;
    pars.startstate(41) = pars.FeIII_s_start ;
    pars.startstate(42) = pars.FeIII_h_start ;
    pars.startstate(43) = pars.FeIII_d_start ;
    pars.startstate(44) = pars.SO4_p_start ;
    pars.startstate(45) = pars.SO4_di_start ;
    pars.startstate(46) = pars.SO4_s_start ;
    pars.startstate(47) = pars.SO4_h_start ;
    pars.startstate(48) = pars.SO4_d_start ;
    pars.startstate(49) = pars.FeII_p_start ;
    pars.startstate(50) = pars.FeII_di_start ;
    pars.startstate(51) = pars.FeII_s_start ;
    pars.startstate(52) = pars.FeII_h_start ;
    pars.startstate(53) = pars.FeII_d_start ;
    pars.startstate(54) = pars.H2S_p_start ;
    pars.startstate(55) = pars.H2S_di_start ;
    pars.startstate(56) = pars.H2S_s_start ;
    pars.startstate(57) = pars.H2S_h_start ;
    pars.startstate(58) = pars.H2S_d_start ;
    pars.startstate(59) = pars.d34s_SO4_p_start * pars.SO4_p_start;
    pars.startstate(60) = pars.d34s_SO4_di_start * pars.SO4_di_start;
    pars.startstate(61) = pars.d34s_SO4_s_start * pars.SO4_s_start;
    pars.startstate(62) = pars.d34s_SO4_h_start * pars.SO4_h_start;
    pars.startstate(63) = pars.d34s_SO4_d_start * pars.SO4_d_start;
    pars.startstate(64) = pars.G_start ;
    pars.startstate(65) = pars.C_start ;
    pars.startstate(66) = pars.PYR_start ;
    pars.startstate(67) = pars.GYP_start ;
    pars.startstate(68) = pars.d13c_G_start ;
    pars.startstate(69) = pars.d13c_C_start ;
    pars.startstate(70) = pars.d34s_PYR_start ;
    pars.startstate(71) = pars.d34s_GYP_start ;
    
    %%% Vectors to store the results

    statevector =1:35

    state.pO2_af = statevector;
    state.Atmospheric_CO2_ppmf = statevector;
    state.DIC_conc_pf = statevector;
    state.DIC_conc_dif = statevector;
    state.DIC_conc_sf = statevector;
    state.DIC_conc_hf = statevector;
    state.DIC_conc_df = statevector;
    state.ALK_conc_pf = statevector;
    state.ALK_conc_dif = statevector;
    state.ALK_conc_sf = statevector;
    state.ALK_conc_hf = statevector;
    state.ALK_conc_df = statevector;
    state.pH_pf = statevector;
    state.pH_dif = statevector;
    state.pH_sf = statevector;
    state.pH_hf = statevector;
    state.pH_df = statevector;
    state.T_sf = statevector;
    state.T_hf = statevector;
    state.T_df = statevector;
    state.T_contf = statevector;
    state.GASTf = statevector;
    state.ccdegf = statevector;
    state.baswf = statevector;
    state.granwf = statevector;
    state.silwf = statevector;
    state.carbwf = statevector;
    state.mccb_pf = statevector;
    state.mccb_dif = statevector;
    state.mccb_df = statevector;
    state.POC_pf = statevector;
    state.POC_dif = statevector;
    state.POC_sf = statevector;
    state.POC_hf = statevector;
    state.POC_df = statevector;
    state.DP_conc_pf = statevector;
    state.DP_conc_dif = statevector;
    state.DP_conc_sf = statevector;
    state.DP_conc_hf = statevector;
    state.DP_conc_df = statevector;
    state.O2_conc_pf = statevector;
    state.O2_conc_dif = statevector;
    state.O2_conc_sf = statevector;
    state.O2_conc_hf = statevector;
    state.O2_conc_df = statevector;
    state.FeIII_conc_pf = statevector;
    state.FeIII_conc_dif = statevector;
    state.FeIII_conc_sf = statevector;
    state.FeIII_conc_hf = statevector;
    state.FeIII_conc_df = statevector;
    state.SO4_conc_pf = statevector;
    state.SO4_conc_dif = statevector;
    state.SO4_conc_sf = statevector;
    state.SO4_conc_hf = statevector;
    state.SO4_conc_df = statevector;
    state.FeII_conc_pf = statevector;
    state.FeII_conc_dif = statevector;
    state.FeII_conc_sf = statevector;
    state.FeII_conc_hf = statevector;
    state.FeII_conc_df = statevector;
    state.H2S_conc_pf = statevector;
    state.H2S_conc_dif = statevector;
    state.H2S_conc_sf = statevector;
    state.H2S_conc_hf = statevector;
    state.H2S_conc_df = statevector;
    state.O2_conc_pf = statevector;
    state.O2_conc_dif = statevector;
    state.O2_conc_sf = statevector;
    state.O2_conc_hf = statevector;
    state.O2_conc_df = statevector;
    state.FeIIIwf = statevector;
    state.FeIIIscavenging_pf = statevector;
    state.FeIIIscavenging_dif = statevector;
    state.FeIIIscavenging_sf = statevector;
    state.FeIIIscavenging_hf = statevector;
    state.FeIIIscavenging_df = statevector;
    state.pyrwf = statevector;
    state.gypwf = statevector;
    state.mgsbf = statevector;
    state.pyF_pf = statevector;
    state.pyF_dif = statevector;
    state.pyF_sf = statevector;
    state.pyF_hf = statevector;
    state.pyF_df = statevector;
    state.ironO_pf = statevector;
    state.ironO_dif = statevector;
    state.ironO_sf = statevector;
    state.ironO_hf = statevector;
    state.ironO_df = statevector;
    state.SironR_pf = statevector;
    state.SironR_dif = statevector;
    state.SironR_sf = statevector;
    state.SironR_hf = statevector;
    state.SironR_df = statevector;
    state.SideP_pf = statevector;
    state.SideP_dif = statevector;
    state.SideP_sf = statevector;
    state.SideP_hf = statevector;
    state.SideP_df = statevector;
    state.pripr_pf = statevector;
    state.pripr_sf = statevector;
    state.pripr_hf = statevector;
    state.remin_pf = statevector;
    state.remin_dif = statevector;
    state.remin_sf = statevector;
    state.remin_hf = statevector;
    state.remin_df = statevector;
    state.mocb_pf = statevector;
    state.mocb_dif = statevector;
    state.mocb_df = statevector;
    state.phoswf = statevector;
    state.sulfatebentic_pf = statevector;
    state.mgsbf = statevector;
    state.sulfR_pf = statevector;
    state.sulfO_pf = statevector;
    state.SironR_pf = statevector;
    state.pyF_pf = statevector;
    state.pyrwf = statevector;
    state.gypwf = statevector;
    state.SO4_lf = statevector;
    state.mocb_p_FeIIIf = statevector;
    state.mocb_di_FeIIIf = statevector;
    state.mocb_d_FeIIIf = statevector;
    state.ocbother_pf = statevector;
    state.ocbother_dif = statevector;
    state.ocbother_df = statevector;
    state.water_sediment_pf = statevector;
    state.water_sediment_dif = statevector;
    state.water_sediment_df = statevector;
    state.BEf_p = statevector;
    state.BEf_di = statevector;
    state.BEf_d = statevector;
    state.AR_pf = statevector;
    state.AR_dif = statevector;
    state.AR_sf = statevector;
    state.AR_hf = statevector;
    state.AR_df = statevector;
    state.ironR_pf = statevector;
    state.ironR_dif = statevector;
    state.ironR_sf = statevector;
    state.ironR_hf = statevector;
    state.ironR_df = statevector;
    state.sulfR_pf = statevector;
    state.sulfR_dif = statevector;
    state.sulfR_sf = statevector;
    state.sulfR_hf = statevector;
    state.sulfR_df = statevector;
    state.oxygenbentic_pf = statevector;
    state.oxygenbentic_dif = statevector;
    state.oxygenbentic_df = statevector;
    state.sulfatebentic_pf = statevector;
    state.sulfatebentic_dif = statevector;
    state.sulfatebentic_df = statevector;
    state.methanogenesis_pf = statevector;
    state.methanogenesis_dif = statevector;
    state.methanogenesis_df = statevector;

    %%%% note model start time
    tic

    %%%%%%% run the system 
    [rawoutput.T,rawoutput.Y] = ode15s(@SCION_equations,[pars.whenstart pars.whenend],pars.startstate,options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% size of output 
    pars.output_length = length(rawoutput.T) ;
        
    if sensanal == 0
        %%%%%%%%%% model finished output to screen
        fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length ) 
        toc
    end

    %%%%%%%%% print final model states using final state for each timepoint
    %%%%%%%%% during integration
    
    if sensanal == 0
    fprintf('assembling state vectors... \t')
    tic
    end

    %%%% trecords is index of shared values between ode15s output T vector and
    %%%% model recorded workingstate t vector
    [sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

    %%%%%% assemble output state vectors
    field_names = fieldnames(workingstate) ;
    for numfields = 1:length(field_names)
        eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
    end

    %%%%%% save state
    run.state = state ;
    run.gridstate = gridstate ;
    run.pars = pars ;
    run.forcings = forcings ;
    
    if sensanal == 0
        %%%%%% done message
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% only plot if no tuning structure exists, only plot fluxes for quick runs
    if isempty(Gtune) == 1
        if plotrun == 1            
            if runcontrol>-1
                SCION_plot_worldgraphic
            end
            SCION_plot_fluxes
        end
    end   
end
