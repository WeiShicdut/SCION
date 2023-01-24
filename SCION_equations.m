function dy = SCION_equations(t,y)

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
%     1----1001111------------0111        Model equations file                                  %
%      1----1101-------------1101         contains flux and reservoir equations                 %
%        1--111----------------1                                                                %
%           1---------------1                                                                   %
%               111011011                                                                       %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% setup dy array
dy = zeros(71,1);  

%%%%%%% set up global parameters
global stepnumber
global pars
global forcings
global workingstate
global gridstate
global INTERPSTACK
global sensanal
global sensparams

%%%%%%%%%%%%% get variables from Y to make working easier
CO2_a  = y(1) ;
DIC_p  = y(2) ;
DIC_di = y(3) ;
DIC_s  = y(4) ;
DIC_h  = y(5) ;
DIC_d  = y(6) ;
ALK_p  = y(7) ;
ALK_di = y(8) ;
ALK_s  = y(9) ;
ALK_h  = y(10) ;
ALK_d  = y(11) ;
d13c_atm    = y(12) / y(1) ;
d13c_DIC_p  = y(13) / y(2) ;
d13c_DIC_di = y(14) / y(3) ;
d13c_DIC_s  = y(15) / y(4) ;
d13c_DIC_h  = y(16) / y(5) ;
d13c_DIC_d  = y(17) / y(6) ;
POC_p  = y(18);
POC_di = y(19);
POC_s  = y(20);
POC_h  = y(21);
POC_d  = y(22);
DP_p  = y(23);
DP_di = y(24);
DP_s  = y(25);
DP_h  = y(26);
DP_d  = y(27);
d13c_POC_p  = y(28) / y(18) ;
d13c_POC_di = y(29) / y(19) ;
d13c_POC_s  = y(30) / y(20) ;
d13c_POC_h  = y(31) / y(21) ;
d13c_POC_d  = y(32) / y(22) ;
O2_a  = y(33);
O2_p  = y(34);
O2_di = y(35);
O2_s  = y(36);
O2_h  = y(37);
O2_d  = y(38);
FeIII_p  = y(39);
FeIII_di = y(40);
FeIII_s  = y(41);
FeIII_h  = y(42);
FeIII_d  = y(43);
SO4_p  = y(44);
SO4_di = y(45);
SO4_s  = y(46);
SO4_h  = y(47);
SO4_d  = y(48);
FeII_p  = y(49);
FeII_di = y(50);
FeII_s  = y(51);
FeII_h  = y(52);
FeII_d  = y(53);
H2S_p  = y(54);
H2S_di = y(55);
H2S_s  = y(56);
H2S_h  = y(57);
H2S_d  = y(58);
d34s_SO4_p  = y(59) / y(44) ;
d34s_SO4_di = y(60) / y(45) ;
d34s_SO4_s  = y(61) / y(46) ;
d34s_SO4_h  = y(62) / y(47) ;
d34s_SO4_d  = y(63) / y(48) ;
G = y(64);
C = y(65);
PYR = y(66);
GYP = y(67);
d13c_G = y(68) / y(64) ;
d13c_C = y(69) / y(65) ;
d34s_PYR = y(70) / y(66) ;
d34s_GYP = y(71) / y(67) ;


%%%% geological time in Ma (S=-2 sets a present day run)
if pars.runcontrol == -2
    t_geol = 0 ;
else
    t_geol = t*(1e-6) ;
end

%%%%%%% atmospheric fraction of total CO2, atfrac(A)
RCO2 = ( CO2_a / pars.CO2_a_0 )^2 ;
%%%% atmospheric CO2 in ppm
CO2ppm = RCO2 * 280 ;
CO2atm = RCO2 * (280e-6) ;
%%%% pCO2 in PAL
pCO2_a = CO2_a / pars.CO2_a_0 ;
%%%% pO2 in PAL
pO2_a = O2_a / pars.O2_a_0 ;

%%%%% mixing ratio of oxygen (not proportional to O reservoir)
mrO2 = ( O2_a / pars.O2_a_0 ) / (( O2_a / pars.O2_a_0 ) + pars.copsek16 ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Interpolate forcings for this timestep   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COPSE Reloaded forcing set
E_reloaded = interp1qr( forcings.t', forcings.E' , t_geol ) ;
W_reloaded = interp1qr( forcings.t', forcings.W' , t_geol ) ;
%%%% Additional forcings
GR_BA = interp1qr( forcings.GR_BA(:,1)./1e6 , forcings.GR_BA(:,2) , t_geol ) ;
newGA = interp1qr( forcings.newGA(:,1)./1e6 , forcings.newGA(:,2) , t_geol ) ;
D_combined_mid = interp1qr( forcings.D_force_x' , forcings.D_force_mid , t_geol) ;
D_combined_min = interp1qr( forcings.D_force_x' , forcings.D_force_min , t_geol) ;
D_combined_max = interp1qr( forcings.D_force_x' , forcings.D_force_max , t_geol) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Choose forcing functions  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%% uplift %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UPLIFT = 1 ;
%%%% degassing or metamorphism %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEGASS = D_sbz_rift ;
% DEGASS = D_merdith ;
DEGASS = D_combined_mid ;
%%%% Plant enhancement of weathering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = W_reloaded ;
%%%% Plant evolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EVO = E_reloaded ;
%%%% Plant stoichiometry, C/P on the land %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CPLAND = 1 ;
% CPLAND = CP_reloaded ;
%%%% pelagic calcification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bforcing = interp1qr([-1000 -150 -100 0]',[0.75 0.75 1 1]',t_geol) ;
%%%% basalt area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BAS_AREA = GR_BA ;
%%%% granite area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRAN_AREA = newGA ;
%%%% organic area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ORG_AREA = 1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PREPLANT = 1/4 ;

%%%% parametersfor C and S fractionation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
capdelS = 27 ; %%% why not 30 or 35 ?
capdelC_land = 27 ;
capdelDIC_g = 8 ; %%% why not 2 or 3 ?
capdelC_marine = 35 ;

%%%% SHORELINE
SHORELINE = interp1qr(forcings.shoreline_time',forcings.shoreline_relative',t_geol) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Sensitivity analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% all sensparams vary between [-1 +1]
if sensanal == 1
    
    %%%% Very degassing between upper and lower bounds
    if sensparams.randminusplus1 > 0
        DEGASS = (1 - sensparams.randminusplus1) * DEGASS + sensparams.randminusplus1*D_combined_max ;
    else
        DEGASS = (1 + sensparams.randminusplus1) * DEGASS - sensparams.randminusplus1*D_combined_min ;
    end
        
    %%%% simple +/- 20% variation
    BAS_AREA = BAS_AREA * (1 + 0.2*sensparams.randminusplus2) ;
    GRAN_AREA = GRAN_AREA * (1 + 0.2*sensparams.randminusplus3) ;
    
    %%%% preplant varies from 1/1 to 1/7
    PREPLANT = 1/ ( 4 + 3*sensparams.randminusplus4) ;
    
    %%%%
    capdelS = 30 + 10 * sensparams.randminusplus5 ;
    capdelC_land = 25 + 5 * sensparams.randminusplus6 ;
    capdelC_marine = 30 + 5 * sensparams.randminusplus7 ;
end

capdelB = capdelC_marine - capdelDIC_g ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Spatial fields from stack   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Fetch keyframe grids   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% find past and future keyframes
key_past_time = max ( INTERPSTACK.time( (INTERPSTACK.time - t_geol) <= 0 ) ) ;
key_future_time = min ( INTERPSTACK.time( (INTERPSTACK.time - t_geol) >= 0 ) ) ;
if isempty(key_past_time) == 1
    key_past_time = key_future_time ;
end

%%%% find keyframe indexes and fractional contribution
key_past_index = find( INTERPSTACK.time == key_past_time ) ;
key_future_index = find( INTERPSTACK.time == key_future_time ) ;
dist_to_past = abs( key_past_time - t_geol ) ;
dist_to_future = abs( key_future_time - t_geol ) ;
%%%% fractional contribution of each keyframe
if dist_to_past + dist_to_future == 0
    contribution_past = 1 ;
    contribution_future = 0 ;
else
    contribution_past = dist_to_future / ( dist_to_past + dist_to_future ) ;
    contribution_future = dist_to_past / ( dist_to_past + dist_to_future ) ;
end

%%%% intrepolate keyframe CO2 concentrations to generate keyframe fields
%%%% find INTERPSTACK keyframes using model CO2
key_upper_CO2 = min ( INTERPSTACK.CO2( (INTERPSTACK.CO2 - CO2ppm) >= 0 ) ) ;
key_lower_CO2 = max ( INTERPSTACK.CO2( (INTERPSTACK.CO2 - CO2ppm) <= 0 ) ) ;
%%%% find keyframe indexes and fractional contribution
key_upper_CO2_index = find( INTERPSTACK.CO2 == key_upper_CO2 )  ;
key_lower_CO2_index = find( INTERPSTACK.CO2 == key_lower_CO2 ) ;
dist_to_upper = abs( key_upper_CO2 - CO2ppm ) ;
dist_to_lower = abs( key_lower_CO2 - CO2ppm ) ;

%%%% fractional contribution of each keyframe
if dist_to_upper + dist_to_lower == 0
    contribution_lower = 1 ;
    contribution_upper = 0 ;
else
    contribution_upper = dist_to_lower / ( dist_to_upper + dist_to_lower ) ;
    contribution_lower = dist_to_upper / ( dist_to_upper + dist_to_lower ) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Create time keyframes using CO2 keyfield contributions   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Runoff
RUNOFF_past = contribution_upper.*INTERPSTACK.runoff(:,:,key_upper_CO2_index,key_past_index) + contribution_lower.*INTERPSTACK.runoff(:,:,key_lower_CO2_index,key_past_index); 
RUNOFF_future = contribution_upper.*INTERPSTACK.runoff(:,:,key_upper_CO2_index,key_future_index) + contribution_lower.*INTERPSTACK.runoff(:,:,key_lower_CO2_index,key_future_index); 
%%%% Tair
Tair_past = contribution_upper.*INTERPSTACK.Tair(:,:,key_upper_CO2_index,key_past_index) + contribution_lower.*INTERPSTACK.Tair(:,:,key_lower_CO2_index,key_past_index); 
Tair_future = contribution_upper.*INTERPSTACK.Tair(:,:,key_upper_CO2_index,key_future_index) + contribution_lower.*INTERPSTACK.Tair(:,:,key_lower_CO2_index,key_future_index); 

%%%% time kayframes that don't depend on CO2
%%%% Topography
TOPO_past = INTERPSTACK.topo(:,:,key_past_index) ; 
TOPO_future = INTERPSTACK.topo(:,:,key_future_index) ; 

%%%% last keyframe land recorded for plot
land_past = INTERPSTACK.land(:,:,key_past_index) ; 
land_future = INTERPSTACK.land(:,:,key_future_index) ; 

%%%% gridbox area
GRID_AREA_km2 = INTERPSTACK.gridarea ;

%%%% topographic slope from interpstack
tslope_past = INTERPSTACK.slope(:,:,key_past_index) ;
tslope_future = INTERPSTACK.slope(:,:,key_future_index) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Spatial silicate weathering   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% West / Maffre weathering approximation

%%%% runoff in mm/yr
Q_past = RUNOFF_past ;
Q_past(Q_past<0) = 0 ;
Q_future = RUNOFF_future ;
Q_future(Q_future<0) = 0 ;

%%%% temp in kelvin
T_past = Tair_past + 273 ;
T_future = Tair_future + 273 ;

%%%%% pierre erosion calculation, t/m2/yr
k_erosion = 3.3e-3 ; %%%% for 16Gt present day erosion in FOAM
EPSILON_past = k_erosion .* (Q_past.^0.31) .* tslope_past .* max(Tair_past,2) ;
EPSILON_future = k_erosion .* (Q_future.^0.31) .* tslope_future .* max(Tair_future,2) ;
%%%%% check total tonnes of erosion - should be ~16Gt
EPSILON_per_gridbox_past = EPSILON_past .* GRID_AREA_km2 .* 1e6 ; %%% t/m2/yr * m2
EPSILON_per_gridbox_future = EPSILON_future .* GRID_AREA_km2 .* 1e6 ; %%% t/m2/yr * m2
erosion_tot_past = sum(sum(EPSILON_per_gridbox_past)) ;
erosion_tot_future = sum(sum(EPSILON_per_gridbox_future)) ;
erosion_tot = erosion_tot_past*contribution_past + erosion_tot_future*contribution_future ;

%%%% Pierre weathering equation params
Xm = 0.1 ;
K = 6e-5 ; 
kw = 1e-3 ;
Ea = 20 ; 
z = 10 ; 
sigplus1 = 0.9 ; 
T0 = 286 ;
R = 8.31e-3 ;

%%%% equations
R_T_past = exp( ( Ea ./ (R.*T0) ) - ( Ea ./ (R.*T_past) ) ) ;
R_T_future = exp( ( Ea ./ (R.*T0) ) - ( Ea ./ (R.*T_future) ) ) ;
R_Q_past = 1 - exp( -1.*kw .* Q_past ) ;
R_Q_future = 1 - exp( -1.*kw .* Q_future ) ;
R_reg_past = ( (z./EPSILON_past).^sigplus1 ) ./ sigplus1 ;
R_reg_future = ( (z./EPSILON_future).^sigplus1 ) ./ sigplus1 ;

%%%% equation for CW per km2 in each box
CW_per_km2_past = 1e6 .* EPSILON_past .* Xm .* ( 1 - exp( -1.* K .* R_Q_past .* R_T_past .* R_reg_past ) ) ; 
CW_per_km2_future = 1e6 .* EPSILON_future .* Xm .* ( 1 - exp( -1.* K .* R_Q_future .* R_T_future .* R_reg_future ) ) ; 
%%%% CW total
CW_past = CW_per_km2_past .* GRID_AREA_km2 ;
CW_future = CW_per_km2_future .* GRID_AREA_km2 ;
%%%% world CW
CW_past(isnan(CW_past)==1) = 0 ;
CW_future(isnan(CW_future)==1) = 0 ;
CW_sum_past = sum(sum(CW_past)) ;
CW_sum_future = sum(sum(CW_future)) ;
CW_tot = CW_sum_past*contribution_past + CW_sum_future*contribution_future ;

%%%% carbonate weathering spatial approximation, linear with runoff
k_carb_scale = 200 ; %%%% scaling parameter to recover present day rate
CWcarb_per_km2_past = k_carb_scale * Q_past ; 
CWcarb_per_km2_future = k_carb_scale * Q_future ; 
%%%% CW total
CWcarb_past = CWcarb_per_km2_past .* GRID_AREA_km2 ;
CWcarb_future = CWcarb_per_km2_future .* GRID_AREA_km2 ;
%%%% world CWcarb
CWcarb_past(isnan(CWcarb_past)==1) = 0 ;
CWcarb_future(isnan(CWcarb_future)==1) = 0 ;
CWcarb_sum_past = sum(sum(CWcarb_past)) ;
CWcarb_sum_future = sum(sum(CWcarb_future)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Grid interpolated variables   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% silicate weathering scale factor by present day rate in cation tonnes
silw_scale = 4.2e8 ; %%%% for k erosion 3.3e-3
%%%% overall spatial weathering
silw_spatial = CW_tot * ( (pars.k_basw + pars.k_granw) / silw_scale) ;
carbw_spatial = ( CWcarb_sum_past*contribution_past + CWcarb_sum_future*contribution_future ) ;

%%%% global average surface temperature
GAST = mean(mean( Tair_past .* pars.rel_contrib )) * contribution_past + mean(mean( Tair_future .* pars.rel_contrib ))*contribution_future  ;
% GAST = 288 + Climate_Sensitivity * ( log( Atmospheric_CO2_ppm / 280 ) / log(2) )  - 7.4*(t_geol/-570)  ;
%%%% set assumed ice temperature
Tcrit = -10 ;
%%%% ice line calculations
Tair_past_ice = Tair_past ;
Tair_past_ice(Tair_past_ice >= Tcrit) = 0 ;
Tair_past_ice(Tair_past_ice < Tcrit) = 1 ;
Tair_future_ice = Tair_future ;
Tair_future_ice(Tair_future_ice >= Tcrit) = 0 ;
Tair_future_ice(Tair_future_ice < Tcrit) = 1 ;
%%%% count only continental ice
Tair_past_ice = Tair_past_ice.* land_past ;
Tair_future_ice = Tair_future_ice.* land_future ; 
%%%% sum into lat bands
latbands_past = sum(Tair_past_ice,2) ;
latbands_future = sum(Tair_future_ice,2) ;
latbands_past(latbands_past>0) = 1 ;
latbands_future(latbands_future>0) = 1 ;
%%%% find appropiate lat
latresults_past =  INTERPSTACK.lat .* latbands_past' ;
latresults_future =  INTERPSTACK.lat .* latbands_future' ;
latresults_past(latresults_past == 0) = 90 ;
latresults_future(latresults_future == 0) = 90 ;
%%%% lowest glacial latitude
iceline_past = min(abs(latresults_past)) ;
iceline_future = min(abs(latresults_future)) ;
iceline = iceline_past * contribution_past +  iceline_future * contribution_future ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% transfer flux calculations   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Thermohaline speed (Sv)
circ_TH_Sv = 10 ;
circ_coastal_Sv = 2 ; 

%%%% Thermohaline speed (m3/yr)
circ_TH_m3_yr = ( circ_TH_Sv * 1e6 ) * 3.15e7 ; 
circ_coastal_m3_yr = ( circ_coastal_Sv * 1e6 ) * 3.15e7 ;

%%%% mixing coefficient (m3/yr)
mixcoeff_dip_m3_yr = pars.mixcoeff_dip * 1e6 * 3.15e7 ;
mixcoeff_ds_m3_yr  = pars.mixcoeff_ds  * 1e6 * 3.15e7 ;
mixcoeff_dh_m3_yr  = pars.mixcoeff_dh  * 1e6 * 3.15e7 ;

%%%% DIC concentration in mol/m3
DIC_conc_p  = DIC_p  / pars.vol_p ;
DIC_conc_di = DIC_di / pars.vol_di ;
DIC_conc_s  = DIC_s  / pars.vol_s ;
DIC_conc_h  = DIC_h  / pars.vol_h ;
DIC_conc_d  = DIC_d  / pars.vol_d ;

ALK_conc_p  = ALK_p  / pars.vol_p ;
ALK_conc_di = ALK_di / pars.vol_di ;
ALK_conc_s  = ALK_s  / pars.vol_s ;
ALK_conc_h  = ALK_h  / pars.vol_h ;
ALK_conc_d  = ALK_d  / pars.vol_d ;

POC_conc_p  = POC_p  / pars.vol_p ;
POC_conc_di = POC_di / pars.vol_di ;
POC_conc_s  = POC_s  / pars.vol_s ;
POC_conc_h  = POC_h  / pars.vol_h ;
POC_conc_d  = POC_d  / pars.vol_d ;

DP_conc_p  = DP_p  / pars.vol_p ;
DP_conc_di = DP_di / pars.vol_di ;
DP_conc_s  = DP_s  / pars.vol_s ;
DP_conc_h  = DP_h  / pars.vol_h ;
DP_conc_d  = DP_d  / pars.vol_d ;

O2_conc_p  = O2_p  / pars.vol_p ;
O2_conc_di = O2_di / pars.vol_di ;
O2_conc_s  = O2_s  / pars.vol_s ;
O2_conc_h  = O2_h  / pars.vol_h ;
O2_conc_d  = O2_d  / pars.vol_d ;

FeIII_conc_p  = FeIII_p  / pars.vol_p ;
FeIII_conc_di = FeIII_di / pars.vol_di ;
FeIII_conc_s  = FeIII_s  / pars.vol_s ;
FeIII_conc_h  = FeIII_h  / pars.vol_h ;
FeIII_conc_d  = FeIII_d  / pars.vol_d ;

SO4_conc_p  = SO4_p  / pars.vol_p ;
SO4_conc_di = SO4_di / pars.vol_di ;
SO4_conc_s  = SO4_s  / pars.vol_s ;
SO4_conc_h  = SO4_h  / pars.vol_h ;
SO4_conc_d  = SO4_d  / pars.vol_d ;

FeII_conc_p  = FeII_p  / pars.vol_p ;
FeII_conc_di = FeII_di / pars.vol_di ;
FeII_conc_s  = FeII_s  / pars.vol_s ;
FeII_conc_h  = FeII_h  / pars.vol_h ;
FeII_conc_d  = FeII_d  / pars.vol_d ;

H2S_conc_p  = H2S_p  / pars.vol_p ;
H2S_conc_di = H2S_di / pars.vol_di ;
H2S_conc_s  = H2S_s  / pars.vol_s ;
H2S_conc_h  = H2S_h  / pars.vol_h ;
H2S_conc_d  = H2S_d  / pars.vol_d ;

%%%% Transport fluxes in mol/yr
Tran_DIC_p_s  = DIC_conc_p  * circ_coastal_m3_yr ;
Tran_DIC_s_h  = DIC_conc_s  * circ_TH_m3_yr ;
Tran_DIC_h_d  = DIC_conc_h  * circ_TH_m3_yr ;
Tran_DIC_d_di = DIC_conc_d  * circ_coastal_m3_yr ;
Tran_DIC_di_p = DIC_conc_di * circ_coastal_m3_yr ;
Tran_DIC_d_s  = DIC_conc_d  * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_ALK_p_s  = ALK_conc_p  * circ_coastal_m3_yr ;
Tran_ALK_s_h  = ALK_conc_s  * circ_TH_m3_yr ;
Tran_ALK_h_d  = ALK_conc_h  * circ_TH_m3_yr ;
Tran_ALK_d_di = ALK_conc_d  * circ_coastal_m3_yr ;
Tran_ALK_di_p = ALK_conc_di * circ_coastal_m3_yr ;
Tran_ALK_d_s  = ALK_conc_d  * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_POC_p_s  = POC_conc_p * circ_coastal_m3_yr ;
Tran_POC_s_h  = POC_conc_s * circ_TH_m3_yr ;
Tran_POC_p_di = 0.103 * POC_p ;
Tran_POC_s_d  = 0.167 * POC_s ;
Tran_POC_h_d  = 0.131 * POC_h ;
Tran_POC_d_di = POC_conc_d * circ_coastal_m3_yr ;

Tran_DP_p_s  = DP_conc_p  * circ_coastal_m3_yr ;
Tran_DP_s_h  = DP_conc_s  * circ_TH_m3_yr ;
Tran_DP_h_d  = DP_conc_h  * circ_TH_m3_yr ;
Tran_DP_d_di = DP_conc_d  * circ_coastal_m3_yr ;
Tran_DP_di_p = DP_conc_di * circ_coastal_m3_yr ;
Tran_DP_d_s  = DP_conc_d  * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_O2_p_s  = O2_conc_p  * circ_coastal_m3_yr ;
Tran_O2_s_h  = O2_conc_s  * circ_TH_m3_yr ;
Tran_O2_h_d  = O2_conc_h  * circ_TH_m3_yr ;
Tran_O2_d_di = O2_conc_d  * circ_coastal_m3_yr ;
Tran_O2_di_p = O2_conc_di * circ_coastal_m3_yr ;
Tran_O2_d_s  = O2_conc_d  * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_FeIII_p_s  = FeIII_conc_p  * circ_coastal_m3_yr ;
Tran_FeIII_s_h  = FeIII_conc_s  * circ_TH_m3_yr ;
Tran_FeIII_h_d  = FeIII_conc_h  * circ_TH_m3_yr ;
Tran_FeIII_d_di = FeIII_conc_d  * circ_coastal_m3_yr ;
Tran_FeIII_di_p = FeIII_conc_di * circ_coastal_m3_yr ;
Tran_FeIII_d_s  = FeIII_conc_d  * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_SO4_p_s  = SO4_conc_p  * circ_coastal_m3_yr ;
Tran_SO4_s_h  = SO4_conc_s  * circ_TH_m3_yr ;
Tran_SO4_h_d  = SO4_conc_h  * circ_TH_m3_yr ;
Tran_SO4_d_di = SO4_conc_d  * circ_coastal_m3_yr ;
Tran_SO4_di_p = SO4_conc_di * circ_coastal_m3_yr ;
Tran_SO4_d_s  = SO4_conc_d  * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_FeII_p_s  = FeII_conc_p  * circ_coastal_m3_yr ;
Tran_FeII_s_h  = FeII_conc_s  * circ_TH_m3_yr ;
Tran_FeII_h_d  = FeII_conc_h  * circ_TH_m3_yr ;
Tran_FeII_d_di = FeII_conc_d  * circ_coastal_m3_yr ;
Tran_FeII_di_p = FeII_conc_di * circ_coastal_m3_yr ;
Tran_FeII_d_s  = FeII_conc_d  * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_H2S_p_s  = H2S_conc_p  * circ_coastal_m3_yr ;
Tran_H2S_s_h  = H2S_conc_s  * circ_TH_m3_yr ;
Tran_H2S_h_d  = H2S_conc_h  * circ_TH_m3_yr ;
Tran_H2S_d_di = H2S_conc_d  * circ_coastal_m3_yr ;
Tran_H2S_di_p = H2S_conc_di * circ_coastal_m3_yr ;
Tran_H2S_d_s  = H2S_conc_d  * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

%%%% Water mixing terms
Mix_DP_di_p =  mixcoeff_ds_m3_yr * (DP_conc_di - DP_conc_p) ;
Mix_DP_d_s  =  mixcoeff_ds_m3_yr * (DP_conc_d  - DP_conc_s) ;
Mix_DP_d_h  =  mixcoeff_dh_m3_yr * (DP_conc_d  - DP_conc_h) ;

Mix_DIC_di_p = mixcoeff_ds_m3_yr * (DIC_conc_di - DIC_conc_p);
Mix_DIC_d_s =  mixcoeff_ds_m3_yr * (DIC_conc_d  - DIC_conc_s);
Mix_DIC_d_h =  mixcoeff_dh_m3_yr * (DIC_conc_d  - DIC_conc_h) ;

Mix_ALK_di_p = mixcoeff_ds_m3_yr * (ALK_conc_di - ALK_conc_p);
Mix_ALK_d_s  = mixcoeff_ds_m3_yr * (ALK_conc_d  - ALK_conc_s);
Mix_ALK_d_h  = mixcoeff_dh_m3_yr * (ALK_conc_d  - ALK_conc_h) ;

Mix_O2_di_p = mixcoeff_ds_m3_yr * (O2_conc_di - O2_conc_p);
Mix_O2_d_s  = mixcoeff_ds_m3_yr * (O2_conc_d  - O2_conc_s);
Mix_O2_d_h  = mixcoeff_dh_m3_yr * (O2_conc_d  - O2_conc_h) ;

Mix_FeIII_di_p = mixcoeff_ds_m3_yr * (FeIII_conc_di - FeIII_conc_p);
Mix_FeIII_d_s  = mixcoeff_ds_m3_yr * (FeIII_conc_d  - FeIII_conc_s);
Mix_FeIII_d_h  = mixcoeff_dh_m3_yr * (FeIII_conc_d  - FeIII_conc_h) ;

Mix_SO4_di_p = mixcoeff_ds_m3_yr * (SO4_conc_di - SO4_conc_p);
Mix_SO4_d_s  = mixcoeff_ds_m3_yr * (SO4_conc_d  - SO4_conc_s);
Mix_SO4_d_h  = mixcoeff_dh_m3_yr * (SO4_conc_d  - SO4_conc_h) ;

Mix_FeII_di_p = mixcoeff_ds_m3_yr * (FeII_conc_di - FeII_conc_p);
Mix_FeII_d_s  = mixcoeff_ds_m3_yr * (FeII_conc_d  - FeII_conc_s);
Mix_FeII_d_h  = mixcoeff_dh_m3_yr * (FeII_conc_d  - FeII_conc_h) ;

Mix_H2S_di_p = mixcoeff_ds_m3_yr * (H2S_conc_di - H2S_conc_p);
Mix_H2S_d_s  = mixcoeff_ds_m3_yr * (H2S_conc_d  - H2S_conc_s);
Mix_H2S_d_h  = mixcoeff_dh_m3_yr * (H2S_conc_d  - H2S_conc_h) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Global variables   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Global average surface temperature
% T_s = 298 + (GAST - 15) * 0.66  ;  
% T_h = max( 275.5 + (GAST - 15) , 271 ) ;
% T_d = max( 275.5 + (GAST - 15) , 271 ) ;
% T_p = T_s ;
% T_di = T_s ;

T_s = GAST ;
T_h = GAST ;
T_d = GAST ;
T_p = GAST ;
T_di = GAST ;

%%%% Oxygen level that is in equilibrium with atmosphere in mol/m3
O2_eq_ap = exp(-173.4292+249.6339*100/T_p+143.3483*log(T_p/100)-21.8492*T_p/100+35*(-0.033096+0.014259*T_p/100- 0.0017000*T_p/100*T_p/100))/22.4*O2_a / pars.O2_a_0 ;
O2_eq_as = exp(-173.4292+249.6339*100/T_s+143.3483*log(T_s/100)-21.8492*T_s/100+35*(-0.033096+0.014259*T_s/100- 0.0017000*T_s/100*T_s/100))/22.4*O2_a / pars.O2_a_0 ;
O2_eq_ah = exp(-173.4292+249.6339*100/T_h+143.3483*log(T_h/100)-21.8492*T_h/100+35*(-0.033096+0.014259*T_h/100- 0.0017000*T_h/100*T_h/100))/22.4*O2_a / pars.O2_a_0 ;

%%%% Atmosphere seawater oxygen exchange (mol/yr)
AirSeaO2_p = 3.572831e+14 * ( O2_eq_ap - O2_conc_p) ;
AirSeaO2_s = 1.476284e+15 * ( O2_eq_as - O2_conc_s) ;
AirSeaO2_h = 2.210798e+16 * ( O2_eq_ah - O2_conc_h) ;

%%%% Carbonate chemistry parameters 
k_2  = 7.4e-10 ;
k_carb_p = 0.000575 + 0.000006 * ( T_p - 278 ) ;
KCO2_p = 0.035 + 0.0019 * ( T_p - 278 ) ;
k_carb_di = 0.000575 + 0.000006 * ( T_di - 278 ) ;
k_carb_s = 0.000575 + 0.000006 * ( T_s - 278 ) ;
KCO2_s = 0.035 + 0.0019 * ( T_s - 278 ) ;
k_carb_h = 0.000575 + 0.000006 * ( T_h - 278 ) ;
KCO2_h = 0.035 + 0.0019 * ( T_h - 278 ) ;
k_carb_d = 0.000575 + 0.000006 * ( T_d - 278 ) ;

%%%% Carbonate speciation proximal
HCO3_p = ( DIC_conc_p - ( DIC_conc_p^2 - ALK_conc_p * ( 2 * DIC_conc_p - ALK_conc_p ) * ( 1 - 4 * k_carb_p ) )^0.5  ) / ( 1 - 4 * k_carb_p ) ;
CO3_p = ( ALK_conc_p - HCO3_p ) / 2 ;
H_p =  k_2 * HCO3_p / CO3_p ;
pH_p = -1 * log10(H_p) ;
%%%% Carbonate speciation distal
HCO3_di = ( DIC_conc_di - ( DIC_conc_di^2 - ALK_conc_di * ( 2 * DIC_conc_di - ALK_conc_di ) * ( 1 - 4 * k_carb_di ) )^0.5  ) / ( 1 - 4 * k_carb_di ) ;
CO3_di = ( ALK_conc_di - HCO3_di ) / 2 ;
H_di =  k_2 * HCO3_di / CO3_di ;
pH_di = -1 * log10(H_di) ;
%%%% Carbonate speciation lowlat
HCO3_s = ( DIC_conc_s - ( DIC_conc_s^2 - ALK_conc_s * ( 2 * DIC_conc_s - ALK_conc_s ) * ( 1 - 4 * k_carb_s ) )^0.5  ) / ( 1 - 4 * k_carb_s ) ;
CO3_s = ( ALK_conc_s - HCO3_s ) / 2 ;
H_s =  k_2 * HCO3_s / CO3_s ;
pH_s = -1 * log10(H_s) ;
%%%% Carbonate speciation hilat
HCO3_h = ( DIC_conc_h - ( DIC_conc_h^2 - ALK_conc_h * ( 2 * DIC_conc_h - ALK_conc_h ) * ( 1 - 4 * k_carb_h ) )^0.5  ) / ( 1 - 4 * k_carb_h ) ;
CO3_h = ( ALK_conc_h - HCO3_h ) / 2 ;
H_h =  k_2 * HCO3_h / CO3_h ;
pH_h = -1 * log10(H_h) ;
%%%% Carbonate speciation deep
HCO3_d = ( DIC_conc_d - ( DIC_conc_d^2 - ALK_conc_d * ( 2 * DIC_conc_d - ALK_conc_d ) * ( 1 - 4 * k_carb_d ) )^0.5 ) / ( 1 - 4 * k_carb_d ) ;
CO3_d = ( ALK_conc_d - HCO3_d ) / 2 ;
H_d = k_2 * HCO3_d / CO3_d ;
pH_d = -1 * log10(H_d) ;

%%%% Air-sea exchange of CO2 (mol/yr)
%%%% pCO2 in PAL
pCO2_p = KCO2_p * ( ( HCO3_p^2 ) / CO3_p ) ;
pCO2_s = KCO2_s * ( ( HCO3_s^2 ) / CO3_s ) ;
pCO2_h = KCO2_h * ( ( HCO3_h^2 ) / CO3_h ) ;
AirSea_p = 5e16 * 0.1 * 0.1 * ( pCO2_a - pCO2_p ) ;
AirSea_s = 5e16 * 0.765 * 0.1 * ( pCO2_a - pCO2_s ) ;
AirSea_h = 5e16 * 0.135 * 0.1 * ( pCO2_a - pCO2_h ) ;

AirSea_ap = 5e16 * 0.1 * 0.1 * ( pCO2_a) ;
AirSea_pa = 5e16 * 0.1 * 0.1 * ( pCO2_p) ;
AirSea_as = 5e16 * 0.765 * 0.1 * ( pCO2_a) ;
AirSea_sa = 5e16 * 0.765 * 0.1 * ( pCO2_s) ;
AirSea_ah = 5e16 * 0.135 * 0.1 * ( pCO2_a ) ;
AirSea_ha = 5e16 * 0.135 * 0.1 * ( pCO2_h ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Continental fluxes  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%55%%%% effect of temp on VEG %%%% fixed %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_T = 1 - (((GAST - 25)/25 )^2) ;
%%%% effect of CO2 on VEG
P_atm = CO2atm*1e6 ;
P_half = 183.6 ;
P_min = 10 ;
V_co2 = (P_atm - P_min) / (P_half + P_atm - P_min) ;
%%%% effect of O2 on VEG
V_o2 = 1.5 - 0.5 *(O2_a / pars.O2_a_0) ; 
%%%% full VEG limitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_npp = 2 * EVO * V_T * V_o2 * V_co2 ;
%%%% COPSE reloaded fire feedback
ignit = min(max(48 * mrO2 - 9.08, 0), 5) ;
firef = pars.kfire /(pars.kfire - 1 + ignit) ;
%%%%% Mass of terrestrial biosphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
VEG = V_npp * firef ;
%%%%%% basalt and granite temp dependency - direct and runoff
Tsurf = GAST + 273 ;
TEMP_gast = Tsurf ;
%%%% COPSE reloaded fbiota %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = VEG ;
f_biota = (1 - min(V * W, 1)) * PREPLANT * (RCO2^0.5) + (V * W) ;

%%%% basalt and granite weathering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basw = pars.k_basw * BAS_AREA * PG * f_biota * f_T_bas ;
basw = silw_spatial * ( pars.basfrac * BAS_AREA / ( pars.basfrac * BAS_AREA + (1 - pars.basfrac) * GRAN_AREA ) ) ;
% granw = pars.k_granw * UPLIFT^silconst * GRAN_AREA * PG * f_biota * f_T_gran ;
granw = silw_spatial * ((1 - pars.basfrac) * GRAN_AREA / (pars.basfrac * BAS_AREA + (1 - pars.basfrac) * GRAN_AREA )) ;
%%% silicate weathering %%%% add fbiota %%%%%%%%%%%%%%%%%%%%%%%
basw = basw * f_biota ;
granw = granw * f_biota ;
silw = basw + granw ;

%%%% carbonate weathering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% carbw = pars.k_carbw * CARB_AREA * UPLIFT^carbconst * PG * f_biota * g_T ;
carbw = carbw_spatial * f_biota ;
carbw_relative = (carbw / pars.k_carbw);

%%%% oxidative weathering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oxidw = pars.k_oxidw * carbw_relative * UPLIFT^pars.silconst * ORG_AREA * sigmf((O2_a/pars.O2_a_0),[5,-7]) * ((O2_a/pars.O2_a_0)^0.5) * (G/pars.G_0);

%%%% FeIII weathering fluxes mol yr-1 %%%%%%%%%%%%%%%%%%%%%%%%%
FeIIIw = pars.k_FeIIIw * silw /(pars.k_basw + pars.k_granw);

%%%% weathering fluxes of S mol yr-1 %%%%%%%%%%%%%%%%%%%%%%%%%%
pyrw = pars.k_pyrw * carbw_relative * 32 / (6e21) * PYR * sigmf(log10(O2_a/pars.O2_a_0),[5 -7]) * (O2_a/pars.O2_a_0)^0.5;
gypw = pars.k_gypw * carbw_relative * 32 / (3e21) * GYP;

% %%%%% seafloor weathering, revised following Brady and Gislason but not directly linking to CO2
% f_T_sfw = exp(0.0608 * (Tsurf - 288)); 
% sfw = pars.k_sfw * f_T_sfw * DEGASS; %%% assume spreading rate follows degassing here

%%%% P weathering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pfrac_silw = 0.8 ;
pfrac_carbw = 0.14 ;
pfrac_oxidw = 0.06 ;
phosw = pars.k_phosw * ((pfrac_silw)*(silw/pars.k_silw) + (pfrac_carbw)*(carbw/pars.k_carbw) + (pfrac_oxidw)*(oxidw/pars.k_oxidw )) ;
pland = pars.k_landfrac * VEG * phosw;
pland0 = pars.k_landfrac * pars.k_phosw;
psea = phosw - pland ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Degassing fluxes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccdeg  = pars.k_ccdeg  * DEGASS * (C/pars.C_0) * Bforcing ;
ocdeg  = pars.k_ocdeg  * DEGASS * (G/pars.G_0) * sigmf((O2_a/pars.O2_a_0),[5,-7]);
pyrdeg = pars.k_pyrdeg * DEGASS * (PYR/pars.PYR_0) * sigmf(log10(O2_a/pars.O2_a_0),[5 -7]);
gypdeg = pars.k_gypdeg * DEGASS * (GYP/pars.GYP_0);

%%%% generic CO2 input forcing to test model
CO2_input = 0;
d13c_CO2_input = -5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Burial fluxes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CaCO3 saturation proximal
Ca_conc_p = 1.397e19 / 1.35e18 ; %%%% assume constant hetrogeneous [Ca]
ksp_p = 0.76 ; %%% in mM^2
sat_p = (Ca_conc_p * CO3_p) / ksp_p ;
satpresent_p = 3 ;

%%%% CaCO3 saturation distal
Ca_conc_di = 1.397e19 / 1.35e18 ; %%%% assume constant hetrogeneous [Ca]
ksp_di = 0.76 ; %%% in mM^2
sat_di = (Ca_conc_di * CO3_di) / ksp_di ;
satpresent_di = 3 ;

%%%% CaCO3 saturation deep
Ca_conc_d = 1.397e19 / 1.35e18 ; %%%% assume constant hetrogeneous [Ca]
ksp_d = 0.76 ; %%% in mM^2
sat_d = ( Ca_conc_d * CO3_d ) / ksp_d ;
satpresent_d = 1 ;

%%%% carbonate burial 
deepfrac_mccb = 0.1 ;
distalfrac_mccb = 0.2 ;
sat_minus_1p = max(sat_p - 1, 0) ;
mccb_p = pars.k_mccb * (1 - deepfrac_mccb - distalfrac_mccb) * (1 / satpresent_p) * ((sat_minus_1p)^1.7) ;
%%%% distal carbonate burial
sat_minus_1di = max(sat_di - 1, 0) ;
mccb_di = pars.k_mccb * distalfrac_mccb * (1 / satpresent_di) * ((sat_minus_1di)^1.7) ;
%%%% deep carbonate burial
sat_minus_1d = max(sat_d - 1, 0) ;
mccb_d = pars.k_mccb * deepfrac_mccb * (1 / satpresent_d) * ((sat_minus_1d)^1.7) ;

%%%% total mccb
mccb = mccb_p + mccb_di + mccb_d ;

%%%% land-derived organic C burial
locb = pars.k_locb * (pland/pland0) * CPLAND ;
%%%% burial of sulfate
% sulfateb = pars.k_sulfateb*SO4_p/pars.SO4_p_0;
mgsb = pars.k_mgsb * SO4_p/pars.SO4_p_0 * (1/SHORELINE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Primary productivity  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mol/yr
pripr_p = 9.142857e+17 * min(DP_conc_p*DP_conc_p/(DP_conc_p+pars.KP),FeIII_conc_p/pars.Red_Fe_P*FeIII_conc_p/(FeIII_conc_p+pars.KFe));
pripr_s = 7.80e+20     * min(DP_conc_s*DP_conc_s/(DP_conc_s+pars.KP),FeIII_conc_s/pars.Red_Fe_P*FeIII_conc_s/(FeIII_conc_s+pars.KFe));
pripr_h = 2.22e+18     * min(DP_conc_h*DP_conc_h/(DP_conc_h+pars.KP),FeIII_conc_h/pars.Red_Fe_P*FeIII_conc_h/(FeIII_conc_h+pars.KFe));

%%% oxygen fluxes to atmosphere due to eutrophication
Oxygenfluxtoatmos = pripr_p /pars.Red_C_O * 0.2 * sigmf(pripr_p/280e12,[2 4]);
pripr_p_O = pripr_p /pars.Red_C_O *( 1- 0.2 * sigmf(pripr_p/280e12,[2 4]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   OC Remineralization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OC Remineralization via oxygen 
AR_p = 0.736  * POC_p  * O2_conc_p  / (O2_conc_p  + pars.KmO2);
AR_di = 0.240 * POC_di * O2_conc_di / (O2_conc_di + pars.KmO2);
AR_s = 0.967  * POC_s  * O2_conc_s  / (O2_conc_s  + pars.KmO2);
AR_h = 0.836  * POC_h  * O2_conc_h  / (O2_conc_h  + pars.KmO2);
AR_d = 0.0098 * POC_d  * O2_conc_d  / (O2_conc_d  + pars.KmO2);
%%% OC Remineralization via FeIII
ironR_p = 0.736  * POC_p  * pars.KmO2 / (O2_conc_p  + pars.KmO2) * FeIII_conc_p  / (FeIII_conc_p  + pars.KmFeIII);
ironR_di = 0.240 * POC_di * pars.KmO2 / (O2_conc_di + pars.KmO2) * FeIII_conc_di / (FeIII_conc_di + pars.KmFeIII);
ironR_s = 0.967  * POC_s  * pars.KmO2 / (O2_conc_s  + pars.KmO2) * FeIII_conc_s  / (FeIII_conc_s  + pars.KmFeIII);
ironR_h = 0.836  * POC_h  * pars.KmO2 / (O2_conc_h  + pars.KmO2) * FeIII_conc_h  / (FeIII_conc_h  + pars.KmFeIII);
ironR_d = 0.0098 * POC_d  * pars.KmO2 / (O2_conc_d  + pars.KmO2) * FeIII_conc_d  / (FeIII_conc_d  + pars.KmFeIII);
%%% OC Remineralization via sulfate
sulfR_p = 0.736  * POC_p  * 0.2 * pars.KmO2 /(O2_conc_p  + pars.KmO2) * pars.KmFeIII/(FeIII_conc_p  + pars.KmFeIII) * SO4_conc_p  / (SO4_conc_p  + pars.KmSO4);
sulfR_di = 0.240 * POC_di * 0.2 * pars.KmO2 /(O2_conc_di + pars.KmO2) * pars.KmFeIII/(FeIII_conc_di + pars.KmFeIII) * SO4_conc_di / (SO4_conc_di + pars.KmSO4);
sulfR_s = 0.967  * POC_s  * 0.2 * pars.KmO2 /(O2_conc_s  + pars.KmO2) * pars.KmFeIII/(FeIII_conc_s  + pars.KmFeIII) * SO4_conc_s  / (SO4_conc_s  + pars.KmSO4);
sulfR_h = 0.836  * POC_h  * 0.2 * pars.KmO2 /(O2_conc_h  + pars.KmO2) * pars.KmFeIII/(FeIII_conc_h  + pars.KmFeIII) * SO4_conc_h  / (SO4_conc_h  + pars.KmSO4);
sulfR_d = 0.0098 * POC_d  * 0.2 * pars.KmO2 /(O2_conc_d  + pars.KmO2) * pars.KmFeIII/(FeIII_conc_d  + pars.KmFeIII) * SO4_conc_d  / (SO4_conc_d  + pars.KmSO4);
%%% total OC Remineralization
remin_p  = AR_p  + ironR_p  + sulfR_p;
remin_di = AR_di + ironR_di + sulfR_di;
remin_s  = AR_s  + ironR_s  + sulfR_s;
remin_h  = AR_h  + ironR_h  + sulfR_h;
remin_d  = AR_d  + ironR_d  + sulfR_d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Other reactions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pyrite burial
pyF_p  = pars.kpy * max((FeII_conc_p  * H2S_conc_p /(pars.KspFeSaq + H2S_conc_p)  - pars.STFeSaq), 0) * H2S_conc_p  * pars.vol_p ;
pyF_di = pars.kpy * max((FeII_conc_di * H2S_conc_di/(pars.KspFeSaq + H2S_conc_di) - pars.STFeSaq), 0) * H2S_conc_di * pars.vol_di ;
pyF_s  = pars.kpy * max((FeII_conc_s  * H2S_conc_s /(pars.KspFeSaq + H2S_conc_s)  - pars.STFeSaq), 0) * H2S_conc_s  * pars.vol_s ;
pyF_h  = pars.kpy * max((FeII_conc_h  * H2S_conc_h /(pars.KspFeSaq + H2S_conc_h)  - pars.STFeSaq), 0) * H2S_conc_h  * pars.vol_h ;
pyF_d  = pars.kpy * max((FeII_conc_d  * H2S_conc_d /(pars.KspFeSaq + H2S_conc_d)  - pars.STFeSaq), 0) * H2S_conc_d  * pars.vol_d ;
%%% FeII oxidation via oxygen
ironO_p  = pars.kironO * FeII_conc_p  * O2_conc_p  * pars.vol_p ;
ironO_di = pars.kironO * FeII_conc_di * O2_conc_di * pars.vol_di ;
ironO_s  = pars.kironO * FeII_conc_s  * O2_conc_s  * pars.vol_s ;
ironO_h  = pars.kironO * FeII_conc_h  * O2_conc_h  * pars.vol_h ;
ironO_d  = pars.kironO * FeII_conc_d  * O2_conc_d  * pars.vol_d ;
%%% H2S oxidation via oxygen
sulfO_p  = pars.ksulfO * H2S_conc_p  * O2_conc_p  * pars.vol_p  ;
sulfO_di = pars.ksulfO * H2S_conc_di * O2_conc_di * pars.vol_di ;
sulfO_s  = pars.ksulfO * H2S_conc_s  * O2_conc_s  * pars.vol_s ;
sulfO_h  = pars.ksulfO * H2S_conc_h  * O2_conc_h  * pars.vol_h ;
sulfO_d  = pars.ksulfO * H2S_conc_d  * O2_conc_d  * pars.vol_d ;
%%% H2S oxidation via FeIII
SironR_p  = pars.kSironR * H2S_conc_p  * FeIII_conc_p  * pars.vol_p ;
SironR_di = pars.kSironR * H2S_conc_di * FeIII_conc_di * pars.vol_di ;
SironR_s  = pars.kSironR * H2S_conc_s  * FeIII_conc_s  * pars.vol_s ;
SironR_h  = pars.kSironR * H2S_conc_h  * FeIII_conc_h  * pars.vol_h ;
SironR_d  = pars.kSironR * H2S_conc_d  * FeIII_conc_d  * pars.vol_d ;
%%% saturation state (饱和状态) of siderite the apparent solubility consiant (表观溶解度常数)
omegaside_p  = FeII_conc_p  * CO3_p  / pars.Kspside; 
omegaside_di = FeII_conc_di * CO3_di / pars.Kspside;
omegaside_s  = FeII_conc_s  * CO3_s  / pars.Kspside;
omegaside_h  = FeII_conc_h  * CO3_h  / pars.Kspside;
omegaside_d  = FeII_conc_d  * CO3_d  / pars.Kspside;
%%% the precipitation (沉淀) rate of siderite
SideP_p  = pars.kSide * max((omegaside_p-1),0)  * pars.vol_p ;
SideP_di = pars.kSide * max((omegaside_di-1),0) * pars.vol_di ;
SideP_s  = pars.kSide * max((omegaside_s-1),0)  * pars.vol_s ;
SideP_h  = pars.kSide * max((omegaside_h-1),0)  * pars.vol_h ;
SideP_d  = pars.kSide * max((omegaside_d-1),0)  * pars.vol_d ;
%%%% FeIII scavenging
FeIIIscavenging_p  = max(3.18e+18*(FeIII_conc_p  - 0.58e-6),0);
FeIIIscavenging_di = max(1.54e+18*(FeIII_conc_di - 0.58e-6),0);
FeIIIscavenging_s  = max(3.63e+18*(FeIII_conc_s  - 0.58e-6),0);
FeIIIscavenging_h  = max(3.63e+18*(FeIII_conc_h  - 0.58e-6),0);
FeIIIscavenging_d  = max(3.63e+18*(FeIII_conc_d  - 0.58e-6),0);

%%%% marine organic C burial, P burial and benthic fluxes
water_sediment_p  = 0.0476  * POC_p;
water_sediment_di = 0.0400  * POC_di;
water_sediment_d  = 1.79e-4 * POC_d;

FeIIbenthic_p  = 62.2e9;
FeIIbenthic_di = 30.5e9;
FeIIbenthic_d  = 57.6e9;

%%% iron bound organic burial
mocb_p_FeIII  = min(0.0476  * POC_p,  4*(pars.sFeIII_p  + FeIIIscavenging_p  - FeIIbenthic_p));
mocb_di_FeIII = min(0.0400  * POC_di, 4*(pars.sFeIII_di + FeIIIscavenging_di - FeIIbenthic_di ));
mocb_d_FeIII  = min(1.79e-4 * POC_d,  4*(pars.sFeIII_d  + FeIIIscavenging_s + FeIIIscavenging_h + FeIIIscavenging_d - FeIIbenthic_d));

BE_p   = 0.128;
BE_di  = 0.125;
BE_d   = 0.076;

ocbother_p  = (0.0476 * POC_p  -  mocb_p_FeIII)  * BE_p;
ocbother_di = (0.04 * POC_di   -  mocb_di_FeIII) * BE_di;
ocbother_d  = (1.79e-4 * POC_d -  mocb_d_FeIII)  * BE_d;

mocb_p  = ocbother_p  + mocb_p_FeIII  ;
mocb_di = ocbother_di + mocb_di_FeIII ;
mocb_d  = ocbother_d  + mocb_d_FeIII ;

oxygenbentic_p  = -(water_sediment_p  - mocb_p)  / pars.Red_C_O * O2_conc_p /(O2_conc_p  + pars.KmO2);
oxygenbentic_di = -(water_sediment_di - mocb_di) / pars.Red_C_O * O2_conc_di/(O2_conc_di + pars.KmO2) ;
oxygenbentic_d  = -(water_sediment_d  - mocb_d)  / pars.Red_C_O * O2_conc_d /(O2_conc_d  + pars.KmO2);

sulfatebentic_p  = -0.5 *(water_sediment_p  - mocb_p) /pars.Red_C_O * pars.KmO2/(O2_conc_p + pars.KmO2) * SO4_conc_p /(SO4_conc_p  + pars.KmSO4); 
sulfatebentic_di = -0.5 *(water_sediment_di - mocb_di)/pars.Red_C_O * pars.KmO2/(O2_conc_di+ pars.KmO2) * SO4_conc_di/(SO4_conc_di + pars.KmSO4); 
sulfatebentic_d  = -0.5 *(water_sediment_d  - mocb_d) /pars.Red_C_O * pars.KmO2/(O2_conc_d + pars.KmO2) * SO4_conc_d /(SO4_conc_d  + pars.KmSO4); 

methanogenesis_p  = 0.5 *(water_sediment_p  - mocb_p) /pars.Red_C_O * pars.KmO2/(O2_conc_p + pars.KmO2) * pars.KmSO4 /(SO4_conc_p  + pars.KmSO4);
methanogenesis_di = 0.5 *(water_sediment_di - mocb_di)/pars.Red_C_O * pars.KmO2/(O2_conc_di+ pars.KmO2) * pars.KmSO4 /(SO4_conc_di + pars.KmSO4);
methanogenesis_d  = 0.5 *(water_sediment_d  - mocb_d) /pars.Red_C_O * pars.KmO2/(O2_conc_d + pars.KmO2) * pars.KmSO4 /(SO4_conc_d  + pars.KmSO4);

DICbenthic_p  = water_sediment_p  - mocb_p;
DICbenthic_di = water_sediment_di - mocb_di;
DICbenthic_d  = water_sediment_d  - mocb_d;

%%%% record total Corg burial
mocb = mocb_p + mocb_di + mocb_d ;

% f_anoxic_p = 1 - 0.9975 * O2_a / pars.O2_a_0 ;
% CPratio_p = CPoxic * CPanoxic / ( (1 - f_anoxic_p )*CPanoxic + f_anoxic_p*CPoxic ) ;
CPratio_p = 250 ;
POPburial_p = mocb_p /CPratio_p ;
PFe_p = 11.8e9 * DP_p/pars.DP_p_0;
Pauth_p = 24e9 * DP_p/pars.DP_p_0 ;

% f_anoxic_di = 1 - 0.9975 * O2_di / pars.O2_di_0 ;
% CPratio_di = CPoxic * CPanoxic / ( (1 - f_anoxic_p )*CPanoxic + f_anoxic_p*CPoxic ) ;
CPratio_di = 250 ;
POPburial_di = mocb_di  / CPratio_di ;
PFe_di = 3.93e9 * DP_di / pars.DP_di_0;
Pauth_di = 8e9  * DP_di / pars.DP_di_0 ; 

% f_anoxic_d = 1 - 0.9975 * O2_d / pars.O2_d_0 ;
% CPratio_d = CPoxic * CPanoxic / ( (1 - f_anoxic_p )*CPanoxic + f_anoxic_p*CPoxic ) ;
CPratio_d = 250 ;
POPburial_d = mocb_d / CPratio_d ;
PFe_d = 7e9  * DP_d / pars.DP_d_0;
Pauth_d = 14e9 * DP_d / pars.DP_d_0;

phosbenthic_p  = water_sediment_p  / pars.Red_C_P - POPburial_p  - PFe_p  - Pauth_p ;
phosbenthic_di = water_sediment_di / pars.Red_C_P - POPburial_di - PFe_di - Pauth_di ;
phosbenthic_d  = water_sediment_d  / pars.Red_C_P - POPburial_d  - PFe_d  - Pauth_d;


%%%% Fe2 hydrothermal input
FeIIhydro = pars.FeIIhydro ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Reservoir calculations   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% CO2_a
% dy(1) = ccdeg + ocdeg + oxidw - AirSea_p - AirSea_s - AirSea_h - locb - carbw - 2*silw + CO2_input ;

%%%% DIC_p
% dy(2) = AirSea_p + Tran_DIC_di_p - Tran_DIC_p_s + Mix_DIC_di_p + 2*carbw + 2*silw - mccb_p - pripr_p + remin_p + DICbenthic_p - SideP_p;

%%%% DIC_di
% dy(3) =  Tran_DIC_d_di - Tran_DIC_di_p - Mix_DIC_di_p - mccb_di  + remin_di + DICbenthic_di - SideP_di;

%%%% DIC_s
% dy(4) = AirSea_s + Tran_DIC_d_s + Tran_DIC_p_s + Mix_DIC_d_s - Tran_DIC_s_h - pripr_s + remin_s - SideP_s;

%%%% DIC_h
% dy(5) = AirSea_h + Tran_DIC_s_h + Mix_DIC_d_h - Tran_DIC_h_d - pripr_h + remin_h - SideP_h;

%%%% DIC_d
% dy(6) =  Tran_DIC_h_d - Tran_DIC_d_s - Tran_DIC_d_di - Mix_DIC_d_s - Mix_DIC_d_h - mccb_d + remin_d + DICbenthic_d - SideP_d;

%%%% ALK_p
% dy(7) =  Tran_ALK_di_p - Tran_ALK_p_s + Mix_ALK_di_p + 2 * carbw + 2 * silw - 2 * mccb_p - 2 * SideP_p + 0.25 * pyF_p - 2 * ironO_p - sulfO_p + 15 * SironR_p  -4*ironR_p + 0.5*sulfR_p - 3*FeIIIscavenging_p;

%%%% ALK_di
% dy(8) =  Tran_ALK_d_di - Tran_ALK_di_p - Mix_ALK_di_p - 2*mccb_di - 2*SideP_di + 0.25*pyF_di - 2*ironO_di - sulfO_di + 15*SironR_di  -4*ironR_di + 0.5*sulfR_di - 3*FeIIIscavenging_di;

%%%% ALK_s
% dy(9) =  Tran_ALK_d_s + Tran_ALK_p_s - Tran_ALK_s_h + Mix_ALK_d_s - 2*SideP_s + 0.25*pyF_s - 2*ironO_s - sulfO_s + 15*SironR_s -4*ironR_s + 0.5*sulfR_s - 3*FeIIIscavenging_s;

%%%% ALK_h
% dy(10) =  Tran_ALK_s_h - Tran_ALK_h_d + Mix_ALK_d_h - 2*SideP_h + 0.25*pyF_h - 2*ironO_h - sulfO_h + 15*SironR_h - 4*ironR_h + 0.5*sulfR_h - 3*FeIIIscavenging_h;

%%%% ALK_d
% dy(11) =  Tran_ALK_h_d - Tran_ALK_d_s - Tran_ALK_d_di - Mix_ALK_d_s - Mix_ALK_d_h - 2*mccb_d - 2*SideP_d + 0.25*pyF_d - 2*ironO_d - sulfO_d + 15*SironR_d - 4*ironR_d + 0.5*sulfR_d - 3*FeIIIscavenging_d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Carbon isotopes   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CO2_a * d13c_a
% dy(12) = - AirSea_ap*d13c_atm - AirSea_as*d13c_atm - AirSea_ah*d13c_atm + AirSea_pa*(d13c_DIC_p - capdelDIC_g) + AirSea_sa*(d13c_DIC_s - capdelDIC_g) + AirSea_ha*(d13c_DIC_h - capdelDIC_g)  + ccdeg*d13c_C + ocdeg*d13c_G + oxidw*d13c_G - locb*(d13c_atm - capdelC_land) - carbw*d13c_atm  - 2*silw*d13c_atm + CO2_input*d13c_CO2_input;

%%%% DIC_p * d13c_p
% dy(13) = AirSea_ap*d13c_atm - AirSea_pa*(d13c_DIC_p - capdelDIC_g) + Tran_DIC_di_p*d13c_DIC_di - Tran_DIC_p_s*d13c_DIC_p + mixcoeff_dip_m3_yr*(DIC_conc_di * d13c_DIC_di - DIC_conc_p*d13c_DIC_p) + carbw*d13c_C + carbw*d13c_atm + 2*silw*d13c_atm - mccb_p*d13c_DIC_p - pripr_p*(d13c_DIC_p - capdelB) + remin_p*d13c_POC_p + DICbenthic_p*d13c_POC_p - SideP_p*d13c_DIC_p ;

%%%% DIC_di * d13c_di
% dy(14) = Tran_DIC_d_di*d13c_DIC_d - Tran_DIC_di_p*d13c_DIC_di - mixcoeff_dip_m3_yr * (DIC_conc_di*d13c_DIC_di - DIC_conc_p*d13c_DIC_p)  - mccb_di*d13c_DIC_di + remin_di*d13c_POC_di + DICbenthic_di*d13c_POC_di - SideP_di*d13c_DIC_di ;

%%%% DIC_s * d13c_s
% dy(15) = AirSea_as*d13c_atm - AirSea_sa*(d13c_DIC_s - capdelDIC_g) + Tran_DIC_p_s * d13c_DIC_p + Tran_DIC_d_s*d13c_DIC_d - Tran_DIC_s_h*d13c_DIC_s + mixcoeff_ds_m3_yr*(DIC_conc_d *d13c_DIC_d - DIC_conc_s*d13c_DIC_s) - pripr_s*(d13c_DIC_s - capdelB) + remin_s*d13c_POC_s - SideP_s*d13c_DIC_s;

%%%% DIC_h * d13c_h
% dy(16) = AirSea_ah*d13c_atm - AirSea_ha*(d13c_DIC_h - capdelDIC_g) + Tran_DIC_s_h*d13c_DIC_s - Tran_DIC_h_d*d13c_DIC_h + mixcoeff_dh_m3_yr*(DIC_conc_d *d13c_DIC_d - DIC_conc_h*d13c_DIC_h) - pripr_h*(d13c_DIC_h - capdelB) + remin_h*d13c_POC_h - SideP_h*d13c_DIC_h;

%%%% DIC_d * d13c_d
% dy(17) = Tran_DIC_h_d*d13c_DIC_h - Tran_DIC_d_s*d13c_DIC_d - Tran_DIC_d_di*d13c_DIC_d - mixcoeff_ds_m3_yr*(DIC_conc_d *d13c_DIC_d - DIC_conc_s*d13c_DIC_s) - mixcoeff_dh_m3_yr*(DIC_conc_d *d13c_DIC_d - DIC_conc_h*d13c_DIC_h) - mccb_d*d13c_DIC_d + (remin_d+ DICbenthic_d)*d13c_POC_d - SideP_d*d13c_DIC_d;

%%%% POC_p
% dy(18) = pripr_p - remin_p - Tran_POC_p_s - water_sediment_p - Tran_POC_p_di;

%%%% POC_di
% dy(19) = - remin_di + Tran_POC_d_di + Tran_POC_p_di - water_sediment_di;

%%%% POC_s
% dy(20) = pripr_s - remin_s + Tran_POC_p_s - Tran_POC_s_h - Tran_POC_s_d ;

%%%% POC_h
% dy(21) = pripr_h - remin_h + Tran_POC_s_h - Tran_POC_h_d ;

%%%% POC_d
% dy(22) = - remin_d + Tran_POC_s_d + Tran_POC_h_d - Tran_POC_d_di - water_sediment_d ;

%%%% DP_p
% dy(23) = psea + phosbenthic_p - Tran_DP_p_s + Tran_DP_di_p + Mix_DP_di_p - pripr_p/pars.Red_C_P + remin_p/pars.Red_C_P ;

%%%% DP_di
% dy(24) = phosbenthic_di + Tran_DP_d_di - Tran_DP_di_p - Mix_DP_di_p + remin_di/pars.Red_C_P ;

%%%% DP_s
% dy(25) = Tran_DP_p_s + Tran_DP_d_s + Mix_DP_d_s + remin_s/pars.Red_C_P - pripr_s/pars.Red_C_P  - Tran_DP_s_h ;

%%%% DP_h
% dy(26) = Tran_DP_s_h - Tran_DP_h_d + Mix_DP_d_h - pripr_h/pars.Red_C_P + remin_h/pars.Red_C_P ;

%%%% DP_d
% dy(27) = phosbenthic_d + Tran_DP_h_d - Tran_DP_d_di - Tran_DP_d_s - Mix_DP_d_s - Mix_DP_d_h  + remin_d/pars.Red_C_P ;

%%%% d13c_POC_p
% dy(28) = pripr_p*(d13c_DIC_p - capdelB) - remin_p*d13c_POC_p - Tran_POC_p_s*d13c_POC_p - water_sediment_p*d13c_POC_p - Tran_POC_p_di*d13c_POC_p;

%%%% d13c_POC_di
% dy(29) =  - remin_di*d13c_POC_di + Tran_POC_d_di*d13c_POC_d + Tran_POC_p_di*d13c_POC_p - water_sediment_di*d13c_POC_di;

%%%% d13c_POC_s
% dy(30) = pripr_s*(d13c_DIC_s - capdelB) - remin_s*d13c_POC_s + Tran_POC_p_s*d13c_POC_p - Tran_POC_s_h*d13c_POC_s - Tran_POC_s_d*d13c_POC_s ;

%%%% d13c_POC_h
% dy(31) = pripr_h*(d13c_DIC_h - capdelB) - remin_h*d13c_POC_h + Tran_POC_s_h*d13c_POC_s - Tran_POC_h_d*d13c_POC_h ;

%%%% d13c_POC_d
% dy(32) = - remin_d*d13c_POC_d + Tran_POC_s_d*d13c_POC_s + Tran_POC_h_d*d13c_POC_h - Tran_POC_d_di*d13c_POC_d - water_sediment_d*d13c_POC_d ;

%%%% O2_a
% dy(33) = - ocdeg/pars.Red_C_O  - oxidw/pars.Red_C_O  + locb/pars.Red_C_O - AirSeaO2_p - AirSeaO2_s - AirSeaO2_h - 2*pyrw  + Oxygenfluxtoatmos - 2*(methanogenesis_p  + methanogenesis_di + methanogenesis_d) ; 

%%%% O2_p
% dy(34) = oxygenbentic_p + AirSeaO2_p + Tran_O2_di_p - Tran_O2_p_s + Mix_O2_di_p  + pripr_p_O - AR_p/pars.Red_C_O - 0.25*ironO_p - 2*sulfO_p;

%%%% O2_di
% dy(35) = oxygenbentic_di + Tran_O2_d_di - Tran_O2_di_p - Mix_O2_di_p - AR_di/pars.Red_C_O - 0.25*ironO_di - 2*sulfO_di;

%%%% O2_s
% dy(36) = AirSeaO2_s + Tran_O2_p_s + Tran_O2_d_s - Tran_O2_s_h + Mix_O2_d_s + pripr_s/pars.Red_C_O - AR_s/pars.Red_C_O - 0.25*ironO_s - 2*sulfO_s;

%%%% O2_h
% dy(37) = AirSeaO2_h + Tran_O2_s_h - Tran_O2_h_d  + Mix_O2_d_h + pripr_h/pars.Red_C_O - AR_h/pars.Red_C_O - 0.25*ironO_h - 2*sulfO_h;

%%%% O2_d
% dy(38) = oxygenbentic_d + Tran_O2_h_d - Tran_O2_d_di - Tran_O2_d_s - Mix_O2_d_s - Mix_O2_d_h - AR_d/pars.Red_C_O - 0.25*ironO_d - 2*sulfO_d ;

%%%% FeIII_p
% dy(39) = FeIIIw + Tran_FeIII_di_p - Tran_FeIII_p_s + Mix_FeIII_di_p - pripr_p/pars.Red_C_Fe + remin_p/pars.Red_C_Fe - FeIIIscavenging_p - 4*ironR_p + ironO_p - 8*SironR_p;

%%%% FeIII_di
% dy(40) = Tran_FeIII_d_di - Tran_FeIII_di_p - Mix_FeIII_di_p + remin_di/pars.Red_C_Fe - FeIIIscavenging_di - 4*ironR_di + ironO_di - 8*SironR_di;

%%%% FeIII_s
% dy(41) = pars.FeIIIa_s + Tran_FeIII_p_s - Tran_FeIII_s_h + Tran_FeIII_d_s + Mix_FeIII_d_s - pripr_s/pars.Red_C_Fe + remin_s/pars.Red_C_Fe - FeIIIscavenging_s - 4*ironR_s  + ironO_s - 8*SironR_s;

%%%% FeIII_h
% dy(42) = pars.FeIIIa_h + Tran_FeIII_s_h - Tran_FeIII_h_d + Mix_FeIII_d_h - pripr_h/pars.Red_C_Fe + remin_h/pars.Red_C_Fe - FeIIIscavenging_h - 4*ironR_h  + ironO_h - 8*SironR_h;

%%%% FeIII_d
% dy(43) = - Tran_FeIII_d_di - Tran_FeIII_d_s + Tran_FeIII_h_d - Mix_FeIII_d_s - Mix_FeIII_d_h + remin_d/pars.Red_C_Fe - FeIIIscavenging_d - 4*ironR_d  + ironO_d - 8*SironR_d;

%%%% SO4_p
% dy(44) = sulfatebentic_p + pyrw + gypw - mgsb + Tran_SO4_di_p - Tran_SO4_p_s + Mix_SO4_di_p - 0.5*sulfR_p/pars.Red_C_O + sulfO_p + SironR_p - 0.25*pyF_p;

%%%% SO4_di
% dy(45) = sulfatebentic_di + Tran_SO4_d_di + sulfO_di + SironR_di - Tran_SO4_di_p - Mix_SO4_di_p - 0.5*sulfR_di/pars.Red_C_O - 0.25*pyF_di;

%%%% SO4_s
% dy(46) = Tran_SO4_p_s - Tran_SO4_s_h + Tran_SO4_d_s + Mix_SO4_d_s - 0.5*sulfR_s/pars.Red_C_O + sulfO_s + SironR_s - 0.25*pyF_s;

%%%% SO4_h
% dy(47) = Tran_SO4_s_h - Tran_SO4_h_d + Mix_SO4_d_h - 0.5*sulfR_h/pars.Red_C_O + sulfO_h + SironR_h - 0.25*pyF_h;

%%%% SO4_d
% dy(48) = sulfatebentic_d + Tran_SO4_h_d - Tran_SO4_d_di - Tran_SO4_d_s - Mix_SO4_d_s - Mix_SO4_d_h - 0.5*sulfR_d/pars.Red_C_O + sulfO_d + SironR_d - 0.25*pyF_d;

%%%% FeII_p
% dy(49) = FeIIbenthic_p + Tran_FeII_di_p - Tran_FeII_p_s + Mix_FeII_di_p + 4*ironR_p - pyF_p - ironO_p + 8*SironR_p - SideP_p;

%%%% FeII_di
% dy(50) = FeIIbenthic_di + Tran_FeII_d_di - Tran_FeII_di_p - Mix_FeII_di_p + 4*ironR_di - pyF_di - ironO_di + 8*SironR_di - SideP_di;

%%%% FeII_s
% dy(51) = Tran_FeII_p_s - Tran_FeII_s_h + Tran_FeII_d_s + Mix_FeII_d_s + 4*ironR_s - pyF_s - ironO_s + 8*SironR_s - SideP_s;

%%%% FeII_h
% dy(52) = Tran_FeII_s_h - Tran_FeII_h_d + Mix_FeII_d_h + 4*ironR_h - pyF_h - ironO_h + 8*SironR_h - SideP_h;

%%%% FeII_d
% dy(53) = FeIIbenthic_d + FeIIhydro + Tran_FeII_h_d - Tran_FeII_d_di - Tran_FeII_d_s - Mix_FeII_d_s - Mix_FeII_d_h + 4*ironR_d - pyF_d - ironO_d + 8*SironR_d - SideP_d;

%%%% H2S_p
% dy(54) = Tran_H2S_di_p - Tran_H2S_p_s + Mix_H2S_di_p + 0.5*sulfR_p/pars.Red_C_O - 1.75*pyF_p - sulfO_p - SironR_p;

%%%% H2S_di
% dy(55) = Tran_H2S_d_di - Tran_H2S_di_p - Mix_H2S_di_p  + 0.5*sulfR_di/pars.Red_C_O - 1.75*pyF_di - sulfO_di - SironR_di;

%%%% H2S_s
% dy(56) = Tran_H2S_p_s - Tran_H2S_s_h + Tran_H2S_d_s + Mix_H2S_d_s + 0.5*sulfR_s/pars.Red_C_O - 1.75*pyF_s - sulfO_s - SironR_s;

%%%% H2S_h
% dy(57) = Tran_H2S_s_h - Tran_H2S_h_d + Mix_H2S_d_h + 0.5*sulfR_h/pars.Red_C_O - 1.75*pyF_h - sulfO_h - SironR_h;

%%%% H2S_d
% dy(58) = Tran_H2S_h_d - Tran_H2S_d_di - Tran_H2S_d_s - Mix_H2S_d_s - Mix_H2S_d_h + 0.5*sulfR_d/pars.Red_C_O - 1.75*pyF_d - sulfO_d - SironR_d;

%%%% d34s_SO4_p
% dy(59) = sulfatebentic_p*(d34s_SO4_p - capdelS) + (pyrw + pyrdeg)*d34s_PYR + (gypw + gypdeg)*d34s_GYP - mgsb*d34s_SO4_p + Tran_SO4_di_p*d34s_SO4_di - Tran_SO4_p_s*d34s_SO4_p + Mix_SO4_di_p*d34s_SO4_di - 0.5*sulfR_p/pars.Red_C_O*(d34s_SO4_p - capdelS) + sulfO_p*(d34s_SO4_p - capdelS) + SironR_p*(d34s_SO4_p - capdelS) - 0.25*pyF_p*(d34s_SO4_p - capdelS);

%%%% d34s_SO4_di
% dy(60) = sulfatebentic_di*(d34s_SO4_di - capdelS) + Tran_SO4_d_di*d34s_SO4_d + sulfO_di*(d34s_SO4_di - capdelS) + SironR_di*(d34s_SO4_di - capdelS) - Tran_SO4_di_p*d34s_SO4_di - Mix_SO4_di_p*d34s_SO4_di - 0.5*sulfR_di/pars.Red_C_O*(d34s_SO4_p - capdelS) - 0.25*pyF_di*(d34s_SO4_p - capdelS);

%%%% d34s_SO4_s
% dy(61) = Tran_SO4_p_s*d34s_SO4_p - Tran_SO4_s_h*d34s_SO4_s + Tran_SO4_d_s*d34s_SO4_d + Mix_SO4_d_s*d34s_SO4_d - 0.5*sulfR_s/pars.Red_C_O*(d34s_SO4_s - capdelS) + sulfO_s*(d34s_SO4_s - capdelS) + SironR_s*(d34s_SO4_s - capdelS) - 0.25*pyF_s*(d34s_SO4_s - capdelS);

%%%% d34s_SO4_h
% dy(62) = Tran_SO4_s_h*d34s_SO4_s - Tran_SO4_h_d*d34s_SO4_h + Mix_SO4_d_h*d34s_SO4_d - 0.5*sulfR_h/pars.Red_C_O*(d34s_SO4_h - capdelS) + sulfO_h*(d34s_SO4_h - capdelS) + SironR_h*(d34s_SO4_h - capdelS) - 0.25*pyF_h*(d34s_SO4_h - capdelS);

%%%% d34s_SO4_d
% dy(63) = sulfatebentic_d * (d34s_SO4_d - capdelS) + Tran_SO4_h_d*d34s_SO4_h - Tran_SO4_d_di*d34s_SO4_d - Tran_SO4_d_s*d34s_SO4_d - Mix_SO4_d_s*d34s_SO4_d - Mix_SO4_d_h*d34s_SO4_d - 0.5*sulfR_d/pars.Red_C_O*(d34s_SO4_d - capdelS) + sulfO_d*(d34s_SO4_d - capdelS) + SironR_d*(d34s_SO4_d - capdelS) - 0.25*pyF_d*(d34s_SO4_d - capdelS);

%%% Buried organic carbon (G)
% dy(64) = locb + mocb_p + mocb_di + mocb_d - oxidw - ocdeg ;

%%% Buried carbnate (C)
% dy(65) = mccb_p + mccb_di + mccb_d  - carbw - ccdeg ;

%%% Buried pyrite (PYR)
% dy(66) = pyF_p + pyF_di + pyF_s + pyF_h + pyF_d - sulfatebentic_p - sulfatebentic_di - sulfatebentic_d - pyrw - pyrdeg;

%%% Buried gypsum (GYP)
% dy(67) = mgsb - gypw - gypdeg ;

%%% delta_G * G (young)
% dy(68) = locb*(d13c_atm - capdelC_land) + mocb_p*(d13c_DIC_p - capdelB) + mocb_di*(d13c_DIC_di - capdelB) + mocb_d*(d13c_DIC_d - capdelB) - oxidw*d13c_G - ocdeg*d13c_G ;

%%% delta_C * C (young)
% dy(69) = mccb_p*d13c_DIC_p + mccb_di*d13c_DIC_di + mccb_d*d13c_DIC_d - carbw*d13c_C - ccdeg*d13c_C ; %%%% How to sea floor weathering

%%% delta_PYR * PYR (young)
% dy(70) = pyF_p*(d34s_SO4_p - capdelS) + pyF_di*(d34s_SO4_di - capdelS) + pyF_s*(d34s_SO4_s - capdelS) + pyF_h*(d34s_SO4_h - capdelS) + pyF_d*(d34s_SO4_d - capdelS) - sulfatebentic_p*(d34s_SO4_p - capdelS) - sulfatebentic_di*(d34s_SO4_di - capdelS) - sulfatebentic_d*(d34s_SO4_d - capdelS) - pyrw*d34s_PYR - pyrdeg*d34s_PYR ;

%%% delta_GYP * GYP (young)
%%% Buried gypsum (GYP)
% dy(71) = mgsb*d34s_SO4_p - gypw*d34s_GYP - gypdeg*d34s_GYP ;

%%%% reductant input
reductant_input = pars.k_reductant_input * DEGASS ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Record full states for single run   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sensanal == 0
    workingstate.TEMP_gast(stepnumber,1) = TEMP_gast ;
    workingstate.tempC(stepnumber,1) = TEMP_gast - 273 ;
    workingstate.Tsurf(stepnumber,1) = Tsurf ;
    workingstate.DEGASS(stepnumber,1) = DEGASS ;
    workingstate.W(stepnumber,1) = W ;
    workingstate.EVO(stepnumber,1) = EVO ;
    workingstate.CPLAND(stepnumber,1) = CPLAND ;
    workingstate.Bforcing(stepnumber,1) = Bforcing ;
    workingstate.BAS_AREA(stepnumber,1) = BAS_AREA ;
    workingstate.GRAN_AREA(stepnumber,1) = GRAN_AREA ;
    workingstate.VEG(stepnumber,1) = VEG ;
    workingstate.G(stepnumber,1) = G ;
    workingstate.C(stepnumber,1) = C ;
    workingstate.PYR(stepnumber,1) = PYR ;
    workingstate.GYP(stepnumber,1) = GYP ;
    workingstate.CO2_input(stepnumber,1) = CO2_input ;
    workingstate.d13c_CO2_input(stepnumber,1) = d13c_CO2_input ;
    workingstate.CO2_a(stepnumber,1) = CO2_a ;
    workingstate.CO2ppm(stepnumber,1) = CO2ppm ;
    workingstate.CO2atm(stepnumber,1) = CO2atm ;
    workingstate.pCO2_a(stepnumber,1) = pCO2_a ;
    workingstate.pO2_a(stepnumber,1) = pO2_a ;
    workingstate.DIC_p(stepnumber,1) = DIC_p ;
    workingstate.DIC_di(stepnumber,1) = DIC_di ;
    workingstate.DIC_s(stepnumber,1) = DIC_s ;
    workingstate.DIC_h(stepnumber,1) = DIC_h ;
    workingstate.DIC_d(stepnumber,1) = DIC_d ;
    workingstate.ALK_p(stepnumber,1) = ALK_p ;
    workingstate.ALK_di(stepnumber,1) = ALK_di ;
    workingstate.ALK_s(stepnumber,1) = ALK_s ;
    workingstate.ALK_h(stepnumber,1) = ALK_h ;
    workingstate.ALK_d(stepnumber,1) = ALK_d ;
    workingstate.DIC_conc_p(stepnumber,1) = DIC_conc_p ;
    workingstate.DIC_conc_di(stepnumber,1) = DIC_conc_di ;
    workingstate.DIC_conc_s(stepnumber,1) = DIC_conc_s ;
    workingstate.DIC_conc_h(stepnumber,1) = DIC_conc_h ;
    workingstate.DIC_conc_d(stepnumber,1) = DIC_conc_d ;
    workingstate.ALK_conc_p(stepnumber,1) = ALK_conc_p ;
    workingstate.ALK_conc_di(stepnumber,1) = ALK_conc_di ;
    workingstate.ALK_conc_s(stepnumber,1) = ALK_conc_s ;
    workingstate.ALK_conc_h(stepnumber,1) = ALK_conc_h ;
    workingstate.ALK_conc_d(stepnumber,1) = ALK_conc_d ;
    workingstate.pH_p(stepnumber,1) = pH_p ;
    workingstate.pH_di(stepnumber,1) = pH_di ;
    workingstate.pH_s(stepnumber,1) = pH_s ;
    workingstate.pH_h(stepnumber,1) = pH_h ;
    workingstate.pH_d(stepnumber,1) = pH_d ;
    workingstate.T_p(stepnumber,1) = T_p ;
    workingstate.T_di(stepnumber,1) = T_di ;
    workingstate.T_s(stepnumber,1) = T_s ;
    workingstate.T_h(stepnumber,1) = T_h ;
    workingstate.T_d(stepnumber,1) = T_d ;
    workingstate.GAST(stepnumber,1) = GAST ;
    workingstate.granw(stepnumber,1) = granw ;
    workingstate.basw(stepnumber,1) = basw ;
    workingstate.silw(stepnumber,1) = silw ;
    workingstate.carbw(stepnumber,1) = carbw ;
    workingstate.phosw(stepnumber,1) = phosw ;
    workingstate.pland(stepnumber,1) = pland ;
    workingstate.psea(stepnumber,1) = psea ;
    workingstate.phosw(stepnumber,1) = phosw ;
    workingstate.ccdeg(stepnumber,1) = ccdeg ;
    workingstate.gypdeg(stepnumber,1) = gypdeg ;
    workingstate.mccb_p(stepnumber,1) = mccb_p ;
    workingstate.mccb_di(stepnumber,1) = mccb_di ;
    workingstate.mccb_d(stepnumber,1) = mccb_d ;
    workingstate.mccb(stepnumber,1) = mccb ;
    workingstate.mocb_p(stepnumber,1) = mocb_p ;
    workingstate.mocb_di(stepnumber,1) = mocb_di ;
    workingstate.mocb_d(stepnumber,1) = mocb_d ;
    workingstate.mocb(stepnumber,1) = mocb ;
    workingstate.HCO3_p(stepnumber,1) = HCO3_p ;
    workingstate.CO3_p(stepnumber,1) = CO3_p ;
    workingstate.H_p(stepnumber,1) = H_p ;
    workingstate.HCO3_di(stepnumber,1) = HCO3_di ;
    workingstate.CO3_di(stepnumber,1) = CO3_di ;
    workingstate.H_di(stepnumber,1) = H_di ;
    workingstate.HCO3_s(stepnumber,1) = HCO3_s ;
    workingstate.CO3_s(stepnumber,1) = CO3_s ;
    workingstate.H_s(stepnumber,1) = H_s ;
    workingstate.HCO3_h(stepnumber,1) = HCO3_h ;
    workingstate.CO3_h(stepnumber,1) = CO3_h ;
    workingstate.H_h(stepnumber,1) = H_h ;
    workingstate.HCO3_d(stepnumber,1) = HCO3_d ;
    workingstate.CO3_d(stepnumber,1) = CO3_d ;
    workingstate.H_d(stepnumber,1) = H_d ;
    workingstate.d13c_atm(stepnumber,1) = d13c_atm ;
    workingstate.d13c_DIC_p(stepnumber,1) = d13c_DIC_p ;
    workingstate.d13c_DIC_di(stepnumber,1) = d13c_DIC_di ;
    workingstate.d13c_DIC_s(stepnumber,1) = d13c_DIC_s ;
    workingstate.d13c_DIC_h(stepnumber,1) = d13c_DIC_h ;
    workingstate.d13c_DIC_d(stepnumber,1) = d13c_DIC_d ;
    workingstate.POC_p(stepnumber,1) = POC_p ;
    workingstate.POC_di(stepnumber,1) = POC_di ;
    workingstate.POC_s(stepnumber,1) = POC_s ;
    workingstate.POC_h(stepnumber,1) = POC_h ;
    workingstate.POC_d(stepnumber,1) = POC_d ;
    workingstate.POC_conc_p(stepnumber,1) = POC_conc_p ;
    workingstate.POC_conc_di(stepnumber,1) = POC_conc_di ;
    workingstate.POC_conc_s(stepnumber,1) = POC_conc_s ;
    workingstate.POC_conc_h(stepnumber,1) = POC_conc_h ;
    workingstate.POC_conc_d(stepnumber,1) = POC_conc_d ;
    workingstate.DP_p(stepnumber,1) = DP_p ;
    workingstate.DP_di(stepnumber,1) = DP_di ;
    workingstate.DP_s(stepnumber,1) = DP_s ;
    workingstate.DP_h(stepnumber,1) = DP_h ;
    workingstate.DP_d(stepnumber,1) = DP_d ;
    workingstate.DP_conc_p(stepnumber,1) = DP_conc_p ;
    workingstate.DP_conc_di(stepnumber,1) = DP_conc_di ;
    workingstate.DP_conc_s(stepnumber,1) = DP_conc_s ;
    workingstate.DP_conc_h(stepnumber,1) = DP_conc_h ;
    workingstate.DP_conc_d(stepnumber,1) = DP_conc_d ;
    workingstate.d13c_POC_p(stepnumber,1) = d13c_POC_p ;
    workingstate.d13c_POC_di(stepnumber,1) = d13c_POC_di ;
    workingstate.d13c_POC_s(stepnumber,1) = d13c_POC_s ;
    workingstate.d13c_POC_h(stepnumber,1) = d13c_POC_h ;
    workingstate.d13c_POC_d(stepnumber,1) = d13c_POC_d ;
    workingstate.pO2_a(stepnumber,1) = pO2_a ;
    workingstate.O2_conc_p(stepnumber,1) = O2_conc_p ;
    workingstate.O2_conc_di(stepnumber,1) = O2_conc_di ;
    workingstate.O2_conc_s(stepnumber,1) = O2_conc_s ;
    workingstate.O2_conc_h(stepnumber,1) = O2_conc_h ;
    workingstate.O2_conc_d(stepnumber,1) = O2_conc_d ;
    workingstate.FeIII_conc_p(stepnumber,1) = FeIII_conc_p ;
    workingstate.FeIII_conc_di(stepnumber,1) = FeIII_conc_di ;
    workingstate.FeIII_conc_s(stepnumber,1) = FeIII_conc_s ;
    workingstate.FeIII_conc_h(stepnumber,1) = FeIII_conc_h ;
    workingstate.FeIII_conc_d(stepnumber,1) = FeIII_conc_d ;
    workingstate.pripr_p(stepnumber,1) = pripr_p ;
    workingstate.pripr_s(stepnumber,1) = pripr_s ;
    workingstate.pripr_h(stepnumber,1) = pripr_h ;
    workingstate.remin_p(stepnumber,1) = remin_p ;
    workingstate.remin_di(stepnumber,1) = remin_di ;
    workingstate.remin_s(stepnumber,1) = remin_s ;
    workingstate.remin_h(stepnumber,1) = remin_h ;
    workingstate.remin_d(stepnumber,1) = remin_d ;
    workingstate.SO4_conc_p(stepnumber,1) = SO4_conc_p ;
    workingstate.SO4_conc_di(stepnumber,1) = SO4_conc_di ;
    workingstate.SO4_conc_s(stepnumber,1) = SO4_conc_s ;
    workingstate.SO4_conc_h(stepnumber,1) = SO4_conc_h ;
    workingstate.SO4_conc_d(stepnumber,1) = SO4_conc_d ;
    workingstate.FeII_conc_p(stepnumber,1) = FeII_conc_p ;
    workingstate.FeII_conc_di(stepnumber,1) = FeII_conc_di ;
    workingstate.FeII_conc_s(stepnumber,1) = FeII_conc_s ;
    workingstate.FeII_conc_h(stepnumber,1) = FeII_conc_h ;
    workingstate.FeII_conc_d(stepnumber,1) = FeII_conc_d ;
    workingstate.H2S_conc_p(stepnumber,1) = H2S_conc_p ;
    workingstate.H2S_conc_di(stepnumber,1) = H2S_conc_di ;
    workingstate.H2S_conc_s(stepnumber,1) = H2S_conc_s ;
    workingstate.H2S_conc_h(stepnumber,1) = H2S_conc_h ;
    workingstate.H2S_conc_d(stepnumber,1) = H2S_conc_d ;
    workingstate.Tran_FeIII_p_s(stepnumber,1) = Tran_FeIII_p_s ;
    workingstate.Tran_FeIII_s_h(stepnumber,1) = Tran_FeIII_s_h ;
    workingstate.Tran_FeIII_d_s(stepnumber,1) = Tran_FeIII_d_s ;
    workingstate.Mix_FeIII_d_s(stepnumber,1) = Mix_FeIII_d_s ;
    workingstate.O2_a(stepnumber,1) = O2_a ;
    workingstate.O2_p(stepnumber,1) = O2_p ;
    workingstate.O2_di(stepnumber,1) = O2_di ;
    workingstate.O2_s(stepnumber,1) = O2_s ;
    workingstate.O2_h(stepnumber,1) = O2_h ;
    workingstate.O2_d(stepnumber,1) = O2_d ;
    workingstate.FeIII_p(stepnumber,1) = FeIII_p ;
    workingstate.FeIII_di(stepnumber,1) = FeIII_di ;
    workingstate.FeIII_s(stepnumber,1) = FeIII_s ;
    workingstate.FeIII_h(stepnumber,1) = FeIII_h ;
    workingstate.FeIII_d(stepnumber,1) = FeIII_d ;
    workingstate.SO4_p(stepnumber,1) = SO4_p ;
    workingstate.SO4_di(stepnumber,1) = SO4_di ;
    workingstate.SO4_s(stepnumber,1) = SO4_s ;
    workingstate.SO4_h(stepnumber,1) = SO4_h ;
    workingstate.SO4_d(stepnumber,1) = SO4_d ;
    workingstate.FeII_p(stepnumber,1) = FeII_p ;
    workingstate.FeII_di(stepnumber,1) = FeII_di ;
    workingstate.FeII_s(stepnumber,1) = FeII_s ;
    workingstate.FeII_h(stepnumber,1) = FeII_h ;
    workingstate.FeII_d(stepnumber,1) = FeII_d ;
    workingstate.H2S_p(stepnumber,1) = H2S_p ;
    workingstate.H2S_di(stepnumber,1) = H2S_di ;
    workingstate.H2S_s(stepnumber,1) = H2S_s ;
    workingstate.H2S_h(stepnumber,1) = H2S_h ;
    workingstate.H2S_d(stepnumber,1) = H2S_d ;
    workingstate.DICbenthic_p(stepnumber,1) = DICbenthic_p ;
    workingstate.DICbenthic_di(stepnumber,1) = DICbenthic_di ;
    workingstate.DICbenthic_d(stepnumber,1) = DICbenthic_d ;
    workingstate.FeIIbenthic_p(stepnumber,1) = FeIIbenthic_p ;
    workingstate.FeIIbenthic_d(stepnumber,1) = FeIIbenthic_d ;
    workingstate.FeIIIw(stepnumber,1) = FeIIIw;
    workingstate.FeIIIscavenging_p(stepnumber,1) = FeIIIscavenging_p;
    workingstate.FeIIIscavenging_di(stepnumber,1) = FeIIIscavenging_di;
    workingstate.FeIIIscavenging_s(stepnumber,1) = FeIIIscavenging_s;
    workingstate.FeIIIscavenging_h(stepnumber,1) = FeIIIscavenging_h;
    workingstate.FeIIIscavenging_d(stepnumber,1) = FeIIIscavenging_d;
    workingstate.pyF_p(stepnumber,1) = pyF_p;
    workingstate.pyF_di(stepnumber,1) = pyF_di;
    workingstate.pyF_s(stepnumber,1) = pyF_s;
    workingstate.pyF_h(stepnumber,1) = pyF_h;
    workingstate.pyF_d(stepnumber,1) = pyF_d;
    workingstate.ironO_p(stepnumber,1) = ironO_p;
    workingstate.ironO_di(stepnumber,1) = ironO_di;
    workingstate.ironO_s(stepnumber,1) = ironO_s;
    workingstate.ironO_h(stepnumber,1) = ironO_h;
    workingstate.ironO_d(stepnumber,1) = ironO_d;
    workingstate.SironR_p(stepnumber,1) = SironR_p;
    workingstate.SironR_di(stepnumber,1) = SironR_di;
    workingstate.SironR_s(stepnumber,1) = SironR_s;
    workingstate.SironR_h(stepnumber,1) = SironR_h;
    workingstate.SironR_d(stepnumber,1) = SironR_d;
    workingstate.SideP_p(stepnumber,1) = SideP_p;
    workingstate.SideP_di(stepnumber,1) = SideP_di;
    workingstate.SideP_s(stepnumber,1) = SideP_s;
    workingstate.SideP_h(stepnumber,1) = SideP_h;
    workingstate.SideP_d(stepnumber,1) = SideP_d;
    workingstate.ocdeg(stepnumber,1) = ocdeg;
    workingstate.oxidw(stepnumber,1) = oxidw;
    workingstate.locb(stepnumber,1) = locb;
    workingstate.AirSeaO2_p(stepnumber,1) = AirSeaO2_p;
    workingstate.AirSeaO2_s(stepnumber,1) = AirSeaO2_s;
    workingstate.AirSeaO2_h(stepnumber,1) = AirSeaO2_h;
    workingstate.Oxygenfluxtoatmos(stepnumber,1) = Oxygenfluxtoatmos;
    workingstate.sulfR_p(stepnumber,1) = sulfR_p;
    workingstate.sulfO_p(stepnumber,1) = sulfO_p;
    workingstate.SironR_p(stepnumber,1) = SironR_p;
    workingstate.sulfatebentic_p(stepnumber,1) = sulfatebentic_p;
    workingstate.sulfatebentic_di(stepnumber,1) = sulfatebentic_di;
    workingstate.sulfatebentic_d(stepnumber,1) = sulfatebentic_d;
    workingstate.pyrw(stepnumber,1) = pyrw;
    workingstate.gypw(stepnumber,1) = gypw;
    workingstate.mgsb(stepnumber,1) = mgsb;
    workingstate.mocb_p_FeIII(stepnumber,1) = mocb_p_FeIII;
    workingstate.mocb_di_FeIII(stepnumber,1) = mocb_di_FeIII;
    workingstate.mocb_d_FeIII(stepnumber,1) = mocb_d_FeIII;
    workingstate.ocbother_p(stepnumber,1) = ocbother_p;
    workingstate.ocbother_di(stepnumber,1) = ocbother_di;
    workingstate.ocbother_d(stepnumber,1) = ocbother_d;
    workingstate.water_sediment_p(stepnumber,1) = water_sediment_p;
    workingstate.water_sediment_di(stepnumber,1) = water_sediment_di;
    workingstate.water_sediment_d(stepnumber,1) = water_sediment_d;
    workingstate.AR_p(stepnumber,1) = AR_p;
    workingstate.AR_di(stepnumber,1) = AR_di;
    workingstate.AR_s(stepnumber,1) = AR_s;
    workingstate.AR_h(stepnumber,1) = AR_h;
    workingstate.AR_d(stepnumber,1) = AR_d;
    workingstate.ironR_p(stepnumber,1) = ironR_p;
    workingstate.ironR_di(stepnumber,1) = ironR_di;
    workingstate.ironR_s(stepnumber,1) = ironR_s;
    workingstate.ironR_h(stepnumber,1) = ironR_h;
    workingstate.ironR_d(stepnumber,1) = ironR_d;
    workingstate.sulfR_p(stepnumber,1) = sulfR_p;
    workingstate.sulfR_di(stepnumber,1) = sulfR_di;
    workingstate.sulfR_s(stepnumber,1) = sulfR_s;
    workingstate.sulfR_h(stepnumber,1) = sulfR_h;
    workingstate.sulfR_d(stepnumber,1) = sulfR_d;
    workingstate.oxygenbentic_p(stepnumber,1) = oxygenbentic_p;
    workingstate.oxygenbentic_di(stepnumber,1) = oxygenbentic_di;
    workingstate.oxygenbentic_d(stepnumber,1) = oxygenbentic_d;
    workingstate.sulfatebentic_p(stepnumber,1) = sulfatebentic_p;
    workingstate.sulfatebentic_di(stepnumber,1) = sulfatebentic_di;
    workingstate.sulfatebentic_d(stepnumber,1) = sulfatebentic_d;
    workingstate.methanogenesis_p(stepnumber,1) = methanogenesis_p;
    workingstate.methanogenesis_di(stepnumber,1) = methanogenesis_di;
    workingstate.methanogenesis_d(stepnumber,1) = methanogenesis_d;
    workingstate.erosion_tot(stepnumber,1) = erosion_tot ;
    workingstate.GAST(stepnumber,1) = GAST ;
    workingstate.iceline(stepnumber,1) = iceline ;
    workingstate.d34s_SO4_p(stepnumber,1) = d34s_SO4_p ;
    workingstate.d34s_SO4_di(stepnumber,1) = d34s_SO4_di ;
    workingstate.d34s_SO4_s(stepnumber,1) = d34s_SO4_s ;
    workingstate.d34s_SO4_h(stepnumber,1) = d34s_SO4_h ;
    workingstate.d34s_SO4_d(stepnumber,1) = d34s_SO4_d ;
    workingstate.d13c_G(stepnumber,1) = d13c_G ;
    workingstate.d13c_C(stepnumber,1) = d13c_C ;
    workingstate.d34s_PYR(stepnumber,1) = d34s_PYR ;
    workingstate.d34s_GYP(stepnumber,1) = d34s_GYP ;
    workingstate.time_myr(stepnumber,1) = t_geol ;
    workingstate.time(stepnumber,1) = t;
    workingstate.reductant_input(stepnumber,1) = reductant_input;

    %%%%%%%% print time
    if pars.runcontrol == -2
        workingstate.time_myr(stepnumber,1) = t.*1e-6 ;
    else
        workingstate.time_myr(stepnumber,1) = t_geol ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Record gridstates for single run   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sensanal == 0
    %%%% print a gridstate when each keytime threshold is crossed, or at model end
    next_stamp = pars.next_gridstamp ;
    if pars.finishgrid == 0
        if t_geol > next_stamp || t_geol == 0

            %%%% write gridstates
            gridstate.time_myr(pars.gridstamp_number,1) = next_stamp ;
            gridstate.land(:,:,pars.gridstamp_number) = land_past ;
            gridstate.Q(:,:,pars.gridstamp_number) = Q_past ;
            gridstate.Tair(:,:,pars.gridstamp_number) = Tair_past ;
            gridstate.TOPO(:,:,pars.gridstamp_number) = TOPO_past ;
            gridstate.CW(:,:,pars.gridstamp_number) = CW_per_km2_past ; %%% t/km2/yr
            gridstate.CWcarb(:,:,pars.gridstamp_number) = CWcarb_past ;
            gridstate.EPSILON(:,:,pars.gridstamp_number) = EPSILON_past * 1e6 ; %%% t/km2/yr

            %%%% set next boundary
            if t_geol < 0
                pars.gridstamp_number = pars.gridstamp_number + 1 ;
                pars.next_gridstamp = pars.runstamps(pars.gridstamp_number) ;
            else
                pars.finishgrid = 1 ;
            end

        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Record plotting states only in sensanal   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sensanal == 1
    workingstate.BAS_AREA(stepnumber,1) = BAS_AREA;
    workingstate.GRAN_AREA(stepnumber,1) = GRAN_AREA;
    workingstate.DEGASS(stepnumber,1) = DEGASS;
    workingstate.delta_mccb(stepnumber,1) = delta_mccb ;
    workingstate.d34s_S(stepnumber,1) = d34s_S ;
    workingstate.delta_OSr(stepnumber,1) = delta_OSr ;
    workingstate.SmM(stepnumber,1) = 28*S/pars.S0 ;
    workingstate.CO2ppm(stepnumber,1) = RCO2*280 ;
    workingstate.mrO2(stepnumber,1) = mrO2 ;
    workingstate.iceline(stepnumber,1) = iceline ;
    workingstate.T_gast(stepnumber,1) = TEMP_gast - 273 ;
    workingstate.ANOX(stepnumber,1) = ANOX ;
    workingstate.P(stepnumber,1) = P/pars.P0 ;
    workingstate.N(stepnumber,1) = N/pars.N0 ;
    workingstate.time_myr(stepnumber,1) = t_geol ;
    workingstate.time(stepnumber,1) = t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Final actions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% output timestep if specified
if sensanal == 0
    if pars.telltime ==1
        if mod(stepnumber,pars.display_resolution) == 0 
            %%%% print model state to screen
            fprintf('Model step: %d \t', stepnumber); fprintf('time: %d \t', t_geol) ; fprintf('next keyframe: %d \n', next_stamp)
        end
    end
end


%%%% record current model step
stepnumber = stepnumber + 1 ;


%%%% option to bail out if model is running aground
if stepnumber > pars.bailnumber
   terminateExecution
end



%%%% end function
end



