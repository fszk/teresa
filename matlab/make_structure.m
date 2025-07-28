
clear mr L2 L3 L4 QC

%%
meta = struct ; 
meta.datetime_start = dt(fu(1)) ; 
meta.datenum_start = dn(fu(1)) ;
meta.datetime_end = dt(fu(end)) ; 
meta.datenum_end = dn(fu(end)) ;

meta.duration_minutes = minutes(dt(fu(end))-dt(fu(1))) ; 
meta.duration_hours = hours(dt(fu(end))-dt(fu(1))) ; 

%meta.dt = dt(fu) ; 
%meta.dn = dn(fu) ; 
%meta.lon = (section_lon) ; 
%meta.lat = (section_lat) ;

meta.dt = dt_eps ; 
meta.dn = dn_eps ; 
meta.lon = lon_eps ; 
meta.lat = lat_eps ; 

meta.navig_dir = navig_dir ;
meta.navig_dir_str = navig_dir_str ;

meta.speed_mean = speed_mean ; 
meta.w_mean = w_mean ; 

meta.lon_start = section_lon(1) ; 
meta.lat_start = section_lat(1) ;
meta.lon_end = section_lon(end) ; 
meta.lat_end = section_lat(end) ;
meta.distance = m_lldist([meta.lon_start,meta.lon_end],[meta.lat_start,meta.lat_end]) ;
meta.lon_min = min(section_lon ,[],'omitnan') ;
meta.lat_min = min(section_lat ,[],'omitnan') ;
meta.lon_max = max(section_lon ,[],'omitnan') ; 
meta.lat_max = max(section_lat,[],'omitnan') ;

meta.lon_mean = mean(section_lon,'omitnan') ; 
meta.lat_mean = mean(section_lat,'omitnan') ;

meta.number_of_sections = lu ; 
meta.section_number = u ; 
meta.section_direct = dirfu ;
meta.section_direct_str = dirfu_str ;
meta.pres_start_end = [x.P_fast(fu(1)),x.P_fast(fu(end))] ; 

meta.pres_min = min(x.P_fast(fu),[],'omitnan') ;
meta.pres_max = max(x.P_fast(fu),[],'omitnan') ;
meta.pres_delta = pres_delta ;

meta.section_indexes = [fu1,fue] ; 
meta.length_section = ls ; 
meta.length_origina = lx ; 
meta.percen_origina = round(100*ls/lx,3) ; 


meta.fullPath = x.fullPath ; 
meta.setupfilestr = x.setupfilestr ; 
meta.odas_version = x.odas_version ; 

meta.profile_dir = x.vehicle_info.profile_dir ; 
meta.speed_algorithm = x.vehicle_info.speed_algorithm ; 
meta.tau = x.vehicle_info.tau ; 
meta.aoa = x.vehicle_info.aoa ; 

meta.speed_cutout = x.params.speed_cutout ; 
meta.speed_tau = x.params.speed_tau ; 
meta.vehicle = x.params.vehicle ; 

meta.thermistor_source = x.params.thermistor_source ; 
meta.constant_temp = x.params.constant_temp ; 
meta.time_offset = x.params.time_offset ; 


%meta.temperature_for_dissipation = temperature_for_dissipation

%meta.speed_source = x.params.speed_source ;
speed_source = x.params.speed_source  ; 
meta.speed_source = x.params.speed_source ;

meta.direction_section = dirfu ; %1 down -1 up

% %%
% meta_simpler = struct ;
% meta_simpler.dt = dt_eps ; 
% meta_simpler.dn = dn_eps ; 
% meta_simpler.lon = lon_eps ; 
% meta_simpler.lat = lat_eps ; 
% %meta_simpler.lat = meta.direction_section  ; 


%%

L2 = struct ; 

L2.sh1 = sh1;  
L2.sh2 = sh2; 

L2.sh1hpa = sh1hpa;  
L2.sh2hpa = sh2hpa;  

L2.sh1hpa_dsp = sh1hpa_dsp ; 
L2.sh2hpa_dsp = sh2hpa_dsp ;

L2.sh1hpa_dsp_hp = sh1hpa_dsp_hp ; 
L2.sh2hpa_dsp_hp = sh2hpa_dsp_hp ; 

%%
L2_params = struct ; 

L2_params.u = u ; 
L2_params.section_id = section_id ; 
L2_params.Pmin = Pmin ; 
L2_params.Wmin = Wmin ;  
L2_params.minDuration = minDuration ;  
%L2_params.speed_method = speed_method ;  
%L2.minSpeed = minSpeed ; 
%L2.minPitch = minPitch ; 
L2_params.Ndown_Nup = [N_down,N_upwd] ;  
L2_params.lu = lu ;  

%%

L2_params.fc_hi = fc_hi ;  
L2_params.fc_hi_order = order ; 

L2_params.fc_lo = fc_lo ;   
L2_params.fc_lo_order = fc_lo_order ;  

%%

L2_params.despike_thresh = despike_thresh ; 
L2_params.despike_fcut = despike_fcut ; 
L2_params.despike_N = despike_N ; 

L2_params.despike_sh1_passcount = despike_sh1_passcount ; 
L2_params.despike_sh1_ratio = despike_sh1_ratio ; 
L2_params.despike_sh1_spike = despike_sh1_spike ; 
L2_params.P_spikes_sh1 = P_spikes_sh1 ; 

L2_params.despike_sh2_passcount = despike_sh2_passcount ; 
L2_params.despike_sh2_ratio = despike_sh2_ratio ; 
L2_params.despike_sh2_spike = despike_sh2_spike ; 
L2_params.P_spikes_sh2 = P_spikes_sh2 ; 

%%

L2_params.vehicle_length = vehicle_length ; 

L2_params.speed_mean = speed_mean ; 
L2_params.w_mean = w_mean ; 
L2_params.speed_for_fft = speed_for_fft ; 

L2_params.tau_to_avoidadv = tau_to_avoidadv ; 
L2_params.tau_to_resolve_05cpm = tau_to_resolve_05cpm ; 
L2_params.tau_to_resolve_10cpm = tau_to_resolve_10cpm ; 
L2_params.tau_to_resolve_20cpm = tau_to_resolve_20cpm ;

L2_params.tau_fft = tau_fft ; 
L2_params.N_fft = N_fft ; 
L2_params.L_fft = L_fft ; 
L2_params.kl = kl ;
L2_params.Ntimes = Ntimes ;
L2_params.Ntimes_N_fft = Ntimes_N_fft ;
L2_params.percent_overlap = percent_overlap ;
L2_params.Fhp = Fhp ;


%%
L2_params.c1 = c1 ; 
L2_params.c2 = c2 ; 

L2_params.Th_source_logic = Th_source_logic ; 
L2_params.Th_source_string = Th_source_string ; 

L2_params.speed_source_logic = speed_source_logic ;
L2_params.speed_source = x.params.speed_source ; 

L2_params.T_source_logic = T_source_logic ;
L2_params.T_source_string = T_source_string ;


%%

L2_params.fs_sampling = fs_sampling ; 
L2_params.fc_hi_shear = fc_hi_shear ; 
L2_params.fc_hi_shear_order = fc_hi_shear_order ; 


%%

L3_params = struct ;  

%L3_params.diss_info = diss_info  ;  TOO HEAVY with speed ....

L3_params.c1 = L2_params.c1 ;
L3_params.c2 = L2_params.c2 ;
L3_params.Th_source_logic = L2_params.Th_source_logic ;
L3_params.Th_source_string =  L2_params.Th_source_string ;
L3_params.T_source_logic = L2_params.T_source_logic ;
L3_params.T_source_string = L2_params.T_source_string ;


L3_params.speed_source_logic = L2_params.speed_source_logic ;
L3_params.speed_source =  L2_params.speed_source ;

L3_params.fft_length = diss_info.fft_length ; 
L3_params.Ntimes = diss_info.Ntimes ; 
L3_params.diss_length = diss_info.diss_length ; 
L3_params.overlap = diss_info.overlap ; 
L3_params.fs_fast = diss_info.fs_fast ; 
L3_params.fs_slow = diss_info.fs_slow ; 
%L3_params.goodman = diss_info.goodman ; 





L3_params.Niw = Niw ;
L3_params.Now = Now ;
L3_params.Nv = Nv ;
L3_params.RNvNow = RNvNow;

L3_params.var_ln_psi = var_ln_psi;
L3_params.sig_ln_psi = sig_ln_psi ;


L3_params.gradT_dis_spec = 'Piccolroaz et al. 2020, Lake Garda' ;


L3_params.k1 = k1 ;
L3_params.f_AA_max = f_AA_max ;
L3_params.kB_S = kB_S ;            
L3_params.W = W ;               
L3_params.sL  = sL ;   
L3_params.sW  = sW ;   
L3_params.sOV = sOV   ;     
L3_params.Tdis = Tdis; 
L3_params.Tdis_logic = Tdis_logic;   
L3_params.q =  q ;                   
L3_params.tau_0 = tau_0 ;    
L3_params.time_corr = time_corr ;      
L3_params.time_corr_logic = time_corr_logic ;      
L3_params.npoles = npoles ; 
L3_params.npoles_logic = npoles_logic ;  
L3_params.int_range =int_range  ;  
L3_params.int_range_logic = int_range_logic ;



%%
L3 = struct ;  

L3 = diss_info ; 

L3.Nasmyth_spec1 = Nasmyth_spec1 ; 
L3.sh_1 = sh_1 ; 
L3.sh_clean1 = sh_clean1 ; 

L3.Nasmyth_spec2 = Nasmyth_spec2 ; 
L3.sh_2 = sh_2 ; 
L3.sh_clean2 = sh_clean2 ; 


%%



L4 = struct ; 
L4.diss_function_name = diss_function_name ; 

L4.speed = speed_e; 
L4.speed_std = speed_std_e; 

L4.speed_source_logic     = speed_source_logic ; 
L4.speed_source   = x.params.speed_source ; 
L4.AOA  = AOA ; 

L4.f_AA = diss.f_AA ; 
L4.f_limit = diss.f_limit ;
L4.dof_spec = diss.dof_spec ; 

L4.Ne = Ne ;  
L4.dz = dz ; 
L4.dt_eps = dt_eps ;  
L4.dn_eps = dn_eps ; 
L4.lon_eps = lon_eps ;  
L4.lat_eps = lat_eps ;  
L4.z = Z_e ; 
L4.P = P_e ;  
L4.T = T_e ;  

L4.C = C_e ;  
L4.SP = SP_e ;  
L4.SA = SA_e ;  
L4.rho = rho_e ;  

L4.nu = diss.nu ;  
L4.T1 = T1_e ;  
L4.T2 = T2_e ;  

L4.pitch_gl = pitch_e ;
L4.roll_gl = roll_e ;

L4.pitch = pitch ;

L4.roll = roll ;


L4.method_sh1 = method_sh1 ;  
L4.dof_e1 = dof_e1 ;  
L4.e1 = e1 ;  
L4.fm1 = fm1 ; 
L4.fm_old1 = fm_old1 ; 
L4.mad1 = mad1 ; 
L4.mad10_1 = mad10_1 ;  

L4.kmax1 = kmax1 ; 


L4.method_sh2 = method_sh2 ;  
L4.dof_e2 = dof_e2 ;  
L4.e2 = e2 ; 
L4.fm2 = fm2 ; 
L4.fm_old2 = fm_old2 ; 
L4.mad2 = mad2 ; 
L4.mad10_2 = mad10_2 ;  

L4.kmax2 = kmax2 ; 



%---------------------- vsr
L4.spectral_model_VSR  = 'VSR Lueck et al. 2024' ;
L4.method_sh1_vsr = method_sh1_vsr ; 
L4.method_sh2_vsr = method_sh2_vsr ; 

L4.le = le ;
L4.sig_ln_ee = sig_ln_ee ;

L4.kolmo1 = kolmo1 ;
L4.k_hat_u1 = k_hat_u1 ;
L4.Vf1 = Vf1 ;
L4.Lf1 = Lf1 ;
L4.var_ln_e1 = var_ln_e1 ;
L4.sig_ln_e1 = sig_ln_e1 ;
L4.upp_e1_vsr = upp_e1_vsr ;
L4.low_e1_vsr = low_e1_vsr ;

L4.kolmo2 = kolmo2 ;
L4.k_hat_u2 = k_hat_u2 ;
L4.Vf2 = Vf2 ;
L4.Lf2 = Lf2 ;
L4.var_ln_e2 = var_ln_e2 ;
L4.sig_ln_e2 = sig_ln_e2 ;
L4.upp_e2_vsr = upp_e2_vsr ;
L4.low_e2_vsr = low_e2_vsr ;

%---------------------- isr
L4.spectral_model_ISR = 'ISR Lueck et al. 2024' ;
L4.method_sh1_isr = method_sh1_isr ; 
L4.method_sh2_isr = method_sh2_isr ; 
L4.sig_ln_psi = sig_ln_psi ; 
L4.sig_ln_psi_on_sqrt_Ns = sig_ln_psi_on_sqrt_Ns ;
L4.Ns1 = Ns1 ;
L4.TM1 = TM1 ;
L4.sig_ln_psi_on_sqrt_Ns1 = sig_ln_psi_on_sqrt_Ns1 ;
L4.upp_e1_isr = upp_e1_isr ;
L4.low_e1_isr = low_e1_isr ;


L4.Ns2 = Ns2 ;
L4.TM2 = TM2 ;
L4.sig_ln_psi_on_sqrt_Ns2 = sig_ln_psi_on_sqrt_Ns2 ;
L4.upp_e2_isr = upp_e2_isr ;
L4.low_e2_isr = low_e2_isr ;

%%

L4.gradT_dis_spec = 'FP07 Piccolroaz et al. 2020, Lake Garda' ; 
L4.thermistor_source = thermistor_source ; 
L4.spectral_model_THERM = L3_params.Tdis; 

L4.therm1_k1 = dissT1.k1  ;
L4.therm1_f_AA_max = dissT1.f_AA_max ;
L4.therm1_TD = dissT1.TD;
L4.therm1_NU = dissT1.NU;
L4.therm1_kB_S = dissT1.kB_S;
L4.therm1_W = dissT1.W;
L4.therm1_sL = dissT1.sL;
L4.therm1_sW = dissT1.sW;
L4.therm1_sOV = dissT1.sOV;
%L4.therm1_Tdis= dissT1.Tdis;

L4.therm1_q = dissT1.q;
L4.therm1_tau_0= dissT1.tau_0;
%L4.therm1_time_corr= dissT1.time_corr;
%L4.therm1_npoles= dissT1.npoles;
L4.therm1_int_range= dissT1.int_range;




L4.NeT1  = NeT1  ; 
L4.P_eT1 = P_eT1 ; 
L4.D_eT1 = D_eT1 ; 
L4.Z_eT1 = Z_eT1 ; 
L4.dzT1  = dzT1  ;
L4.dt_eT1 = dt_eT1 ; 
L4.dn_eT1 = dn_eT1 ; 
L4.t_eT1  = t_eT1  ; 
L4.lon_eT1 = lon_eT1 ;   
L4.lat_eT1 = lat_eT1 ;  
L4.T_eT1   = T_eT1   ; 
L4.nu_eT1  = nu_eT1  ; 
L4.td_eT1  = td_eT1  ; 

L4.eT1   = dissT1.eT   ; 
L4.xiv1  = dissT1.xiv  ; 
L4.xit1  = dissT1.xit  ; 
L4.xist1 = dissT1.xist ; 

L4.likehratio1 = dissT1.likehratio ; 
L4.qc_flag_T1  = dissT1.qc_flag_T  ; 
L4.qc_mad_T1  = dissT1.qc_mad_T  ; 
L4.qc_mad_ST1 = dissT1.qc_mad_ST ; 
L4.qc_mad_c1 = dissT1.qc_mad_c ; 


L4.kB_T1  = dissT1.kB_T ;
L4.kL_T1  = dissT1.kL_T ;
L4.kU_T1  = dissT1.kU_T ;
L4.kP_T1  = dissT1.kP_T ;
L4.krange_T1 = dissT1.krange_T ;


%----------------------
L4.therm2_k1 = dissT2.k1  ;
L4.therm2_f_AA_max = dissT2.f_AA_max ;
L4.therm2_TD = dissT2.TD;
L4.therm2_NU = dissT2.NU;
L4.therm2_kB_S = dissT2.kB_S;
L4.therm2_W = dissT2.W;
L4.therm2_sL = dissT2.sL;
L4.therm2_sW = dissT2.sW;
L4.therm2_sOV = dissT2.sOV;
L4.therm2_Tdis= dissT2.Tdis;
%L4.therm2_Tdis_logic = dissT2.Tdis_logic;
L4.therm2_q = dissT2.q;
L4.therm2_tau_0= dissT2.tau_0;
%L4.therm2_time_corr= dissT2.time_corr;
%L4.therm2_time_corr_logic = dissT2.time_corr_logic ;
%L4.therm2_npoles= dissT2.npoles;
L4.therm2_int_range= dissT2.int_range;


L4.NeT2  = NeT2  ; 
L4.P_eT2 = P_eT2 ; 
L4.D_eT2 = D_eT2 ; 
L4.Z_eT2 = Z_eT2 ; 
L4.dzT2 = dzT2  ;
L4.dt_eT2 = dt_eT2 ; 
L4.dn_eT2 = dn_eT2 ; 
L4.t_eT2 = t_eT2  ; 
L4.lon_eT2 = lon_eT2 ;   
L4.lat_eT2 = lat_eT2 ;  
L4.T_eT2   = T_eT2   ; 
L4.nu_eT2  = nu_eT2  ; 
L4.td_eT2  = td_eT2  ; 

L4.eT2   = dissT2.eT   ; 
L4.xiv2  = dissT2.xiv  ; 
L4.xit2  = dissT2.xit  ; 
L4.xist2 = dissT2.xist ; 

L4.likehratio2 = dissT2.likehratio ; 
L4.qc_flag_T2  = dissT2.qc_flag_T  ; 
L4.qc_mad_T2  = dissT2.qc_mad_T  ; 
L4.qc_mad_ST2 = dissT2.qc_mad_ST ; 
L4.qc_mad_c2 = dissT2.qc_mad_c ; 


L4.kB_T2  = dissT2.kB_T ;
L4.kL_T2  = dissT2.kL_T ;
L4.kU_T2  = dissT2.kU_T ;
L4.kP_T2  = dissT2.kP_T ;
L4.krange_T2 = dissT2.krange_T ;






%%
QC = struct ; 



%%

QC.despike_thresh = despike_thresh ; 
QC.despike_fcut = despike_fcut ; 
QC.despike_N = despike_N ; 

QC.despike_sh1_passcount = despike_sh1_passcount ; 
QC.despike_sh1_ratio = despike_sh1_ratio ; 
QC.dspk_frac2 = dspk_frac2 ; 
QC.despike_sh2_passcount = despike_sh2_passcount ; 
QC.despike_sh2_ratio = despike_sh2_ratio ; 
QC.dspk_frac1 = dspk_frac1 ; 

QC.figofmer_fail = figofmer_fail ;

QC.spikfrac_fail = spikfrac_fail ;
QC.spikfrac_failr = spikfrac_failr ;
QC.spikpass_fail = spikpass_fail ;

QC.epsratio_fail_vsr = epsratio_fail_vsr ;
QC.epsratio_fail_isr = epsratio_fail_isr ;
QC.varresol_fail = varresol_fail  ;

QC.aoa_thresh1 = aoa_thresh1 ; 
QC.aoa_thresh2 = aoa_thresh2 ; 



%%

QC.Q_figofmer1 = Q_figofmer1 ; 
QC.Q_spikfrac1 = Q_spikfrac1 ; 
QC.Q_spikpass1 = Q_spikpass1 ; 
QC.Q_varresol1 = Q_varresol1 ; 
QC.Q_e1_eratio = Q_e1_eratio ; 
QC.Q_aoa = Q_aoa ; 

QC.Q1 = Q1 ; 

QC.fok1 = fok1 ; 
QC.lok1 = lok1 ; 
QC.pok1 = pok1 ; 


QC.Q_figofmer2 = Q_figofmer2 ; 
QC.Q_spikfrac2 = Q_spikfrac2 ; 
QC.Q_spikpass2 = Q_spikpass2 ; 
QC.Q_varresol2 = Q_varresol2 ; 
QC.Q_e2_eratio = Q_e2_eratio ; 

QC.Q2 = Q2 ; 
QC.fok2 = fok2 ; 
QC.lok2 = lok2 ; 
QC.pok2 = pok2 ; 


QC.qc_flag_T1  = dissT1.qc_flag_T  ; 
QC.likehratio1 = dissT1.likehratio ; 
QC.qc_mad_T1  = dissT1.qc_mad_T  ; 
QC.qc_mad_ST1 = dissT1.qc_mad_ST ; 
QC.qc_mad_c1 = dissT1.qc_mad_c ; 

QC.qc_flag_T2  = dissT2.qc_flag_T  ; 
QC.likehratio2 = dissT2.likehratio ; 
QC.qc_mad_T2  = dissT2.qc_mad_T  ; 
QC.qc_mad_ST2 = dissT2.qc_mad_ST ; 
QC.qc_mad_c2 = dissT2.qc_mad_c ; 

QC.fvsr = fvsr ; 
QC.ee_vsr = ee_vsr ; 
QC.ok_vsr = ok_vsr ; 
QC.ee_vsr_qc1 = ee_vsr_qc1 ;  
QC.ok_vsr_qc1 = ok_vsr_qc1 ;  


QC.fisr = fisr ; 
QC.ee_isr = ee_isr ;  
QC.ok_isr = ok_isr ;  
QC.ee_isr_qc1 = ee_isr_qc1 ;  
QC.ok_isr_qc1 = ok_isr_qc1 ;  

QC.Q_T1 = Q_T1 ; 
QC.Q_T2 = Q_T2 ;  
QC.Q_T1T2_SPREAD = Q_T1T2_SPREAD ;  
QC.QC_FP07_1 = QC_FP07_1 ;  
QC.QC_FP07_2 = QC_FP07_2 ;  


% 
% QC.ee_fp07 = ee_fp07 ; 
% QC.xi_fp07 = xi_fp07 ;  
% QC.ok_fp07 = ok_fp07 ;  
% 
% QC.ee_fp07_i = ee_fp07_i;
% QC.xi_fp07_i = xi_fp07_i;
%    



%%
MR              = struct  ; 
MR.meta         = meta    ;
MR.L2_params    = L2_params ;
MR.L2           = L2      ;
MR.L3_params    = L3_params   ;
MR.L3           = L3   ;
MR.L4           = L4   ;
MR.QC           = QC   ;
%mr.glider       = glider  ;

MR_Mb = sprintf('%04d',round(0.7*whos('MR').bytes / (1024 * 1024)) );
disp(MR_Mb)


%%

mr           = struct  ; 
mr.meta      = meta    ;
mr.L2_params = L2_params ; 
mr.L3_params = L3_params   ;
mr.L4        = L4   ;
mr.QC        = QC   ;
%mr_qc.glider    = glider  ;

mr_Mb = sprintf('%04d',round(0.7*whos('mr').bytes / (1024 * 1024)) );
disp(mr_Mb)



