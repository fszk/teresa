%% get_diss
% Input Structure
% get_diss_odas_FK25
%
diss_info               = struct ;
diss_info.fft_length    = N_fft  ;
diss_info.Ntimes        = Ntimes ;  
diss_info.diss_length   = Ntimes_N_fft  ;
diss_info.overlap       = round(percent_overlap*diss_info.diss_length) ;
diss_info.fs_fast       = x.fs_fast ;
diss_info.fs_slow       = x.fs_slow ;
diss_info.t_fast        = x.t_fast(fu) ;

diss_info.T             = temperature_for_dissipation(fu) ; 
diss_info.P             = x.P_fast(fu) ;
diss_info.A             = [x.Ax(fu) x.Ay(fu)]    ;
%diss_info.goodman       = logical(1)    ; %% by default : 'true'
diss_info.speed         = glider_speed(fu) ;

diss_info.AOA = glider_aoa(fu) ;

%diss_info.fit_2_isr =
%* [fit_2_isr] Value of dissipation rate at which the estimation
%       will be derived from a fit to the inertial subrange rather than by
%       spectral integration. Default fit_2_isr = 1.5e-5 W/kg.



%% DISS

try 
   diss = get_diss_odas_FK25([sh1hpa_dsp_hp sh2hpa_dsp_hp],[x.Ax(fu) x.Ay(fu)],diss_info) ; 
   diss_function_name = 'get_diss_odas_FK25()' ; 

   %diss = get_diss_odas([sh1hpa_dsp_hp sh2hpa_dsp_hp],[x.Ax(fu) x.Ay(fu)],diss_info) ; 
   % diss_function_name = 'get_diss_odas()' ; 
   get_diss_odas_worked = 1 ; 
catch ME
    get_diss_odas_worked = 0 ; 
    disp(ME.message)
end






