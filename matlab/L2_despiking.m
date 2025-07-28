


%% DESPIKING


thresh  = 8 ;  %---------------- METADATA
fcut   = 0.5 ; %---------------- METADATA
N = round(0.04*Fs) ; %---------------- METADATA

despike_thresh = thresh ;
despike_fcut = fcut ; 
despike_N = N ; 



%
[sh1hpa_dsp, spike, pass_count, ratio] = despike(sh1hpa, thresh, fcut, Fs, N);
% fraction of data affected ?
% number of iterations ?
despike_sh1_passcount = pass_count;  %---------------- METADATA
despike_sh1_ratio = ratio ;      %---------------- METADATA
despike_sh1_spike = spike ;      %---------------- METADATA
sh1hpa_dsp ;

Pfast = x.P_fast(fu) ; 
P_spikes_sh1 = Pfast(despike_sh1_spike) ; 

%
[sh2hpa_dsp, spike, pass_count, ratio] = despike(sh2hpa, thresh, fcut, Fs, N);
% fraction of data affected ?
% number of iterations ?
despike_sh2_passcount = pass_count;  %---------------- METADATA
despike_sh2_ratio = ratio ;      %---------------- METADATA
despike_sh2_spike = spike ;      %---------------- METADATA
sh2hpa_dsp ;


Pfast = x.P_fast(fu) ; 
P_spikes_sh2 = Pfast(despike_sh2_spike) ; 

%
%figure; 
%plot(dt(fu),sh1hpa)
%hold on
%plot(dt(fu),sh1hpa_dsp)
%
%figure; 
%plot(dt(fu),sh2hpa)
%hold on
%plot(dt(fu),sh2hpa_dsp)