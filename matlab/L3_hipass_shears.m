%% Hi pass shears


Fs = x.fs_fast    ;  % Sampling frequency (Hz)
Fc = Fhp          ;
order = 1         ;  % Filter order
% Design Butterworth filter
[b, a] = butter(order, Fc / (Fs / 2), 'high');
sh1hpa_dsp_hp = filtfilt(b, a, sh1hpa_dsp) ; 
sh2hpa_dsp_hp = filtfilt(b, a, sh2hpa_dsp) ;

fs_sampling = Fs ; 
fc_hi_shear = Fc ; 
fc_hi_shear_order = order ; 


