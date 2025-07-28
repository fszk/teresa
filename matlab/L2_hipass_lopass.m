%% HP


sh1 = x.sh1(fu) ;
sh2 = x.sh2(fu) ; 



Fs = x.fs_fast    ;  % Sampling frequency (Hz)
Fc = 0.1          ;  % Cutoff frequency (Hz) 
order = 1         ;  % Filter order          
fc_hi = Fc ; 
fc_hi_order = order ; 

% Design Butterworth filter
[b, a] = butter(order, Fc / (Fs / 2), 'high');
sh1hp = filtfilt(b, a, sh1) ; 
sh2hp = filtfilt(b, a, sh2) ; 
sh1hpa = abs(sh1hp) ; 
sh2hpa = abs(sh2hp) ; 


sh1hpa = sh1 ; 
sh2hpa = sh2 ; 


Fs = x.fs_fast    ;  
Fc = 1            ;  %------------------- FEAT - QC
order = 1         ;  %------------------- FEAT - QC
fc_lo = Fc ;
fc_lo_order = order ; 





% Design Butterworth filter
[b, a] = butter(order, Fc / (Fs / 2), 'low');
sh1hpalp = filtfilt(b, a, sh1hpa) ; 
sh2hpalp = filtfilt(b, a, sh2hpa) ; 