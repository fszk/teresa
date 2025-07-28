%% PARAMETERS FOR FFT
% in function of vehicle_length (1.5 m), and mean speed 
% tau_fft = round(min([tau_to_avoidadv,tau_to_resolve_05cpm])) ;
%%
vehicle_length = 1.5 ; 

speed_mean = abs(round(mean(glider_speed(fu),'omitnan'),2)) ;  
w_mean = abs(round(mean(x.W_fast(fu)),2)) ; 
w_mean_speed = w_mean ; 

%%
speed_for_fft = speed_mean ; 
tau_to_avoidadv = round(vehicle_length / speed_for_fft) ; % tau to not exceed vehicle length

tau_to_resolve_05cpm = round((1/ 0.5) / speed_for_fft) ; 
tau_to_resolve_10cpm = round((1/ 1.0) / speed_for_fft) ; 
tau_to_resolve_20cpm = round((1/ 2.0) / speed_for_fft) ; 
%tau_fft = round(mean([tau_to_avoidadv,tau_to_resolve_05cpm])) ;
tau_fft = round(min([tau_to_avoidadv,tau_to_resolve_05cpm])) ;

%%
N_fft = round(tau_fft/(1/x.fs_fast)) ;       %---------------- METADATA
L_fft = tau_fft*speed_for_fft ;                     %---------------- METADATA
kl = 1/(L_fft) ;                        %---------------- METADATA
Fhp =  (1/tau_fft)/2  ;                %---------------- METADATA % Cutoff frequency (Hz)
Ntimes = 4 ;       
if ls / N_fft < Ntimes
    Ntimes = floor(ls / N_fft) ; 
end

%---------------- METADATA NL*diss_info.fft_length
Ntimes_N_fft = Ntimes .*  N_fft ;

percent_overlap = 0.5 ;                       %---------------- METADATA % overlap
%disp([lfft,k1,Fhp]) ;  

Niw = Ntimes ;
Now = 2*Niw -1  ;
Nv = 2 ;
RNvNow = 1 - 1.02 * Nv/Now;

var_ln_psi = (5/4)*(Now-Nv).^(-7/9);
sig_ln_psi = sqrt(var_ln_psi);




