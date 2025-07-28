%%


%% gradT_dis_spec

%[Xiv,Xi_ST,Xi_T,kB,eps_T,MAD_ST,MAD_T,MADc,LR,kL,kU,krange,kP,fit_flag_T]
% =
% gradT_dis_spec(pres,x0,k1,fn,kB_S,W,sL,sOV,Tdis,q,tau_0,time_corr,npoles,int_range,D,visco,T_dT,T_string,setupstr,plt,presplot,Tplot)

% GOAL
% Estimate of TKE and temperature variance dissipation rates eps and Xi by 
% integration of the measured temperature gradient wavenumber spectrum

% INPUT:
% pres: pressure vector (db)
% x0: grad temperature vector (°C/m)
% k1: minimium wavenumber for integration (small number, e.g., 0.1 cpm)
% fn: maximum frequency for calculations (90% of anti-aliasing filter f_AA. E.g., if f_AA=98 Hz, fn=0.9*98=88.2 Hz)
% kB_S: Batchelor wavenumber determined from shear probe used to caluclate Xi_ST (if 0 does not calculate. When two sh are available and accepted, the average of kB_S is used.)
% W: mean profiling speed (m/s)
% sL: length of segments for fft (scans)
% sOV: overlap for fft (scans)
% Tdis: type of theoretical spectrum (B: Batchelor, K: Kraichnan)
% q: turbulent parameter
% tau_0: nominal response time (s)
% time_corr: time correction approach: KOC, RSI, NAS, SOM (see below for details)
% npoles: transfer function for time response correction 'single' or 'double' pole pole 
% plt: if !=0 shows the spectra
% int_range: lower limit of the integration range according to Steinbuck et al 2009 ('S') or Luketina and Imberger 2001 ('L')
% D: thermal diffusivity (m^2/s)
% visco: kinematic viscosity of water (m^2/s)
% T_dT (*): raw pre-emphasized temperature signal
% T_string (*): name of the FP07 channel 'T1_dT1' or 'T2_dT2'
% setupstr (*): setupfilestr configuration file provided by ODAS libraries
% plt: flag for making the figure (0 no figure, ~=0 make figure)
% presplt: pressure vector for making the plot
% Tplt: temperature vector for making the plot

% (*) these entries are required for caluclating the noise spectrum according 
% to the function profided by RSI in the ODAS libraries.
% They are obtained reading the .p files with the ODAS libraries. If the
% ODAS libraries are not accessible by the user, an alternative for
% caluclating the noise spectrum is suggested below.

% OUTPUT
% Xiv: Xi from spectral integration in the well resolved part (°C^2/s)
% Xi_ST: Xi corrected with kB_S obtained from shear probes (°C^2/s)
% Xi_T: Xi after MLE spectral fitting  (°C^2/s)
% kB: Batchelor wavenumber after MLE spectral fitting (cpm) 
% eps_T: TKE dissipation rate after MLE spectral fitting (m^2/s^3)
% MAD_ST: Mean Absolute Deviation between observed and empirical spectra using kB_S from shear
% MAD_T: Mean Absolute Deviation between observed and empirical spectra using MLE fitting
% MADc: Threshold for the Mean Absolute Deviation between observed and empirical spectra
% LR: likelihood ratio
% kL: Lower integration wavenumber (cpm)
% kU: Upper integration wavenumber (cpm)
% krange: wavenumber range used for spectral integration 
% kP: wavenumber corresponding to the peak of the fitted theoretical spectrum (cpm)
% fit_flag_ST: acceptance flag of Xi_ST according to quality metrics: 0 rejected, 1 accepted
% fit_flag_T: acceptance flag of eps_T and Xi_T according to quality metrics: 0 rejected, 1 accepted

dn_section = dn(fu)    ;
pres_in = x.P_fast(fu) ;                                    % pres: pressure vector (db)

% x0: grad temperature vector (°C/m)
% T_dT (*): raw pre-emphasized temperature signal





k1 = kl ;                       % 0.1 ; k1: minimium wavenumber for integration (small number, e.g., 0.1 cpm)
f_AA_max = 0.9 * diss.f_AA ;    % fn: maximum frequency for calculations (90% of anti-aliasing filter f_AA. E.g., if f_AA=98 Hz, fn=0.9*98=88.2 Hz)

TD = mean(TD_fast(fu),'omitnan');   % D: thermal diffusivity (m^2/s)
NU = mean(NU_fast(fu),'omitnan');   % visco: kinematic viscosity of water (m^2/s)

kB_S = 0 ;                      % kB_S: Batchelor wavenumber determined from shear probe used to caluclate Xi_ST 
                                % (if 0 does not calculate. When two sh are available and accepted, the average of kB_S is used.)
                                % kB = (eps / (nu * DT^2))^(1/4)
                                % nu Kinematic viscosity of the fluid
                                % DT Thermal diffusivity 

W = speed_mean ;                % or w_mean % W: mean profiling speed (m/s)

sL  = diss.fft_length  ;        % sL: length of segments for fft (scans)
sW  = diss.diss_length/2 ;        % sL: length of segments for fft (scans)
sOV = diss.overlap ;            % sOV: overlap for fft (scans)
sOV = sL/2 ; 

Tdis = 'K' ;                    % Tdis: type of theoretical spectrum (B: Batchelor, K: Kraichnan)
Tdis_logic = 2 ; 

q =  5.26  ;                    % q: turbulent parameter
tau_0 = 7e-3 ;                  % tau_0: nominal response time (s)
time_corr = 'RSI' ;             % time_corr: time correction approach: KOC, RSI, NAS, SOM (see below for details)
time_corr_logic = 2 ;
npoles = 'single' ;  % npoles: transfer function for time response correction 'single' or 'double' pole pole 
npoles_logic = 1 ; 

int_range = 'S'  ;              % int_range: lower limit of the integration range according to Steinbuck et al 2009 ('S') or Luketina and Imberger 2001 ('L')
int_range_logic = 1 ; 

setupstr = x.setupfilestr ;     % setupstr (*): setupfilestr configuration file provided by ODAS libraries
pltYN = 0 ;                     % plt: flag for making the figure (0 no figure, ~=0 make figure)
plt = 0 ;                       % plt: if !=0 shows the spectra
P_for_plot = diss.P ;            % presplt: pressure vector for making the plot
T_for_plot = diss.T ;         % Tplt: temperature vector for making the plot




segments = 1:sW:length(pres_in);
%
eT   = [] ;
xiv  = [] ;
xit  = [] ;
xist = [] ;
%
P_eT   = []   ;
dn_eT   = []   ;

likehratio = []   ;
qc_flag_T  = []   ;
qc_mad_T   = []   ;
qc_mad_ST  = []   ;
qc_mad_c  = []   ;

kB_T  = []   ;
kL_T  = []   ;
kU_T  = []   ;
kP_T  = []   ;
krange_T = [];
%
for s = 1:(length(segments)-1)


    ss = segments(s):segments(s+1);

try 
       [Xiv,Xi_ST,Xi_T,kB,eps_T,MAD_ST,MAD_T,MAD_c,LR,kL,kU,krange,kP,fit_flag_T] = gradT_dis_spec( ...
 pres_in(ss),x0(ss),k1,f_AA_max,kB_S,W,sL,sOV,Tdis,q,tau_0,time_corr,npoles,int_range,TD,NU,T_dT(ss), ...
 T_string,setupstr,pltYN,P_for_plot,T_for_plot);


eT   = [eT eps_T] ; 
xiv = [xiv Xiv] ;
xit = [xit Xi_T] ;
xist = [xit Xi_ST] ;

kB_T  = [kB_T kB]   ;
kL_T  = [kL_T kL]   ;
kU_T  = [kU_T kU]   ;
kP_T  = [kP_T kP]   ;
krange_T  = [krange_T krange]   ;

%
P_eT = [P_eT mean(pres_in(ss))] ;
dn_eT = [dn_eT mean(dn_section(ss))] ;
likehratio = [likehratio LR] ;
qc_flag_T  = [qc_flag_T fit_flag_T];
%
qc_mad_T  = [qc_mad_T MAD_T];
qc_mad_ST = [qc_mad_ST MAD_ST];
qc_mad_c  = [qc_mad_c MAD_c];

catch ME
  
    disp(ME.message)

end






end
%

%clc



dissT = struct ;

dissT.k1 = kl ;                       
dissT.f_AA_max = f_AA_max ;  
dissT.TD = TD ; 
dissT.NU = NU ; 
dissT.kB_S = kB_S ; 
dissT.W = W ;              
dissT.sL  = sL ; 
dissT.sW  = sW ; 
dissT.sOV = sOV ; 
dissT.Tdis = Tdis ;                    
dissT.q =  q  ;                   
dissT.tau_0 =tau_0 ;                 
dissT.time_corr = time_corr ;             
dissT.npoles = npoles ;             
dissT.int_range = int_range ;         



dissT.eT   = eT ;
dissT.xiv  = xiv ;
dissT.xit  = xit ;
dissT.xist = xist ;
%
dissT.P_eT   = P_eT   ;
dissT.dn_eT   = dn_eT   ;

dissT.likehratio = likehratio   ;
dissT.qc_flag_T  = qc_flag_T  ;
dissT.qc_mad_T   = qc_mad_T   ;
dissT.qc_mad_ST  = qc_mad_ST   ;
dissT.qc_mad_c  = qc_mad_c  ;

dissT.kB_T  = kB_T  ;
dissT.kL_T  = kL_T  ;
dissT.kU_T  = kU_T    ;
dissT.kP_T  = kP_T   ;
dissT.krange_T = krange_T;

%mean(diff(dissT.P_eT))