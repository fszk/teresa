# Teresa

## Step Alpha
Files are retrieved and organized in a systematic way by folders:
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2015_agosto/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2015_july/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2017_marzo_aprile/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2017_MISC/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2018_maggio/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2022/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2023/data/

glider_microstructure_data_conversion.m
Then .p files are converted to .mat files with the ODAS Toolbox from RSI.
The .mat files contain all the measurements converted from voltage to physically readable values.  

## Step Beta
teresa_microrider.ipynb
Files are listed:
liste_P.txt : list of .P 
liste_M.txt : list of .mat 
liste_F.txt : list of failed conversions 

: 657 x .mat files
total of  601 Go of .mat

: 922 x .P files
total of  99  Go of .P files

A metadata file (META.csv) is created from the scan of each .mat file, where each line contains
meta = [fn,fnshort,date,time,N_slow,N_fast,delta_t,
        P_min, P_max, P_mean,P_std,
        T2_min,T2_max,T2_mean,T2_std,
        T1_min,T1_max,T1_mean,T1_std,
        Tfast_min,Tfast_max,Tfast_mean,Tfast_std,
        W_min,W_max,W_mean,W_std,
        speed_min,speed_max,speed_mean,speed_std,
        sh1_min,sh1_max,sh1_mean,sh1_std,
        sh2_min,sh2_max,sh2_mean,sh2_std,
        V_Bat_min,V_Bat_max,V_Bat_mean,V_Bat_std,
        Incl_Y_min,Incl_Y_max,Incl_Y_mean,Incl_Y_std,
        Incl_X_min,Incl_X_max,Incl_X_mean,Incl_X_std,
        N_profiles_downward_fast,N_profiles_upward_fast,
        peaks_fast_i,peaks_fast_e,difpres_fast,
        N_profiles_downward_slow,N_profiles_upward_slow,
        peaks_slow_i,peaks_slow_e,difpres_slow,
        setupfilestr]


## Step Gamma
quick_epsilon.m
calculate_epsilon.m

The liste_M.txt is used to loop over each .mat file and produce values of epsilon, the turbulent kinetic energy dissipation rates (W/kg). 

General parameters are applied
%% 
% spikes
Ntresh  = 8 ;
fcut    = 0.5 ;
fs      = 512 ;
cm_cut  = 0.05 ; 
%
% Filter parameters
freq_HP_cut      = 1.5 ;  %   % default_HP_cut % frequency [Hz] of high-pass filter for shear probe data
freq_LP_cut      = 0.1 ;  %   % default_HP_cut % frequency [Hz] of low-pass  filter for shear probe data 
%
segment_fft = [] ;   % default will be 1m % Ndata1m_fast, around 1400 data
N_segment_fft = 2 ;
%N_window_fft = 2 ; % N x (N_segment_fft x segment_fft)
N_window_fft = 3 ;
Percent_overlap_fft = 50./100 ; % round(Percent_overlap_fft.* (diss_info.fft_length) )

Values are conserved through a matlab structure QUICK_EPSILON for each file.
PROJECT_DATAFILENAME_DATE_TIME.mat
e.g. TERESA_microrider_DAT_008_2018-04-23_11-17.mat

QUICK_EPSILON.profileNumber  = N ; %------------------ to be discarded, yoyo are now processed as a unique vector
QUICK_EPSILON.totalProfiles  = N_profiles ; %--------  to be discarded, yoyo are now processed as a unique vector
QUICK_EPSILON.date_profile = date_profile ;
QUICK_EPSILON.time_profile = time_profile ;
QUICK_EPSILON.lon_profile = lon_profile ;
QUICK_EPSILON.lat_profile = lat_profile ;
QUICK_EPSILON.yyyy = yyyy ;
QUICK_EPSILON.mm = mm ;
QUICK_EPSILON.dd = dd ;
QUICK_EPSILON.HH = HH ;
QUICK_EPSILON.MM = MM ;
QUICK_EPSILON.SS = SS;
QUICK_EPSILON.mSS = mSS;
QUICK_EPSILON.dn = dn;
QUICK_EPSILON.PROFILE_ID = PROFILE_ID;

QUICK_EPSILON.odas_version        = x.odas_version;
QUICK_EPSILON.vehicle_info        = x.vehicle_info;
QUICK_EPSILON.fullPath = fullPath ;
QUICK_EPSILON.setupfilestr = setupfilestr ;

QUICK_EPSILON.Ntresh  = Ntresh ;
QUICK_EPSILON.fcut    = fcut ;
QUICK_EPSILON.fs      = fs ;
QUICK_EPSILON.cm_cut  = cm_cut ;
QUICK_EPSILON.freq_HP_cut      = freq_HP_cut ;  
QUICK_EPSILON.freq_LP_cut      = freq_LP_cut ;  
QUICK_EPSILON.segment_fft = segment_fft; 
QUICK_EPSILON.N_segment_fft = N_segment_fft;
QUICK_EPSILON.N_window_fft = N_window_fft;
QUICK_EPSILON.Percent_overlap_fft = Percent_overlap_fft;
QUICK_EPSILON.fft_length = diss_info.fft_length ;
QUICK_EPSILON.diss_length = diss_info.diss_length ;
QUICK_EPSILON.overlap = diss_info.overlap ;

QUICK_EPSILON.PE = PE;
QUICK_EPSILON.FM = FM;
QUICK_EPSILON.DE = DE;
QUICK_EPSILON.ee = ee;
QUICK_EPSILON.e1 = e1;
QUICK_EPSILON.e2 = e2;
QUICK_EPSILON.fm1 = fm1;
QUICK_EPSILON.fm2 = fm2;
QUICK_EPSILON.ze = ze;
QUICK_EPSILON.te = te;
        
