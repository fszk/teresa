# Teresa

## Step Alpha
## organizing and sorting files / data conversion 
Files are retrieved and organized in a systematic way by folders:<br>
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2015_agosto/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2015_july/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2017_marzo_aprile/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2017_MISC/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2018_maggio/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2022/data/
/Volumes/DENISE/data/glider/teresa/data/teresa_microrider_2023/data/
<br><br>
glider_microstructure_data_conversion.m<br>
Then .p files are converted to .mat files with the ODAS Toolbox from RSI.<br>
The .mat files contain all the measurements converted from voltage to physically readable values.  

## Step Beta
## metadata
teresa_microrider.ipynb<br><br>
Files are listed:<br>
liste_P.txt : list of .P <br>
liste_M.txt : list of .mat <br>
liste_F.txt : list of failed conversions <br>
<br><br>
: 657 x .mat files<br>
total of  601 Go of .mat<br>
<br><br>
: 922 x .P files<br>
total of  99  Go of .P files<br>
<br><br>
A metadata file (META.csv) is created from the scan of each .mat file, where each line contains<br>
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
## a quick estimation of epsilon values
quick_epsilon.m<br>
calculate_epsilon.m<br>
<br><br>
The liste_M.txt is used to loop over each .mat file and produce values of epsilon, the turbulent kinetic energy dissipation rates (W/kg). <br>
<br><br>
General parameters are applied<br>
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
N_window_fft = 2 ;
Percent_overlap_fft = 50./100 ; % round(Percent_overlap_fft.* (diss_info.fft_length) )

<br><br>
Values are conserved through a matlab structure QUICK_EPSILON for each file.<br>
PROJECT_DATAFILENAME_DATE_TIME.mat<br>
e.g. TERESA_microrider_DAT_008_2018-04-23_11-17.mat<br>
<br><br>
QUICK_EPSILON.profileNumber  = N ; %------------------ to be discarded, yoyo are now processed as a unique vector<br>
QUICK_EPSILON.totalProfiles  = N_profiles ; %--------  to be discarded, yoyo are now processed as a unique vector<br>
QUICK_EPSILON.date_profile = date_profile ;<br>
QUICK_EPSILON.time_profile = time_profile ;<br>
QUICK_EPSILON.lon_profile = lon_profile ;<br>
QUICK_EPSILON.lat_profile = lat_profile ;<br>
QUICK_EPSILON.yyyy = yyyy ;<br>
QUICK_EPSILON.mm = mm ;<br>
QUICK_EPSILON.dd = dd ;<br>
QUICK_EPSILON.HH = HH ;<br>
QUICK_EPSILON.MM = MM ;<br>
QUICK_EPSILON.SS = SS;<br>
QUICK_EPSILON.mSS = mSS;<br>
QUICK_EPSILON.dn = dn;<br>
QUICK_EPSILON.PROFILE_ID = PROFILE_ID;<br>
<br>
QUICK_EPSILON.odas_version        = x.odas_version;<br>
QUICK_EPSILON.vehicle_info        = x.vehicle_info;<br>
QUICK_EPSILON.fullPath = fullPath ;<br>
QUICK_EPSILON.setupfilestr = setupfilestr ;<br>
<br>
QUICK_EPSILON.Ntresh  = Ntresh ;<br>
QUICK_EPSILON.fcut    = fcut ;<br>
QUICK_EPSILON.fs      = fs ;<br>
QUICK_EPSILON.cm_cut  = cm_cut ;<br>
QUICK_EPSILON.freq_HP_cut      = freq_HP_cut ;  <br>
QUICK_EPSILON.freq_LP_cut      = freq_LP_cut ;  <br>
QUICK_EPSILON.segment_fft = segment_fft; <br>
QUICK_EPSILON.N_segment_fft = N_segment_fft;<br>
QUICK_EPSILON.N_window_fft = N_window_fft;<br>
QUICK_EPSILON.Percent_overlap_fft = Percent_overlap_fft;<br>
QUICK_EPSILON.fft_length = diss_info.fft_length ;<br>
QUICK_EPSILON.diss_length = diss_info.diss_length ;<br>
QUICK_EPSILON.overlap = diss_info.overlap ;<br>
<br>
QUICK_EPSILON.PE = PE;<br>
QUICK_EPSILON.FM = FM;<br>
QUICK_EPSILON.DE = DE;<br>
QUICK_EPSILON.ee = ee;<br>
QUICK_EPSILON.e1 = e1;<br>
QUICK_EPSILON.e2 = e2;<br>
QUICK_EPSILON.fm1 = fm1;<br>
QUICK_EPSILON.fm2 = fm2;<br>
QUICK_EPSILON.ze = ze;<br>
QUICK_EPSILON.te = te;<br>
        
