%% controlboard
tic
% startup 
startup_microrider

% %%
% volume_path = '/Volumes/DENISE/' ;
% data_pth = '/Volumes/DENISE/data/glider/teresa/microrider/data/' ;
% config_pth = '/Volumes/DENISE/data/glider/teresa/microrider/config/' ;
% 
% mird_pth = '/Volumes/DENISE/data/glider/teresa/microrider/processing/L0/' ;
% outp_pth = '/Volumes/DENISE/data/glider/teresa/microrider/processing/L1/' ;
% outm_pth = '/Volumes/DENISE/data/glider/teresa/microrider/processing/L1/mat/' ;
% 
% hotel_file_MR = '/Volumes/DENISE/data/glider/teresa/glider/HOTEL/slocum_glider_teresa_hotel_for_MR.mat' ;
% 
% %hotel = load(hotel_file_MR) ;

%%
volume_path = '/Volumes/ORSO/' ;
data_pth = '/Volumes/ORSO/TERESA/microrider/data/' ;
config_pth = '/Volumes/ORSO/TERESA/microrider/config/' ;

mird_pth = '/Volumes/ORSO/TERESA/microrider/processing/L0/' ;
outp_pth = '/Volumes/ORSO/TERESA/microrider/processing/L1/' ;
outm_pth = '/Volumes/ORSO/TERESA/microrider/processing/L1/mat/' ;

hotel_file_MR = '/Volumes/ORSO/TERESA/glider/HOTEL/slocum_glider_teresa_hotel_for_MR.mat' ;

%hotel = load(hotel_file_MR) ;

%%
% TEMRa = readtable([mird_pth,'TEMR_L0_organize.csv']);
% TEMRb = readtable([mird_pth,'TEMR_L0_2024_organize.csv']);
% TEMR = [TEMRa; TEMRb];


TEMR = readtable([mird_pth,'TEMR_L0_2015_2024_organize.csv']);

s1 = size(TEMR,1);
s2 = size(TEMR,2);
disp(TEMR.Properties.VariableNames)
list_L0 = TEMR.filepaths ;
list_L0 = cellstr(TEMR.filepaths);

status_list = zeros(s1, 1); 
fid_mat = strings(s1, 1); 
out_pth = strings(s1, 1); 
errormsg = strings(s1, 1); 
toffset = zeros(s1, 1); 
repatched = zeros(s1, 1); 
dn = zeros(s1, 1); 
lo = zeros(s1, 1); 
la = zeros(s1, 1); 

% 0 not processed
% 1 ok
% 2 erreur
% 3 erreur

str1 = 'ok' ; 
str2 = 'Index exceeds the number of array elements. Index must not exceed 0.' ; 
str3 = 'still exotic' ; 

%%
p = 400;
p = 450;

for p = 1:s1
    disp([num2str(p),'/',num2str(s1),'  ',num2str(round(100.*p./s1)),'%'])
    pause(2)
    fpath = list_L0{p};
    disp(fpath)
    
    tokens = regexp(fpath, 'microrider_(\d{4})', 'tokens');
    YYYY_file = str2double(tokens{1}{1});
    parts = split(fpath, '/');
    last_part = split(parts{end},'.');
    last_part = last_part{1} ; 
    
    % Extract setup_cfg and save it
    try
        extract_setupstr(fpath, [config_pth,'before/p',num2str(p),'_',num2str(YYYY_file),'_',last_part,'.cfg']);
    catch ME
    end
    
    % repatch of early years (Anneke's update)
    cfg_file = [config_pth,'update/setup_216_corrected.cfg'] ; 
    if YYYY_file <= 2022 
        disp([num2str(YYYY_file),'  repatched'])
        try
            patch_setupstr(fpath,cfg_file)
            repatched(p) = 1 ;
            cannotpatch = 0 ; 
        catch ME
            status_list(p) = 3;
            errormsg(p)= ME.message ;
            cannotpatch = 1 ; 
            repatched(p) = -1 ;
        end

    end
    
    if cannotpatch == 1
        continue
    end
    
    
    % Conversion infos
    conversion_info = odas_p2mat();
    %
    conversion_info.constant_temp = 'T_gl' ; 
    conversion_info.hotel_file = [hotel_file_MR] ;
    conversion_info.vehicle = 'slocum_glider' ;
    conversion_info.time_offset = 0 ;
    %conversion_info.time_offset = 2*toffset ;
    conversion_info;
    pause(4)
    clc
    
    %% odas
    clear x
    try
        x = odas_p2mat(fpath,conversion_info);
        status_list(p) = 1;
        errormsg(p)= str1 ;
    
        %
        dt0 = datetime(x.Year,x.Month,x.Day,x.Hour,x.Minute,x.Second,x.Milli) ;
        dt = dt0 + seconds(x.t_slow) ;
        dt00 = datetime(x.Year,x.Month,x.Day,x.Hour,x.Minute,x.Second,x.Milli) ;
        yyyy = year(dt00);
        dst_start = datetime(yyyy,3,31) - days(weekday(datetime(yyyy,3,31))-1); % dernier dimanche de mars
        dst_end = datetime(yyyy,10,31) - days(weekday(datetime(yyyy,10,31))-1); % dernier dimanche d’octobre
        if dt00 >= dst_start && dt00 < dst_end
            offset = 7200; % UTC+2
        else
            offset = 3600; % UTC+1
        end
        disp(offset)
        
        toffset(p) = offset  ;
        dn(p) = datenum(dt0) ;
        lo(p) = mean(x.lon_gl_slow,'omitnan') ;
        la(p) = mean(x.lat_gl_slow,'omitnan') ;
    
    catch ME
    
            if strcmp(ME.message,str2)
                status_list(p) = 2;
                errormsg(p)= str2 ;
    
            else 
                status_list(p) = 3;
                errormsg(p)= ME.message ;
                end
    
    end
    pause(4)
    clc

clear x
end

TEMR.L1_dataconversion = status_list;
TEMR.L1_errormsg = errormsg;
%TEMR.L1_matfile = fid_mat;
%TEMR.L1_matfilepath = out_pth;
TEMR.L1_repatched = repatched;
TEMR.L1_time_offset = toffset;
TEMR.L1_dn = dn;
TEMR.L1_lo = lo;
TEMR.L1_la = la;

% Exporter la table mise à jour vers un fichier CSV
output_csv = [outp_pth,'TEMR_L1_overview.csv'];  % Nom du fichier CSV de sortie
output_csv = [outp_pth,'TEMR_L1_2015_2024_overview.csv'];  % Nom du fichier CSV de sortie

writetable((TEMR), output_csv);
disp('OVERVIEW DONE !')

toc

% 
% 
% OVERVIEW DONE !
% Elapsed time is 23671.407402 seconds.
% 23671/60
% 
% ans =
% 
%   394.5167
% 
% 23671/3600
% 
% ans =
% 
%     6.5753



