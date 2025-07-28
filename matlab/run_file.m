%%

%ppp = 500 ;

%% filepath
disp([num2str(ppp),'/',num2str(s1),'  ',num2str(round(100.*ppp./s1)),'%']); pause(2)
fpath = list_L0{ppp};
disp(fpath); pause(2)

%% get filename and year
tokens = regexp(fpath, 'microrider_(\d{4})', 'tokens');
YYYY_file = str2double(tokens{1}{1});
parts = split(fpath, '/');
last_part = split(parts{end},'.');
last_part = last_part{1} ; 
    

%% Extract setup_cfg and save it
try
    extract_setupstr(fpath, [config_pth,'before/p',num2str(ppp),'_',num2str(YYYY_file),'_',last_part,'.cfg']);
catch ME
    disp('original setup_cfg already extracted and archived')
end

%% repatch of early years (Anneke's update)
cfg_file = [config_pth,'update/setup_216_corrected.cfg'] ; 
if YYYY_file <= 2022 
    disp([num2str(YYYY_file),'  repatched with setup_216_corrected.cfg'])
    try
        patch_setupstr(fpath,cfg_file)
        disp('done')
    catch ME
        disp('/ not patched')
    end
end

%% get offset
toffset = TEMR_OK.L1_time_offset(ppp);
coffset = 0 ; %-5

%% conversion structure
conversion_info = odas_p2mat();
conversion_info.constant_temp = 'T_gl' ; 
conversion_info.hotel_file = [hotel_file_MR]  ; %%%% NOT UTC, alignment is done after
conversion_info.vehicle = 'slocum_glider' ;
conversion_info.time_offset = -toffset-coffset ;
conversion_info;

%% log
write_log([processing_out_path,'log/'], 'boucle_log.txt', ppp);

%% odas_p2mat

try
    disp('odas_p2mat')
    x = odas_p2mat(fpath,conversion_info);
    odas_ok = 1 ; 
    clc
catch ME
    disp('odas_p2mat')
    disp(ME)
    odas_ok = 0 ; 
    clc
end

if odas_ok == 1
%% make file ID
disp('L1_make_FID')
L1_make_FID; disp(FID)

fname = FID ; 
fpath = [fname,'.mat'];  % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fparts = strsplit(fpath, filesep);

%% get time vectors
disp('L1_time') 

lx = size(x.t_fast,1); disp(['length:  ',num2str(lx)])

dt = datetime(x.Year,x.Month,x.Day,x.Hour,x.Minute,x.Second,x.Milli) ; 
dt = dt + seconds(x.t_fast) ; 
dn = datenum(dt);
disp(dt(1))

YEAR = x.Year ; 

dt_slow = datetime(x.Year,x.Month,x.Day,x.Hour,x.Minute,x.Second,x.Milli) ; 
dt_slow = dt_slow + seconds(x.t_slow) ; 
dn_slow = datenum(dt_slow);
%% get glider fields nested from hotel
disp('L1_glider') 
L1_glider

%% get continuous sections/profiles

disp('L2_sections') 
L2_sections
disp([num2str(lu),' sections'])


end

