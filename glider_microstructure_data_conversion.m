%% controlboard

% startup VMP

matlab_path_1 = ['/Users/floriankokoszka/Desktop/matlab/general/functions'] ;
matlab_path_2 = ['/Users/floriankokoszka/Desktop/matlab/general/toolboxes'] ;
addpath(matlab_path_1)
addpath(matlab_path_2)


matlab_path_mmap = ['/Users/floriankokoszka/Desktop/matlab/general/toolboxes/m_map'] ;
addpath(matlab_path_mmap)

%matlabtools_path4 = ['/Users/floriankokoszka/Desktop/matlab/general/toolboxes/gsw_matlab_v3_04'] ; % GIBBS SEAWATER
gsw_path_1 = ['/Users/floriankokoszka/Desktop/matlab/general/toolboxes/gsw_matlab_v3_06_11'] ; % GIBBS SEAWATER
gsw_path_2 = ['/Users/floriankokoszka/Desktop/matlab/general/toolboxes/gsw_matlab_v3_06_11/html'] ; % GIBBS SEAWATER
gsw_path_3 = ['/Users/floriankokoszka/Desktop/matlab/general/toolboxes/gsw_matlab_v3_06_11/library'] ; % GIBBS SEAWATER
gsw_path_4 = ['/Users/floriankokoszka/Desktop/matlab/general/toolboxes/gsw_matlab_v3_06_11/pdf'] ; % GIBBS SEAWATER
addpath(gsw_path_1)
addpath(gsw_path_2)
addpath(gsw_path_3)
addpath(gsw_path_4)



odas_path = ['/Users/floriankokoszka/Desktop/matlab/general/toolboxes/odas_v4.5.1/']  ;
addpath(odas_path)




%
volume_path = '/Volumes/DENISE/' ;

% 
TERESA_path = [volume_path,'data/glider/teresa/'] ;
TERESA_path
addpath(TERESA_path)


%%

dpath = [TERESA_path,'data/teresa_microrider_2023/data/'] ;
cd(dpath) ;
lf = dir('DAT*.P') 
sl = size(lf) 


dpath = [TERESA_path,'data/teresa_microrider_2022/data/'] ;
cd(dpath) ;
lf = dir('DAT*.P') 
sl = size(lf) 


dpath = [TERESA_path,'data/teresa_microrider_2018_maggio/data/'] ;
cd(dpath) ;
lf = dir('DAT*.P') 
sl = size(lf) 


dpath = [TERESA_path,'data/teresa_microrider_2017_marzo_aprile/data/'] ;
cd(dpath) ;
lf = dir('DAT*.P') 
sl = size(lf) 


dpath = [TERESA_path,'data/teresa_microrider_2017_MISC/data/'] ;
cd(dpath) ;
lf = dir('DAT*.P') 
sl = size(lf) 


dpath = [TERESA_path,'data/teresa_microrider_2015_july/data/'] ;
cd(dpath) ;
lf = dir('DAT*.P') 
sl = size(lf) 


dpath = [TERESA_path,'data/teresa_microrider_2015_agosto/data/'] ;
cd(dpath) ;
lf = dir('DAT*.P') 
sl = size(lf) 




for k = 1:sl(1)
fname = lf(k).name;
display(fname)
fpath = [dpath,fname] ;
ql_info = quick_look ;
x = quick_look(fpath,1,1000,ql_info) 
pause(3)
close all
clc
end

