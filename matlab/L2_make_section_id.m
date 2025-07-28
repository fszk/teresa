

toto = split(fname,'/');
toto = toto{end} ;
toto = split(toto,'_');



if YEAR == 2024 
    sfx1 = strjoin({toto{1},toto{2},'L2L3L4QC',toto{4},toto{5},toto{6},toto{7},toto{8}}, '_');
    sfx11 = strjoin({toto{1},toto{2},'L4QC',toto{4},toto{5},toto{6},toto{7},toto{8}}, '_');
else
    sfx1 = strjoin({toto{1},toto{2},'L2L3L4QC',toto{4},toto{5},toto{6},toto{7}}, '_');
    sfx11 = strjoin({toto{1},toto{2},'L4QC',toto{4},toto{5},toto{6},toto{7}}, '_');

end

sfx2 = datestr(dt(fu(1)), 'yyyy_mm_dd_HH_MM_SS');


lon_str = sprintf('%02d_%04d', floor(mean(section_lon,'omitnan')), round(mod(mean(section_lon,'omitnan'), 1) * 10000)); % Lon
lat_str = sprintf('%02d_%04d', floor(mean(section_lat,'omitnan')), round(mod(mean(section_lat,'omitnan'), 1) * 10000)); % Lat
sfx3 = strjoin({'lat',lat_str,'lon',lon_str}, '_');

sfx4 = strjoin({'nav',navig_dir_str(1)}, '_'); 

pmin_str = sprintf('%04d',  round(min(x.P_fast(fu))) ); 
pmax_str = sprintf('%04d',  round(max(x.P_fast(fu))) ); 
sfx5 = strjoin({'pmin',pmin_str,'pmax',pmax_str}, '_');

sec_str = sprintf('%03d',  u ); 
tot_str = sprintf('%03d',  lu ); 
sfx6 = strjoin({'sec',sec_str,'on',tot_str}, '_');

sfx7 = strjoin({'glid',dirfu_str}, '_');







FID =  strjoin({sfx1,sfx2,sfx3,sfx4,sfx5,sfx6,sfx7}, '_');
section_id = FID ;

FID_QC =  strjoin({sfx11,sfx2,sfx3,sfx4,sfx5,sfx6,sfx7}, '_');
section_id_qc = FID_QC ;

%  ['TERESA_MR_L2L3L4_converted_file_0012_' ...
%      '2024_05_22_10_20_40_l' ...
%      'on_03_5373_lat_39_5809_' ...
%      'E_' ...
%      'pmin_003_pmax_058_' ...
%      'section_001_on_040_' ...
%      'dir_down']






