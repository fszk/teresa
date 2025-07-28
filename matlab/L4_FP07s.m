%%
%% from thermistor

QT1 = dissT1.qc_flag_T ; 
NeT1 = length(dissT1.P_eT)  ; 
P_eT1 = dissT1.P_eT  ;
D_eT1 = - gsw_z_from_p(P_eT1,mean(glider_lat(fu)).*ones(NeT1,1)); 
Z_eT1 = D_eT1 ; 
dzT1 = mean(abs(diff(Z_eT1)));
dn_eT1 = dissT1.dn_eT ;
dt_eT1 = datetime(dn_eT1, 'ConvertFrom', 'datenum') ;
if isempty(dt_eT1)
    t_eT1 = [] ; 
else
    t_eT1 = seconds(dt_eT1-dt_eT1(1)) ;
end
lon_eT1 = interp1(dn(fu), glider_lon(fu),dn_eT1,'linear');  
lat_eT1 = interp1(dn(fu), glider_lat(fu),dn_eT1,'linear');  
T_eT1 = interp1(dn(fu),glider_temp(fu),dn_eT1,'linear'); 
nu_eT1 = interp1(dn(fu),glider_nu(fu),dn_eT1,'linear'); 
td_eT1 = interp1(dn(fu),glider_td(fu),dn_eT1,'linear'); 



%% from thermistor

QT2 = dissT2.qc_flag_T ; 
NeT2 = length(dissT2.P_eT)  ; 
P_eT2 = dissT2.P_eT  ;
D_eT2 = - gsw_z_from_p(P_eT2,mean(glider_lat(fu)).*ones(NeT2,1)); 
Z_eT2 = D_eT2 ; 
dzT2 = mean(abs(diff(Z_eT2)));
dn_eT2 = dissT2.dn_eT ; 
dt_eT2 = datetime(dn_eT2, 'ConvertFrom', 'datenum') ;
if isempty(dt_eT2)
    t_eT2 = [] ; 
else
    t_eT2 = seconds(dt_eT2-dt_eT2(1)) ;
end
lon_eT2 = interp1(dn(fu), glider_lon(fu),dn_eT2,'linear');  
lat_eT2 = interp1(dn(fu), glider_lat(fu),dn_eT2,'linear');  
T_eT2 = interp1(dn(fu),glider_temp(fu),dn_eT2,'linear'); 
nu_eT2 = interp1(dn(fu),glider_nu(fu),dn_eT2,'linear'); 
td_eT2 = interp1(dn(fu),glider_td(fu),dn_eT2,'linear'); 




