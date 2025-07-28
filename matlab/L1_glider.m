%% get glider variables 
%% eventually align/sync check
%
checkplot = 0 ;
if checkplot == 1
figure; plot(dt_slow,x.P_slow); hold on; plot(dt_slow,x.P_gl_slow.*10)
figure; plot(dt_slow,-x.Incl_Y ); hold on; plot(dt_slow,x.pitch_gl_slow)
figure; plot(dt_slow,x.speed_slow ); hold on; 

figure; plot(dt_slow,x.T1_slow ); hold on; plot(dt_slow,x.T2_slow)

end
%%

glider_speed = x.speed_fast ; 
speed_source_logic = 1 ;
%%
Incl_Y = x.Incl_Y ; 
Incl_X = x.Incl_X ; 
Incl_Y_fast = interp1(dn_slow, Incl_Y, dn, 'linear') ; 
Incl_X_fast = interp1(dn_slow, Incl_X, dn, 'linear') ; 

%%
T_gl_slow = x.T_gl_slow ; 
T_gl_fast = interp1(dn_slow, T_gl_slow, dn, 'linear') ; 
glider_temp = T_gl_fast ; 



%%
lat_gl_slow = x.lat_gl_slow ; 
lon_gl_slow = x.lon_gl_slow ; 

lon_gl_fast = interp1(dn_slow, lon_gl_slow, dn, 'linear') ; 
lat_gl_fast = interp1(dn_slow, lat_gl_slow, dn, 'linear') ; 

glider_lon = lon_gl_fast ; 
glider_lat = lat_gl_fast ; 

%%
aoa_slow = x.aoa_slow ; 
aoa_fast = interp1(dn_slow, aoa_slow, dn, 'linear') ;
glider_aoa = aoa_fast ; 

%%
SP_gl_slow = x.SP_gl_slow ; 
SP_gl_fast = interp1(dn_slow, SP_gl_slow, dn, 'linear') ;
glider_sali = SP_gl_fast ; 

SA_gl_slow = x.SA_gl_slow ; 
SA_gl_fast = interp1(dn_slow, SA_gl_slow, dn, 'linear') ;
glider_asal = SA_gl_fast ;

%%
C_gl_slow = x.C_gl_slow ; 
C_gl_fast = interp1(dn_slow, C_gl_slow, dn, 'linear') ;
glider_cond = C_gl_fast ; 


%%
SP_fast = SP_gl_fast ; 
NU_fast = SW_Kviscosity(T_gl_fast,'C', SP_fast,'ppt') ;
TD_fast = SW_Diffusivity(T_gl_fast,'C', SP_fast,'ppt') ;
glider_nu = NU_fast ; 
glider_td = TD_fast ; 

%%
rho_gl_slow = x.rho_gl_slow ; 
rho_gl_fast = interp1(dn_slow, rho_gl_slow, dn, 'linear') ;
glider_rho = rho_gl_fast ; 


%%
pitch_gl_slow = x.pitch_gl_slow ; 
pitch_gl_fast = interp1(dn_slow, pitch_gl_slow, dn, 'linear') ;
glider_pitch = pitch_gl_fast ; 

roll_gl_slow = x.roll_gl_slow ; 
roll_gl_fast = interp1(dn_slow, roll_gl_slow, dn, 'linear') ;
glider_roll = roll_gl_fast ; 


buoyancy_change_gl_slow = x.buoyancy_change_gl_slow ; 
buoyancy_change_gl_fast = interp1(dn_slow, buoyancy_change_gl_slow, dn, 'linear') ;
glider_buoy = buoyancy_change_gl_fast ; 

%%
% lon_gl_slow: [2188736×1 double]
%                 lat_gl_slow: [2188736×1 double]
%                    aoa_slow: [2188736×1 double]
%                   C_gl_slow: [2188736×1 double]
%                   T_gl_slow: [2188736×1 double]
%                   P_gl_slow: [2188736×1 double]
%                  SP_gl_slow: [2188736×1 double]
%                  SA_gl_slow: [2188736×1 double]
%                 rho_gl_slow: [2188736×1 double]
%               pitch_gl_slow: [2188736×1 double]
%                roll_gl_slow: [2188736×1 double]
%     buoyancy_change_gl_slow: [2188736×1 double]
