%% from shear probes

P_e = diss.P  ;
D_e = - gsw_z_from_p(diss.P,mean(glider_lat(fu)).*ones(length(diss.P),1)); % ADD ON
Z_e = D_e ;
t_e = diss.t ;
T_e = diss.T ;
speed_e = diss.speed ; 
speed_std_e = diss.speed_std ; 

dt_eps = datetime(x.Year,x.Month,x.Day,x.Hour,x.Minute,x.Second,x.Milli) ; 
dt_eps = dt_eps + seconds(t_e) ; 
dn_eps = datenum(dt_eps);

AOA = interp1(dn(fu), glider_aoa(fu),dn_eps,'linear');  %---------------- METADATA

lon_eps = interp1(dn(fu), glider_lon(fu),dn_eps,'linear');  %---------------- METADATA
lat_eps = interp1(dn(fu), glider_lat(fu),dn_eps,'linear');  %---------------- METADATA

T1_e = interp1(dn(fu), x.T1_fast(fu),dn_eps,'linear');  %---------------- METADATA
T2_e = interp1(dn(fu), x.T2_fast(fu),dn_eps,'linear'); 

dz = mean(abs(diff(Z_e))); %---------------- METADATA
Ne = length(Z_e)  ;  %---------------- METADATA

SP_e = interp1(dn(fu), glider_sali(fu),dn_eps,'linear');  %---------------- METADATA
SA_e = interp1(dn(fu), glider_asal(fu),dn_eps,'linear');  %---------------- METADATA
C_e = interp1(dn(fu), glider_cond(fu),dn_eps,'linear');  %---------------- METADATA
rho_e = interp1(dn(fu), glider_rho(fu),dn_eps,'linear');  %---------------- METADATA

pitch_e = interp1(dn(fu), glider_pitch(fu),dn_eps,'linear');  %---------------- METADATA
roll_e = interp1(dn(fu), glider_roll(fu),dn_eps,'linear');  %---------------- METADATA
pitch = interp1(dn(fu), -Incl_Y_fast(fu),dn_eps,'linear');  %---------------- METADATA
roll = interp1(dn(fu), Incl_X_fast(fu),dn_eps,'linear');  %---------------- METADATA

%%

dof_spec = diss.dof_spec ; %------------------
k = diss.K(:,1) ;  %------------------
nu = diss.nu' ;

diss.Nasmyth_spec = diss.Ref_spec ; 

% SH1_auto = diss.sh_clean(:,2,1,s) ; 
% SH2_auto = diss.sh_clean(:,1,2,s) ;

% [M N N P],
% [M 1 1 P],
% [M 1 2 P], auto
% [M 2 1 P], auto
% [M 2 2 P],

method_sh1 = diss.method(1,:) ;
dof_e1 = diss.dof_e(1,:) ;
e1 = diss.e(1,:) ;
fm1 = diss.FM(1,:) ;
fm_old1 = diss.FM_old(1,:) ;
mad10_1 = diss.mad10(1,:) ;

mad1 = diss.mad(1,:) ;
kmax1 = diss.K_max(1,:) ;
Nasmyth_spec1 = diss.Nasmyth_spec(:,1,:) ;
sh_1 = diss.sh(:,1,1,:) ;
sh_clean1 = diss.sh_clean(:,1,1,:) ;

method_sh2 = diss.method(2,:) ;
dof_e2 = diss.dof_e(2,:) ;
e2 = diss.e(2,:) ;
fm2 = diss.FM(2,:)   ;
fm_old2 = diss.FM_old(2,:) ;

mad2 = diss.mad(2,:) ;
mad10_2 = diss.mad10(2,:) ;

kmax2 = diss.K_max(2,:) ;
Nasmyth_spec2 = diss.Nasmyth_spec(:,2,:) ;
sh_2 = diss.sh(:,2,2,:) ;
sh_clean2 = diss.sh_clean(:,2,2,:) ;

%% spectral integration vsr


method_sh1_vsr = zeros(size(method_sh1));
method_sh2_vsr = zeros(size(method_sh2));
method_sh1_vsr(method_sh1==0) = 1 ;
method_sh2_vsr(method_sh2==0) = 1 ;
fvsr = method_sh1_vsr.*method_sh2_vsr ; 

le = speed_for_fft*Ntimes_N_fft./x.fs_fast ;
%le = w_mean_speed*N_fft./x.fs_fast ;

kolmo1 = (nu.^3./e1).^(1/4) ;
kolmo2 = (nu.^3./e2).^(1/4) ;

k_hat_u1 = kmax1 .* kolmo1 ;
k_hat_u2 = kmax2 .* kolmo2 ;

I_N0 = @(k) tanh(48 * k.^(4/3)) - 2.9 * k.^(4/3) .* exp(-22.3 * k.^(4/3));

Vf1 = I_N0(k_hat_u1);
Vf2 = I_N0(k_hat_u2);
%----------
Vf1 = diss.VR(1,:) ; 
Vf2 = diss.VR(2,:) ; 
%----------

Lf1 = (le./kolmo1).*Vf1.^(3/4);
Lf2 = (le./kolmo2).*Vf2.^(3/4);

var_ln_e1 = 5.5./(1 + (Lf1/4).^(7/9));
var_ln_e2 = 5.5./(1 + (Lf2/4).^(7/9));

% figure; plot(sig_ln_e1); hold on; plot(diss.sigma_e(1,:))
sig_ln_e1 = sqrt(var_ln_e1);
sig_ln_e2 = sqrt(var_ln_e2);
%--------
sig_ln_e1 = diss.sigma_e(1,:) ; 
sig_ln_e2 = diss.sigma_e(2,:) ; 
%--------

sig_ln_ee = (sig_ln_e1 + sig_ln_e2)/2;

upp_e1_vsr = e1.* exp(sqrt(2).*1.96.*sig_ln_e1);
low_e1_vsr = e1.* exp(sqrt(2).*-1.96.*sig_ln_e1);
upp_e2_vsr = e2.* exp(sqrt(2).*1.96.*sig_ln_e2);
low_e2_vsr = e2.* exp(sqrt(2).*-1.96.*sig_ln_e2);

%% isr

method_sh1_isr = zeros(size(method_sh1));
method_sh2_isr = zeros(size(method_sh2));
method_sh1_isr(method_sh1>0) = 1 ;
method_sh2_isr(method_sh2>0) = 1 ;

fisr = method_sh1_isr.*method_sh2_isr ; 

Ns1 = dof_e1 ;
Ns2 = dof_e2 ;
%--------
Ns1 = diss.N_s(1,:) ; 
Ns2 = diss.N_s(2,:) ; 
%--------

TM1 = 0.8 + 1.25./sqrt(Ns1) ;
TM2 = 0.8 + 1.25./sqrt(Ns2) ;

sig_ln_psi_on_sqrt_Ns1 = sig_ln_psi./sqrt(Ns1);
sig_ln_psi_on_sqrt_Ns2 = sig_ln_psi./sqrt(Ns2);
% figure; plot(sig_ln_psi_on_sqrt_Ns1); hold on; 

sig_ln_psi_on_sqrt_Ns = (sig_ln_psi_on_sqrt_Ns1 + sig_ln_psi_on_sqrt_Ns2)/2;

upp_e1_isr = e1.* exp(sqrt(2).*1.96.*(3/2).*sig_ln_psi_on_sqrt_Ns1);
low_e1_isr = e1.* exp(sqrt(2).*1.96.*(3/2).*sig_ln_psi_on_sqrt_Ns1);
upp_e2_isr = e2.* exp(sqrt(2).*1.96.*(3/2).*sig_ln_psi_on_sqrt_Ns2);
low_e2_isr = e2.* exp(sqrt(2).*1.96.*(3/2).*sig_ln_psi_on_sqrt_Ns2);
