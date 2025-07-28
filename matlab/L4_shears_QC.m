%% We define QC-fail conditions
%

%%
figofmer_fail = 1.4; 
spikfrac_fail = 5/100; 
spikfrac_failr = 15/100; 
spikpass_fail = 9; 
epsratio_fail_vsr = 2.77.* sig_ln_ee ;
epsratio_fail_isr = 4.2 .* sig_ln_psi_on_sqrt_Ns ;
varresol_fail = 0.6; 
aoa_thresh1 = 1.5 ;
aoa_thresh2 = 4.5 ;

%%

Q_figofmer1 = zeros(size(e1));
Q_spikfrac1 = zeros(size(e1));
Q_spikpass1 = zeros(size(e1));
Q_varresol1 = zeros(size(e1));
Q_e1_eratio = zeros(size(e1));

Q_figofmer2 = zeros(size(e1));
Q_spikfrac2 = zeros(size(e1));
Q_spikpass2 = zeros(size(e1));
Q_varresol2 = zeros(size(e1));
Q_e2_eratio = zeros(size(e1));

Q_aoa = zeros(size(e1));



%% FOM
Q_figofmer1(fm1>figofmer_fail) = 1 ;
Q_figofmer2(fm2>figofmer_fail) = 1 ;

%% VARIANCE RESOLVED

Q_varresol1(Vf1<varresol_fail) = 16 ;
Q_varresol2(Vf2<varresol_fail) = 16 ;

%% SPIKES
%

dspk_frac1 = ones(size(e1)) .* despike_sh1_ratio ; 
dspk_frac2 = ones(size(e1)) .* despike_sh2_ratio ; 

if despike_sh1_ratio > spikfrac_failr
    Q_spikfrac1 = Q_spikfrac1 + 2 ;
elseif  (despike_sh1_ratio > spikfrac_fail) && (despike_sh1_ratio <= spikfrac_failr)
    Q_spikfrac1 = Q_spikfrac1 + 32 ;
end


%
if despike_sh1_passcount > spikpass_fail
    Q_spikpass1 = Q_spikpass1 + 8 ;
end


%

%
if despike_sh2_ratio > spikfrac_failr
    Q_spikfrac2 = Q_spikfrac2 + 2 ;
elseif (despike_sh2_ratio > spikfrac_fail) && (despike_sh2_ratio <= spikfrac_failr)
    Q_spikfrac2 = Q_spikfrac2 + 32 ;
end

%
if despike_sh2_passcount > spikpass_fail
    Q_spikpass2 = Q_spikpass2 + 8 ;
end

%% AOA




Q_aoa(abs(AOA)<aoa_thresh1) = 128 ;
Q_aoa(abs(AOA)>aoa_thresh2) = 128 ;

%% RATIO (and METHODS)

%
for m =1:length(method_sh1)
    e12 = abs(log(e1(m)) - log(e2(m))) ;

    if (method_sh1(m)+method_sh2(m))==0
        if  e12 > epsratio_fail_vsr(m)
            if e1(m) > e2(m)
                Q_e1_eratio(m) = 4 ;
            elseif e2(m) > e1(m)
                Q_e2_eratio(m) = 4 ;
            end
        end
        
    elseif (method_sh1(m)+method_sh2(m))==2
         if  e12 > epsratio_fail_isr(m)
            if e1(m) > e2(m)
                Q_e1_eratio(m) = 4 ;
            elseif e2(m) > e1(m)
                Q_e2_eratio(m) = 4 ;
            end
         end

    elseif (method_sh1(m)+method_sh2(m))==1
        Q_e1_eratio(m) = 64 ;
        Q_e2_eratio(m) = 64 ;
  
    end
end



%%

Q1 = Q_figofmer1 + Q_spikfrac1 + Q_spikpass1 + Q_varresol1 + Q_e1_eratio  + Q_aoa;
Q2 = Q_figofmer2 + Q_spikfrac2 + Q_spikpass2 + Q_varresol2 + Q_e2_eratio  + Q_aoa;

%%
fok1 = find( (Q1==0)  | (Q1==32) );
lok1 = length(fok1) ;
pok1 = round(100*lok1/Ne,3);

fok2 = find( (Q2==0) | (Q2==32) );
lok2 = length(fok2) ;
pok2 = round(100*lok2/Ne,3) ;

%%
checkplotcheckplot = 0 ;
if checkplotcheckplot == 1
figure; 
plot(dt_eps,log10(e1),'bo' ); 
hold on;
plot(dt_eps,log10(e2),'rs' );

figure; 
plot(dt_eps(fok1),log10(e1(fok1)),'bo' ); 
hold on;
plot(dt_eps(fok2),log10(e2(fok2)),'rs' );
end

%% OLD




%%


%------
ee_vsr = nan(size(e1)) ; 
ok_vsr = [] ; 
f = find(Q1 == 0 & Q2 == 0 & fvsr == 1) ; 
if isempty(f)==0
    ee_vsr(f) = (e1(f)+e2(f))/2 ;
    ok_vsr = f ; 
end

ee_isr = nan(size(e1)) ; 
ok_isr = [] ; 
f = find(Q1 == 0 & Q2 == 0 & fisr == 1) ; 
if isempty(f)==0
    ee_isr(f) = (e1(f)+e2(f))/2 ;
    ok_isr = f ; 
end

%------
ee_vsr_qc1 = nan(size(e1)) ; 
ok_vsr_qc1 = [] ; 
f = find(Q1 <= 1 & Q2 <= 1 & fvsr == 1) ; 
if isempty(f)==0
    ee_vsr_qc1(f) = (e1(f)+e2(f))/2 ;
    ok_vsr_qc1 = f ; 
end

ee_isr_qc1 = nan(size(e1)) ;
ok_isr_qc1 = [] ; 
f = find(Q1 <= 1 & Q2 <= 1 & fisr == 1) ; 
if isempty(f)==0
    ee_isr(f) = (e1(f)+e2(f))/2 ;
    ok_isr_qc1 = f ; 
end


