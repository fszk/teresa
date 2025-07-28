%% Fs Fc from velocity

P_fast = x.P_fast ; 
W_fast = x.W_fast ;
speed_fast = x.speed_fast ; 
%speed_method = 'speed_fast' ;

speed_mean_0 = abs(round(mean(speed_fast,'omitnan'),2)) ;

Fs = x.fs_fast    ;  
Fc = round(mean(speed_fast),2)          ;  %------------------- FEAT - QC
order = 1         ;  %------------------- FEAT - QC
fc_lo = Fc ;
fc_lo_order = order ; 
% Design Butterworth filter
[b, a] = butter(order, Fc / (Fs / 2), 'low');
W_fast_lp = filtfilt(b, a, W_fast) ; 

%% minima P W Duration
Pmin = 3 ; 
Wmin = 0.01; 
minDuration = 60 ;


%% get profile / section
W_fast_in = W_fast_lp ;
%W_fast_in = W_fast ;

ip_down_fast = get_profile( P_fast, W_fast_in, Pmin, Wmin, 'down', minDuration, x.fs_fast ) ;
sp = size(ip_down_fast) ;
N_down = sp(2) ;
% figure ;
% hold on
%   for p = 1:sp(2)
%     plot(dt(ip_down_fast(1,p):ip_down_fast(2,p)),...
%         x.P_fast(ip_down_fast(1,p):ip_down_fast(2,p)),'-','LineWidth',2,'Color',rand(1, 3))
% end


ip_upwd_fast = get_profile( P_fast, W_fast_in, Pmin, -Wmin, 'up', minDuration, x.fs_fast ) ;

sp = size(ip_upwd_fast) ;
N_upwd = sp(2) ;
%   for p = 1:sp(2)
%      plot(dt(ip_upwd_fast(1,p):ip_upwd_fast(2,p)),...
%          x.P_fast(ip_upwd_fast(1,p):ip_upwd_fast(2,p)),'-','LineWidth',2,'Color',rand(1, 3))
%   end


%% direction 
sections_down = zeros(length(x.P_fast),1) ;
sections_upwd = zeros(length(x.P_fast),1) ;
sections_dire = zeros(length(x.P_fast),1) ;

pp = 0 ;
for p = 1:size(ip_down_fast,2)
ip1 = ip_down_fast(1,p) ;
ip2 = ip_down_fast(2,p) ;
sections_down(ip1:ip2) = pp+2*p ;
sections_dire(ip1:ip2) = 1 ;
%disp(pp+2*p)
end

pp = 1 ;
for p = 1:size(ip_upwd_fast,2)
ip1 = ip_upwd_fast(1,p) ;
ip2 = ip_upwd_fast(2,p) ;
sections_upwd(ip1:ip2) = pp+2*p ;
sections_dire(ip1:ip2) = -1 ;
%disp(pp+2*p)
end

sections_down ; %---------------- KEEP
sections_upwd ; %---------------- KEEP
sections_dire ; %---------------- KEEP
sections = sections_down + sections_upwd ;  %---------------- KEEP
%
mask_sections = nan(length(P_fast),1) ;
f = find(sections>0) ;
mask_sections(f) = 1 ;

%%
[us, is] = unique(sections) ; us = us(2:end); is = is(2:end) ; 
lu = length(us) ; %---------------- METADATA
sections ; % %---------------- METADATA


