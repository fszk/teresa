%% 
%% QC
% 
Q_T1 = double(~QT1) ; % 0 = pass, 1 = fail
Q_T2 = double(~QT2) ;
% 
Q_T1(Q_T1 == 1) = 2;
Q_T2(Q_T2 == 1) = 2;

Q_T1T2_SPREAD = zeros(size(Q_T1)); 
f = find(abs(log10((dissT2.xit)./(dissT1.xit))) > 1)  ; 
Q_T1T2_SPREAD(f) = 1 ; 

QC_FP07_1 = Q_T1 + Q_T1T2_SPREAD ;
QC_FP07_2 = Q_T2 + Q_T1T2_SPREAD ; 


%%
f1 = find( (QC_FP07_1==1)   );
l1 = length(f1) ;
p1 = round(100*l1/length(QC_FP07_1),3);

f2 = find( (QC_FP07_2==1)   );
l2 = length(f2) ;
p2 = round(100*l2/length(QC_FP07_2),3);



%%
checkplotcheckplot = 0 ;
if checkplotcheckplot == 1
figure; 
plot(dt_eT1,log10(dissT1.xit),'bo' ); 
hold on;
plot(dt_eT2,log10(dissT2.xit),'rs' );

figure; 
plot(dt_eps(f1),log10(dissT1.xit(f1)),'bo' ); 
hold on;
plot(dt_eps(f2),log10(dissT1.xit(f2)),'rs' );
end
