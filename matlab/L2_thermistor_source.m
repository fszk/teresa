
%% FP07 master identification, 1 or 2, 
% first: more similar to T_glider 
% then from variance

%%
T1 = x.T1_slow ; 
T2 = x.T2_slow ; 
T3 = T_gl_slow ; 



var_thresh = 1e-5;
cor_thresh = 0.3;

%%
c1 = corr(T1(:), T3(:));
c2 = corr(T2(:), T3(:));

disp(['FP07_1/T_glider: ', num2str(c1),'  FP07_2/T_glider: ',num2str(c2)])
switch true

      case c1 > cor_thresh && c2 > cor_thresh
        %disp('at least one correlation significant');
        if c1 > c2
            %disp('T1 est le plus similaire à T3');
            selected = 'T1_fast';
            Th_source_logic = 11 ; 
            Th_source_string = selected ; 
            disp('both FP07 correlates with T | FP07_1 is chosen')
        else
            %disp('T2 est le plus similaire à T3');
            selected = 'T2_fast';
            Th_source_logic = 22 ;
            Th_source_string = selected ; 
            disp('both FP07 correlates with T | FP07_2 is chosen')
        end


    case c1 > cor_thresh || c2 > cor_thresh
        %disp('at least one correlation significant');
        if c1 > c2
            %disp('T1 est le plus similaire à T3');
            selected = 'T1_fast';
            Th_source_logic = 1 ; 
            Th_source_string = selected ; 
            disp('FP07_1 correlates with T | FP07_1 is chosen')
        else
            %disp('T2 est le plus similaire à T3');
            selected = 'T2_fast';
            Th_source_logic = 2 ;
            Th_source_string = selected ; 
             disp('FP07_2 correlates with T | FP07_2 is chosen')
        end


    case c1 <= cor_thresh && c2 <= cor_thresh
        %disp('none  significant');
        
        var1 = var(x.T1_fast, 'omitnan');
        var2 = var(x.T2_fast, 'omitnan');
        
        if var1 < var_thresh && var2 < var_thresh   %---------- both below thresh
                % both low variance = Functioning error
                % error('both potentially noisy');
                %selected = x.params.temperature_source ;
                selected = 'WARNING_potential_thermistors_malfunctions' ; 
                Th_source_logic = 0 ; 
                Th_source_string = selected ; 
                disp('none of FP07 correlates with T | both have strange variance')
        elseif var1 < var_thresh
                selected = 'T2_fast';
                Th_source_logic = 200 ;
                Th_source_string = selected ; 
                disp('none of FP07 correlates with T | 1 is strange | 2 is chosen')
        elseif var2 < var_thresh
                selected = 'T1_fast';
                Th_source_logic = 100 ; 
                Th_source_string = selected ; 
                disp('none of FP07 correlates with T | 2 is strange | 1 is chosen')
        else
                                        % ---------- both above thresh, we pick the one with lowest var
                if var1 < var2
                    selected = 'T1_fast';
                    Th_source_logic = 10 ; 
                    Th_source_string = selected ; 
                    disp('none of FP07 correlates with T | both seem OK | 1 is chosen')
                else
                    selected = 'T2_fast';
                    Th_source_logic = 20 ;
                    Th_source_string = selected ; 
                    disp('none of FP07 correlates with T | both seem OK | 2 is chosen')
                    
                end
        end


end




thermistor_source = Th_source_string;

x.params.thermistor_source = thermistor_source ;
x.params.Th_source_string = thermistor_source ;
x.params.Th_source_logic = Th_source_logic ;


%%
checkplot = 0 ;
if checkplot == 1
figure; 
plot(x.t_fast,x.T1_fast)
hold on; 
plot(x.t_fast,x.T2_fast)
plot(x.t_fast,T_gl_fast)
end



