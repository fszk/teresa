%% temperature source for dissipation
%
% 11 or 22 : default T is T_glider, and both Thermistor are suitable, best
% is chosen
%
% 1 or 2  : default T is T_glider, and one Thermistor is suitable 1 or 2
%
%  0 : no correlation  with T_glider and possibly both thermistor have
%  problems, default T is still T_glider
%
%  100,200 : no correlation  with T_glider, FP07 are below thresh, lower variance is
%  chosen, default T is still T_glider
%
%  10,20 : no correlation  with T_glider, FP07 are above thresh, lower variance is
%  chosen, default T is chosen FP07


if Th_source_logic == 11 || Th_source_logic == 22
    temperature_for_dissipation = T_gl_fast ; 
    temperature_source = 'T_gl_fast' ; 
    T_source_string = temperature_source ;
    T_source_logic = Th_source_logic ; 
    

elseif Th_source_logic == 1 || Th_source_logic == 2
    temperature_for_dissipation = T_gl_fast ; 
    temperature_source = 'T_gl_fast' ; 
    T_source_string = temperature_source ;
    T_source_logic = Th_source_logic ; 

elseif  Th_source_logic == 0
    temperature_for_dissipation = T_gl_fast ; 
    temperature_source = 'T_gl_fast' ; 
    T_source_string = temperature_source ;
    T_source_logic = Th_source_logic ; 

elseif Th_source_logic == 100 || Th_source_logic == 200
    temperature_for_dissipation = T_gl_fast ; 
    temperature_source = 'T_gl_fast' ; 
    T_source_string = temperature_source ;
    T_source_logic = Th_source_logic ; 

elseif Th_source_logic == 10 
    temperature_source = thermistor_source ; 
    temperature_for_dissipation = x.T1_fast ; 
    T_source_string = temperature_source ;
    T_source_logic = Th_source_logic ; 

elseif Th_source_logic == 20 
    temperature_source = thermistor_source ; 
    temperature_for_dissipation = x.T2_fast ; 
    T_source_string = temperature_source ;
    T_source_logic = Th_source_logic ; 

end

x.params.temperature_source = temperature_source ;
x.params.T_source_string = T_source_string ;
x.params.T_source_logic = T_source_logic ;


% 
% %% NEW
% 
% % --- Sélection de la sonde de température basée sur la variance ---
% var_thresh = 1e-5;
% 
% % Calcul des variances (omit NaNs)
% var1 = var(x.T1_fast, 'omitnan');
% var2 = var(x.T2_fast, 'omitnan');
% 
% % Logique de sélection
% if var1 < var_thresh && var2 < var_thresh
%     % Les deux variances sont trop faibles : échec ou variable de repli
%     %error('Les variances de T1_fast et T2_fast sont toutes deux inférieures au seuil ; aucun choix possible.');
%     selected = x.params.temperature_source ;
% elseif var1 < var_thresh
%     selected = 'T2_fast';
%      T_source_logic = 2 ; 
% elseif var2 < var_thresh
%     selected = 'T1_fast';
%      T_source_logic = 1 ; 
% else
%     % Les deux sont au-dessus du seuil : on choisi la moins bruitée (variance min)
%     if var1 < var2
%         selected = 'T1_fast';
%         T_source_logic = 1 ; 
%     else
%         selected = 'T2_fast';
%         T_source_logic = 2 ; 
%     end
% end
% 
% % Enregistrement de la source choisie
% x.params.temperature_source = selected;
% 
% % Assignation de la série de température retenue
% switch selected
%     case 'T1_fast'
%         temperature_for_dissipation = x.T1_fast;
%     case 'T2_fast'
%         temperature_for_dissipation = x.T2_fast;
% end
% 
% 
% T_source_string = x.params.temperature_source ;
% 
% T_source_string
% T_source_logic
