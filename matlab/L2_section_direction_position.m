%% Get indexes _fast from the individual section
% identifies if down/or/up-ward 
% get lon lat
% get navigation direction (simpler than heading)
%%
fu = find(sections == us(u)) ;


fu1 = fu(1)      ;                  
fue = fu(end)    ;                 
dirfu = mean(sections_dire(fu)) ;
if dirfu > 0 
    dirfu_str = 'down';
else
    dirfu_str = 'upwd' ;
end
ls = length(fu) ;                   

pres_delta = max(x.P_fast(fu),[],'omitnan') - min(x.P_fast(fu),[],'omitnan') ;

%%
% %% OLD, before nesting glider data at conversion
% section_pres = glider_pres_mr(fu) ;
% fpresfin = find(isfinite(section_pres)==1) ; 
% 
% section_lat = glider_lat_mr(fu) ;
% section_lon = glider_lon_mr(fu) ;
% 
% glider_found_section = length(fpresfin) ; 

%%

section_pres = x.P_fast(fu) ; 
fpresfin = find(isfinite(section_pres)==1) ; 

section_lat = lat_gl_fast ; 
section_lon = lon_gl_fast ; 

glider_found_section = length(fpresfin) ; 


%%
navig_dir = sign(diff([section_lon(1),section_lon(end)]));
if navig_dir > 0 
    navig_dir_str = 'Eastward' ;
else
    navig_dir_str = 'Westward' ;
end
disp(navig_dir_str)

