
%u = 1 ; 

%%
disp('L2_section_direction_position') 
L2_section_direction_position

%
disp('L2_make_section_id') 
L2_make_section_id
disp(section_id) 
if dirfu  < 0
    disp([num2str(length(list_L0)),'  ',num2str(ppp),'  ',num2str(u),'/',num2str(lu),' traits  /'])
else
    disp([num2str(length(list_L0)),'  ',num2str(ppp),'  ',num2str(u),'/',num2str(lu),' traits  \'])
end

%
disp('L2_thermistor_source') 
L2_thermistor_source
disp(['master for thermistor:  ',thermistor_source,'   ',num2str(Th_source_logic)])
pause(2)

%
disp('L2_temperature_source') 
L2_temperature_source
disp(['source for dissipation: ',temperature_source,'   ',num2str(T_source_logic)])
pause(2)

%
disp('L2_hipass_lopass') 
L2_hipass_lopass

%
disp('L2_despiking') 
L2_despiking

%
disp('L2_fft_parameters') 
L2_fft_parameters %L3_fft_parameters

%%
    if Ntimes_N_fft/ls <= 0.5
    
        disp('L2_hipass_shears') 
        L2_hipass_shears %L3_hipass_shears
        
        disp('L3_spectra_sh') 
        L3_spectra_sh_FK25
    
            if get_diss_odas_worked == 1             
                %
                disp('L3_spectra_th 1') 
                x0   = x.gradT1(fu) ;   
                T_dT = x.T1_dT1(fu) ; 
                T_string = 'T1_dT1' ; % T_string (*): name of the FP07 channel 'T1_dT1' or 'T2_dT2'
                L3_spectra_th
                dissT1 = dissT ; 
                
                %
                disp('L3_spectra_th 2') 
                x0 = x.gradT2(fu) ;   
                T_dT = x.T2_dT2(fu) ;
                T_string = 'T2_dT2' ; % T_string (*): name of the FP07 channel 'T1_dT1' or 'T2_dT2'
                L3_spectra_th
                dissT2 = dissT ; 
                %
                %structT1 = dissT1 ; 
                %structT2 = dissT2 ; 
                
                        
                champs = fieldnames(dissT1);
                a_tronquer = false;             
                for i = 1:length(champs)
                    champ = champs{i};                  
                    if isfield(dissT2, champ)
                        a = dissT1.(champ);
                        b = dissT2.(champ);                        
                        if isvector(a) && isvector(b) && ~isscalar(a) && ~isscalar(b) ...
                                && ~isempty(a) && ~isempty(b)
                            a_tronquer = true;
                            n = min(length(a), length(b));
                            dissT1.(champ) = a(1:n);
                            dissT2.(champ) = b(1:n);
                        end
                    end
                end
                
                if ~a_tronquer
                    return
                end
                
                
                %%
                disp('L4_shears') 
                L4_shears
                %
                disp('L4_shears_QC') 
                L4_shears_QC
                
                %%
                disp('L4_FP07s') 
                L4_FP07s
                %
                disp('L4_FP07s_QC') 
                L4_FP07s_QC
                
                %%
                disp('make structure') 
                make_structure

                %---------
                outm_pth = [processing_out_path,'L4/mat/'] ; %
                disp(['saving in L4/mat   ',FID_QC])
                save([outm_pth,FID_QC,'.mat'], 'mr','-v7.3'); 
                %---------
                pause(3)
                clc

            
            else
                disp('*** get_diss_odas() skipped due to error ***')
            end
    
    else
        disp('*** too short for FFT ***')
    end