function [simulation_results,settings] = seq_loop(settings)
% Control of loops through T1/TR/Ratio values etc.

%% Preallocate variables
if settings.HR_TR == 1
    % Only re-run simulations if variable TR is on, otherwise Monte Carlo
    NRepeats = settings.Repeats;
else
    NRepeats = 1;
end
N_imaging_RF = 0; N_prep_RF = 0; N_imaging_RF2 = 0;
IT1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s),length(settings.Velocities),length(settings.Diff_coeffs),1,NRepeats);
IT2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s),length(settings.Velocities),length(settings.Diff_coeffs),1,NRepeats);
if settings.T1Corr == 1
    IT3 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s),length(settings.Velocities),length(settings.Diff_coeffs),1,NRepeats);
end
if isfield(settings,'N_TRs')
    % This is used when a variable TR is simulated (emulating a variable heartrate)
    seq_TRs = zeros(settings.N_TRs,length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s),NRepeats);
end
FAs = zeros(size(settings.Tx_FA_map,3),settings.Scan_Size(2));
if strcmpi(settings.Scheme,'SA2RAGE')
    FA2s = zeros(size(settings.Tx_FA_map,3),settings.Scan_Size(2));
end
if strcmpi(settings.Scheme,'SA2RAGE') || strcmpi(settings.Scheme,'SatTFL') || strcmpi(settings.Scheme,'Sandwich') || strcmpi(settings.Scheme,'DREAM')
    PP_FAs = zeros(size(settings.Tx_FA_map,3),settings.Scan_Size(2));
end
Total_Energy = zeros(1,length(settings.Dynamic_Range));

%% Run loops
for Dynamic_Range_n = 1:size(settings.Dynamic_Range,2)
    if settings.verbose == 1 && any(Dynamic_Range_n == round(length(settings.Dynamic_Range).*(0.1:0.1:1)))
        disp(['Simulations ',num2str(100*Dynamic_Range_n./length(settings.Dynamic_Range),'%.2f'),'% complete.']);
    end
    
    % If synthetic data is being used then don't bother simulating values not in the mask
    if settings.UseSyntheticData == 0 || (settings.UseSyntheticData == 1 && settings.Long_Synthetic_Mask(Dynamic_Range_n))
        for T1_n = 1:length(settings.T1s)
            if settings.UseSyntheticData == 1
                T1 = settings.Long_Synthetic_T1s(Dynamic_Range_n);
            else
                T1 = settings.T1s(T1_n);
            end
            
            for B0_n = 1:length(settings.B0_Range_Hz)
                if settings.UseSyntheticData == 1 && all(settings.Dynamic_Range(:,Dynamic_Range_n) == 0)
                    break
                end
                
                for Mode_n = 1:settings.Modes
                    FAs(Mode_n,:) = Simulate_RF_Pulse(settings.Hz_per_Volt*abs(settings.Dynamic_Range(Mode_n,Dynamic_Range_n))*settings.RF_Pulse,settings.RF_Time,settings.B0_Range_Hz(B0_n)+settings.Slice_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                    if strcmpi(settings.Scheme,'SA2RAGE')
                        FA2s(Mode_n,:) = Simulate_RF_Pulse(settings.Hz_per_Volt*abs(settings.Dynamic_Range(Mode_n,Dynamic_Range_n))*settings.RF_Pulse2,settings.RF_Time,settings.B0_Range_Hz(B0_n)+settings.Slice_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                    end
                    if strcmpi(settings.Scheme,'SA2RAGE') || strcmpi(settings.Scheme,'SatTFL') || strcmpi(settings.Scheme,'Sandwich') || strcmpi(settings.Scheme,'DREAM')
                        PP_FAs(Mode_n,:) = Simulate_RF_Pulse(settings.Hz_per_Volt*abs(settings.Dynamic_Range(Mode_n,Dynamic_Range_n))*settings.PP_RF_Pulse,settings.PP_RF_Time,settings.B0_Range_Hz(B0_n)+settings.PP_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                    end
                end
                if settings.UseSyntheticData == 1
                    % Set RF phase if synthetic data is being used
                    RF_Phase = angle(settings.Dynamic_Range(:,Dynamic_Range_n));
                else
                    RF_Phase = 0;
                end
                for Flow_n = 1:length(settings.Velocities)
                    for Diff_n = 1:length(settings.Diff_coeffs)
                        Diff_co = settings.Diff_coeffs(Diff_n);
                        for Repeat_n = 1:NRepeats
                            
                            % Is this a voxel in the blood pool?
                            if settings.UseSyntheticData == 1 && settings.Long_Synthetic_T1s(Dynamic_Range_n) == settings.T1_blood
                                Ejection_Fraction = 0;
                                Velocity = settings.Velocities(Flow_n);
                                Angle = settings.Angles(Flow_n);
                            elseif settings.UseSyntheticData == 1
                                Ejection_Fraction = settings.Ejection_Fraction;
                                Velocity = 0;
                                Angle = 0;
                            else % else not synthetic body simulation
                                Ejection_Fraction = settings.Ejection_Fraction;
                                Velocity = settings.Velocities(Flow_n);
                                Angle = settings.Angles(Flow_n);
                            end
                            
                            % Define seq_TRs(:,Dynamic_Range_n,B0_n,T1_n,Repeat_n)
                            if settings.HR_TR == 1 && isfield(settings,'N_TRs') % Variable TRs simulating gating to a heart rate
                                HR_TR = zeros(1,settings.N_TRs);
                                for TR_n = 1:settings.N_TRs
                                    HR_TR(TR_n) = settings.TR + settings.HR_SD*settings.TR*randn; % Variable TRs
                                    if HR_TR(TR_n) < settings.HR_minTR % Enforce a minimum TR by skipping a TR if too short
                                        HR_TR(TR_n) = HR_TR(TR_n) + settings.TR + settings.HR_SD*settings.TR*randn;
                                    end
                                end
                                seq_TRs(:,Dynamic_Range_n,B0_n,T1_n,Repeat_n) = HR_TR;
                            elseif isfield(settings,'N_TRs')
                                seq_TRs(:,Dynamic_Range_n,B0_n,T1_n,Repeat_n) = settings.TR.*ones(1,settings.N_TRs);
                            end
                            
                            % Simulate two compartments for ejection fraction
                            if settings.Ejection_Fraction ~= 0
                                N_Compartments = 2;
                            else
                                N_Compartments = 1;
                            end
                            
                            settings = Check_Mag_Track_Flag(Dynamic_Range_n,T1_n,B0_n,Flow_n,Diff_n,settings);
                            Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),N_Compartments);
                            Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),N_Compartments);
                            Train3 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),N_Compartments);
                            for Compartment_n = 1:N_Compartments
                                if strcmpi(settings.Scheme,'SatTFL')
                                    [Train1(:,:,:,Compartment_n),Train2(:,:,:,Compartment_n),Cumulative_Time,N_imaging_RF,N_prep_RF,Mag_Track] = epg_sattfl(T1,FAs,PP_FAs,RF_Phase,Velocity,Angle,Diff_co,settings);
                                elseif strcmpi(settings.Scheme,'Sandwich')
                                    [Train1(:,:,:,Compartment_n),Train2(:,:,:,Compartment_n),Train3(:,:,:,Compartment_n),Cumulative_Time,N_imaging_RF,N_prep_RF,Mag_Track] = epg_sandwich(T1,FAs,PP_FAs,RF_Phase,Velocity,Angle,Diff_co,seq_TRs(:,Dynamic_Range_n,B0_n,T1_n,Repeat_n),Compartment_n,settings); % Simulate sequence
                                elseif strcmpi(settings.Scheme,'SA2RAGE')
                                    [Train1(:,:,:,Compartment_n),Train2(:,:,:,Compartment_n),Cumulative_Time,N_imaging_RF,N_imaging_RF2,N_prep_RF,Mag_Track] = epg_sa2rage(T1,FAs,FA2s,PP_FAs,RF_Phase,Velocity,Angle,Diff_co,settings);
                                elseif strcmpi(settings.Scheme,'AFI')
                                    [Train1(:,:,:,Compartment_n),Train2(:,:,:,Compartment_n),Cumulative_Time,N_imaging_RF,Mag_Track] = epg_afi(T1,FAs,RF_Phase,Velocity,Angle,Diff_co,settings);
                                elseif strcmpi(settings.Scheme,'DREAM')
                                    [Train1(:,:,:,Compartment_n),Train2(:,:,:,Compartment_n),Cumulative_Time,N_imaging_RF,N_prep_RF,Mag_Track] = epg_dream(T1,FAs,PP_FAs,RF_Phase,Velocity,Angle,Diff_co,settings);
                                elseif strcmpi(settings.Scheme,'GRE')
                                    [Train1(:,:,:,Compartment_n),Cumulative_Time,N_imaging_RF,seq_TRs(:,Dynamic_Range_n,B0_n,T1_n,Repeat_n),Mag_Track] = epg_gre(T1,FAs,RF_Phase,Velocity,Angle,Diff_co,settings); % Simulate sequence
                                else
                                    error('ABORTED: Scheme not recognised, please input either ''SatTFL'', ''Sandwich'', ''DREAM'', ''AFI'', ''SA2RAGE'' OR ''ALL''.')
                                end
                            end % End of dual compartment simulation for EF
                            
                            if settings.Ejection_Fraction ~= 0
                                % Mix compartments according to requested EF
                                IT1(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,1,Repeat_n) = ((1-Ejection_Fraction).*Train1(:,:,:,1)) + (Ejection_Fraction.*Train1(:,:,:,2));
                                IT2(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,1,Repeat_n) = ((1-Ejection_Fraction).*Train2(:,:,:,1)) + (Ejection_Fraction.*Train2(:,:,:,2));
                                IT3(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,1,Repeat_n) = ((1-Ejection_Fraction).*Train3(:,:,:,1)) + (Ejection_Fraction.*Train3(:,:,:,2));
                            else
                                IT1(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,1,Repeat_n) = Train1(:,:,:,1);
                                IT2(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,1,Repeat_n) = Train2(:,:,:,1);
                                IT3(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,1,Repeat_n) = Train3(:,:,:,1);
                            end
                        end
                        
                        if any(settings.Mag_Track_Flags == 1)
                            simulation_results.Mag_Track{settings.Mag_Track_Flags == 1} = Mag_Track;
                        end
                    end
                end
            end
        end
    end
    
    % Calculate energy
    Imaging_RF_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.RF_Pulse,settings.RF_Time); % Calculate Imaging RF pulse energy
    if strcmpi(settings.Scheme,'SA2RAGE')
        Imaging_RF2_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.RF_Pulse2,settings.RF_Time); % Calculate Imaging RF pulse energy
    else
        Imaging_RF2_Energy = 0;
    end
    if strcmpi(settings.Scheme,'SA2RAGE') || strcmpi(settings.Scheme,'SatTFL') || strcmpi(settings.Scheme,'Sandwich') || strcmpi(settings.Scheme,'DREAM')
        Preparation_RF_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.PP_RF_Pulse,settings.PP_RF_Time); % Calculate Preparation RF pulse energy
    else
        Preparation_RF_Energy = 0;
    end
    Total_Energy(Dynamic_Range_n) = N_imaging_RF*Imaging_RF_Energy + N_imaging_RF2*Imaging_RF2_Energy + N_prep_RF*Preparation_RF_Energy;
end
Average_10s_Power = 10*Total_Energy./Cumulative_Time;


%% Store results
simulation_results.IT1 = IT1;
simulation_results.IT2 = IT2;
if settings.T1Corr == 1
    simulation_results.IT3 = IT3;
end
simulation_results.N_imaging_RF = N_imaging_RF;
simulation_results.N_imaging_RF2 = N_imaging_RF2;
simulation_results.N_prep_RF = N_prep_RF;
simulation_results.Cumulative_Time = Cumulative_Time;
simulation_results.Total_Energy = Total_Energy;
if isfield(settings,'N_TRs')
simulation_results.seq_TRs = seq_TRs;
end
simulation_results.Average_10s_Power = Average_10s_Power;
end
