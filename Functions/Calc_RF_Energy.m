function [RF_Energy,peak_RF] = Calc_RF_Energy(RF_Pulse,RF_Pulse_Time)
% Calculate RF pulse energy.
% Sinc and HS8 have matched bandwidths of 3kHz.

        Ohms = 50; 
        RF_Sample_Time = RF_Pulse_Time./size(RF_Pulse,2);
        %peak_RF = 
        RF_Energy = sum((RF_Pulse.*conj(RF_Pulse)).*RF_Sample_Time,2)/Ohms;
end

