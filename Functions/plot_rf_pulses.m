function plot_rf_pulses(settings)
figure('color','w');
str = {};
if isfield(settings,'IT_RF_Pulse') && any(settings.IT_RF_Pulse ~= 0)
RF_Pulse_Size = size(settings.IT_RF_Pulse,2);  
RF_Pulse_Time = settings.IT_RF_Time;   
RF_Sample_Time = RF_Pulse_Time./RF_Pulse_Size;
% Pad out for plotting
Samples_To_Pad = 6e-3./RF_Sample_Time - RF_Pulse_Size;
RF_Pulse1 = abs(cat(2,zeros(1,ceil(Samples_To_Pad/2)),settings.IT_RF_Pulse,zeros(1,floor(Samples_To_Pad/2))));
plot((1:size(RF_Pulse1,2))*RF_Sample_Time*1e3,RF_Pulse1,'b'); hold on;
ylabel('Pulse Amplitude [Volts]');
xlabel('Pulse Time [ms]');
str = {'Nominal Excitation Pulse'};

[IT_mz] = Simulate_RF_Pulse_for_plotting(settings.IT_RF_Pulse,settings.IT_RF_Time,settings.Dynamic_Range,settings.Gamma,settings.Ref_Voltage);
else 
    RF_Pulse1 = 0;
end

if isfield(settings,'IT1_RF_Pulse') && any(settings.IT1_RF_Pulse ~= 0)
RF_Pulse_Size = size(settings.IT1_RF_Pulse,2);  
RF_Pulse_Time = settings.IT_RF_Time;   
RF_Sample_Time = RF_Pulse_Time./RF_Pulse_Size;
% Pad out for plotting
Samples_To_Pad = 6e-3./RF_Sample_Time - RF_Pulse_Size;
RF_Pulse2 = abs(cat(2,zeros(1,ceil(Samples_To_Pad/2)),settings.IT1_RF_Pulse,zeros(1,floor(Samples_To_Pad/2))));
plot((1:size(RF_Pulse2,2))*RF_Sample_Time*1e3,RF_Pulse2,'b'); hold on;
ylabel('Pulse Amplitude [Volts]');
xlabel('Pulse Time [ms]');
str = {'Nominal Excitation Pulse 1'};

[IT1_mz] = Simulate_RF_Pulse_for_plotting(settings.IT1_RF_Pulse,settings.IT_RF_Time,settings.Dynamic_Range,settings.Gamma,settings.Ref_Voltage);
else 
    RF_Pulse2 = 0;
end

if isfield(settings,'IT2_RF_Pulse') && any(settings.IT2_RF_Pulse ~= 0)
RF_Pulse_Size = size(settings.IT2_RF_Pulse,2);  
RF_Pulse_Time = settings.IT_RF_Time;   
RF_Sample_Time = RF_Pulse_Time./RF_Pulse_Size;
% Pad out for plotting
Samples_To_Pad = 6e-3./RF_Sample_Time - RF_Pulse_Size;
RF_Pulse3 = abs(cat(2,zeros(1,ceil(Samples_To_Pad/2)),settings.IT2_RF_Pulse,zeros(1,floor(Samples_To_Pad/2))));
plot((1:size(RF_Pulse3,2))*RF_Sample_Time*1e3,RF_Pulse3,'b'); hold on;
ylabel('Pulse Amplitude [Volts]');
xlabel('Pulse Time [ms]');
str = {'Nominal Excitation Pulse 2'};

[IT2_mz] = Simulate_RF_Pulse_for_plotting(settings.IT2_RF_Pulse,settings.IT_RF_Time,settings.Dynamic_Range,settings.Gamma,settings.Ref_Voltage);
else 
    RF_Pulse3 = 0;
end

if isfield(settings,'PP_RF_Pulse') && any(settings.PP_RF_Pulse ~= 0)
RF_Pulse_Size = size(settings.PP_RF_Pulse,2);  
RF_Pulse_Time = settings.PP_RF_Time;   
RF_Sample_Time = RF_Pulse_Time./RF_Pulse_Size;
% Pad out for plotting
Samples_To_Pad = 6e-3./RF_Sample_Time - RF_Pulse_Size;
RF_Pulse4 = abs(cat(2,zeros(1,ceil(Samples_To_Pad/2)),settings.PP_RF_Pulse,zeros(1,floor(Samples_To_Pad/2))));
plot((1:size(RF_Pulse4,2))*RF_Sample_Time*1e3,RF_Pulse4,'r'); hold on;
ylabel('Pulse Amplitude [Volts]');
xlabel('Pulse Time [ms]');
str = [str,'Nominal Preparation Pulse'];

[PP_mz] = Simulate_RF_Pulse_for_plotting(settings.PP_RF_Pulse,settings.PP_RF_Time,settings.Dynamic_Range,settings.Gamma,settings.Ref_Voltage);
else 
    RF_Pulse4 = 0;
end

if isfield(settings,'RF_Pulse') && any(settings.RF_Pulse ~= 0)
RF_Pulse_Size = size(settings.RF_Pulse,2);  
RF_Pulse_Time = settings.RF_Time;   
RF_Sample_Time = RF_Pulse_Time./RF_Pulse_Size;
% Pad out for plotting
Samples_To_Pad = 6e-3./RF_Sample_Time - RF_Pulse_Size;
RF_Pulse5 = abs(cat(2,zeros(1,ceil(Samples_To_Pad/2)),settings.RF_Pulse,zeros(1,floor(Samples_To_Pad/2))));
plot((1:size(RF_Pulse5,2))*RF_Sample_Time*1e3,RF_Pulse5,'r'); hold on;
ylabel('Pulse Amplitude [Volts]');
xlabel('Pulse Time [ms]');
str = [str,'Nominal RF Pulse'];

[RF_mz] = Simulate_RF_Pulse_for_plotting(settings.RF_Pulse,settings.RF_Time,settings.Dynamic_Range,settings.Gamma,settings.Ref_Voltage);
else 
    RF_Pulse5 = 0;
end

Vpeak = max([RF_Pulse1,RF_Pulse2,RF_Pulse3,RF_Pulse4,RF_Pulse5]);
disp(['Peak voltage = ',num2str(Vpeak),' V.']);

ylim([min([RF_Pulse1,RF_Pulse2,RF_Pulse3,RF_Pulse4,RF_Pulse5]),  Vpeak + 10])
title('RF Pulse Envelope');
legend(str)
hold off;

% Plot simulation of RF pulse
figure('color','w'); tiledlayout('flow','TileSpacing','compact','Padding','none')
B0_Hz = -5000:10:5000;
if exist('IT_mz','var')
nexttile; imagesc(settings.Dynamic_Range,B0_Hz,IT_mz)
ylabel('B_0 (kHz)'); 
xlabel('B_1 Dynamic Range'); 
title('Imaging Excitation');
end
if exist('IT1_mz','var')
nexttile; imagesc(settings.Dynamic_Range,B0_Hz,IT1_mz)
ylabel('B_0 (kHz)'); 
xlabel('B_1 Dynamic Range'); 
title('Imaging 1 Excitation');
end
if exist('IT2_mz','var')
nexttile; imagesc(settings.Dynamic_Range,B0_Hz,IT2_mz)
ylabel('B_0 (kHz)'); 
xlabel('B_1 Dynamic Range'); 
title('Imaging 2 Excitation');
end
if exist('PP_mz','var')
nexttile; imagesc(settings.Dynamic_Range,B0_Hz,PP_mz)
ylabel('B_0 (kHz)'); 
xlabel('B_1 Dynamic Range'); 
title('Preparation Pulse');
end
if exist('RF_mz','var')
nexttile; imagesc(settings.Dynamic_Range,B0_Hz,RF_mz)
ylabel('B_0 (kHz)'); 
xlabel('B_1 Dynamic Range'); 
title('RF Pulse');
end

end

