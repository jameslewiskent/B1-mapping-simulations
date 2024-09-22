function [Dynamic_Range,Tx_FA_map,Enc_Mat] = Calc_Tx(RF_Voltage,settings,Enc_Mat)
% Based on nominal RF pulse voltage, now calculate FA for tx

if strcmpi(settings.UseSyntheticData,'Duke')
SyntheticDuke = load(['Data',filesep,'SyntheticDuke.mat']); % Reads in if current folder is Masterscript
% Centre body in FOV
B1Tx = circshift(SyntheticDuke.B1Tx,15,1); % Transmit field 139 x 178 x 124 slices x 8 channels
B1Rx = circshift(SyntheticDuke.B1Rx,15,1);
PDimage = circshift(SyntheticDuke.PDimage,15,1);
mask = circshift(SyntheticDuke.mask,15,1);
rho = circshift(SyntheticDuke.rho,15,1);
sigma = circshift(SyntheticDuke.sigma,15,1);
clear SyntheticDuke
elseif strcmpi(settings.UseSyntheticData,'Phantom')
WaterPhantom = load(['Data',filesep,'WaterPhantom.mat']); % Reads in if current folder is Masterscript
B1Tx = WaterPhantom.PhantomMaps; settings.Syn_Slice = 1;
mask = WaterPhantom.Mask;
CP = [1.00000000000000 + 2.44929359829471e-16i,0.707106781186548 - 0.707106781186548i,6.12323399573677e-17 - 1.00000000000000i,-0.707106781186548 - 0.707106781186548i,-1.00000000000000 - 1.22464679914735e-16i,-0.707106781186548 + 0.707106781186548i,-1.83697019872103e-16 + 1.00000000000000i,0.707106781186547 + 0.707106781186548i];
B1Tx = bsxfun(@times,mask.*B1Tx,permute(conj(CP).*16.45e-8,[1 3 4 2])); % Apply CP mode, arbitary scaling factor
end

if ~exist('Enc_Mat','var')
[Enc_Mat] = Calc_Enc_Mat(settings.Enc_Scheme,settings.Modes); % Calculate Encoding Matrix
end

if settings.Modes == 1
    %B1Tx= sum(bsxfun(@times,B1Tx(:,:,settings.Syn_Slice,:),exp(1i*angle(conj(B1Tx(76,100,settings.Syn_Slice,:))))),4); % Aligned fields to give max B1 field in centre of heart voxel (i=73,j=97) 53 103
    B1Tx = sum(bsxfun(@times,B1Tx(:,:,settings.Syn_Slice,:),permute(Enc_Mat,[1,3,4,2])),4).*RF_Voltage/sqrt(50); % B1Tx in (T/V)
    Tx_FA_map = B1Tx.*settings.Gamma.*1e-3.*360; % Convert B1 map to FA (degrees) (this is complex)
else
    B1Tx_modes = zeros(size(B1Tx,1),size(B1Tx,2),settings.Modes); Tx_FA_map = zeros(size(B1Tx_modes));
    for mode = 1:settings.Modes
        B1Tx_modes(:,:,mode) = sum(bsxfun(@times,B1Tx(:,:,settings.Syn_Slice,:),permute(Enc_Mat(mode,:),[1,3,4,2])),4).*RF_Voltage/sqrt(50); % Sum of channels for Slice in (T)
        Tx_FA_map(:,:,mode) = B1Tx_modes(:,:,mode).*settings.Gamma.*1e-3.*360; % Convert B1 map to FA (degrees) (this is complex)
    end
end
    Dynamic_Range = permute(squeeze(reshape(Tx_FA_map./(settings.nom_FA.*(180/pi)),[],1,settings.Modes)),[2 1]);
end

