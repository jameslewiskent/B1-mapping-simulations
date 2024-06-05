function [Dynamic_Range,Tx_FA_map,Enc_Mat] = Calc_Tx(RF_Voltage,settings)
% Based on nominal RF pulse voltage, now calculate FA for tx

SyntheticDuke = load(['Data',filesep,'SyntheticDuke.mat']); % Reads in if current folder is Masterscript
[Enc_Mat] = Calc_Enc_Mat(settings.Enc_Scheme,settings.Modes); % Calculate Encoding Matrix

if settings.Modes == 1
    % Set up B1 Tx Field
    B1Tx = SyntheticDuke.B1Tx; % Transmit field 139 x 178 x 124 slices x 8 channels
    %B1Tx_sumv = sum(bsxfun(@times,B1Tx(:,:,settings.Syn_Slice,:),exp(1i*angle(conj(B1Tx(76,100,settings.Syn_Slice,:))))),4); % Aligned fields to give max B1 field in centre of heart voxel (i=73,j=97) 53 103
    B1Tx_CP = sum(bsxfun(@times,B1Tx(:,:,settings.Syn_Slice,:),permute(Enc_Mat,[1,3,4,2])),4).*RF_Voltage/sqrt(50); % B1Tx in (T/V)
    Tx_FA_map = B1Tx_CP.*settings.Gamma.*1e-3.*360; % Convert B1 map to FA (degrees) (this is complex)
else
    % Set up B1 Tx Field
    B1Tx = SyntheticDuke.B1Tx; % 139 x 178 x 124 slices x 8 channels
    B1Tx_modes = zeros(size(B1Tx,1),size(B1Tx,2),settings.Modes); Tx_FA_map = zeros(size(B1Tx_modes));
    for mode = 1:settings.Modes
        B1Tx_modes(:,:,mode) = sum(bsxfun(@times,B1Tx(:,:,settings.Syn_Slice,:),permute(Enc_Mat(mode,:),[1,3,4,2])),4).*RF_Voltage/sqrt(50); % Sum of channels for Slice in (T)
        Tx_FA_map(:,:,mode) = B1Tx_modes(:,:,mode).*settings.Gamma.*1e-3.*360; % Convert B1 map to FA (degrees) (this is complex)
    end

    %Tx_Channel_FA_map = squeeze(abs(B1Tx(:,:,Syn_Slice,:)).*RF_Voltage.*Gamma.*1e-3.*360)/sqrt(50); % B1Tx in (T/V); % Convert B1 map to FA (degrees)
    %save('Tx_Channel_FA_map.mat','Tx_Channel_FA_map');
end
    Dynamic_Range = permute(squeeze(reshape(Tx_FA_map./(settings.nom_FA.*(180/pi)),[],1,settings.Modes)),[2 1]);
end

