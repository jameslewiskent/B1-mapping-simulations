function [Dynamic_Range,Tx_FA_map,Enc_Mat,W_Mat] = Calc_Tx(RF_Voltage,settings)
% Based on nominal RF pulse voltage, now calculate FA for tx

SyntheticDuke = load('Data\SyntheticDuke.mat'); % Reads in if current folder is Masterscript

if settings.MTx == 0
    settings.Modes = 'NA';
    % Set up B1 Tx Field
    B1Tx = SyntheticDuke.B1Tx; % Transmit field 139 x 178 x 124 slices x 8 channels
    B1Tx_sumv = sum(bsxfun(@times,B1Tx(:,:,settings.Syn_Slice,:),exp(1i*angle(conj(B1Tx(76,100,Syn_Slice,:))))),4); % Aligned fields to give max B1 field in centre of heart voxel (i=73,j=97) 53 103
    %B1Tx_sumv =sum(bsxfun(@times,B1Tx(:,:,Syn_Slice,:),permute(exp(1i.*[0,45,90,135,180,225,270,315]),[1,3,4,2])),4); % Sum of channels for Slice 60
    B1Tx_v = B1Tx_sumv/sqrt(50); % B1Tx in (T/V)
    B1Tx = B1Tx_v.*RF_Voltage.*1e4; % Complex Magnitude of B1Tx in (G)
    Tx_FA_map = abs(B1Tx).*1e-4.*settings.Gamma.*1e-3.*360; % Convert B1 map to nominal FA
    
    imagesc(abs(Tx_FA_map(1:114,:))) % Plot B1Tx FA maps for each mode
    title('Simulated Flip Angle Map for a Synthetic Body Model')
    axis image off
    cb = colorbar;
    cb.Label.String = ['Flip angle, (',char(176),')'];
elseif settings.MTx == 1
    % Calculate Encoding Matrix
    [Enc_Mat,W_Mat] = Calc_Enc_Mat(settings.Enc_Scheme,settings.Modes);
    
    % Set up B1 Tx Field
    B1Tx = SyntheticDuke.B1Tx; % 139 x 178 x 124 slices x 8 channels
    B1Tx_modes = zeros(size(B1Tx,1),size(B1Tx,2),settings.Modes); Tx_FA_map = zeros(size(B1Tx_modes));
    for mode = 1:settings.Modes
        B1Tx_modes(:,:,mode) = sum(bsxfun(@times,B1Tx(:,:,settings.Syn_Slice,:),permute(Enc_Mat(mode,:),[1,3,4,2])),4).*RF_Voltage/sqrt(50); % Sum of channels for Slice in (T)
        Tx_FA_map(:,:,mode) = (B1Tx_modes(:,:,mode)).*settings.Gamma.*1e-3.*360; % Convert B1 map to FA (degrees) (this is complex)
    end
    
    figure();
    imagesc(imtile(abs(Tx_FA_map),'Gridsize',[1 settings.Modes])) % Plot B1Tx FA maps for each mode
    title('Simulated Transmit Modes for Synthetic Body Model')
    axis image
    set(gca,'YTick',[]); set(gca,'XTick',[]);
    xlabel('Transmit Mode')
    xticks(((1:settings.Modes)).*size(Tx_FA_map,2) - size(Tx_FA_map,2)/2);
    xticklabels(1:settings.Modes);

    Dynamic_Range = permute(squeeze(reshape(Tx_FA_map./(settings.nomIT_FA.*(180/pi)),[],1,8)),[2 1]);
    
    
    %Tx_Channel_FA_map = squeeze(abs(B1Tx(:,:,Syn_Slice,:)).*RF_Voltage.*Gamma.*1e-3.*360)/sqrt(50); % B1Tx in (T/V); % Convert B1 map to FA (degrees)
    %save('Tx_Channel_FA_map.mat','Tx_Channel_FA_map');
end
end

