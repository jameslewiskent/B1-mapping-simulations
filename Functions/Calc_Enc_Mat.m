function [Enc_Mat,W_Mat] = Calc_Enc_Mat(Enc_Scheme,Modes)

Channels = 8; % Channels
Enc_Mat = zeros(Modes,Channels); % Pre-allocate Encoding Matrix
for mode = 1:Modes
    for channel = 1:Channels
        
        if strcmp(Enc_Scheme,'FE')
            Enc_Mat(mode,channel) = exp((2*pi*1i*(mode-1)*channel)/Modes);
            
        elseif strcmp(Enc_Scheme,'Invert')
            if channel == mode
                Enc_Mat(mode,channel) = -exp((2*pi*1i*(channel-1))/Channels);
            elseif channel ~= mode
                Enc_Mat(mode,channel) = exp((2*pi*1i*(channel-1))/Channels);
            end
            
        elseif strcmp(Enc_Scheme,'Indiv')
            Enc_Mat = eye(Modes,Channels); % Execute one transmit channel at a time
            
            if channel == mode
                Enc_Mat(mode,channel) =  exp((2*pi*1i*(channel-1))/Channels);
            end
            
        elseif strcmp(Enc_Scheme,'OneOFF')
            if channel ~= mode
                Enc_Mat(mode,channel) = exp((2*pi*1i*(channel-1))/Channels);
            end
        else
            error(['Encoding Matrix for ',Enc_Scheme,' not calculated succesfully. Check spelling of supplied Enc_Scheme is one of the available options.']);
            
        end
    end
end
W_Mat = eye(Modes,Modes); % Weighting matrix set to indentity for now

end

