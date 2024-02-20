function [Enc_Mat] = Calc_Enc_Mat(Enc_Scheme,Modes)

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
                     
        elseif strcmp(Enc_Scheme,'OneOFF')
            if channel ~= mode
                Enc_Mat(mode,channel) = exp((2*pi*1i*(channel-1))/Channels);
            end
            
        elseif strcmp(Enc_Scheme(1:3),'One')
            if channel ~= mode
                Enc_Mat(mode,channel) = exp((2*pi*1i*(channel-1))/Channels);
            else
                Enc_Mat(mode,channel) = sscanf(string(Enc_Scheme),'One%d')/100;
            end
        else
            error(['Encoding Matrix for ',Enc_Scheme,' not calculated succesfully. Check spelling of supplied Enc_Scheme is one of the available options.']);
            
        end
    end
end

end

