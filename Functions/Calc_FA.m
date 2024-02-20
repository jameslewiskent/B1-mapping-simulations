function [Measured_FA] = Calc_FA(settings,Max_Val_IT1,Max_Val_IT2,Max_Val_IT3)
% Calculate flip angle for various schemes
if strcmpi(settings.Scheme,'DREAM')
    %if strcmp(settings.Echo_Order,'STEFirst')
    Measured_FA = atan(sqrt(abs(2.*(Max_Val_IT1./Max_Val_IT2)))); % Eq. [5] Nehrke, K. et al. MRM. 2012. 68:1517-1526.
    %elseif strcmp(settings.Echo_Order,'FIDFirst')
    %Measured_FA = atan(sqrt(abs(2.*(Max_Val_IT1./Max_Val_IT2)))); % Eq. [5] Nehrke, K. et al. MRM. 2012. 68:1517-1526.
    %end
    
elseif strcmpi(settings.Scheme,'AFI')
    n = settings.TR2/settings.TR1; % TR2/TR1 for AFI Scheme (E.g. n = 100 ms / 20 ms = 5)
    r = (Max_Val_IT2./Max_Val_IT1);
    Measured_FA = acos( ((r.*n) - 1) ./ (n - r) ); % Eq. [6] Yarnykh VL. AFI. Magn Reson Med 2007;57:192?200. https://doi.org/10.1002/mrm.21120.
    
elseif strcmpi(settings.Scheme,'SA2RAGE')
    if settings.Use_Previous_Lookup == 0
        Generate_Lookup_Table(Max_Val_IT1./Max_Val_IT2,settings);
    end
    Measured_FA = Apply_Lookup_Table(Max_Val_IT1./Max_Val_IT2,settings);
    
elseif strcmpi(settings.Scheme,'SatTFL')
    if settings.Lookup_T1 ~= 0
        if settings.Use_Previous_Lookup == 0
            Generate_Lookup_Table(Max_Val_IT2./Max_Val_IT1,settings);
        end
        Measured_FA = Apply_Lookup_Table(Max_Val_IT2./Max_Val_IT1,settings);
    else
        disp('Calculating flip angle using arcosine instead of lookup table.')
        Measured_FA = acos(Max_Val_IT2./Max_Val_IT1); % Eq. 2 from Chung, S. et al. MRM. 2010. 64(2):439-446.
    end
elseif  strcmpi(settings.Scheme,'Sandwich')
    Image_Ratio = Max_Val_IT2./Max_Val_IT1;
    if settings.T1Corr == 1
    Image_Ratio = (Max_Val_IT2 + 0.5.*(Max_Val_IT2 - Max_Val_IT3))./Max_Val_IT1;  
    end
    if settings.Lookup_T1 ~= 0
        if settings.Use_Previous_Lookup == 0
            Generate_Lookup_Table(Image_Ratio,settings);
        end
        Measured_FA = Apply_Lookup_Table(Image_Ratio,settings);
    else
        disp('Calculating flip angle using arcosine instead of lookup table.')
        Measured_FA = acos(Image_Ratio); % Eq. 2 from Chung, S. et al. MRM. 2010. 64(2):439-446.
    end
elseif strcmpi(settings.Scheme,'GRE')
        Measured_FA = Max_Val_IT1;
end

if ~strcmpi(settings.Scheme,'GRE')
    % Don't convert relative maps to FA (and keep complex)
    Measured_FA = (180/pi)*real(Measured_FA);
end


% Image_Maps = Image_Maps.*exp(1i*angle(conj(Image_Maps(indi,indj,Mode_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n))))

end

