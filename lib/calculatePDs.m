function [tuning] = calculatePDs(velInput,FR,doBoot,bootreps)
% CALCULATETUNING Calculates glm given velocity of the hand and the firing
% rate. Bootstrapping is optional, but necessary to calculate confidence
% intervals (doBoot=true for bootstrapping). Bootreps gives number of
% bootstrap samples (irrelevant if not bootstrapping).

bootfunc = @(X,y) glmBootfunc(X,y,'Distribution','poisson');

if doBoot
    % get bootstrapped coefficients
    bootCoeffs = bootstrp(bootreps,@(X,y) bootfunc(X,y),velInput,FR);
    
    % get mean coefficients
    mean_coeffs = mean(bootCoeffs);
    
    % get PD
    velDir = atan2(mean_coeffs(3),mean_coeffs(2));
    
    % get bootstrapped directions
    bootDirs = atan2(bootCoeffs(:,3),bootCoeffs(:,2));
    
    % center bootDirs on mean dir
    centeredDirs = bootDirs-velDir;
    while(sum(centeredDirs<-pi))
        centeredDirs(centeredDirs<-pi) = centeredDirs(centeredDirs<-pi)+2*pi;
    end
    while(sum(centeredDirs>pi))
        centeredDirs(centeredDirs>pi) = centeredDirs(centeredDirs>pi)-2*pi;
    end
    
    velDirCI = prctile(centeredDirs,[2.5 97.5])+velDir;
    
    % could add moddepth here, but don't care
else
    mean_coeffs = bootfunc(velInput,FR);
    velDir = atan2(mean_coeffs(3),mean_coeffs(2));
    velDirCI = [NaN NaN];
end

tuning = table(velDir,velDirCI,'VariableNames',{'velDir','velDirCI'});

end

