function spectraSmooth = SmoothSpectrum (spectra, parameters)
% function spectraSmooth = SmoothSpectrum (spectra, parameters)
% 
% The function smoothes 'spectra' according to the smoothing span and
% method in 'parameters'.
% Smoothing span of 0 denotes no smoothing. Span of 1 is effectively no
% smoothing (use every sample to smooth...)
if ((parameters.nSmoothingSpan)&&(parameters.nSmoothingSpan~=1))
    for k=size(spectra,2):-1:1
        spectraSmooth(:,k)=smooth(spectra(:,k),...
            parameters.nSmoothingSpan,parameters.sSmoothingMethod);
    end
    % A hack to make sure we do not get <0 Fsm, resulting in complex logs
    spectraSmooth(spectraSmooth<0)=min(spectraSmooth(spectraSmooth>0));
end
