function [resModes,varargout]  = detect_resonant_modes(resLambda,lambda,tol)
% DETECT_RESONANT_MODES This function selects modes for resonant master subpace
% for a given eigenvalue in an eigen-spectrum. The selection is based on
% the closeness of the ratio of eigenvalues to integers
% 
% [RESMODES,VARARGOUT] = DETECT_RESONANT_MODES(RESLAMBDA, LAMBDA, TOL)
%
% resLambda: master eigenvalue
% lambda:    spectrum of eigenvalues
% resModes:  subspace resonant with master eigenvalue
% varargout: internal resonance relation vector

absResLambda = abs(imag(resLambda));
absLambda    = abs(imag(lambda));
numLambda = numel(lambda);
ratio     = zeros(numLambda,1);
mFreqs    = ones(numLambda,1);

for k=1:numLambda
    if absLambda(k)>absResLambda
        ratio(k)  = absLambda(k)/absResLambda;
        mFreqs(k) = round(ratio(k)); 
    else
        ratio(k)  = absResLambda/absLambda(k);
        mFreqs(k) = 1/round(ratio(k));
    end
end

% kill ones larger than 5
ratio(ratio>5) = inf;

flag  = abs(ratio-round(ratio))<tol;
modes = 1:numLambda;
resModes = modes(flag);
varargout{1} = mFreqs(flag);

if numel(resModes)>2
    disp(['The master subspace has internal resonances: [',num2str(mFreqs(flag).'),']']);
end
end