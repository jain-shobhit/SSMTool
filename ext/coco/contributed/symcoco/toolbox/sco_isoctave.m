function isoct=sco_isoctave()
%% check if octave is running instead of Matlab
%
% $Id: dde_isoctave.m 212 2017-07-09 22:31:00Z jansieber $
%
if exist('OCTAVE_VERSION','builtin')
    isoct=true;
else
    isoct=false;
end
end
