function SSM_TR2tor(obj,oid,run,lab,parName,parRange,outdof,varargin)
% SSM_TR2tor This function performs continuation of 2D tori of slow
% dynamics. Each 2D tori corresponds to a 3D torus (quasi-periodic)
% response in regular time dynamics. The continuation here starts by
% switching to branch of tori at torus (Neimark-Sacker) bifurcation point
% found in the continuation of periodic orbits. 
%
% FRCIRS = SSM_TR2TOR(OBJ,OID,RUN,LAB,PARNAME,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of current continuation
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the
%           label of a torus point
% parName:  freq-amp/amp/freq continuation parameter. In the case of
%           freq-amp, both two parameters are free to change and hence
%           varrho is fixed. In contrast, varrho is fixed if there is only
%           one continuation parameter
% parRange: continuation domain of parameter, which is given in the form of
%           {[om1,om2],[f1,f2]} for freq-amp, [f1,f2] for amp and [om1,om2]
%           for freq. It is noted that the domain [om1,om2] should contain
%           the natural frequency with index 1 in the mFreq vector
% outdof:   output for dof in physical domain
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_PO2TR, SSM_TOR2PTOR

obj.SSM_cont_tor('TR',oid,run,lab,parName,parRange,outdof,varargin{:});
end
