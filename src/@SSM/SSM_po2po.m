function SSM_po2po(obj,oid,run,lab,parName,parRange,outdof,varargin)
% SSM_PO2PO This function performs continuation of periodic orbits of slow
% dynamics. Each periodic orbit corresponds to a torus (quasi-periodic)
% response in regular time dynamics. The continuation here starts a saved
% solution
%
% FRCIRS = SSM_PO2PO(OBJ,OID,RUN,LAB,PARNAME,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of current continuation
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution
% parName:  amp/freq continuation parameter
% parRange: continuation domain of parameter, which should be near the
%           value of natural frequency with index 1 in the mFreq if continuation
%           parameter is freq
% outdof:   output for dof in physical domain
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_HB2PO, SSM_BP2PO

obj.SSM_cont_po('po',oid,run,lab,[],parName,parRange,outdof,varargin{:});
end
