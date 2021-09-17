function SSM_BP2po(obj,oid,run,lab,parName,parRange,outdof,varargin)
% SSM_BP2PO This function performs continuation of periodic orbits of slow
% dynamics. Each periodic orbit corresponds to a torus (quasi-periodic)
% response in regular time dynamics. The continuation here follows the
% secondary branch passing through a branch point.
%
% FRCIRS = SSM_BP2PO(OBJ,OID,RUN,LAB,PARNAME,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of current continuation
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the
%           label of a branch point
% parName:  amp/freq continuation parameter
% parRange: continuation domain of parameter, which should be near the
%           value of natural frequency with index 1 in the mFreq if continuation
%           parameter is freq
% outdof:   output for dof in physical domain
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_HB2PO, SSM_PO2PO

obj.SSM_cont_po('BP',oid,run,lab,[],parName,parRange,outdof,varargin{:});
end
