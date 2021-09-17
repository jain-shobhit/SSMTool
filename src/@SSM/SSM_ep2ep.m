function varargout = SSM_ep2ep(obj,oid,run,lab,parName,parRange,outdof,varargin)
% SSM_EP2EP This function performs continuation of equilibrium points of
% slow dynamics. Each equilibrium point corresponds to a periodic orbit in
% the regular time dynamics. The continuation here starts from a saved
% solution.
%
% FRCIRS = SSM_EP2EP(OBJ,OID,RUN,LAB,PARNAME,PARRANGE,OUTDOF,VARARGIN)
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
% See also: SSM_ISOL2EP SSM_BP2EP SSM_CONT_EP

FRC = obj.SSM_cont_ep('ep',oid,run,lab,parName,parRange,outdof,varargin{:});
varargout{1} = FRC;
end