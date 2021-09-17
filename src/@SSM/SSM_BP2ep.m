function varargout = SSM_BP2ep(obj,oid,run,lab,parName,parRange,outdof,varargin)
% SSM_BP2EP This function performs continuation of equilibrium points of
% slow dynamics. Each equilibrium point corresponds to a periodic orbit in
% the regular time dynamics. The continuation here starts from a saved
% solution, which is a BRANCH point. The continuation follows the secondary
% branch passing through this point
%
% FRC = ODE_BP2EP(OBJ, OID, RUN, LAB, PARNAME, PARRANGE, OUTDOF, VARARGIN)
%
% oid:      runid of current continuation
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the label of
%           a branch point
% parName:  amp/freq continuation parameter
% parRange: continuation domain of parameter, which should be near the
%           value of natural frequency with index 1 in the mFreq if continuation
%           parameter is freq
% outdof:   dofs for output in physical domain
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_ISOL2EP SSM_EP2EP SSM_cont_ep

FRC = obj.SSM_cont_ep('BP',oid,run,lab,parName,parRange,outdof,varargin{:});
varargout{1} = FRC;
end