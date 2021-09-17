function varargout = SSM_ep2HB(obj,oid,run,lab,parRange,outdof,varargin)
% SSM_EP2HB This function performs continuation of hopf bifurcation (HB)
% equilibirium points of slow dynamics. HB bifurcation is of codimension
% one and hence two parameters are free to vary to yield an one-dimensional
% manifold of HB points. Each HB point corresponds to a TR bifurcation
% periodic orbit in the regular time dynamics. The continuation here starts
% from a saved solution, which is a HB point.
%
% FRCIRS = SSM_EP2SN(OBJ,OID,RUN,LAB,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of current continuation
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the
%           label of a saddle-node point
% parRange: continuation domain of parameters. It is of the form
%           {[om1,om2],[f1,f2]}, where [om1,om2] and [f1,f2] specify the
%           continuation domain of excitation frequency and amplitude
%           respectively. You can give empty array and then no domain is
%           specified, e.g., {[],[f1,f2]} only presents the domain of
%           forcing amplitude
% outdof:   output for dof in physical domain
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_ISOL2EP, SSM_EP2EP, SSM_BP2EP, SSM_EP2SN

FRC = obj.SSM_cont_ep('HB',oid,run,lab,[],parRange,outdof,varargin{:});
varargout{1} = FRC;
end