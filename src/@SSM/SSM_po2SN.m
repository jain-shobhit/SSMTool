function SSM_po2SN(obj,oid,run,lab,parRange,outdof,varargin)
% SSM_PO2SN This function performs continuation of saddle-node (SN)
% periodic orbits of slow dynamics. Each periodic orbit corresponds to a
% torus (quasi-periodic) response in regular time dynamics. The
% continuation here starts from a saved SN solution.
%
% FRCIRS = SSM_PO2SN(OBJ,OID,RUN,LAB,PARRANGE,OUTDOF,VARARGIN)
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
% See also: SSM_HB2PO, SSM_PO2PO, SSM_BP2PO

obj.SSM_cont_po('SN',oid,run,lab,[],[],parRange,outdof,varargin{:});
end
