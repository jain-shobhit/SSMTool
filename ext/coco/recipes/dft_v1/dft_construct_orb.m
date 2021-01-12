function prob = dft_construct_orb(prob, tbid, data, sol)
%DFT_CONSTRUCT_ORB   Append an instance of 'dft' to problem.
%
% Add collocation and continuity conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters.
%
% PROB = DFT_CONSTRUCT_ORB(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: dft_construct_orb.m 2839 2015-03-05 17:09:01Z fschild $

data.tbid = tbid;
data = coco_func_data(data); % Convert to func_data class for shared access
prob = coco_add_func(prob, tbid, @dft_F, @dft_DFDU, data, ...
  'zero', 'u0', sol.u0);
uidx = coco_get_func_data(prob, tbid, 'uidx');
fid  = coco_get_id(tbid, 'period');
prob = coco_add_pars(prob, fid, uidx(data.T_idx), fid, 'active'); % Monitor period
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  prob = coco_add_pars(prob, fid, uidx(data.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');
prob = coco_add_slot(prob, tbid, @dft_update, data, 'update'); % Update number of active modes

efid = coco_get_id(tbid, {'err' 'err_TF' 'NMOD'});
prob = coco_add_func(prob, efid{1}, @dft_error, data, ...
  'regular', efid, 'uidx', uidx); % Monitor discretization error estimate and number of active modes
prob = coco_add_event(prob, 'MXCL', 'MX', efid{2}, '>', 1); % MXCL - event type

end

function [data y] = dft_F(prob, data, u)
%DFT_F   Implement 'dft' toolbox zero problem.
%
% Encode the spectral zero problem and the correponding phase condition.

xf = u(data.xf_idx); % Extract active Fourier modes
xd = u(data.xd_idx); % Extract dummy Fourier modes
T  = u(data.T_idx);  % Extract period
p  = u(data.p_idx);  % Extract problem parameters

dim  = data.dim;

pp  = repmat(p, data.p_rep);
xp = reshape(data.FinvsW*xf, data.x_shp); % Values at mesh points
f  = data.fhan(xp, pp);
f  = T*data.Fs*f(:)-data.Wp*xf;
f  = [real(f(1:dim)); real(f(dim+1:end)); imag(f(dim+1:end))]; % Spectral zero problem
y  = [f; xd; data.phs*xf]; % data.phs*xf = 0 is phase condition

end

function [data J] = dft_DFDU(prob, data, u)
%DFT_DFDU   Linearization of 'dft' zero problem.
%
% Encode linearization of spectral zero problem as well as phase condition.

xf = u(data.xf_idx); % Extract active Fourier modes
T  = u(data.T_idx);  % Extract period
p  = u(data.p_idx);  % Extract period

dim  = data.dim;
pdim = data.pdim;
mdim = data.mdim;
ddim = data.ddim;

pp  = repmat(p, data.p_rep);
xp = reshape(data.FinvsW*xf, data.x_shp); % Values at mesh points
if isempty(data.dfdxhan)
  df = coco_ezDFDX('f(x,p)v', data.fhan, xp, pp);
else
  df = data.dfdxhan(xp, pp);
end
df = sparse(data.dxrows, data.dxcols, df(:));
df = T*data.Fs*df*data.FinvsW - data.Wp; % W.r.t. Fourier modes
df = [real(df(1:dim,:)); real(df(dim+1:end,:)); imag(df(dim+1:end,:))];

f = data.fhan(xp, pp);
f = data.Fs*f(:); % W.r.t. to period
f = [real(f(1:dim)); real(f(dim+1:end)); imag(f(dim+1:end))];

J1 = [df, sparse(mdim, ddim), f; data.jac];

if isempty(data.dfdphan)
  df = coco_ezDFDP('f(x,p)v', data.fhan, xp, pp);
else
  df = data.dfdphan(xp, pp);
end
df = sparse(data.dprows, data.dpcols, df(:));
df = T*data.Fs*df; % W.r.t. problem parameters
df = [real(df(1:dim,:)); real(df(dim+1:end,:)); imag(df(dim+1:end,:))];

if pdim>0
  dfcont = sparse(ddim+1, pdim);
else
  dfcont = [];
end

J2 = [df; dfcont];

J = sparse([J1 J2]);

end

function data = dft_update(prob, data, cseg, varargin)
%DFT_UPDATE   Update spectral discretization.
%
% Use information about current solution to update the number of active and
% dummy modes, provided that within given bounds.

uidx = coco_get_func_data(prob, data.tbid, 'uidx'); % 'dft' context-dependent index array
u    = cseg.src_chart.x(uidx); % Current chart
xf   = u(data.xf_idx); % Extract active Fourier modes
p    = u(data.p_idx);  % Extract problem parameters

pp   = repmat(p, data.p_rep);
xp   = reshape(data.FinvsW*xf, data.x_shp); % Values at mesh points
err  = data.fhan(xp, pp)*data.F; % Last columns provide neglected residual
err  = real(err.*conj(err));
err  = sqrt(sum(err(:))); % Discretization error estimate

dft  = data.dft;
NMOD = dft.NMOD;
NMAX = dft.NMAX;
NMIN = dft.NMIN;
dim  = data.dim;
sol.x0 = reshape(xf, [dim 2*NMOD+1]);
NMODi = min(ceil(NMOD*1.1+1),   NMAX);
NMODd = max(ceil(NMOD/1.025-1), NMIN);
if err>dft.TOLINC && NMOD~=NMODi % If above upper bound of adaptation window
  data.dft.NMOD = NMODi;
  sol.x0 = [sol.x0 zeros(dim, 2*(NMODi-NMOD))];
  data = dft_init_modes(data, sol); % Reconstruct data
elseif err<dft.TOLDEC && NMOD~=NMODd % If below lower bound of adaptation window
  data.dft.NMOD = NMODd;
  sol.x0 = sol.x0(:,1:2*NMODd+1);
  data = dft_init_modes(data, sol); % Reconstruct data
end

end

function [data y] = dft_error(prob, data, u)
%DFT_ERROR   Evaluate estimate of approximation error.
%
% Estimate error in discretization by examining the neglected residuals in
% the spectral formulation.

xf   = u(data.xf_idx); % Extract active Fourier modes
p    = u(data.p_idx);  % Extract problem parameters

pp   = repmat(p, data.p_rep);
xp   = reshape(data.FinvsW*xf, data.x_shp); % Values at mesh points
err  = data.fhan(xp, pp)*data.F; % Last columns provide neglected residual
err  = real(err.*conj(err));
err  = sqrt(sum(err(:))); % Discretization error estimate

dft = data.dft;
y   = [err; err/dft.TOL; dft.NMOD];

end
