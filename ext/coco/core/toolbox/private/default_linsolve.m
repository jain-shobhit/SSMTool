function [opts x] = default_linsolve(opts, A, b)
%DEFAULT_LINSOLVE  Default implementation for solving Ax=b.
%
%   [OPTS X] = COCO_DEFAULT_LINSOLVE(OPTS, A, B) solves the system of
%   linear equations A*X=B using Matlab's backslash operator.
%
if nargout<2
	error('%s: too few output arguments', mfilename);
end

[m n] = size(A);
if m~=n
  fprintf(2, '%s: warning, matrix A not square, size(A) = [%d %d]\n', ...
    mfilename, m, n);
end

% current_spparms = spparms();
% spparms('umfpack', 0);

% fprintf(2, '%s: condest(A) = %.4e\n', mfilename, condest(A));

% mldivide uses:
% thresh = [0.1, 0.001];
% [L,U,P,Q,R] = lu(A, thresh); % we have P*(R\A)*Q = L*U
%                              % with this we get X = Q*(U\L\(P*(R\B)))
%                              % or, use X = linsolve(A,B,opts)

x = A\b;

% spparms(current_spparms);
