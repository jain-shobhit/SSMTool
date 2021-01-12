function opts = coco_set_xperm(opts, xidx, perm) % no longer in use?
%COCO_SET_XINFO   Set permutation vectors.
%
%   OPTS = COCO_SET_XINFO([OPTS], XIDX, PIDX) allows toolbox developers to
%   define permutation vectors for variables and parameters. This is useful
%   if one wants to define an algorithm as a map F:(x,p)->y, where x and y
%   have the same dimension, x plays the role of variables and p are free
%   parameters. The usual ordering of x and p in the combined vector u is
%   u=[x;p]. Here, XIDX = 1:n, PIDX = n+(1:m), for x n-by-1 and p m-by-1.
%   However, for many algorithms it is advantageous to reorder variables
%   and parameters. XIDX and PIDX are defined as follows:
%
%   x = u(XIDX),
%   p = u(PIDX),
%
%   that is, XIDX and PIDX contain the indices of the positions of the x-
%   and p- variables in the combined vector u.

%% affected fields in opts
%
%    opts.efunc.xidx    - order of variables

%% set permutation info fields of efunc

if max(abs(sort(xidx)-sort(perm)))>0
  error('%s: arguments do not form a valid permutation', mfilename);
end

opts.efunc.xidx(xidx) = perm;
