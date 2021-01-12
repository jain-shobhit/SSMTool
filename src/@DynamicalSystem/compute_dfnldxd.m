function dfnl = compute_dfnldxd(obj,x,xd)
% COMPUTE_DFNLDXD This function computes the Jacobian of  nonlinear internal 
% force f with respect to the velocities xd in a second-order
% mechanical system. Currently, we do not treat velocity dependent
% nonlinearity.

assert(obj.order == 2, ' dfnldxd can only be computed for second-order systems')

dfnl = sparse(obj.n,obj.n);
% for j = 1:length(obj.fnl)
%     dfnl = dfnl + expand_tensor_derivative(obj.fnl{j},x);
% end