function dfnl = compute_dfnldx(obj,x,xd)
% COMPUTE_DFNLDX This function computes the Jacobian of  nonlinear internal 
% force with respect to the displacements x in a second-order
% mechanical system. Currently, we do not treat velocity dependent
% nonlinearity.

assert(obj.order == 2, ' dfnldx can only be computed for second-order systems')

dfnl = sparse(obj.n,obj.n);
for j = 1:length(obj.fnl)
    dfnl = dfnl + expand_tensor_derivative(obj.fnl{j},x);
end