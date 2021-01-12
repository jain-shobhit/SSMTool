function [W_0_out,R_0_out] = unified_output(W_0,R_0,order)

% Puts struct in each cell entry that contains the multi-indices the
% coefficients correspond to

W_0_out = cell(1,order);
R_0_out = cell(1,order);
l = size(R_0{1},1);
for i = 1:order
    %create multi-indices
    K = sparse(sortrows(nsumk(l,i,'nonnegative')));

    % SSM coefficients with multi-indices
    idx_W_0     = all(W_0{i}==0);
    W_0i.coeffs = W_0{i}(:,~idx_W_0);
    W_0i.ind    = K(~idx_W_0,:);
    W_0_out{i}  = W_0i;
    
    % Reduced dynamics coefficients with multi-indices
    [~,idx_R_0] = find(R_0{i});
    idx_R_0     = unique(idx_R_0);
    R_0i.coeffs = R_0{i}(:,idx_R_0);
    R_0i.ind    = K(idx_R_0,:);
    R_0_out{i}  = R_0i;    

end
end