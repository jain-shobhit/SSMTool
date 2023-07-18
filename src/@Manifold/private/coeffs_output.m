function [W_0_out,R_0_out] = coeffs_output(W_0,R_0,order)
%%
% Storing the coefficients in suitable format for function output.
% They are stored in lexicographic ordering of the multi-indices
%%
W_0_out = repmat(struct('coeffs',[],'ind',[]),1,order);
R_0_out = repmat(struct('coeffs',[],'ind',[]),1,order);

l = size(R_0(1).coeffs,1);
for i = 1:order
    %create multi-indices
    if l >1
    K = sparse(sortrows(nsumk(l,i,'nonnegative')));
    else
        K=i;
    end
    % SSM coefficients with multi-indices
    idx_W_0     = all(W_0(i).coeffs==0);
    W_0_out(i).coeffs = W_0(i).coeffs(:,~idx_W_0);
    W_0_out(i).ind    =  K(~idx_W_0,:);
    
    % Reduced dynamics coefficients with multi-indices
    [~,idx_R_0] = find(R_0(i).coeffs);
    idx_R_0     = unique(idx_R_0);
    R_0_out(i).coeffs = R_0(i).coeffs(:,idx_R_0);
    R_0_out(i).ind   = K(idx_R_0,:);    

end
end