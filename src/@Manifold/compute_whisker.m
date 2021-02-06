function [W_0, R_0] = compute_whisker(obj, order)
% Invariant manifold in the autonomous system limit

% Initialize
W_0 = cell(1,order); R_0 = cell(1,order);

% Leading-order terms
Lambda_E = obj.E.spectrum;
W_01     = obj.E.basis;

switch obj.Options.notation
    
    case 'tensor'
        W_0{1} = sptensor(W_01);
        R_0{1} = sptensor(diag(Lambda_E));
        multi_input = [];
    case 'multiindex'
        % Set up system for Multi-index notation version
        [W_0{1},R_0{1},multi_input] = multi_index_setup(obj,order);
end

%% 
% *Outer resonance case* : issue warning
if obj.resonance.outer.occurs
    prompt = 'Due to (near) outer resonance, the exisitence of the manifold is questionable and the underlying computation may suffer.';
    disp(prompt)
    disp('Attempting manifold computation')
end


for j = 2:order
    %recursively approximate at j-th order
    startOrderj = tic;
    [W_0{j},R_0{j},multi_input] = cohomological_solution(obj,j,W_0,R_0,multi_input);
    obj.solInfo.timeEstimate(j) = toc(startOrderj);
    disp(['Manifold computation time at order ' num2str(j) ' = ' datestr(datenum(0,0,0,0,0,obj.solInfo.timeEstimate(j)),'HH:MM:SS')])
    fprintf('Estimated memory usage at order % 2i = %05.2E MB\n', j, obj.solInfo.memoryEstimate(j))
end

switch obj.Options.notation
    case 'tensor'
        for j = 1:order
            W_0{j} = tensor_to_multi_index(W_0{j});
            R_0{j} = tensor_to_multi_index(R_0{j});
        end
        
    case 'multiindex'
        [W_0,R_0] = conjugate_to_lexicographic_coefficients(multi_input,order,W_0,R_0);
        
        %%
        % Add multi-indices field to the coefficients field at every order
        [W_0,R_0] = unified_output(W_0,R_0,order);
        
end