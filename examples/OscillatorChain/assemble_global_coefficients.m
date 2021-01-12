function [K,C,F2,F3] = assemble_global_coefficients(k,kappa2,kappa3,c,n)
% This function assembles sparse nonlinear coefficients
%% This function provides the nonlinear (cubic) spring and damping force
% coefficients at an element level in a oscillator chain

% Linear level
subs = [1 1;
        1 2;
        2 1;
        2 2];
vals = [ 1;
        -1;
        -1;
         1];


%% nonlinear coefficients in tensor format for the element
% %quadratic (x_2 - x_1)^2
% using sum over repeated to ensure uniquness of tensor representation
subs_2 = [1, 1, 1;
    1, 1, 2;
    1, 2, 1
    1, 2, 2;
    2, 1, 1;
    2, 1, 2;
    2, 2, 1;
    2, 2, 2];
vals_2 = [ -1; 1; 1;-1; 1; -1; -1; 1];

%cubic (x_2 - x_1)^3
subs_3 = [ones(8,1), subs_2;
        2*ones(8,1), subs_2];
vals_3 = [- vals_2;
            vals_2];

% scale values
vals_2 = kappa2 * vals_2;
vals_3 = kappa3 * vals_3;



%% Perform sparse assembly
n_el = n+1;
N = n+2; % total number of DOFs including Dirichlet DOFs
SUBS = zeros(n_el*4, 2);
VALS = zeros(n_el*4,1);
SUBS2 = zeros(n_el*8, 3);
VALS2 = zeros(n_el*8, 1);
SUBS3 = zeros(n_el*16, 4);
VALS3 = zeros(n_el*16, 1);

for j = 1:n_el
    % second order matrices
    index = j-1; % zero index for the j-th element
    SUBS(4*(j-1)+1:4*j, :) = repmat(index,[4,2]) + subs;
    VALS(4*(j-1)+1:4*j, :) = vals;
    % second-order nonlinearity     
    SUBS2(8*(j-1)+1:8*j, :) = repmat(index,[8,3]) + subs_2;
    VALS2(8*(j-1)+1:8*j, :) = vals_2;
    
    SUBS3(16*(j-1)+1:16*j, :) = repmat(index,[16,4]) + subs_3;
    VALS3(16*(j-1)+1:16*j, :) = vals_3;
end

F2 = sptensor(SUBS2,VALS2,[N,N,N]);
F3 = sptensor(SUBS3,VALS3,[N,N,N,N]);
K = k*sparse(SUBS(:,1),SUBS(:,2),VALS,N,N);
C = K*c/k;

