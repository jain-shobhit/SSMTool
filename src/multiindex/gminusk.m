function [S,I_k,I_g] = gminusk(g,k)
%% GMINUSK gminusk
k_l = size(k,2);    
g_l = size(g,2);
S   = cell(1,g_l);
I_k = cell(1,g_l);
I_g = cell(1,g_l);
%% 
% The output is a matrix, containing all columns of $\texttt{k}$ subtracted 
% all columns of $\texttt{g}$ that do not lead to columns with negative entries. 
% The matrix $\texttt{I\_g}$ contains in position $j$ the column index of $\texttt{g}$ 
% corresponding to the column $j$ of $\texttt{S}$. The same holds for $\texttt{I\_k}$.
% 
% Example:
% 
% $\texttt{I\_g(3) = 1}$, $\texttt{I\_k(3) = 5}$ implies that  $\texttt{S(:,3)}  
% = \texttt{k(:,5)} - \texttt{g(:,1)}$.
for i = 1:g_l
    S{i}   = k - g(:,i);
    I_k{i} = 1:k_l;
    I_g{i} = ones(1,k_l)*i;
end
S          = cell2mat(S);
I_k        = cell2mat(I_k);
I_g        = cell2mat(I_g);
[~,neg]    = find(S<0);
S(:,neg)   = []; %delete columns with negative entries
I_g(:,neg) = [];
I_k(:,neg) = [];
end