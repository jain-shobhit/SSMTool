function [g, g_comb,nv_un,nv_ic] = nsumkvector(nv,v)
%  NSUMKVECTOR This function finds all combinations of $nv$ positive integer vectors adding
% 
% up to another positive integer vector $\mathbf{v}$:
% 
% Let $\mathbf{v} \in \mathbb{N}^l$. It finds all combinations consisting of 
% $n$ positive, nonzero integer vectors such that they sum up to
% 
% $\mathbf{v}$: $\sum_{i=1}^n \mathbf{g}_i = \mathbf{v}$. All of those combinations 
% are stored as 2D arrays, where the individual columns are filled with thevectors, 
% in column $i$ the vector $\mathbf{g}_i$ is stored. The various combinations 
% are stacked upon each other in a 3D array, which contains a third dimension 
% with size equal to amount of combinations of vectors that sum up to $\mathbf{v}$.
% 
% The result is then stored in a cell array, that contains these combinations 
% for all tuples 
% 
% $i,j$ , corresponding to $\texttt{v(:,i), nv(j)}$.
%there can be multiple vectors in v
z_k = size(v,2);
l = size(v,1);
S = cell(1,l);
[nv_un,~,nv_ic] = unique(nv);
g = cell(z_k,size(nv_un,2));
g_comb = zeros(z_k,size(nv_un,2));
n_run = 1;
for n = nv_un
    for f = 1:z_k
        
        for i = 1:l
            %special case
            if v(i,f) == 1 && n == 1
                S{i} = 1;
            else
                S{i} = nsumk(n,v(i,f),'nonnegative').' ;
                %finds all combinations for each index in v
            end
        end
%% 
% $\texttt{combvec}$ finds all possible combinations of the vectors in $\texttt{S}$.
        g_fn = combvec(S{:});
        g_fn = reshape(g_fn,n,l,[]);
        g_fn = permute(g_fn,[2 1 3]);
        
        
        g{f,n_run}  = g_fn;
        g_comb(f,n_run) = size(g_fn,3);
    end
    
    n_run = n_run +1;
end
end