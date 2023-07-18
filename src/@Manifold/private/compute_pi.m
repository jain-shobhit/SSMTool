function [pis] = compute_pi(n_mat,k_mat, field)
%%
% Computes pi_n_k for every combination of multi-indices in the two input
% matrices
% Multi index arrays n_mat and k_mat need to contain multi indices of same
% order |n| and |k| resepectively
% Pi are the forcing contribution coefficients as given in $\textit{Computing_spectral_submanifolds_using_multi_index_notation}$
%%
n_pos_idx = any(n_mat);
k_pos_idx = any(k_mat);

% No calculations for zero multi-indices needed
n_zero_idx = find(~n_pos_idx);
k_zero_idx = find(~k_pos_idx);

switch field.ordering
    case 'revlex'

        pis_pos = pink(n_mat(:,n_pos_idx),k_mat(:,k_pos_idx), field);
       
    case 'conjugate'
        %Needs all multi-indices in n_mat and k_mat to have same order
        pis_pos = pink_conj(n_mat(:,n_pos_idx),k_mat(:,k_pos_idx),field);
end
%% Read results into output array
%insert zero columns
if ~isempty(k_zero_idx)
    pis_cols = zeros(size(pis_pos,1),size(k_mat,2));
    pis_cols(:,k_pos_idx) = pis_pos;
else 
    pis_cols = pis_pos;
end

%insert zero rows
if ~isempty(n_zero_idx)
    pis = zeros(size(n_mat,2),size(k_mat,2));
    pis(n_pos_idx,:) = pis_cols;
else
    pis = pis_cols;
end
%% Special case
% fill pis with ones where k and n_mat both have a zero multiindex
for n_i = n_zero_idx
    for k_i = k_zero_idx
        pis(n_i,k_i) = 1;
    end
end
end

function [pis] = pink(n_mat,K,field)
assert(nnz(n_mat<0)==0,'Negative multiindices')

% Conjugate center index at all orders, for
% imaginary system this is empty since then no symmetry is
% inherent and we only need the number of multi-indices that
% exist at order k, given by nchoosek(k+l-1,l-1)
l           = size(K,1);        % SSM dimension
z_k         = size(K,2);
% Phase space dimension
N           = size(n_mat,1);
% nonzero entries in the multi-index array
F_multi_pos   = (n_mat ~= 0);
sz          = size(n_mat,2);
col_pos     = 1:sz;
%%
% Initialise space for the $\pi_{\mathbf{m},\mathbf{k}_f}$.
pis = zeros(size(col_pos,2),z_k);

%%
% The array $\texttt{sz\_pos}$ is a row vector containing in each row the number
% of nonzero entries of the multi-index the row corresponds to. This is useful,
% since $H_{i,0,\mathbf{h}} = 0$ for any nonzero $\mathbf{h}$. Therefore for each
% $\mathbf{m}$ in $\texttt{F\_multi}$ we only have to find the composition coefficients
% corresponding to the nonzero entries of $\mathbf{m}$. Furthermore then condition
% 1 looks like $\sum_{v,m_v \neq 0} \mathbf{h}_v = \mathbf{k}_f$. For each $f$
% and $\mathbf{m}$ all such combinations have to be found.
sz_pos   = sum(F_multi_pos,1);

%%
% This function computes all those combinations for all order $m$ multi-indices
% $\mathbf{m}$ and all $f$. It does so by considering unique numbers in $\texttt{sz\_pos}$
% and using them as the number of multi-indices $\mathbf{h}_i$ that have to sum
% up to $\mathbf{k}_f$ for all $f$.
%
% $\texttt{h}$ is a cell array, that contains in its first dimension the combinations
% for all $f$ and its second dimension corresponds to unique numbers in $\texttt{sz\_pos}$.
% They are stored in $\texttt{sz\_un}$. For each entry in $\texttt{sz\_pos}$ the
% corresponding entry in $\texttt{sz\_un\_ic}$ gives the column that it corresponds
% to in $\texttt{h}$.
%
% $\texttt{h\{f,i\}}$ contains a 3 dimensional array. The first dimension has
% size $l$, the second hase size $\texttt{sz\_un(i)}$. Its third dimension has
% size $\texttt{combnos(i)}$.
%
% This array contains in every matrix a set $\{\mathbf{h}_{v}\}_{v}$ that fulfills
% condition 1 for $\mathbf{k}_f$ that is stored in $\texttt{K2(:,f)}$. Each set
% consists of $\texttt{sz\_un(i)}$ elements.
%
% $\texttt{combnos}$ contains the amount of combinations that exist for every
% tuple $(\texttt{f,i})$.
[h, combnos,sz_un,sz_un_ic] = multi_nsumk(sz_pos,K);
%stores the order of all the multi-indices in h
h_abs = cellfun(@(x) reshape(sum(x,1),size(x,2),[]), h, 'UniformOutput',false);

%loop over all the order k multi-indices with index below center
%index in conjugate ordering
for f = 1:z_k
    
    %loop over unique amounts of nonzero entries in the
    %multi-indices m
    for i = 1:size(sz_un,2)
        
        %position of multi-indices in F_multi that have sz_un(i) nonzero entries
        idx    = sz_un_ic == i;
        
        %number of multi-indices in F_multi that have i nonzero entries
        sumidx = sum(idx);        
        
        if sz_un(i) == 0
            if nnz(K(:,f)) == 0
               pis(idx,f) = ones(sumidx,1);
            end
            continue
        end
        
        %columns of Fm that contain i nonzero elements
        n_idx         = n_mat(:,idx.');
        n_idx_find    = find(n_idx);
        [n_idx_row,~] = ind2sub(size(n_idx),n_idx_find);
        %%
        % $\texttt{Fm\_idx\_pos}$ contains the entries of $\texttt{F\_multi}$ that correspond
        % to the multi-indices that have $i$ nonzero elements. It is therefore a  $\texttt{sumidx
        % *sz\_un(i)}$ by 1 array.
        n_idx_pos     = n_idx(n_idx_find);
        %%
        % This step checks where condition 2 is fulfilled. $\texttt{h\_abs\{f,i\}}$
        % contains a $\texttt{sz\_un(i)}$ by $\texttt{combnos(i)}$ dimensionaly array.
        % In every column it contains the absolute values of all multi_indices of one
        % combination $\{\mathbf{h}_v\}_v$ that fulfills condition 1 for $\texttt{f}$
        % and contains $\texttt{sz\_un(i)}$ elements. Condition 2 is checked simultaneously
        % for all $\mathbf{m}$ of order $\texttt{order}$ that contain $\texttt{sz\_un(i)}$
        % nonzero entries.
        %
        % $$$\texttt{h\_abs\{f,i\}} = \pmatrix{ \mathbf{h_1^1} & \cdots & \mathbf{h_1^\texttt{combnos(i)}}\cr
        % \mathbf{h_2^1} & \cdots & \mathbf{h_2^\texttt{combnos(i)}}\cr \vdots & \cdots
        % & \vdots\cr \mathbf{h_{\texttt{sz\_un(i)}}^1} & \cdots & \mathbf{h_{\texttt{sz\_un(i)}}^\texttt{combnos(i)}}}
        % $.$$
        %
        % And therefore if the size $2n$ multi indices that have $\texttt{sz\_un(i)}$
        % nonzero elements are $\mathbf{m}_1, ... , \mathbf{m}_\texttt{sumidx}$, and we
        % define $\mathbf{m}_^1^+$ as the size $\texttt{sz\_un(i)}$ multi-index containing
        % the nonzero elements of $\mathbf{m}_1$,then condition 2 is checked as follows:
        %
        % $$\texttt{Cond\_dum} = \pmatrix{\texttt{h\_abs\{f,i\}} ==  \mathbf{m}_1^+\cr\texttt{h\_abs\{f,i\}}
        % ==  \mathbf{m}_2^+ \cr \vdots \cr\texttt{h\_abs\{f,i\}}==  \mathbf{m}_\texttt{sumidx}^+}
        % \in \mathbb{N} \text{mod 2}^{\texttt{sumidx}\ \times \ \texttt{combnos(i)}$$
        
        Cond_dum = repmat(h_abs{f,i},sumidx,1) >= full(n_idx_pos);
        
        %%
        % Now we would like to check which combinations fulfill condition 2 for all
        % of their elements. This corresponds to columns of $\texttt{h\_abs\{f,i\}} ==
        % \mathbf{m}_j^+$ where all elements are one. So summing over the rows only the
        % columns that have $\texttt{sz\_un(i)}$ as their entry correspond to valid combinations
        % that fulfill condition 2to for $\mathbf{m}_j^+$, and this has to be done for
        % $j = 1,...,\texttt{sumidx}$.
        %
        % To do this we transpose $\texttt{Cond\_dum}$ and permute the first two dimensions.
        % Furthermore, now a third dimension is added which each corresponds to one $\mathbf{m}_i^+$.
        %reshape, such that every slice of 3rd dim corresponds to
        %conditions for one m, first dimension
        Cond_dum = reshape(Cond_dum.',[],sz_un(i),sumidx);
        
        %permute, every column now contains conds for one combo of
        %h_is, sum all those cond values, if they are all one in a
        %column, they contribute.
        Cond_dum =  sum(permute( Cond_dum,[2,1,3]),1);
        %%
        % Now a slice of the third dimension looks like $\texttt{Cond\_dum(:,:,j)} =
        % \pmatrix{\texttt{h\_abs\{f,i\}} ==  \mathbf{m}_j^+} $. The sum over the rows
        % reveals which combos contribute. The resulting vector is recast, such that every
        % column corresponds to one $\mathbf{m}_j^+$ and the rows correspond to the combinations
        % that fulfill condition 1. The ones that also fulfill condition 2 are logical
        % ones, the others are zero.
        %check if every vector of a combo fulfills cond, if not
        %discard it
        Cond     = sparse(reshape(Cond_dum == full(sz_un(i)),combnos(f,i),[] ));

        if nnz(Cond)==0
            continue
        end
        %%
        % The array indices reads out the index of the valid combinations. Specifically
        % if the different combinations are numbered with numbers 1 to $\texttt{combnos(f,i)}$,
        % then it assigns every valid combination the number it corresponds to.
        %contains indices to read out the ones fulfilling cond2
        %from h
        indices  = repmat([1:combnos(f,i)].',1,sumidx);
        Cond_ind = indices(Cond);
        %%
        % Consequently all valid combinations of multi-indices are explicitly read out
        % into a array that contains all of the multi-indices in its second dimension.
        % All valid combinations for $\mathbf{m}_1^+$ are in it as a sequence, followed
        % by all valid combinations for $\mathbf{m}_2^+$ and so forth.
        %All combos in first dim. each combo has nns(ind)
        %contributions - to get back, reshape to l, nns(ind), []
        h_i      = reshape( h{f,i}(:,:,reshape(Cond_ind,[],1)) ,l,[]);
        %%
        % In order to read out the composition coefficients corresponding to all those
        % multi-indices we have to know their position in a set of reverse lexicographical
        % ordering of all multi-indices of their respective order and size. In order to
        % do this efficient, unique multi-indices are extracted from $\texttt{h\_i}$,
        % and their position is calculated, in either conjugate ordering (real system)
        % or in reverse lex. ordering (non real system). Then the subindices of all multi-indices
        % are stored in the array $\texttt{h\_i\_multi}$.
        [h_i_un,~,h_i_ic] = unique(h_i.','rows');
        h_i_idx_temp    = multi_index_2_ordering(h_i_un.','revlex',[]);
        h_i_idx         = h_i_idx_temp(h_i_ic.');
        %%
        % To read out the composition coefficients at the right order the absolute values
        % of the multi-indices also has to be known. They are stored in $\texttt{h\_i\_abs}$.
        if l > 1
            h_i_abs    = sum(h_i_un.');
        else
            h_i_abs = h_i_un.';
        end
        [h_i_abs_un,~,h_i_abs_ic] = unique(h_i_abs);
        %%
        % $\texttt{list\_1\_tmp}$ contains the positions of the nonzero elements of
        % all the multi-indices $\mathbf{m}_j^+$ in column $j$. Then each column is replicated
        % by the amount of combinations that fulfill both condition 1 and condition 2,
        % which is given by $\texttt{sum(Cond,1)}$, since $\texttt{Cond}$ contains in
        % column $j$ and row $r$ logical indices that indicate wheter combination $r$
        % fulfills both conditions for $\mathbf{m}_j^+$. That is then reshaped, such that
        % the entry $u$ in $\texttt{list\_1}$ contains the index $q$ of the composition
        % coefficient $H_{q,.,.}$ that the multi-index in position $u$ in $\texttt{h\_i}$
        % corresponds to.
        %contains the row index for all of the valid combos for
        %all ms in the same ordering as in HV
        list_1_tmp = reshape(n_idx_row,sz_un(i),[]);
        list_1     = reshape(repelem(list_1_tmp,1,sum(Cond,1)),1,[] );
        %%
        % The same thing is done for the entries (not the positions) of all the multi-indices
        % $\mathbf{m}_j^+$. The entry $u$ in $\texttt{list\_2}$ contains the value $m_q$
        % of the composition coefficient $H_{.,m_q,.}$ that the multi-index in position
        % $u$ in $\texttt{h\_i}$ corresponds to.
        %contains the entries of m for the valid combos for
        %all ms in the same ordering as in HV
        list_2_tmp = reshape(n_idx_pos,sz_un(i),[]);
        list_2     = reshape(repelem(list_2_tmp,1,sum(Cond,1)),1,[]);
        %%
        % What is left now is to multiply and add the indices accordingly to get the
        % right force contribution. Firstly all the composition coefficients are ready
        % out into $\texttt{pi\_dum}$.
        pi_dum  = ones(size(h_i,2),1);
        run_ord = 1;

        for ord = h_i_abs_un
            %%
            % The logical index $\texttt{idx\_dum}$ contains the positions of all multi-indices
            % in $\texttt{h\_i}$ that have order $\texttt{ord}$. The explicit positions of
            % them are stored in $\texttt{I\_ord\_pos}$.
            idx_dum = 1:size(h_i_abs_ic.',2);
            idx_dum = idx_dum(h_i_abs_ic.' == run_ord);
            idx_dum = (h_i_ic == idx_dum);
            idx_dum = sum(idx_dum,2);
            I_ord =  logical(idx_dum);

            %%
            % The values corresponding to those multi_indices of order $\texttt{ord}$ are
            % read out of the arrays that specify the position of the composition coefficients
            % for each of the multi-indices, being the phase space direction they correspond
            % to, the multi-index subindex in the reverse lexicographically ordered set they
            % have and the multi-index entry of the multi-index $\mathbf{m}_j^+$ they correspond
            % to.
            I_1   = list_1(I_ord);
            I_2   = h_i_idx(I_ord);
            I_3   = list_2(I_ord);
            
            %%
            % In order to make use of the symmetry of the SSM-coefficients, all the multi-indices
            % that in conjugate ordering have subindex bigger than the conjugate center index
            % are changed to their conjugate counterpart, and for them then conjugate composition
            % coefficients corresponding to this conjugate multi-index are read out.
            z_ord = nchoosek(ord+l-1,l-1);
            
            %split in the two index parts about conjugate
            %center index at order ord

            lin_idx  = sub2ind([N,z_ord,ord],I_1,I_2,I_3);
            pi_dum(I_ord) =  field.H{ord}(lin_idx);

            
            run_ord = run_ord +1;
        end
        
        %%
        % The following step multiplies all composition coefficients that correspond
        % to one combination of multi-indices fulfilling condition 1 and 2.
        %multiplication
        pi_dum = prod(reshape(pi_dum, sz_un(i),[]),1);
        %%
        % Next all the values of $\texttt{pi\_dum}$ corresponding to the same multi-index
        % $\mathbf{m}_j^+$ have to be added. This is done by creating an array $\texttt{list\_3}$
        % that contains for each entry in $\texttt{pi\_dum}$ the index $j$ of the size
        % $2n$ multi-index it corresponds to. Using $\texttt{accumarray}$ the values that
        % correspond to the same index $j$ are added up and then stored in the array $\texttt{pi}$
        % that contains all the composition terms.
        %sum over combos and put into position
        
        list_3 = repelem([1:sumidx].',sum(Cond,1));

        % Indices to read out the contributions into correct position
        n_contrib = sum(Cond,1)>0;
        idx_pos   = find(idx); 
        
        accum_pi  = accumarray(reshape(list_3,[],1),pi_dum);        
        pis(idx_pos(n_contrib),f) = accum_pi(find(n_contrib));
    end
end
end

function [pis] = pink_conj(F_multi,K,field)

% Conjugate center index at all orders, for
% imaginary system this is empty since then no symmetry is
% inherent and we only need the number of multi-indices that
% exist at order k, given by nchoosek(k+l-1,l-1)
k           = sum(K(:,1),1);    % Order of computation
Z_cci       = field.Z_cci;
z_k         = Z_cci(k);
l           = size(K,1);        % SSM dimension
% Phase space dimension
N           = field.N;

% if system is not real, then full set of order k multi-indices in reverse lex. ordering
revlex2conj = field.revlex2conj;   % Set to convert indices from rev. lex. to conjugate ordering

% nonzero entries in the multi-index array
F_multi_pos   = (F_multi ~= 0);
sz          = size(F_multi,2);
col_pos     = 1:sz;
%%
% Initialise space for the $\pi_{\mathbf{m},\mathbf{k}_f}$.
pis = zeros(size(col_pos,2),z_k);

%%
% The array $\texttt{sz\_pos}$ is a row vector containing in each row the number
% of nonzero entries of the multi-index the row corresponds to. This is useful,
% since $H_{i,0,\mathbf{h}} = 0$ for any nonzero $\mathbf{h}$. Therefore for each
% $\mathbf{m}$ in $\texttt{F\_multi}$ we only have to find the composition coefficients
% corresponding to the nonzero entries of $\mathbf{m}$. Furthermore then condition
% 1 looks like $\sum_{v,m_v \neq 0} \mathbf{h}_v = \mathbf{k}_f$. For each $f$
% and $\mathbf{m}$ all such combinations have to be found.
sz_pos   = sum(F_multi_pos,1);

%%
% This function computes all those combinations for all order $m$ multi-indices
% $\mathbf{m}$ and all $f$. It does so by considering unique numbers in $\texttt{sz\_pos}$
% and using them as the number of multi-indices $\mathbf{h}_i$ that have to sum
% up to $\mathbf{k}_f$ for all $f$.
%
% $\texttt{h}$ is a cell array, that contains in its first dimension the combinations
% for all $f$ and its second dimension corresponds to unique numbers in $\texttt{sz\_pos}$.
% They are stored in $\texttt{sz\_un}$. For each entry in $\texttt{sz\_pos}$ the
% corresponding entry in $\texttt{sz\_un\_ic}$ gives the column that it corresponds
% to in $\texttt{h}$.
%
% $\texttt{h\{f,i\}}$ contains a 3 dimensional array. The first dimension has
% size $l$, the second hase size $\texttt{sz\_un(i)}$. Its third dimension has
% size $\texttt{combnos(i)}$.
%
% This array contains in every matrix a set $\{\mathbf{h}_{v}\}_{v}$ that fulfills
% condition 1 for $\mathbf{k}_f$ that is stored in $\texttt{K2(:,f)}$. Each set
% consists of $\texttt{sz\_un(i)}$ elements.
%
% $\texttt{combnos}$ contains the amount of combinations that exist for every
% tuple $(\texttt{f,i})$.
[h, combnos,sz_un,sz_un_ic] = multi_nsumk(sz_pos,K);

%stores the order of all the multi-indices in h
h_abs = cellfun(@(x) reshape(sum(x,1),size(x,2),[]), h, 'UniformOutput',false);

%loop over all the order k multi-indices with index below center
%index in conjugate ordering
for f = 1:z_k
    
    %loop over unique amounts of nonzero entries in the
    %multi-indices m
    for i = 1:size(sz_un,2)
        
        %position of multi-indices in F_multi that have sz_un(i) nonzero entries
        idx    = sz_un_ic == i;
        
        %number of multi-indices in F_multi that have i nonzero entries
        sumidx = sum(idx);
        
        %columns of Fm that contain i nonzero elements
        Fm_idx         = F_multi(:,idx.');
        Fm_idx_find    = find(Fm_idx);
        [Fm_idx_row,~] = ind2sub(size(Fm_idx),Fm_idx_find);
        %%
        % $\texttt{Fm\_idx\_pos}$ contains the entries of $\texttt{F\_multi}$ that correspond
        % to the multi-indices that have $i$ nonzero elements. It is therefore a  $\texttt{sumidx
        % *sz\_un(i)}$ by 1 array.
        Fm_idx_pos     = Fm_idx(Fm_idx_find);
        %%
        % This step checks where condition 2 is fulfilled. $\texttt{h\_abs\{f,i\}}$
        % contains a $\texttt{sz\_un(i)}$ by $\texttt{combnos(i)}$ dimensionaly array.
        % In every column it contains the absolute values of all multi_indices of one
        % combination $\{\mathbf{h}_v\}_v$ that fulfills condition 1 for $\texttt{f}$
        % and contains $\texttt{sz\_un(i)}$ elements. Condition 2 is checked simultaneously
        % for all $\mathbf{m}$ of order $\texttt{order}$ that contain $\texttt{sz\_un(i)}$
        % nonzero entries.
        %
        % $$$\texttt{h\_abs\{f,i\}} = \pmatrix{ \mathbf{h_1^1} & \cdots & \mathbf{h_1^\texttt{combnos(i)}}\cr
        % \mathbf{h_2^1} & \cdots & \mathbf{h_2^\texttt{combnos(i)}}\cr \vdots & \cdots
        % & \vdots\cr \mathbf{h_{\texttt{sz\_un(i)}}^1} & \cdots & \mathbf{h_{\texttt{sz\_un(i)}}^\texttt{combnos(i)}}}
        % $.$$
        %
        % And therefore if the size $2n$ multi indices that have $\texttt{sz\_un(i)}$
        % nonzero elements are $\mathbf{m}_1, ... , \mathbf{m}_\texttt{sumidx}$, and we
        % define $\mathbf{m}_^1^+$ as the size $\texttt{sz\_un(i)}$ multi-index containing
        % the nonzero elements of $\mathbf{m}_1$,then condition 2 is checked as follows:
        %
        % $$\texttt{Cond\_dum} = \pmatrix{\texttt{h\_abs\{f,i\}} ==  \mathbf{m}_1^+\cr\texttt{h\_abs\{f,i\}}
        % ==  \mathbf{m}_2^+ \cr \vdots \cr\texttt{h\_abs\{f,i\}}==  \mathbf{m}_\texttt{sumidx}^+}
        % \in \mathbb{N} \text{mod 2}^{\texttt{sumidx}\ \times \ \texttt{combnos(i)}$$
        
        Cond_dum = repmat(h_abs{f,i},sumidx,1) >= full(Fm_idx_pos);
        %%
        % Now we would like to check which combinations fulfill condition 2 for all
        % of their elements. This corresponds to columns of $\texttt{h\_abs\{f,i\}} ==
        % \mathbf{m}_j^+$ where all elements are one. So summing over the rows only the
        % columns that have $\texttt{sz\_un(i)}$ as their entry correspond to valid combinations
        % that fulfill condition 2to for $\mathbf{m}_j^+$, and this has to be done for
        % $j = 1,...,\texttt{sumidx}$.
        %
        % To do this we transpose $\texttt{Cond\_dum}$ and permute the first two dimensions.
        % Furthermore, now a third dimension is added which each corresponds to one $\mathbf{m}_i^+$.
        %reshape, such that every slice of 3rd dim corresponds to
        %conditions for one m, first dimension
        Cond_dum = reshape(Cond_dum.',[],sz_un(i),sumidx);
        %permute, every column now contains conds for one combo of
        %h_is, sum all those cond values, if they are all one in a
        %column, they contribute.
        Cond_dum =  sum(permute( Cond_dum,[2,1,3]),1);
        %%
        % Now a slice of the third dimension looks like $\texttt{Cond\_dum(:,:,j)} =
        % \pmatrix{\texttt{h\_abs\{f,i\}} ==  \mathbf{m}_j^+} $. The sum over the rows
        % reveals which combos contribute. The resulting vector is recast, such that every
        % column corresponds to one $\mathbf{m}_j^+$ and the rows correspond to the combinations
        % that fulfill condition 1. The ones that also fulfill condition 2 are logical
        % ones, the others are zero.
        %check if every vector of a combo fulfills cond, if not
        %discard it
        Cond     = reshape(Cond_dum == full(sz_un(i)),combnos(f,i),[] );
        if nnz(Cond) == 0
            continue
        end
        %%
        % The array indices reads out the index of the valid combinations. Specifically
        % if the different combinations are numbered with numbers 1 to $\texttt{combnos(f,i)}$,
        % then it assigns every valid combination the number it corresponds to.
        %contains indices to read out the ones fulfilling cond2
        %from h
        indices  = repmat([1:combnos(f,i)].',1,sumidx);
        Cond_ind = indices(Cond);
        %%
        % Consequently all valid combinations of multi-indices are explicitly read out
        % into a array that contains all of the multi-indices in its second dimension.
        % All valid combinations for $\mathbf{m}_1^+$ are in it as a sequence, followed
        % by all valid combinations for $\mathbf{m}_2^+$ and so forth.
        %All combos in first dim. each combo has nns(ind)
        %contributions - to get back, reshape to l, nns(ind), []
        h_i      = reshape( h{f,i}(:,:,reshape(Cond_ind,[],1)) ,l,[]);
        %%
        % In order to read out the composition coefficients corresponding to all those
        % multi-indices we have to know their position in a set of reverse lexicographical
        % ordering of all multi-indices of their respective order and size. In order to
        % do this efficient, unique multi-indices are extracted from $\texttt{h\_i}$,
        % and their position is calculated, in either conjugate ordering (real system)
        % or in reverse lex. ordering (non real system). Then the subindices of all multi-indices
        % are stored in the array $\texttt{h\_i\_multi}$.
        [h_i_un,~,h_i_ic] = unique(h_i.','rows');
        h_i_idx_temp    = multi_index_2_ordering(h_i_un.',field.ordering,revlex2conj);
        h_i_idx         = h_i_idx_temp(h_i_ic.');
        %%
        % To read out the composition coefficients at the right order the absolute values
        % of the multi-indices also has to be known. They are stored in $\texttt{h\_i\_abs}$.
        h_i_abs    = sum(h_i_un.');
        [h_i_abs_un,~,h_i_abs_ic] = unique(h_i_abs);
        %%
        % $\texttt{list\_1\_tmp}$ contains the positions of the nonzero elements of
        % all the multi-indices $\mathbf{m}_j^+$ in column $j$. Then each column is replicated
        % by the amount of combinations that fulfill both condition 1 and condition 2,
        % which is given by $\texttt{sum(Cond,1)}$, since $\texttt{Cond}$ contains in
        % column $j$ and row $r$ logical indices that indicate wheter combination $r$
        % fulfills both conditions for $\mathbf{m}_j^+$. That is then reshaped, such that
        % the entry $u$ in $\texttt{list\_1}$ contains the index $q$ of the composition
        % coefficient $H_{q,.,.}$ that the multi-index in position $u$ in $\texttt{h\_i}$
        % corresponds to.
        %contains the row index for all of the valid combos for
        %all ms in the same ordering as in HV
        list_1_tmp = reshape(Fm_idx_row,sz_un(i),[]);
        list_1     = reshape(repelem(list_1_tmp,1,sum(Cond,1)),1,[] );
        %%
        % The same thing is done for the entries (not the positions) of all the multi-indices
        % $\mathbf{m}_j^+$. The entry $u$ in $\texttt{list\_2}$ contains the value $m_q$
        % of the composition coefficient $H_{.,m_q,.}$ that the multi-index in position
        % $u$ in $\texttt{h\_i}$ corresponds to.
        %contains the entries of m for the valid combos for
        %all ms in the same ordering as in HV
        list_2_tmp = reshape(Fm_idx_pos,sz_un(i),[]);
        list_2     = reshape(repelem(list_2_tmp,1,sum(Cond,1)),1,[]);
        %%
        % What is left now is to multiply and add the indices accordingly to get the
        % right force contribution. Firstly all the composition coefficients are ready
        % out into $\texttt{pi\_dum}$.
        pi_dum  = ones(size(h_i,2),1);
        run_ord = 1;
        for ord = h_i_abs_un
            %%
            % The logical index $\texttt{idx\_dum}$ contains the positions of all multi-indices
            % in $\texttt{h\_i}$ that have order $\texttt{ord}$. The explicit positions of
            % them are stored in $\texttt{I\_ord\_pos}$.
            idx_dum = 1:size(h_i_abs_ic.',2);
            idx_dum = idx_dum(h_i_abs_ic.' == run_ord);
            idx_dum = (h_i_ic == idx_dum);
            idx_dum = sum(idx_dum,2);
            I_ord =  logical(idx_dum);
            I_ord_pos = find(I_ord);
            %%
            % The values corresponding to those multi_indices of order $\texttt{ord}$ are
            % read out of the arrays that specify the position of the composition coefficients
            % for each of the multi-indices, being the phase space direction they correspond
            % to, the multi-index subindex in the reverse lexicographically ordered set they
            % have and the multi-index entry of the multi-index $\mathbf{m}_j^+$ they correspond
            % to.
            I_1   = list_1(I_ord);
            I_2   = h_i_idx(I_ord);
            I_3   = list_2(I_ord);
            
            %%
            % In order to make use of the symmetry of the SSM-coefficients, all the multi-indices
            % that in conjugate ordering have subindex bigger than the conjugate center index
            % are changed to their conjugate counterpart, and for them then conjugate composition
            % coefficients corresponding to this conjugate multi-index are read out.
            z_ord = nchoosek(ord+l-1,l-1);
            z_cci_ord = Z_cci(ord);
            
            %split in the two index parts about conjugate
            %center index at order ord
            idx_2   = I_2 <= z_cci_ord;
            
            lin_idx_a  = sub2ind([N,Z_cci(ord),ord],I_1(idx_2),I_2(idx_2),I_3(idx_2));
            %[N,Z_cci(ord),ord],I_1(idx_2),I_2(idx_2),I_3(idx_2)
            %field.H{ord},(lin_idx_a)
            
            pi_dum(I_ord_pos(idx_2)) =  field.H{ord}(lin_idx_a);
            
            lin_idx_b  = sub2ind([N,Z_cci(ord),ord],I_1(~idx_2),z_ord-I_2(~idx_2)+1,I_3(~idx_2));
            pi_dum(I_ord_pos(~idx_2)) = conj(field.H{ord}(lin_idx_b));
            
            
            run_ord = run_ord +1;
        end
        %%
        % The following step multiplies all composition coefficients that correspond
        % to one combination of multi-indices fulfilling condition 1 and 2.
        %multiplication
        pi_dum = prod(reshape(pi_dum, sz_un(i),[]),1);
        %%
        % Next all the values of $\texttt{pi\_dum}$ corresponding to the same multi-index
        % $\mathbf{m}_j^+$ have to be added. This is done by creating an array $\texttt{list\_3}$
        % that contains for each entry in $\texttt{pi\_dum}$ the index $j$ of the size
        % $2n$ multi-index it corresponds to. Using $\texttt{accumarray}$ the values that
        % correspond to the same index $j$ are added up and then stored in the array $\texttt{pi}$
        % that contains all the composition terms.
        %sum over combos and put into position
        
        list_3 = repelem([1:sumidx].',sum(Cond,1));

        % Indices to read out the contributions into correct position
        n_contrib = sum(Cond,1)>0;
        idx_pos   = find(idx); 
        
        accum_pi  = accumarray(reshape(list_3,[],1),pi_dum);        
        pis(idx_pos(n_contrib),f) = accum_pi(find(n_contrib));

    end
end
end
