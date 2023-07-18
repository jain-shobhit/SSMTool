function [RES,I_sum1,I_sum2] = multi_addition(SUMMAND1,SUMMAND2)
%% 
% The output is a matrix, containing all columns of $\texttt{SUMMAND1}$
% added to all columns of $\texttt{SUMMAND2}$ that do not lead to columns with negative entries. 
% The matrix $\texttt{I\_g}$ contains in position $j$ the column index of $\texttt{SUMMAND2}$ 
% corresponding to the column $j$ of $\texttt{RES}$. The same holds for $\texttt{I_sum1}$.
% 
% Example:
% 
% $\texttt{I\_sum2(3) = 1}$, $\texttt{I\_sum1(3) = 5}$ implies that  $\texttt{RES(:,3)}  
% = \texttt{SUMMAND1(:,5)} + \texttt{SUMMAND2(:,1)}$.

%%

sz_sum1 = size(SUMMAND1,2);    
sz_sum2 = size(SUMMAND2,2);

% Index of the first summand
I_sum1 = repmat(1:sz_sum1,1,sz_sum2);
% Index of the second summand
I_sum2 = reshape(repmat(1:sz_sum2,sz_sum1,1),1,[]);

RES = SUMMAND1(:,I_sum1)+SUMMAND2(:,I_sum2);

end
