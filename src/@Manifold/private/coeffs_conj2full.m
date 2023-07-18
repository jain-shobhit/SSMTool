function [Coeff_full] = coeffs_conj2full(Coeff_part,row_idx,col_idx, ordering, type)
%%
%This function computes the full coefficients with ordering specified
%in the input array with the same name, given the coefficients in conjugate
%ordering up to the conjugate center index. It is used in a cellfun of
%conj2lex as well for the autonomous computation.

%% Input
% Coeff_part - coefficients in conjugate ordering up to conjugate center
%              index
% x_idx      - option to change coordinate directions (for symmetry of
%              reduced dynamics coefficients
% y_idx      - conjugate center index at the order of the coefficients
% ordering   - array that converts conjugate to other ordering


%% Ouput
% Coeff_full - full coefficients in ordering specified by the array
% ordering

%%
switch type
    case 'TaylorCoeffs'
        if isempty(row_idx)
            Coeff_full.coeffs = [Coeff_part.coeffs, flip(conj(Coeff_part.coeffs(:,1:col_idx,:)),2)];
            Coeff_full.coeffs = Coeff_full.coeffs(:,ordering,:);
        else
            Coeff_full.coeffs = [Coeff_part.coeffs, flip(conj(Coeff_part.coeffs(row_idx,1:col_idx)),2)];
            Coeff_full.coeffs = Coeff_full.coeffs(:,ordering);
        end
    case 'CompCoeffs'   
            z_k = numel(ordering);
            Coeff_full = zeros(size(Coeff_part,1),z_k,size(Coeff_part,3));
            for i = 1:size(Coeff_part,3)
                Coeff_full(:,:,i) = [Coeff_part(:,:,i), flip(conj(Coeff_part(:,1:col_idx,i)),2)];
                Coeff_full(:,:,i) = Coeff_full(:,ordering,i);
            end
end
end
