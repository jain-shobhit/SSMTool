function [Coeff_full] = restore_full_coeff(Coeff_part,x_idx,y_idx, ordering)

%This function computes the full coefficients with ordering specified
%in the input array with the same name, given the coefficients in conjugate
%ordering up to the conjugate center index.

%% Input
% Coeff_part - coefficients in conjugate ordering up to conjugate center
%              index
% x_idx      - option to change coordinate directions (for symmetry of
%              reduced dynamics coefficients
% y_idx      - conjugate center index at the order of the coefficients
% ordering   - array that converts conjugate to other ordering


%% Ouput
% Coeff_full - full coefficients in ordering specified by the array
%               ordering

        if isempty(x_idx)
            Coeff_full = [Coeff_part, flip(conj(Coeff_part(:,1:y_idx,:)),2)];
            Coeff_full = Coeff_full(:,ordering,:);
        else
            Coeff_full = [Coeff_part, flip(conj(Coeff_part(x_idx,1:y_idx)),2)];
            Coeff_full = Coeff_full(:,ordering);
        end
        

end