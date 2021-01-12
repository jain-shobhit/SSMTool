function [x,y,n] = get_contour_xy(C)
%% 
% This function extracts and returns all the x and the y values in the form 
% of row vectors from C.
% 
% The input matrix C is in the following form (see contour documentation)
% 
% % 
% 
j = 1;
x = [];
y = [];
while j<size(C,2)
    j_new = C(2,j) + j;
    x = [x, nan, C(1,j+1:j_new)];
    y = [y, nan, C(2,j+1:j_new)];
    j = j_new + 1;
end

% number of connected components
n = sum(isnan(x));
end