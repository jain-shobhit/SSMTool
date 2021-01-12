function [idx] = conjugate_flip(l_i,l_r)

% This function computes and index array that flips the conjugate coordinate directions  

idx = reshape(1:2*l_i,2,[]);
idx = reshape(flip(idx),1,[]);
idx = [idx,2*l_i+1:(2*l_i+l_r)];
end