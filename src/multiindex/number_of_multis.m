function [Z] = number_of_multis(l,max_order)

%Returns the number of size l multi-indices for all orders up to order max_order

Z = zeros(1,max_order);

for j = 1:max_order
    Z(j) = nchoosek(j+l-1,l-1);
end
end