function [nodes,elements,bnodes] = RectangularMesh(L,H,n,m,w)
%% Description
% L = length of plate
% H = height of plate
% n = number of divisions in L direction
% m = number of divisions in H direction
% w = height of curvature in z direction 

% the plate lies in the xy plane
% z direction is the out-of-plane direction

% output
% nodes : [xcoord ycoord zcoord] 
% elements : [ NodeID1 NodeID2 NodeID3 ] 

nodes = zeros((n+1)*(m+1),3); % Initialize

seeds=1:(n+1)*(m+1);
seeds=reshape(seeds,n+1,m+1)'; 
count=1;
for j=1:m+1
    for i=1:n+1
        nodes(count,:)=[((i-1)/n)*L H-(j-1)/m*H 0];
        count=count+1;
    end
end

if w~=0
    R = 1/2*(w^2+(H/2)^2)/w;
    theta0 = asin(H/2/R);    
    for i = 1:size(nodes,1)
        th = asin((nodes(i,2)-H/2)/R);
        nodes(i,3) = R*cos(th)-cos(theta0);
    end
end

elements = zeros(2*n*m,3); % Initialize

count=1;
for j=1:m
    for i=1:n
    
      elements(count:count+1,:)=[ seeds(j,i)   seeds(j+1,i+1)   seeds(j,i+1);
                                  seeds(j,i) seeds(j+1,i)     seeds(j+1,i+1)];
  
      count=count+2;
       
   end
end

bnodes = cell(1,4);
bnodes{1}=1:n+1:m*(n+1)+1;
bnodes{2}=n+1:n+1:m*(n+1)+n+1;
bnodes{3}=1:n+1;
bnodes{4}=m*(n+1)+1:m*(n+1)+n+1;


