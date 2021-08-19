function y = net(t, x, p)
%VDP   'coll'-compatible encoding of Van der pol vector field.

u  = x(1:4,:); % disp
v  = x(5:8,:); % vel
om = p(1,:);
a  = p(2,:);
u1 = u(1,:);
v1 = v(1,:);
ep = 0.05;

C = diag([-1; 0.8 ;0.8 ;0.8]);
K = [3 -1 0 -1;-1 4 -1 -1;0 -1 2 0;-1 -1 0 3];

y = zeros(8,numel(t));
y(1:4,:) = v;
y(5:8,:) = -ep*C*v-K*u;
y(5,:)   = y(5,:)+ep*a.*cos(om.*t)-ep*u1.^2.*v1;
% y(5:8,:) = -0.075*v-K*u;
% y(5,:)   = y(5,:)+ep*a.*cos(om.*t);

end
