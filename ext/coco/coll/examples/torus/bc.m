function [data, y, J] = bc(prob, data, u)

T  = u(data.T_idx);
x0 = u(data.x0_idx);
x1 = u(data.x1_idx);
p  = u(data.p_idx);

y = [T-2*pi/p(1);data.F*x1-data.RF*x0; x0(2)];

if nargout>=3
  nt = numel(T);
  nx = numel(x0);
  np = numel(p);
  
  J1 = zeros(1,nt+2*nx+np);
  J1(1,nt+2) = 1;
  
  J = [eye(nt), zeros(nt,2*nx),  2*pi/p(1)^2*ones(nt,1), zeros(nt,1);
    zeros(nx,nt), -data.RF, data.F, zeros(nx,np);
    J1];
end

end
