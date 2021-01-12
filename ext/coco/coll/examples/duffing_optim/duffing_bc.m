function y = duffing_bc(~,T,x0,x1,p)  %#ok<INUSD,INUSL>

y = [x1(1)-x0(1); x1(2)-x0(2); x1(3)-x0(3)-2*pi; x0(2)];
  
end
