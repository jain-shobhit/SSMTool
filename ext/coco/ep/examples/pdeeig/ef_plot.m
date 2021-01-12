function thm = ef_plot(thm, ~, ~, bd, ~)

persistent X Y P Q N

if isempty(X) || ~(N==thm.N)
  N = thm.N;
  P = 3*N+1;
  Q = 2*N+1;
  [X, Y] = meshgrid(linspace(0,2,Q), linspace(0,3,P));
end
lab = coco_bd_labs(bd, 'UZ');
U   = coco_bd_val(bd, lab, 'x');
U   = reshape(U, P, Q);
mesh(X,Y,U)

thm.xlab = 'x'; thm.ylab = 'y'; thm.zlab = 'u';

end
