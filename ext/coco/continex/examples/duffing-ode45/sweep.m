% p = [A om la al eps]
p = [1 2.5 0.2 1 1]';

x        = zeros(7,1);
[data c] = exp_F([], x, p);

clf
oms = linspace(0.1, 3, 50);
fws = [];

for om=oms
  p(2)     = om;
  [data c] = exp_F(data, x, p);
  fws      = [fws norm(c)]; %#ok<AGROW>
  
  hold on
  plot(om, norm(c), 'b*')
  drawnow
  hold off
end

bws = [];
for om=fliplr(oms)
  p(2)     = om;
  [data c] = exp_F(data, x, p);
  bws      = [norm(c) bws]; %#ok<AGROW>
  
  hold on
  plot(om, norm(c), 'ro')
  drawnow
  hold off
end

save('sweep', 'oms', 'fws', 'bws');
