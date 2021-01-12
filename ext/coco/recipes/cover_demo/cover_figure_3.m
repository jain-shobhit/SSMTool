function cover_figure_3

% Generate data
pd = cd('../atlas1d_v7_demo');
if ~coco_exist('run2', 'run')
  run demo_atlas1d_v7
end

% Extract data
bd  = coco_bd_read('run2');

% Plot data
figure(1)
clf
hold on

plot_bd(bd)

axis tight
axis off
view([5 45])

camproj perspective
lighting flat
set(gca, 'AmbientLightColor', [1 1 1])
material dull

light('Position', [ 2  2 4], 'Color', 0.15*[1 1 1])
light('Position', [-4 -1 1], 'Color', 0.5*[1 1 1])
light('Position', [-1 -1 1], 'Color', 0.25*[1 1 1])
light('Position', [ 0 -1 1], 'Color', 0.25*[1 1 1])
light('Position', [ 1 -1 1], 'Color', 0.25*[1 1 1])
light('Position', [ 4 -1 1], 'Color', 0.5*[1 1 1])

hold off

cd(pd);

end

function plot_bd(bd)
A  = coco_bd_col(bd, 'A');
x0 = coco_bd_col(bd, 'X0');
st = coco_bd_col(bd,'hspo.test.stab');

[X Y Z] = sphere(20);
R = 0.015;

x = x0(1,:);
y = x0(2,:);

XX    = zeros(0,3);
o     = [];
scale = [80.5 7.6 44.7];
RR    = [3 1 1];

cols = autumn(256);
cols = cols(112:240,:);
for idx = find(st==0)
  xx = [A(idx) x(idx) y(idx)]./scale;
  if isempty(XX) || min( sqrt(sum( ((XX-xx(o,:)).*RR(o,:)).^2 ,2)) )>0.5*R
    c = round(size(cols,1)*(2+x(idx))/8);
    surf(xx(1)+(R/RR(1))*X,xx(2)+(R/RR(2))*Y,xx(3)+(R/RR(3))*Z, ...
      'FaceColor', cols(c,:), 'LineStyle', 'none');
    XX = [ XX ; xx ]; %#ok<*AGROW>
    o  = [ o ; true];
  end
end

for idx = find(st>0)
  xx = [A(idx) x(idx) y(idx)]./scale;
  if isempty(XX) || min( sqrt(sum( ((XX-xx(o,:)).*RR(o,:)).^2 ,2)) )>0.5*R
    c = round(size(cols,1)*(2+x(idx))/8);
    surf(xx(1)+(R/RR(1))*X,xx(2)+(R/RR(2))*Y,xx(3)+(R/RR(3))*Z, ...
      'FaceColor', cols(c,:), 'LineStyle', 'none');
    XX = [ XX ; xx ];
    o  = [ o ; true];
  end
end

end
