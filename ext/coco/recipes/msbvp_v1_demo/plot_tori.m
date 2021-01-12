function plot_tori(run, fac)

if nargin<2
  fac=0.75;
end

bd = coco_bd_read(run);  % Extract bifurcation data
labs = coco_bd_labs(bd); % Extract solution labels
for lab = labs
  fprintf('plotting solution label %d\n', lab);
  [sol data] = msbvp_read_solution('', run, lab); % Extract trajectory segments
  N = data.nsegs;
  M = ceil(fac*size(sol{1}.x,1));
  clf
  hold on
  x0 = zeros(N+1,3);
  x1 = zeros(N+1,3);
  XX = zeros(M,N+1);
  YY = XX;
  ZZ = XX;
  for i=1:N+1
    n       = mod(i-1,N)+1;
    XX(:,i) = sol{n}.x(1:M,1);
    YY(:,i) = sol{n}.x(1:M,2);
    ZZ(:,i) = sol{n}.x(1:M,3);
    x0(i,:) = sol{n}.x(1,:);
    x1(i,:) = sol{n}.x(M,:);
    % plot3(sol{i}.x(1:25,1), sol{i}.x(1:25,2), sol{i}.x(1:25,3),'g-');
  end
  surf(XX, YY, ZZ, 'FaceColor', 0.9*[1 1 1], 'MeshStyle', 'column');
  plot3(x0(:,1), x0(:,2), x0(:,3), 'k.-', 'LineWidth', 2);
  plot3(x1(:,1), x1(:,2), x1(:,3), 'k.-', 'LineWidth', 0.5);
  hold off
  grid on
  axis([-1.5 1.5 -2 2 -0.5 2]);
  view([50 15]);
  drawnow
  pause(0.1)
end

end
