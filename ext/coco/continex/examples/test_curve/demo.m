fprintf('************************************************************\n\n');

% a = 0.001;
% f = @(x,p) (p(1,:)+a*randn(size(x))).^2 - ...
%   25*(x(1,:)+a*randn(size(x)))^2 + ...
%   24*(x(1,:)+a*randn(size(x))).^4-15;

% a = 0.0;
a = 0.2;
f = @(x,p) p(1,:).^2-25*x(1,:)^2+24*x(1,:).^4-15 + a*randn(size(x));

figure(1)
clf
prob = coco_prob();

% prob = coco_add_event(prob, 'UZ', 'H',  0);

prob = coco_set(prob, 'all', 'TOL', max(0.5*a,1.0e-3));

% prob = coco_set(prob, 'cont', 'ela', '-');
% prob = coco_set(prob, 'cont', 'shape', 'parabola');
prob = coco_set(prob, 'cont', 'shape', 'egg');

% test large step size with default values
prob = coco_set(prob, 'cont', 'ItMX', 150);
prob = coco_set(prob, 'cont', 'h_min', 0.075);
prob = coco_set(prob, 'cont', 'h0'   , 0.15 );
prob = coco_set(prob, 'cont', 'h_max', 0.20 );
prob = coco_set(prob, 'cont', 'TrMX', 3);
% prob = coco_set(prob, 'cont', 'GAPMX', 2);
% prob = coco_set(prob, 'cont', 'corr2.GAP', true);

% test medium step size with 3 sample points
% prob = coco_set(prob, 'cont', 'ItMX', 200);
% prob = coco_set(prob, 'cont', 'h_min', 0.05);
% prob = coco_set(prob, 'cont', 'h0',    0.1 );
% prob = coco_set(prob, 'cont', 'h_max', 0.15);
% prob = coco_set(prob, 'cont', 'TrMX', 3);

% test small step size with 5 sample points
% prob = coco_set(prob, 'cont', 'ItMX', 300);
% prob = coco_set(prob, 'cont', 'h_min', 0.025);
% prob = coco_set(prob, 'cont', 'h0'   , 0.05 );
% prob = coco_set(prob, 'cont', 'h_max', 0.1  );
% prob = coco_set(prob, 'cont', 'TrMX', 5);
% prob = coco_set(prob, 'cont', 'NPtAv', 3);

prob = coco_set(prob, 'cont', 'LogLevel', [3 3]);

prob = coco_set(prob, 'corr', 'NPt', 3);
prob = coco_set(prob, 'corr', 'ItMX', 25);
prob = coco_set(prob, 'corr', 'ItMN', 5);
prob = coco_set(prob, 'corr', 'SubItMX', 3);
prob = coco_set(prob, 'cont', 'corr.ItMX', 12);
prob = coco_set(prob, 'cont', 'corr.ItMN', 3);
prob = coco_set(prob, 'cont', 'corr.MADSteps', 4);

prob = coco_set(prob, 'corr', 'phan', [subplot(3,2,5) subplot(3,2,6)]);

prob = coco_set(prob, 'cont', 'xidx', 3);
prob = coco_set(prob, 'cont', 'yidx', 1);
prob = coco_set(prob, 'cont', 'phan', subplot(3,2, [1 3]));

prob = coco_set(prob, 'continex', 'phan', subplot(3,2, [2 4]));
prob = coco_set(prob, 'continex', 'xpar', 'a');
prob = coco_set(prob, 'continex', 'ypar', 'FC');
prob = coco_set(prob, 'continex', 'pfunc', @(x) x);

prob = coco_set(prob, 'continex', 'NJac', 50);
prob = coco_set(prob, 'continex', 'JacH', [0.01 0.001]);

bd1 = coco(prob, 'test1', 'continex', 'isol', 'sol', f, ...
	[1.2], 'a', [1], 'a', [-6, 6]); %#ok<NBRAK>

y = coco_bd_col(bd1, 'FC');
x = coco_bd_col(bd1, 'a');

clf
plot(x, y, 'b.-', 1,1,'r*');
grid on
drawnow
axis equal
