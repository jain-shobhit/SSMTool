if ~exist('sweep.mat', 'file')
  run sweep % make sure sweep.mat exists
end

%     A om  la al eps
p0 = [1  3 0.2  1   1]';

[data x0]      = exp_F([], zeros(7,1), p0);
data.cdata.k1  = -1;
data.cdata.k2  = -1;
data.TransTime = 8;

prob = coco_prob();

prob = coco_set(prob, 'cont', 'ItMX', -100);
prob = coco_set(prob, 'cont', 'NPR', 1, 'NSV', 10);
prob = coco_set(prob, 'cont', 'LogLevel', [3 3]);

prob = coco_set(prob, 'cont', 'h_min', 0.05);
prob = coco_set(prob, 'cont', 'h0'   , 0.10);
prob = coco_set(prob, 'cont', 'h_max', 0.15);
prob = coco_set(prob, 'cont', 'TrMX' , 3);

prob = coco_set(prob, 'cont', 'yidx', 2);
prob = coco_set(prob, 'cont', 'xidx', 13);
prob = coco_set(prob, 'cont', 'phan', subplot(1,2,1));

prob = coco_set(prob, 'corr', 'TOL', 1.0e-3);
prob = coco_set(prob, 'corr', 'ResTOL', 2.0e-4);
prob = coco_set(prob, 'corr', 'NPt', 3);
prob = coco_set(prob, 'corr', 'ItMN', 10);
prob = coco_set(prob, 'corr', 'ItMX', 20);
prob = coco_set(prob, 'corr', 'SubItMX', 3);
prob = coco_set(prob, 'corr', 'DampResMX', 1.0e-2);
prob = coco_set(prob, 'cont', 'corr.ItMN', 3);
prob = coco_set(prob, 'cont', 'corr.ItMX', 8);

prob = coco_set(prob, 'continex', 'NJac', 0);
prob = coco_set(prob, 'continex', 'JacH', [0.05 0.005]);

% prob = coco_set(prob, 'cont', 'xidx', 2);
% prob = coco_set(prob, 'cont', 'yidx', 5);
% prob = coco_set(prob, 'cont', 'phan', subplot(1,2,1));

prob = coco_set(prob, 'continex', 'phase', 1:7);
prob = coco_set(prob, 'continex', 'xpar', 'om');
prob = coco_set(prob, 'continex', 'ypar', '||C||');
prob = coco_set(prob, 'continex', 'phan', subplot(1,2,2));

prob = continex_add_recording(prob);

load sweep oms fws bws;
plot(oms, fws, 'k+', oms, bws, 'ko')
grid on
drawnow

% prob = coco_add_func(prob, 'tanh_ramp', 'continex', @tanh_ramp, [], ...
%   'active', 'ramp', 'vectorised', 'on');

bd1 = coco(prob, 'duff1', 'continex', 'isol', 'sol', ...
  @exp_F, 'fpar', data, x0, {'A' 'om' 'la' 'al' 'eps'}, p0, ...
  {'om' 'continex.phase'}, [0.2 3.1]);

bd = coco_bd_read('duff1');
om = coco_bd_col(bd, 'om');
A  = coco_bd_col(bd, '||C||');
clf
plot(om, A, 'b.-', oms, fws, 'k+', oms, bws, 'ko')
grid on
