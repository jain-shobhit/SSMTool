function [data res] = add_IP(prob, data, command, varargin)
%ADD_IP   Slot function: Add to bifurcation data.
%
% Store trajectory end point at t=0 to bifurcation data

res = {};
switch command
  case 'init'
    res   = 'X0';
  case 'data'
    chart = varargin{1};
    [data uidx] = coco_get_func_data(prob, 'bvp.seg.coll', 'data', 'uidx');
    res  = chart.x(uidx(data.x0_idx));
end

end
