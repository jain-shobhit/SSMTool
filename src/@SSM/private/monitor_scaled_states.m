function [prob, args1, args2] = monitor_scaled_states(prob, ispolar, m, scales)
% MONITOR_STATES This function add state of reduced dynamics as
% continuation parameters and define string array for the state. The prob
% here is a continuation problem. (args1,args2)=(rhoargs,thargs) or 
% (Reargs,Imargs) depending on the value of ispolar

args1 = cell(m,1);
args2 = cell(m,1);
if ispolar
    for k=1:m
        args1{k} = strcat('rho',num2str(k));
        args2{k}  = strcat('th',num2str(k));
    end
    parsRData = struct(); parsTHData = struct();
    parsRData.scale  = scales(1:2:end);
    parsTHData.scale = scales(2:2:end);
    prob = coco_add_func(prob, 'radius', @scale_pars, parsRData,...
        'active', args1(:)', 'uidx', 1:2:2*m-1);
    prob = coco_add_func(prob, 'angle', @scale_pars, parsTHData,...
        'active', args2(:)', 'uidx', 2:2:2*m);    
else
    for k=1:m
        args1{k} = strcat('Rez',num2str(k));
        args2{k} = strcat('Imz',num2str(k));
    end
    parsReData = struct(); parsImData = struct();
    parsReData.scale = scales(1:2:end);
    parsImData.scale = scales(2:2:end);
    prob = coco_add_func(prob, 'realParts', @scale_pars, parsReData,...
        'active', args1(:)', 'uidx', 1:2:2*m-1);
    prob = coco_add_func(prob, 'imagParts', @scale_pars, parsImData,...
        'active', args2(:)', 'uidx', 2:2:2*m); 
end 

end

function [data, y] = scale_pars(prob, data, u)
% CAL_RHOS: This function converts cartesian coordinates to polar
% coordinats and return radius
%
% [DATA, Y] = CAL_RHOS(PROB, DATA, U)
%

u = u(:);
k = data.scale(:);
y = k.*u;

end