function stab_plot3(x,y,z,S,order,varargin)
% This function is adapted from stab_plot in coco_plot_bd in coco

II  = 1;
EI  = numel(x);
I   = II;
figs = [];
stab = [];
ST = cell(2,1);
ST{1} = {'b--','LineWidth',1.5};
ST{2} = {'b-','LineWidth',1.5};
while true
  if II>=EI; break; end
  [I, I0] = next_index(EI, S, I, II);
  fig = plot3(x(II:I0), y(II:I0), z(II:I0), ST{S(II)+1}{:});
  figs = [figs fig];
  stab = [stab S(II)+1];
  II = I0;
end

if isempty(varargin)
    Leg = cell(2,1);
    Leg{1} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable');
    Leg{2} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable');
    numSegs = numel(figs);
    if numSegs==1
        % add legend 
        legend(Leg{stab(1)},'Interpreter','latex');
    else
        legend(Leg{stab(1)},Leg{stab(2)},'Interpreter','latex');
        if numSegs>2            
            % remove repetitive legends from option
            for j=3:numSegs
                set(get(get(figs(j),'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off');
            end
        end     
    end
else
    numSegs = numel(figs);
    for j=1:numSegs
        set(get(get(figs(j),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end    
end
end


function [I, I0] = next_index(EI, S, I, II)
while (I<EI && S(II)==S(I)); I=I+1; end
I0 = I;
end