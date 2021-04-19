function plot_stab_lines(x,y,S,ST,varargin)
% PLOT_STAB_LINES This function is adapted from stab_plot in coco_plot_bd
% in coco. It plots continuation results with stability info.
%
% PLOT_STAB_LINES(X,Y,S,ST,VARARGIN)
%
% X:        Data for x axis
% Y:        Data for y axis
% S:        Data for stability of solution
% ST:       A cell array for line style, for instance, we may have ST = cell(2,1);
%             ST{1} = {'r--','LineWidth',1.0}; % unstable
%             ST{2} = {'r-','LineWidth',1.0};  % stable
% VARARGIN: Strings for legends of unstable/stable solution branches
%

oldlegs = get(gca, 'legend');
if ~isempty(oldlegs)
    legsold = oldlegs.String;
else
    legsold = {};
end
    
II  = 1;
EI  = numel(x);
I   = II;
figs = [];
stab = [];
while true
  if II>=EI; break; end
  [I, I0] = next_index(EI, S, I, II);
  if isnan(S(II))
    fig = plot(x(II:I0), y(II:I0), ST{S(II+1)+1}{:});
    stab = [stab S(II+1)+1];
  else
    fig = plot(x(II:I0), y(II:I0), ST{S(II)+1}{:});
    stab = [stab S(II)+1];
  end
  figs = [figs fig];
  
  II = I0;
end

if ~isempty(varargin)
    Leg = cell(2,1);
    Leg{1} = varargin{1};
    Leg{2} = varargin{2};
    numSegs = numel(figs);

    if numSegs==1
        % add legend 
        legend(legsold{:},Leg{stab(1)},'Interpreter','latex');
    else
        leg2 = setdiff(unique(stab),stab(1));
        legend(legsold{:},Leg{stab(1)},Leg{leg2},'Interpreter','latex');
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
if isnan(S(II)) % BUG: does not support two adjacent nan
    II=II+1;
    I=I+1;
end
while (I<EI && S(II)==S(I)); I=I+1; end
I0 = I;
end