function PBplotSD(n,cs)

order = 3;


figure(); hold on


%% Double Plot Setting
% {
axis_size = 20; % Set legend and Axis size
xwidth = 600; % window size
ywidth = 500;

ssm_sz = 2;
coco_sz = 5;
%}

% Axis limits
xmin = 10.3;
xmax = 11;
ymin = 0;
ymax = 1;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

strings = {[]};
for i = 1:numel(cs)
    strings{i} = strcat('SD_damp',num2str(cs(i)),'n',num2str(n),'.mat');
end

fullom= [];fulleps = [];

colors = [[0.9290    0.6940    0.1250]; [0.4940    0.1840    0.5560]; [0.4660    0.6740    0.1880]; [0.3010    0.7450    0.9330]];

%% Plot results
for i =1:numel(cs)
    load(strings{i})
    h = plot(SD_full.Omega,SD_full.Epsilon,'ok', 'MarkerSize', coco_sz,'DisplayName','Full System');
    
    if i < numel(cs)
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
for i = 1:numel(cs)
    load(strings{numel(cs)-i+1})
    % Reduced results
    plot(SD.omega,SD.epsilon,'-','Linewidth',ssm_sz,'color',colors(i,:),'DisplayName',strcat('SSM-$$\mathcal{O}(',num2str(order),') ,c = ' , num2str(cs(end+1-i)) , '$$'))
    
    % Full results
    
end




%{
%% Shading of the Figure
a = area(SD_full.Omega,[SD_full.Epsilon ; ones(1,numel(SD_full.Omega))*ymax*1.1 - SD_full.Epsilon].' ,'FaceAlpha',.3);

a(1).FaceColor = [0.9 0.9 0.9];
a(1).LineWidth = 0.001;
a(1).DisplayName = 'Stable';

a(2).FaceColor = [0.5 0.5 0.5];
a(2).LineWidth = 0.001;
a(2).DisplayName = 'Unstable';


%}
%% Set axis limits
add_labels('$\Omega$','$\epsilon \mu$')
add_legends();
xlim([xmin,xmax]);
ylim([ymin,ymax]);

%% Parameters for figure size, label size, legend size
set(gca,'fontsize',axis_size)

%% Size of window
set(gcf, 'Position', [100 100 xwidth ywidth])


end

function add_legends()
legend('show');
legend boxon;
end

function add_labels(xlab,ylab)
xlabel(xlab,'Interpreter','latex');
ylabel(ylab,'Interpreter','latex');

grid on, axis tight; legend boxoff;
end