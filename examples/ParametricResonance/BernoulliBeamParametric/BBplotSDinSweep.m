function BBplotSDinSweep(cf,SD_full)


figure(cf)
xlim('manual')
ylim('manual')
% Reduced results
%plot3(SD.omega,SD.epsilon,zeros(1,numel(SD.epsilon)),'ob','MarkerSize',10,'DisplayName',strcat('SSM-$$\mathcal{O}(',num2str(SD.order),')$$'))

% Full results
%plot3(SD_full.Omega,SD_full.Epsilon,zeros(1,numel(SD_full.Epsilon)),'-k', 'LineWidth', 3,'DisplayName','Full System')


%% Fill region inside tongue

curve1 = [SD_full.Omega;SD_full.Epsilon;zeros(1,numel(SD_full.Epsilon))];
y2u     = max(SD_full.Epsilon) * ones(1,numel(SD_full.Omega));
curve2 = [SD_full.Omega; y2u;  zeros(1,numel(SD_full.Epsilon))];

Xu=[curve1(1,:),fliplr(curve2(1,:))];
Yu=[curve1(2,:),fliplr(curve2(2,:))];
Zu=[curve1(3,:),fliplr(curve2(3,:))];

colorSpec=[0.5,0.5,0.5];
fill3(Xu,Yu,Zu,colorSpec,'DisplayName','Unstable','FaceAlpha',.3)

%% Fill region that is stable
y2s     =  zeros(1,numel(SD_full.Omega));

curve2 = [SD_full.Omega; y2s;  zeros(1,numel(SD_full.Epsilon))];


Xs=[curve1(1,:),fliplr(curve2(1,:))];
Ys=[curve1(2,:),fliplr(curve2(2,:))];
Zs=[curve1(3,:),fliplr(curve2(3,:))];

colorSpec=[0.9,0.9,0.9];
fill3(Xs,Ys,Zs,colorSpec,'DisplayName','Stable','FaceAlpha',.3)
end