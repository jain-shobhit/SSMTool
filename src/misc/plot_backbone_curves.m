function plot_backbone_curves(coeffs,exp_idx,rhosamp,orders)
% This function plot the backbone curves at reduced coordinates with
% various orders. It only supports 2D SSM

if exp_idx(1)>exp_idx(2)
    exp_idx = flip(exp_idx);
    coeffs  = flip(coeffs);
end
figure; hold on
y = 0;
for k=1:numel(exp_idx)
    y = y+coeffs(k)*rhosamp.^(exp_idx(k));
    orderk = exp_idx(k)+1;
    if any(orderk==orders)
        plot(y,rhosamp,'LineWidth',1.5,'DisplayName',['SSM-O(',num2str(orderk),')']);
    end
end
legend('show','Location', 'Best'); 
legend boxoff
set(gca,'FontSize',14);
set(gca,'LineWidth',1.2);
grid on
xlabel('$\Omega$','Interpreter',"latex",'FontSize',16);
ylabel('$\rho$','Interpreter',"latex",'FontSize',16);
end
