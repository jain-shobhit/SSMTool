function [rho] = compute_analyticity_domain(obj,appr_order)
% This function  computes the roots of the function a(rho) and returns the
% radius of convergence as computed at approximation order appr_order.
% Detailed information S. Ponsioen et al. 'Analytic prediction of isolated
% forced response curves from spectral submanifolds' and sources given
% therein.

[~,R_0] = obj.compute_whisker(appr_order);

% Define colormap - grey shades increasing in darkness with higher order approximations.
color = flip(gray(appr_order-floor(appr_order/2)));

run_idx = 1;
for order = floor(appr_order/2):appr_order
    p = zeros(1,1+order);
    for j = 1:order
        if any(any(R_0{j}.coeffs))
            p(j+1) =   real(nonzeros(R_0{j}.coeffs(1,:)));
        end
    end
    
    % Compute roots
    p_roots = roots(flip(p));
    
    % Create plot
    if order < appr_order
        plot(p_roots, 'm.','LineWidth',2,'color',color(run_idx,:),'MarkerSize',10,'HandleVisibility','off');
    else
        plot(p_roots, 'm.','LineWidth',2,'DisplayName','Roots of $$a(\rho)$$','color','m','MarkerSize',10);
        
    end
    hold on;
    run_idx = run_idx+1;
end

title("Approximation order " + appr_order, 'Fontsize',20);
legend('show','Fontsize',20)
ylabel('Imaginary part','Fontsize',20);xlabel('Real part','Fontsize',20);

% Domain of analyticity as given by mean of the roots on the circle
% approximating radius of convergence

% We assume there are more non persistent zeros that approximate the radius
% of convergence than spurious zeros -> therefore the median distance of all
% roots will be bigger than the distance of the non spurious zero that is
% furthest away from the origin.

% Take the median distance from origin of all roots
rho_mean = median(abs(p_roots));

% Take the mean of all the roots that have bigger distance from zero than
% the median distance
rho = mean(abs(p_roots(abs(p_roots)>rho_mean)));


if ~isnan(rho)
    xlim([-1.5*rho 1.5*rho])
    ylim([-1.5*rho 1.5*rho])
end

hold off;
end
