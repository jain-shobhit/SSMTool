function [response,Znorm,Aout] = linear_response(DS,omegas)
% LINEAR ANALYSIS This function calcualtes forced response curve of linear
% part of dynamical system DS. The forced responses are calcuated at a
% collection of sampled omegas.
%
% [response,Znorm,Aout] = linear_response(DS,omegas)
%
% response: a cell array of periodic orbits
% Znorm: L2 norm of trajectory (for the whole state)
% Aout: amplitude at outdof

nOmega = length(omegas);
response = cell(nOmega,1);
Znorm = zeros(nOmega,1);
Aout = zeros(nOmega,1);

nt = 128;
phi = linspace(0,2*pi,nt);
outdof = DS.Options.outDOF;

for j = 1:nOmega
    if DS.order==2
        x_kappa = zeros(DS.n,DS.nKappa);
        for k = 1:DS.nKappa
            kOmega = DS.fext.kappas(k) * omegas(j);
            x_kappa(:,k) = (-kOmega^2 * DS.M + 1i * kOmega * DS.C + DS.K)\DS.fext.coeffs(:,k);
        end
        response{j} = DS.fext.epsilon * real( x_kappa * exp(1i * DS.fext.kappas * phi));   
    else
        z_kappa = zeros(DS.N,DS.nKappa);
        for k = 1:DS.nKappa
            kOmega = DS.Fext.kappas(k) * omegas(j);
            z_kappa(:,k) = (DS.A-1i*kOmega*DS.B)\DS.Fext.coeffs(:,k);
        end
        response{j} = DS.Fext.epsilon * real( z_kappa * exp(1i * DS.Fext.kappas * phi));    
    end   
%% 
% $$\|\mathbf{z}\|_{L_2}=\sqrt{\frac{1}{T}\int_0^T\mathbf\|\mathbf{z}(t)\|^2dt}\approx 
% \frac{\|\mathbf{Z}\|_{F}}{\sqrt{n_t-1}}$$
    Znorm(j) = norm(response{j},'fro')/sqrt(nt-1);    
    if ~isempty(outdof)
        Aout(j) = norm( response{j}(outdof,:),'inf');
    end
end


