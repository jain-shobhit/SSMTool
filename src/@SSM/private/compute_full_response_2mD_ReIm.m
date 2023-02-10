function [Aout,Zout,Znorm,Zic] = compute_full_response_2mD_ReIm(W_0, W_1, state, epsilon, nt, mFreqs, outdof)
phi = linspace(0,2*pi,nt)/min(mFreqs);
%reduced coordinate
z = transpose(state).*exp(1i*(mFreqs'.*phi));
m = numel(mFreqs);
p = zeros(2*m,nt);
p(1:2:end-1,:) = z;
p(2:2:end,:)   = conj(z);
% 
% Full system response: $\mathbf{z}(t) = \mathbf{S}(\mathbf{p}(\Omega t), \Omega    
% t)$
z = reduced_to_full_traj(phi,p,W_0,W_1,epsilon,1);
%
% Timeseries at outdof
if isnumeric(outdof)
    Zout = z(outdof,:); noutdof = numel(outdof);
else
    Zout = outdof(z); noutdof = size(Zout,1); % function handle for observables
end
%
% $$\|\mathbf{z}\|_{L_2}=\sqrt{\frac{1}{T}\int_0^T\mathbf\|\mathbf{z}(t)\|^2dt}\approx 
% \frac{\|\mathbf{Z}\|_{F}}{\sqrt{n_t-1}}$$
Znorm = norm(z,'fro')/sqrt(nt-1);   
% Aout = norm(Zout,'inf');
Aout = [];
for k=1:noutdof
    amp = norm(Zout(k,:),'inf');
    Aout = [Aout amp];
end
Zic = z(:,1)';
end