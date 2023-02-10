function [Aout, Zout,Znorm,varargout] = compute_full_response_traj(W_0, W_1, epsilon, t, p, om, outdof)
% Full system response: $\mathbf{z}(t) = \mathbf{S}(\mathbf{p}(\Omega t), \Omega    
% t)$
z = reduced_to_full_traj(t,p,W_0,W_1,epsilon,om);
%
% Timeseries at outdof
if isnumeric(outdof)
    Zout = z(outdof,:); noutdof = numel(outdof);
else
    Zout = outdof(z); noutdof = size(Zout,1);
end
%
% $$\|\mathbf{z}\|_{L_2}=\sqrt{\frac{1}{T}\int_0^T\mathbf\|\mathbf{z}(t)\|^2dt}$$
Znorm = sqrt(trapz(t, sum(z.^2))/(t(end)-t(1)));
% Aout = norm(Zout,'inf');
Aout = [];
for k=1:noutdof
    amp = norm(Zout(k,:),'inf');
    Aout = [Aout amp];
end
varargout{1} = z(:,1);
end