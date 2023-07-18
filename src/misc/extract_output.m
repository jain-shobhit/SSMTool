function [Zout, Aout, Znorm, varargout] = extract_output(z, outdof)
% Aout : Max z at output DOFs
% Zout : Time history at output DOFs
% Znorm: L_2 norm estimate of the periodic response

nt = size(z,2);
% $$\|\mathbf{z}\|_{L_2}=\sqrt{\frac{1}{T}\int_0^T\mathbf\|\mathbf{z}(t)\|^2dt}\approx
% \frac{\|\mathbf{Z}\|_{F}}{\sqrt{n_t-1}}$$
Znorm = norm(z,'fro')/sqrt(nt-1);

if ~isempty(outdof)
    Zout = z(outdof,:);    
    noutdof = numel(outdof);
    Aout = zeros(1, noutdof);
    for k=1:numel(outdof)
        Aout(k) = norm(Zout(k,:),'inf');
    end
    ZoutNorm = norm(z(outdof,:),'fro')/sqrt(nt-1);
    varargout{1} = ZoutNorm;
else
    Zout = [];
    Aout = [];
end
end
