function [Aout,Znorm,Omegas] =  read_num_int_sol(labels,outdof)
% this function reads the numerical integration solution
nOmega = length(labels);
Omegas = nan(nOmega,1);
Znorm = zeros(nOmega,1);
Aout = zeros(nOmega,1);
for j = labels
    load(['po' num2str(j) '.mat'],'x0','Omega','t0');
    Omegas(j) = Omega;
    nt = length(t0);
    Znorm(j) = norm(x0,'fro')/sqrt(nt-1);
    if ~isempty(outdof)
        Aout(j) = norm( x0(:,outdof),'inf');
    end
end

