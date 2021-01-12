function [zInit,Omega] = read_po_ssm_init(labels)
% this function reads the full initial conditions on the periodic orbit 
% from the stored data during SSM FRC extraction
zInit = [];
Omega = [];
for j = labels
    load(['po_ssm' num2str(j) '.mat'],'Z','Omega_0');
    nRho = length(Omega_0);
    for k = 1:nRho
        zInit = [zInit Z{k}(:,1)];
        Omega = [Omega Omega_0(k)];
    end
end

