function [fdata,data_dir] = create_reduced_dynamics_data(beta,kappa,lambda,mFreqs,Nonauto,W_0,W_1,order,resonant_modes)
% CREATE_REDUCED_DYNAMICS_DATA This function create a data structure for
% the vector field of reduced dynamics (leading-order nonautonomous
% approximation)

fdata = struct();
fdata.beta  = beta;
fdata.kappa = kappa;
fdata.lamdRe = lambda.lambdaRe(1:2:end-1);
fdata.lamdIm = lambda.lambdaIm(1:2:end-1);
fdata.mFreqs = mFreqs;
fdata.iNonauto = Nonauto.iNonauto;
fdata.rNonauto = Nonauto.rNonauto;
fdata.kNonauto = Nonauto.kNonauto;
% put W_0 and W_1 in fdata is a bad idea because it will be stored in disk
% for each saved continuation solution. As an alternative, we save W_0 and
% W_1 in disk here under folder data. When needed, they will be loaded into
% memory.
% fdata.W_0   = W_0;
% fdata.W_1   = W_1;
data_dir = fullfile(pwd,'data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
wdir = fullfile(data_dir,'SSM.mat');
SSMcoeffs = struct();
SSMcoeffs.W_0 = W_0;
SSMcoeffs.W_1 = W_1;
save(wdir, 'SSMcoeffs');
fdata.order = order;
fdata.modes = resonant_modes;

end