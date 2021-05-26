function activate_parallel(varargin)
%ACTIVATE_PARALLEL This function starts parallel pool for the first time
%and specifies/detects the number of cores avaliable
%
% varargin - the number of cores used for parallel computing


try    
    pp1 = gcp('nocreate');
    if isempty(pp1)
        h = helpdlg('Starting parallel pool for the first time and detecting number of available cores.', 'Info');
        disp('----------------------------------------------------------------')
        disp('Starting parallel pool for the first time and detecting number')
        disp('of available cores.')
        disp('----------------------------------------------------------------')
        defaultProfile = parallel.defaultClusterProfile;
        myCluster = parcluster(defaultProfile);
        if numel(varargin)>1 && isnumeric(varargin{2})
            parpool(myCluster,varargin{2});
        else
            parpool(myCluster);
        end
        pp2 = gcp('nocreate');
        cpuNum =pp2.NumWorkers;
        save('cluster_info.mat','cpuNum')

        if isvalid(h)
            close(h)
        end
    end
catch ME
    rethrow(ME)
end


end