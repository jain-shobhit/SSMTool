function deactivate_parallel(~)
%DEACTIVATE_PARALLEL This function deletes parallel pool if it exists


try    
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        disp('----------------------------------------------------------------')
        disp('Delete parallel pool')
        disp('----------------------------------------------------------------')
        delete(poolobj)
    end
catch ME
    rethrow(ME)
end


end