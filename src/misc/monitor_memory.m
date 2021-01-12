function [ memory_in_use ] = monitor_memory(workspace)
%MONITOR_MEMORY uses the WHOS command to evaluate the total memory used by
% all the variables in the specified workspace ('caller' or 'base')
% The output is returned in MB.

mem_elements = evalin(workspace,'whos');
if size(mem_elements,1) > 0
    
    for i = 1:size(mem_elements,1)
        memory_array(i) = mem_elements(i).bytes;
    end
    
    memory_in_use = sum(memory_array);
    memory_in_use = memory_in_use/1048576;
else
    memory_in_use = 0;
end
end
