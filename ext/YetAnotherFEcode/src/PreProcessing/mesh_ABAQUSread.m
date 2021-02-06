% MESH_ABAQUSREAD
%
% [nodes, elements, nset, elset, nset_name, elset_name] = mesh_ABAQUSread(filename,savemfile)
%
% Description: read an ABAQUS .inp file to extract nodes and elements.
% The function first reads the input file, then writes and execute a 
% function (named "fun_abaqus2matlab_mesh.m") which produces the outputs.
% This support function is then deleted.
% INPUTS    filename: name of the .inp file to convert (with or without
%                     extension
%           savemat: (optional) if set to 1, outputs are stored in a .mat 
%                     file. Since for large structures the conversion 
%                     process could take a while, at times it may result 
%                     convenient to save and run the .mat file. If 
%                     "filename.mat" is already available in the folder,
%                     this is directly loaded and conversion is aborted.
% OUTPUTS   nodes --> [nodeID, x, y, z]
%           elements --> [elementID, node1_ID ... nodeN_ID]
%           nset --> struct containing node sets (IDs)
%           elset --> struct containing element sets (IDs)
%           nset/elset_name --> struct with set names as defined in ABAQUS
%
% Additional notes: at the time being, this function handles only
% mono-element meshes, for "nodes" and "elements" are matrices and cannot
% have different row sizes. One could fix this problem by using cell arrays
%
% Author: Jacopo Marconi
% Last modified: 12/AUG/2019

function [nodes, elements, nset, elset, nset_name, elset_name] = mesh_ABAQUSread(filename,savemat)

tic
fprintf(' Converting input file...')

if nargin<2
    savemat=0;
end

% check the input file extension (only .inp is supported). If filename does
% not contain the extension, it is assumed it is .inp (this may lead to
% errors)
ind = strfind(filename,'.');
if length(ind)>1
    error(' filename contains too many "."')
end
if contains(filename,'.')
    Ext = filename(ind:end);
    filename(ind:end) = [];
    if ~strcmp(Ext,'.inp')
        error(' Unsupported file extension')
    end
    fid = fopen([filename Ext]);
else
    Ext = '.inp';
    fid = fopen([filename Ext]);
end

% if a .mat file dubbed filename already exist, load it and return
if exist([filename '.mat'],'file')
    load([filename '.mat'],'nodes','elements','nset','elset','nset_name','elset_name')
    fprintf(' %.2f s \n',toc)
    fprintf([' ' filename '.mat found and loaded.\n\n'])
    return
end

% create output file
filename_output = 'fun_abaqus2matlab_mesh.m';
fid2 = fopen(filename_output,'w');

fprintf(fid2,'function [nodes, elements, nset, elset, nset_name, elset_name] = fun_abaqus2matlab_mesh() \n');
fprintf(fid2,' nodes=[]; elements=[]; nset=[]; elset=[]; nset_name=[]; elset_name=[]; \n\n');

flag_nodes = 0; % flag: when set to 1, looks for the nodes endline
flag_elem = 0;	% flag: when set to 1, looks for the elements endline
flag_nset = 0;  % flag: when set to 1, looks for the nset endline
count_nset = 0; % node sets counter
flag_elset = 0; % flag: when set to 1, looks for the elset endline
flag_end = 0;   % when = 1, only remaining n/el-sets will be written
flag_multline = 0;
count_elset =0; % element sets counter
a = 'a';        % initialize a
while ischar(a)
    a = fgetl(fid);         % read (successive) line
    
    % check Abaqus headers
    if strcmp(a(1),'*')
        
        a1 = ['% ' a]; % modify headers to matlab comments
        
        % check nodes
        if strcmp(a,'*Node')
            flag_nodes = 1;
            fprintf(fid2,'%s\n',a1);
            fprintf(fid2,'nodes = [...\n');
            continue
        end
        if flag_nodes
            flag_nodes = 0;
            fprintf(fid2,'];\n');
        end
        
        % check elements
        if strcmp(a(1:min(8,length(a))),'*Element')
            flag_elem = 1;
            fprintf(fid2,'%s\n',a1);
            fprintf(fid2,'elements = [...\n');
            
            % chek multiple lines
            aa1 = fgetl(fid);
            aa2 = fgetl(fid);
            if length(aa1)>length(aa2)
                flag_multline = 1;
                fprintf(fid2,'%s\n',[aa1 aa2]);
            else
                fprintf(fid2,'%s\n',aa1);
                fprintf(fid2,'%s\n',aa2);
            end
            
            continue
        end
        if flag_elem
            flag_elem = 0;
            fprintf(fid2,'];\n');
            flag_end = 1;
            flag_multline = 0;
        end
        
        % check nsets
        if strcmp(a(1:min(5,length(a))),'*Nset')
            % in case a nset starts just after an elset:
            if flag_elset
                fprintf(fid2,'];\n\n'); % close elset
                flag_elset = 0;
            end
            % in case a nset starts just after another nset:
            if flag_nset
                fprintf(fid2,'];\n\n'); % close nset
            end
            % example of header: '*Nset, nset=Set-50, generate'
            flag_nset = 1;
            count_nset = count_nset+1;
            fprintf(fid2,'%s\n',a1);
            nname = replace(a,'*Nset, nset=','');
            nname = replace(nname,', generate','');
            fprintf(fid2,'nset_name{%d} = ''%s'';\n',count_nset,nname);
            fprintf(fid2,'nset{%d} = [...\n',count_nset);
            continue
        end
        if flag_nset
            flag_nset = 0;
            fprintf(fid2,'];\n\n');
        end
        
        % check elsets
        if strcmp(a(1:min(6,length(a))),'*Elset')
            % in case an elset starts just after a nset:
            if flag_nset
                fprintf(fid2,'];\n\n'); % close nset
                flag_nset = 0;
            end
            % in case a nset starts just after another nset:
            if flag_elset
                fprintf(fid2,'];\n\n'); % close elset
            end
            % example of header: '*Elset, elset=Set-50, generate'
            flag_elset = 1;
            count_elset = count_elset+1;
            fprintf(fid2,'%s\n',a1);
            elname = replace(a,'*Nset, nset=','');
            elname = replace(elname,', generate','');
            fprintf(fid2,'elset_name{%d} = ''%s'';\n',count_elset,elname);
            fprintf(fid2,'elset{%d} = [...\n',count_elset);
            continue
        end
        if flag_elset
            fprintf(fid2,'];\n\n');
            flag_elset = 0;
        end
        
        fprintf(fid2,'%s\n',a1);
        continue
    end
    
    if flag_nset
        fprintf(fid2,'%s ... \n',a);
        continue
    end
    
    if flag_elset
        fprintf(fid2,'%s ... \n',a);
        continue
    end
    
    if flag_end
        if ~(flag_nset || flag_elset)
            a1 = ['% ' a]; % modify headers to matlab comments
            fprintf(fid2,'%s\n',a1);
            continue
        end
    end
    
    % write the line read from input file
    if flag_multline == 0
        fprintf(fid2,'%s\n',a);
    else
        aa1 = a;
        aa2 = fgetl(fid);
        fprintf(fid2,'%s\n',[aa1 aa2]);
    end
end

fclose(fid);
fclose(fid2);
fclose('all');

[nodes, elements, nset, elset, nset_name, elset_name] = fun_abaqus2matlab_mesh();
delete('fun_abaqus2matlab_mesh.m')

nodes = nodes(:,2:end);
elements = elements(:,2:end);

if savemat==1
    save(filename,'nodes','elements','nset','elset','nset_name','elset_name')
end

fprintf(' %.2f s \n',toc)
