function [dmn_ids,ebc_ids,nbc_ids]=n_pid(meshFile)
%function reads the mesh file and returns 
%the different physical domain ids dmn_ids , (these are numbered from 1-100) 
%the different natural BC ids - nbc_ids,(these are numbered from 101-200)
%the different essential BC ids - ebc_ids (these are numbered from 201-300)
% open the file


fid=fopen(meshFile,'r');
if ( fid<0 )
  disp(['Error could not open file ',meshFile]);
  dmn_ids=[];
    ebc_ids=[];
    nbc_ids=[];
  return
end

while 1
  
  line=fgetl(fid);            % %  fgetl(FID) returns the next line of a file associated with file identifier FID as a MATLAB string
  if ~ischar(line), break, end %  check End of file - if line is not string.. break and end loop
  
  switch line                 % find section
    
  case '$Elements'          %if Line has "$Elements, elements list has started
    n=str2num(fgetl(fid));  %read no. of elements from next line annd convert to interger
    dmn_ids=[];
    ebc_ids=[];
    nbc_ids=[];
    for i=1:n
      temp=str2num(fgetl(fid));  % get element
 
      pid=temp(4); % pick physical Id of Element
      if pid<=100                                 % if pid<=100, then it is a domain ID
          if  (nnz( ismember(dmn_ids,pid) ) == 0) %if pid not already present in list of domain ids, then add it to the list
              dmn_ids=[dmn_ids;pid];
              continue
          end
          
      else if (pid>100 && pid<=200)                        % if 100<pid<=200, then it is an  EBC ID
          if  (nnz( ismember(ebc_ids,pid) ) == 0) % if pid not already present in list of ebc ids, then add it to the list
              ebc_ids=[ebc_ids;pid];
              continue
          end
          
          else if (pid>200 && pid<=300)                  % if 200<pid<=300, then it is an  NBC ID
          if  (nnz( ismember(nbc_ids,pid) ) == 0)% if pid not already present in list of nbc ids, then add it to the list
              nbc_ids=[nbc_ids;pid];
              continue
          end
          
              else
                   disp(['Error : PIDS not in range: Check MESHFILE ',meshFile]);
                   dmn_ids=[];
                   ebc_ids=[];
                   nbc_ids=[];
                   return
              end
          end
      end
    end
    otherwise
    % skip line
  end
  
end

fclose(fid);
                  
      
      
      
      
      
      
      
