function [element,eid,type]=readelements(meshFile,pids)


%       Reads the element connectivity matrix from a Gmsh ver 2.0 ASCII file
%       Reads the element connectivity matrix from a Gmsh file with 
%       physical ids in pid.
%


% open the file
fid=fopen(meshFile,'r');
if ( fid<0 )
  disp(['Error could not open file ',meshFile]);
  element=[];
  eid=[];
  type=[];
  return
end


while 1
  
  line=fgetl(fid);            % %  fgetl(FID) returns the next line of a file associated with file identifier FID as a MATLAB string
  if ~ischar(line), break, end %  check End of file - if line is not string.. break and end loop
  
  switch line                 % find seciton
    
  case '$Elements'          %if Line has "$Elements, elements list has started
    n=str2num(fgetl(fid));  %read no. of elements from next line annd convert to interger
    
    ne=0; % no. of elements counter
    for i=1:n
      temp=str2num(fgetl(fid));  % get element
 
      pid=temp(4); % pick physical Id of Element
      
      if ( nargin>1 ) %    nargin Number of function input arguments.
                          % Inside the body of a user-defined function, nargin returns
                          %   the number of input arguments that were used to call thefunction.
        if ( nnz( ismember(pids,pid) ) == 0 )  % If PID od element doesnt match the required PID, skip the line
                                              % nnz(S) is the number of nonzero elements in S
                                              % ismember(A,B) for arrays A and B returns an array of the same
                                              % size as A containing true where the elements of A are in B and false otherwise.
                                             
          continue  % skip this element
        end
      end
      % found an element with matching pid
      ne=ne+1; %increase element count with matching PID
      
      eid(ne)=temp(1); % fill element ID
      etype=temp(2);   % fill element type
      [~,nn]=etypestr(etype); %returns Type(string) and no. of nodes per element.
%       type{ne}=estr;
      type{ne}=etype;
      ntags=temp(3);   %No. of tags associated with the element
      element(ne,1:nn)=temp(4+ntags:3+ntags+nn); % pick nodes from the element ID (which start after the tags of the element)
      
    end
    
  otherwise
    % skip line
  end
  
end

fclose(fid);

