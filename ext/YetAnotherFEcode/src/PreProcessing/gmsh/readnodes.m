function [nids,node]=readnodes(meshFile)

% open the file
fid=fopen(meshFile,'r'); % 'r'     open file for reading
if ( fid<0 )
  disp(['Error could not open file ',meshFile]);
  node=[];
  nids=[];
  return
end


while 1
  
  line=fgetl(fid);            %  fgetl(FID) returns the next line of a file associated with file identifier FID as a MATLAB string
  if ~ischar(line), break, end % check if End of file - if line is not string.. break and end loop
  
  switch line                 % find section

    
  case '$Nodes'               % if Line has "$Nodes"
    numnode=str2num(fgetl(fid)); % read no. of nodes from next line annd convert to interger
    node=zeros(numnode,3);     % % Node coordinates matrix initialization
    nids=zeros(numnode,1);     % Node id Vector
    
    for i=1:numnode 
      nodeline=str2num(fgetl(fid)); %get lines of nodes
      node(i,:)=nodeline(2:4); % fill node coordinates in node matrix
      nids(i)=nodeline(1);     % fill node id in nid vector
    end

  otherwise
    % skip line
  end
  
end

fclose(fid);
