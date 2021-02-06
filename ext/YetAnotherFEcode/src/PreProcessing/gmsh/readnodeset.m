function nids=readnodeset(meshFile,pids)


% function nids=readnodeset(meshFile,pids)
%
%       Reads the nodes ids belonging to physical id pid from a Gmsh file
%


% open the file
fid=fopen(meshFile,'r');
if ( fid<0 )
  disp(['Error could not open file ',meshFile]);
  nids=[];
  return
end


 
% read sections

conn=readelements(meshFile,pids);
nids=unique(conn); %C=unique(A) for the array A returns the same values as in A but with no repetitions. C will be sorted.

fclose(fid);
