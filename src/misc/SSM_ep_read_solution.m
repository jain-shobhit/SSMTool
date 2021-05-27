function sol = SSM_ep_read_solution(runid)
% SSM_EP_READ_SOLUTION This function reads contination solutions stored in
% disk. The solutions here are obtained from continuation of equilibrium
% points in run specified by runid
%
% SOL = SSM_EP_READ_SOLUTION(RUN,VARARGIN)
%
% run:      runid of continuation
% 
% See also: SSM_PO_READ_SOLUTION

sfname = 'SSMep.mat';
runid  = [runid,'.ep'];
sfname = coco_fname(runid, sfname);
FRC = load(sfname);
FRC = FRC.FRC;
sol = FRC;

end
