function sol = SSM_po_read_solution(runid, varargin)
% SSM_PO_READ_SOLUTION This function reads contination solutions stored in
% disk. The solutions here are obtained from continuation of periodic
% orbits in run with runid. Each periodic orbit in the slow time scale
% corresponds to a torus or quasi-periodic orbit in the regular time
% dynamics. By specifying the label of solution in varargin, only this
% solution will be included in sol. Otherwise, all solutions will be
% included in sol.
%
% SOL = SSM_PO_READ_SOLUTION(RUN,VARARGIN)
%
% runid:    runid of continuation of periodic orbits
% varargin: label - label of solution
%
% See also: SSM_EP_READ_SOLUTION

sfname = 'SSMpo.mat';
runid  = [runid,'.po'];
sfname = coco_fname(runid, sfname);
FRC = load(sfname);
FRC = FRC.FRC;
if isempty(varargin)
    sol = struct();
    sol.om = FRC.om;
    sol.ep = FRC.ep;
    sol.st = FRC.st;
    sol.lab = FRC.lab;
else
    label = varargin{1};
    idx = find(FRC.lab==label);
    sol = struct();
    sol.om  = FRC.om(idx);
    sol.ep  = FRC.ep(idx);
    sol.st  = FRC.st(idx);
    sol.tTr = FRC.tTr{idx};
    sol.xTr = FRC.zTr{idx};
    sol.z0  = FRC.Z0_frc(:,idx);
    posol   = po_read_solution(runid,label);
    assert(abs(posol.p(1)-sol.om)<=1e-5*abs(sol.om),...
        'the solution read from po is inconsistent with the one from SSM');
    sol.tPo = posol.tbp;
    sol.xPo = posol.xbp;
    sol.Tpo = posol.T;
end

info = struct();
info.SSMorder   = FRC.SSMorder;
info.SSMnonAuto = FRC.SSMnonAuto;
info.SSMispolar = FRC.SSMispolar;
sol.info = info;

end
