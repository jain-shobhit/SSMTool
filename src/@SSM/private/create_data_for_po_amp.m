function optData = create_data_for_po_amp(cind,dind,mFreqs,wcoeffs,optdof,obj,Flead,varargin)
% construct data structure for po_amp

m = numel(mFreqs);
rhat = (cind-dind)*mFreqs(:);
[val,~,ic] = unique(rhat,'stable');
pos = cell(numel(val),1);
for k=1:numel(val)
    pos{k} = find(ic==k);
end
hatr = struct();
hatr.val = val;
hatr.pos = pos;
optData = struct();
optData.hatr = hatr;
optData.cind = cind;
optData.dind = dind;
optData.coeffs = wcoeffs;
optData.optdof = optdof;
isl2norm = false;
if ~isempty(varargin)
    isl2norm = varargin{1}; optData.isl2norm = isl2norm;
end
if isempty(obj.FRCOptions.DBCobjweight)
    optData.Q    = eye(numel(optdof));
    optData.Qbar = optData.Q;
else
    optData.Q    = obj.FRCOptions.DBCobjweight;
    optData.Qbar = 0.5*(optData.Q+optData.Q');
end
uidxpo = (1:2*m)';
isnonauto = obj.Options.contribNonAuto;
optData.isnonauto   = isnonauto;
optData.coordinates = obj.FRCOptions.coordinates;
if isnonauto
    idk = obj.System.kappas==1;
    optData.A = obj.System.A;
    optData.B = obj.System.B;
    optData.kappa = kappa_pos;
    optData.Flead = Flead(:,idk);
    optData.isbaseForce = obj.System.Options.BaseExcitation;
    uidxpo = [uidxpo; 2*m+1; 2*m+2]; %[omega,epsilon]
else
    if ~isl2norm; uidxpo = [uidxpo; 2*m+1]; end% omega
end
optData.uidxpo = uidxpo;

end