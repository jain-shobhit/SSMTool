function efunc_chkDeriv(opts, func, DFDX, DFDX2)
err = max(max(abs(DFDX-DFDX2)));
if err>=func.chkTOL
  msg = sprintf('<a href="matlab:open(''%s'')">%s</a>', ...
    mfilename('fullpath'), mfilename);
  msg = sprintf('%s:\nderivative of function ''%s'' might be incorrect', ...
    msg, func.identifyer);
  msg = sprintf('%s:\nerr [%.2e] >= TOL [%.2e]\n', msg, err, func.chkTOL);
  if func.chkTOL>0
    [er ec] = find(abs(DFDX-DFDX2)>=func.chkTOL);
    coords = sprintf('  %3d  %3d\n', [er(:) ec(:)]');
    msg = sprintf('%slarge errors found at:\n  row  col\n%s', msg, coords);
  end
  coco_warn(opts, 1, 1, '%s', msg);
end
end
