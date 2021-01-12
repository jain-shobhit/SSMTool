function thm = hspo_plot_sol(thm, args, dim, bd, HH)

  function [var, var_f, thm, args] = coordinate(thm, args, bd, col, idx, lab)
    solvars = { 't'   'x'   'tn'  'p' };
    solflds = { 'tbp' 'xbp' 'tau' 'p' };
    cls = args.(col);
    id  = args.(idx);
    var   = cell(1,numel(cls));
    var_f = cell(1,numel(cls));
    for cl=1:numel(cls)
      mtc = strcmp(cls{cl}, solvars);
      if any(mtc(1:3))
        fld = solflds{mtc};
        if ~isa(id, 'function_handle')
          var_f{cl} = @(sol, j, o) sol.(fld)';
        else
          var_f{cl} = @(sol, j, o) sol.(fld);
        end
      elseif mtc(4)
        fld       = solflds{mtc};
        var_f{cl} = @(sol, j, o) sol.(fld);
      else
        output    = coco_bd_col(bd, cls{cl});
        var_f{cl} = @(sol, j, o) output(:,j)*o;
      end
    end
    if ~isa(id, 'function_handle')
      if any(strcmp(fld, {'xbp', 'p'}))
        thm.(lab) = sprintf('%s_%d', cls{1}, id);
      else
        thm.(lab) = cls{1};
      end
      args.(idx) = @(x) x(id,:);
    end
  end

[x, x_f, thm, args] = coordinate(thm, args, bd, 'col1', 'idx1', 'xlab');
[y, y_f, thm, args] = coordinate(thm, args, bd, 'col2', 'idx2', 'ylab');
if dim==3
  [z, z_f, thm, args] = coordinate(thm, args, bd, 'col3', 'idx3', 'zlab');
end

idx  = coco_bd_lab2idx(bd, args.labs);
type = coco_bd_col(bd, 'TYPE', false);

for i=1:numel(idx)
  row = idx(i);
  lab = args.labs(i);
  if ~isempty(args.oidx)
    lspec = thm.sol.RO;
  else
    lspec = thm.sol.RO;
    SP    = type{row}; % bug: use SP to select line style
    if ismember(SP, thm.special) && isfield(thm.sol, SP)
      lspec = thm.sol.(SP);
    end
  end
  for j=1:numel(args.oid)
    sol_hspo = hspo_read_solution(args.oid{j}, args.run, lab);
    o = ones(size(sol_hspo.tbp{1}'));
    for l=1:numel(sol_hspo.tbp)
      sol.tbp = sol_hspo.tbp{l};
      sol.xbp = sol_hspo.xbp{l};
      sol.tau = sol_hspo.tbp{l}/sol_hspo.T{l};
      sol.p   = sol_hspo.p;
      for k=1:numel(x)
        x{k} = x_f{k}(sol,row,o);
      end
      for k=1:numel(y)
        y{k} = y_f{k}(sol,row,o);
      end
      if dim==2
        plot(args.idx1(x{:}), args.idx2(y{:}), lspec{:});
      else
        for k=1:numel(z)
          z{k} = z_f{k}(sol,row,o);
        end
        plot3(args.idx1(x{:}), args.idx2(y{:}), args.idx3(z{:}), lspec{:});
      end
    end
    if HH; hold on; end
  end
end

end
