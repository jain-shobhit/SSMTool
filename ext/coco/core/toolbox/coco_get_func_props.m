function varargout = coco_get_func_props(opts, prname, varargin)

varargout = {};
funcs     = opts.efunc.funcs;

table = cell(0,nargout);
ii = 1;
for i=1:numel(funcs)
  func = funcs(i);
  if isfield(func.props, prname)
    table{ii,1} = func.identifyer;
    table{ii,2} = func.props.(prname);
    for j = 3:nargin
      switch lower(varargin{j-2})
        case 'cdata' % chart data
          table{ii,j} = func.data.cdata;
        case 'data'  % toolbox data
          if isa(func.data, 'coco_func_data')
            table{ii,j} = func.data.protect();
          else
            table{ii,j} = func.data;
          end
        case {'xidx' 'uidx'} % position of x0 and t0 in full solution vector
          table{ii,j} = func.x_idx(:);
        case 'fidx'  % position of F(x) in zero function
          table{ii,j} = func.f_idx(:);
        case 'midx'  % position of F(x) in monitor function
          table{ii,j} = func.m_idx(:);
        case {'x0' 'u0'}  % initial solution point of toolbox
          table{ii,j} = efunc.x0(func.x_idx);
        case 't0'    % initial tangent of toolbox
          table{ii,j} = efunc.tx(func.x_idx);
        case {'v0' 'vecs'}  % initial solution point of toolbox
          table{ii,j} = efunc.V0(func.x_idx,:);
        otherwise
          error('%s: unknown function data field ''%s''', ...
            mfilename, varargin{j-2});
      end
    end
    ii = ii+1;
  end
end

for i=1:nargout
  varargout{i} = table(:,i);
end

end
