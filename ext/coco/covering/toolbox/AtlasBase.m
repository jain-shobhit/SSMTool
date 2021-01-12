classdef AtlasBase
  % Core interface class to atlas classes.
  
  % Copyright (C) Frank Schilder, Harry Dankowicz
  % $Id: AtlasBase.m 2839 2015-03-05 17:09:01Z fschild $
  
  properties
    PrintHeadLine = true; % Switch printing headline on/off.
  end
  
  properties (Constant = true)
    IsAdmissible    = 0; % Continuation direction points inside computational domain.
    IsTangent       = 1; % Continuation direction tangent to boundary of computational domain.
    IsNotAdmissible = 2; % Continuation direction points outside computational domain.
  end
  
  methods
    function atlas = AtlasBase(opts) %#ok<INUSD>
      % Class constructor.
    end
  end

  methods (Abstract=true,Static=true) % construction method
    [opts cont atlas] = create(opts, cont, dim)
  end
  
  methods (Abstract=true) % interface methods to covering algorithm
    [opts atlas cseg correct] = init_prcond (atlas, opts, chart)
    [opts atlas cseg flush  ] = init_atlas  (atlas, opts, cseg)
    [opts atlas cseg fluah  ] = add_chart   (atlas, opts, cseg)
    [opts atlas cseg correct] = predict     (atlas, opts, cseg)
  end
  
  methods % predefined interface methods

%     function [opts atlas cseg correct] = init_prcond(atlas, opts, cseg, varargin)
%       % construct initial tangent space chart.TS
%       [opts chart accept] = CurveSegmentBase.init_tangent_space(opts, chart, varargin{:});
%       
%       % Initialize initial chart
%       chart.R       = 0;
%       chart.pt      = -1;
%       chart.pt_type = 'IP';
%       chart.ep_flag = 1;
%       chart.accept  = accept;
%     end
    
    function [opts atlas cseg] = init_chart(atlas, opts, cseg)
      % Initialize chart after correction of initial point.
      %
      % This function initializes fields required by the subsequent
      % execution of the state init_admissible.
      
      % default: do nothing
      % this function is relevant for atlas codes that can start on
      % boundaries of the computational domain
    end
    
    function [opts atlas cseg] = init_admissible(atlas, opts, cseg, S)
      % Initialize admissible continuation directions.
      %
      % Reject point not inside computational domain.
      
      cseg.curr_chart.pt_type = 'EX';
      cseg.Status = cseg.CurveSegmentCorrupted;
      if ~isempty(cseg.ptlist)
        cseg.ptlist{1} = cseg.curr_chart;
      end
      if S.ep_flag == 1
        coco_warn(opts, 1, opts.cont.LogLevel, ...
          '%s: %s\n%s: %s\n', mfilename, ...
          'initial point is not inside computational domain', ...
          'active boundary or terminal constraints were', S.msg);
      else
        coco_warn(opts, 1, opts.cont.LogLevel, ...
          '%s: %s\n%s: %s\n', mfilename, ...
          'initial point is outside computational domain', ...
          'active boundary or terminal constraints were', S.msg);
      end
    end
    
    function [opts atlas cseg] = flush(atlas, opts, cseg, varargin)
      % Write point list to screen and disk.
      
      if nargin<=3
        mode = 'unique';
      else
        if ischar(varargin{1})
          mode = varargin{1};
        else
          mode = 'individual';
        end
      end
      
      opts.atlas = atlas;
      HLFlag = atlas.PrintHeadLine | cseg.isInitialSegment;
      opts = AtlasBase.bddat_print_headline(opts, HLFlag);
      
      if cseg.Status == cseg.CorrectionFailed
        
        % correction failed, save initial point of segment
        if isempty(cseg.ptlist)
          chart = cseg.src_chart;
        else
          chart = cseg.ptlist{1};
        end
        chart.pt_type = 'MX';
        chart.ep_flag = 2;
        opts = AtlasBase.bddat_add(opts, chart);
        PT   = chart.pt;
        
      elseif any(cseg.Status == [cseg.TerminalEventDetected cseg.CurveSegmentCorrupted])
        
        % MX or problem event ocurred, save end points of current curve
        % segment to allow restart algorithms to access both points.
        if numel(cseg.ptlist) >= 2
          charts = cseg.ptlist([1 end]);
        elseif numel(cseg.ptlist) == 1
          charts = cseg.ptlist(1);
        else
          charts = { cseg.curr_chart };
        end
        opts = AtlasBase.bddat_add(opts, charts{:});
        PT   = charts{end}.pt;
        
      else
        
        % compute point flags
        switch mode
          case 'all'
            pflags = true(1,numel(cseg.ptlist));
          case 'unique'
            pflags = true(1,numel(cseg.ptlist));
            pflags(1) = HLFlag;
          case 'individual'
            pflags = varargin{1};
        end
        
        % Flush new points of point list.
        for i=find(pflags)
          opts = AtlasBase.bddat_add(opts, cseg.ptlist{i});
        end
        PT = cseg.ptlist{end}.pt;
        
        atlas.PrintHeadLine = false;
        chart  = cseg.ptlist{end};
        if chart.pt>0 && chart.ep_flag
          cseg.Status = cseg.BoundaryPoint;
        else
          cseg.Status = cseg.CurveSegmentOK;
        end
      
      end
      coco_log(opts, 1, 1, '\nPOINT %d: computation finished\n', PT);
      coco_log(opts, 1, 1, ...
        '*********************************************************************\n');
    end
    
    function [opts atlas cseg predict] = refine_step(atlas, opts, cseg)
      % Handle failure of convergence of corrector.
      %
      % Stop FSM.
      
      predict = false;
    end
    
  end
  
  methods (Static=true) % gateway functions to bifurcation data
    
    function opts = bddat_init(opts)
      % Initialize interface for writing to BD output array.
      
      opts = coco_add_signal(opts, 'bddat', 'AtlasBase');
      opts = coco_add_signal(opts, 'cont_print', 'AtlasBase');
      opts = coco_add_signal(opts, 'save_full', 'AtlasBase');
      opts = coco_add_signal(opts, 'save_reduced', 'AtlasBase');
      
      op_idx         = opts.efunc.op_idx;
      bddat.op_idx   = opts.efunc.pidx2midx(op_idx);
      bdp_idx        = opts.efunc.bdp_idx;
      bddat.bdp_idx  = opts.efunc.pidx2midx(bdp_idx);
      bddat.insert   = [];
      bddat.next_lab = 1;
      
      % if isempty(op_idx)
      %   bd = { 'PT' 'StepSize' 'TYPE' 'SLAB' 'LAB' '||U||' };
      % else
        op_names   = coco_idx2par(opts, op_idx);
        op_len     = max(cellfun(@numel, op_names), 12);
        op_fmt     = sprintf(' %%%d.4e', op_len);
        op_sfmt    = sprintf(' %%%ds', op_len);
        
        bddat.op_names  = op_names;
        bddat.bdp_names = coco_idx2par(opts, bdp_idx);
        bddat.op_fmt    = op_fmt;
        bddat.op_sfmt   = op_sfmt;
        
        bd = [ {} 'PT' 'StepSize' 'TIME' '||U||' 'SLAB' 'LAB' 'TYPE' ...
          bddat.bdp_names ];
      % end
      
      [opts fids res] = coco_emit(opts, 'bddat', 'init'); %#ok<ASGLU>
      for i=1:size(res,1)
        bd = [bd res{i}]; %#ok<AGROW>
      end
      
      opts.bddat = bddat;
      opts.bd    = bd;
    end
    
    function opts = bddat_set(opts, cmd, varargin)
      % Tune writing to BD output array.
      %
      % PROB = BDDAT_SET(PROB, CMD, VARARGIN)
      %
      % CMD == 'ins_mode', VARARGIN = { 'append' | 'prepend' }
      %   Set insert mode to append data at end of BD or to prepend data at
      %   beginning of BD.
      
      bddat = opts.bddat;
      
      switch cmd
        case 'ins_mode'
          
          switch varargin{1}
            case 'append'
              bddat.insert = @AtlasBase.bddat_append;
            case 'prepend'
              bddat.insert = @AtlasBase.bddat_prepend;
            otherwise
              error('%s: %s: unrecognised mode ''%s''', ...
                mfilename, cmd, varargin{1});
          end
          
        otherwise
          error('%s: unrecognised command ''%s''', mfilename, cmd);
      end
      
      opts.bddat = bddat;
    end
    
    function opts = bddat_print_headline(opts, flag)
      % Print headline to screen.
      if flag
        opts = print_headline(opts, 1);
      else
        opts = print_headline(opts, 2);
      end
    end
    
    function opts = bddat_add(opts, varargin)
      % Add data to storage.
      %
      % Calling this function triggers adding data to the DB output array
      % and saving labelled solution points to disk.
      [opts sol] = save_data(opts, varargin{:});
      opts       = print_data(opts, sol);
      if nargin>2
        opts     = print_data(opts, varargin{2});
      end
    end
    
    function opts = bddat_append(opts, sol)
      % Append data at end of BD.
      bd              = opts.bd;
      [opts new_data] = AtlasBase.bddat_data(opts, sol);
      bd              = [ bd ; new_data ];
      opts.bd         = bd;
    end
    
    function opts = bddat_prepend(opts, sol)
      % Insert data at beginning of BD.
      bd              = opts.bd;
      [opts new_data] = AtlasBase.bddat_data(opts, sol);
      bd              = [ bd(1,:) ; new_data ; bd(2:end,:) ];
      opts.bd         = bd;
    end
    
  end
  
  methods (Static=true,Access=private) % local functions
    
    function [opts new_data] = bddat_data(opts, sol)
      bddat = opts.bddat;
      
      if strcmp(sol.pt_type, 'ROS')
        lab = [];
      else
        lab = sol.lab;
      end
      
      if  isfield(sol, 'p')
        pars   = num2cell(sol.p(bddat.bdp_idx));
      else
        pars   = cell(1, numel(bddat.bdp_idx));
      end
      new_data = { sol.pt sol.R etime(clock, opts.cont.tm) norm(sol.x) ...
        sol.lab lab sol.pt_type pars{1:end} };
      
      [opts fids res] = coco_emit(opts, 'bddat', 'data', sol); %#ok<ASGLU>
      for i=1:size(res,1)
        new_data = [new_data res{i}]; %#ok<AGROW>
      end
    end
    
  end
  
end
