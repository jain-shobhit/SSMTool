classdef coco_func_data
  % COCO_FUNC_DATA function data with access rights management.
  
  % bug: Implementation of pr is wrong. Should be inside coco_func_data.
  % bug: Implementation of protecting fields is wrong.
  % Protect(...) should make given or all fields read-only.
  % on protect: - include all fields in list of protected vars -> read only
  %             - [new] private fields [become]/remain read-write
  %             - new shared fields become rwad-write until protected
  
  properties ( Access = public, Dependent )
    data, pr, sh
  end
  
  properties ( Access = private )
    idx = 0, is_protected = false, S = struct()
    share = true, protected_fields = {}, secret_fields = {}
  end
  
  methods (Static)
    
    function out = pointers(cmd, in)
      persistent ptr
      switch cmd
        case 'get'
          out = ptr(in);
        case 'new'
          out = numel(ptr)+1;
          if out==1
            ptr = coco_func_data_ptr();
          else
            ptr(out) = coco_func_data_ptr();
          end
        case 'set'
          for i=numel(ptr):-1:1
            delete(ptr(i));
          end
          ptr = in;
        case 'restore'
          for i=numel(ptr):-1:1
            delete(ptr(i));
          end
          for i=1:numel(in)
            p      = in(i);
            ptr(i) = coco_func_data_ptr(p.pr,p.sh);
          end
        case 'copy'
          out = ptr;
          for i=1:numel(ptr)
            p      = ptr(i);
            ptr(i) = coco_func_data_ptr(p.pr,p.sh);
          end
        case 'query'
          if nargout==0
            fprintf('%d %s object(s) allocated\n', numel(ptr), mfilename);
          else
            out.count = numel(ptr);
          end
      end
    end
    
    function p = loadobj(S)
      % do not create pointer but copy
      p = coco_func_data([], [], S);
    end
    
  end
  
  methods
    
    function obj = coco_func_data(prot, shrd, S)
      switch nargin
        case 0
          obj.idx = obj.pointers('new');
          
        case 1
          if isa(prot, 'coco_func_data')
            if prot.idx
              obj.idx = prot.idx;
              obj.is_protected = prot.is_protected;
            else
              % create from loaded object
              obj.idx = obj.pointers('new');
              ptr     = obj.pointers('get', obj.idx);
              ptr.pr  = prot.S.pr;
              ptr.sh  = prot.S.sh;
            end
          else
            shmode = true;
            if ischar(prot) && ~isempty(prot)
              sidx = find(strcmpi(prot, { 'share' , 'protect'}), 1, 'first');
              assert(~isempty(sidx), ...
                '%s: second argument must be either ''share'' ot ''protect''', mfilename);
              shmode = (sidx == 1);
              prot = [];
            end
            if isempty(prot)
              prot=struct();
            end
            assert(isstruct(prot), '%s: argument ''prot'' must be a struct', mfilename);
            obj.idx = obj.pointers('new');
            ptr     = obj.pointers('get', obj.idx);
            ptr.pr  = prot;
            obj.share = shmode;
          end
          
        case 2
          if isa(prot, 'coco_func_data')
            assert(ischar(shrd), ...
              '%s: second argument must be either ''share'' ot ''protect''', mfilename);
            sidx = find(strcmpi(shrd, { 'share' , 'protect'}), 1, 'first');
            assert(~isempty(sidx), ...
              '%s: second argument must be either ''share'' ot ''protect''', mfilename);
            if prot.idx
              obj.idx = prot.idx;
              obj.is_protected = prot.is_protected;
            else
              % create from loaded object
              obj.idx = obj.pointers('new');
              ptr     = obj.pointers('get', obj.idx);
              ptr.pr  = prot.S.pr;
              ptr.sh  = prot.S.sh;
            end
            obj.share = (sidx == 1);
            
          else
            if isempty(prot)
              prot=struct();
            end
            assert(isstruct(prot), '%s: argument ''prot'' must be a struct', mfilename);
            shmode = true;
            if ischar(shrd) && ~isempty(shrd)
              sidx = find(strcmpi(shrd, { 'share' , 'protect'}), 1, 'first');
              assert(~isempty(sidx), ...
                '%s: second argument must be either ''share'' ot ''protect''', mfilename);
              shmode = (sidx == 1);
              shrd = [];
            end
            if isempty(shrd)
              shrd=struct();
            end
            assert(isstruct(shrd), '%s: argument ''shrd'' must be a struct', mfilename);
            assert(isempty(intersect(fieldnames(prot), fieldnames(shrd))), ...
              '%s: fields of protected and shared data must be unique', mfilename);
            obj.idx = obj.pointers('new');
            ptr     = obj.pointers('get', obj.idx);
            ptr.pr  = prot;
            ptr.sh  = shrd;
            obj.share = shmode;
          end
          
        case 3
          % create copy, not pointer; used by loadobj
          assert(isempty(prot) && isempty(shrd), ...
            '%s: ''prot'' and ''shrd'' must be empty when assigning copy of data', ...
            mfilename);
          obj.S = S;
      end
    end
    
    function S = saveobj(p)
      if p.idx
        ptr  = p.pointers('get', p.idx);
        S.pr = ptr.pr;
        S.sh = ptr.sh;
      else
        S = p.S;
      end
    end
    
    function fields = fieldnames(p)
      if p.idx
        ptr    = p.pointers('get', p.idx);
        fields = union(fieldnames(ptr.pr), fieldnames(ptr.sh));
      else
        fields = union(fieldnames(p.S.pr), fieldnames(p.S.sh));
      end
    end
    
    function flag = isfield(p, name)
      if p.idx
        ptr  = p.pointers('get', p.idx);
        flag = isfield(ptr.pr, name) || isfield(ptr.sh, name);
      else
        flag = isfield(p.S.pr, name) || isfield(p.S.sh, name);
      end
    end
    
    function p = rmfield(p, fields)
      if p.idx
        ptr    = p.pointers('get', p.idx);
        f      = intersect(fieldnames(ptr.pr), fields);
        ptr.pr = rmfield(ptr.pr,f);
        f      = intersect(fieldnames(ptr.sh), fields);
        ptr.sh = rmfield(ptr.sh,f);
      else
        f      = intersect(fieldnames(p.S.pr), fields);
        p.S.pr = rmfield(p.S.pr,f);
        f      = intersect(fieldnames(p.S.sh), fields);
        p.S.sh = rmfield(p.S.sh,f);
      end
    end
    
    function p = protect(p, varargin)
      if nargin==1
        p.is_protected = true;
        % bug: remove secret fields
      else
        assert(iscellstr(varargin), ...
          '%s: illegal input for field list', mfilename);
        p.protected_fields = union(p.protected_fields, varargin);
      end
    end
    
    function p = secret(p, varargin)
      assert(iscellstr(varargin), ...
        '%s: illegal input for field list', mfilename);
      p.secret_fields = union(p.secret_fields, varargin);
    end
    
    function varargout = subsref(p, S)
      if p.idx
        ptr = p.pointers('get', p.idx);
      else
        ptr = p.S;
      end
      switch S(1).type
        case '.'
          nrefs = numel(S);
          switch S(1).subs
            
            case 'data'
              if nrefs==1
                varargout{1} = p.data;
              else
                if isfield(ptr.pr, S(2).subs)
                  [varargout{1:nargout}] = subsref(ptr.pr, S(2:end));
                elseif isfield(ptr.sh, S(2).subs)
                  [varargout{1:nargout}] = subsref(ptr.sh, S(2:end));
                else
                  error('%s: reference to non-existing field ''%s''', mfilename, S(1).subs);
                end
              end
              
            case 'pr'
              if nrefs==1
                varargout{1} = ptr.pr;
              else
                [varargout{1:nargout}] = subsref(ptr.pr, S(2:end));
              end
              
            case 'sh'
              if nrefs==1
                varargout{1} = ptr.sh;
              else
                [varargout{1:nargout}] = subsref(ptr.sh, S(2:end));
              end
              
            otherwise
              if isfield(ptr.pr, S(1).subs)
                [varargout{1:nargout}] = subsref(ptr.pr, S);
              elseif isfield(ptr.sh, S(1).subs)
                [varargout{1:nargout}] = subsref(ptr.sh, S);
              elseif ismethod(p, S(1).subs)
                if nrefs==1
                  [varargout{1:nargout}] = p.(S(1).subs);
                elseif nrefs==2
                  if strcmp(S(2).type, '()')
                    [varargout{1:nargout}] = p.(S(1).subs)(S(2).subs{:});
                  else
                    [varargout{1:nargout}] = subsref(p.(S(1).subs), S(2));
                  end
                else
                  if strcmp(S(2).type, '()')
                    d = p.(S(1).subs)(S(2).subs{:});
                    [varargout{1:nargout}] = subsref(d, S(3:end));
                  else
                    [varargout{1:nargout}] = subsref(p.(S(1).subs), S(2:end));
                  end
                end
              else
                error('%s: reference to non-existing field ''%s''', mfilename, S(1).subs);
              end
          end
          
        case '()'
          [varargout{1:nargout}] = builtin('subsref', p, S);
        case '{}'
          [varargout{1:nargout}] = builtin('subsref', p, S);
          
        otherwise
          error('%s: access type ''%s'' not supported', mfilename, S(1).type);
      end
    end
    
    function p = subsasgn(p, S, b)
      if ~isa(p, 'coco_func_data')
        p = builtin('subsasgn', p, S, b);
        return
      end
      
      if p.idx
        ptr = p.pointers('get', p.idx);
      else
        ptr = p.S;
      end
      switch S(1).type
        case '.'
          nsubs = numel(S);
          switch S(1).subs
            
            case 'data'
              if nsubs==1
                p.data = b;
              else
                assert(S(2).type=='.', ...
                  '%s access type ''%s'' not supported on data', S(2).type);
                if isfield(ptr.pr, S(2).subs)
                  S(1).subs = 'pr';
                elseif p.share
                  S(1).subs = 'sh';
                else
                  S(1).subs = 'pr';
                end
                p = subsasgn(p, S, b);
              end
              
            case 'pr'
              if nsubs==1
                p.pr = b;
              else
                assert(S(2).type=='.', ...
                  '%s access type ''%s'' not supported on prot', S(2).type);
                assert(~p.is_protected, ...
                  '%s: attempt to assign to protected field(s)', mfilename);
                assert(~any(strcmp(S(2).subs, fieldnames(ptr.sh))), ...
                  '%s: attempt to assign a shared field as protected', mfilename);
                ptr.pr = subsasgn(ptr.pr, S(2:end), b);
              end
              
            case 'sh'
              if nsubs==1
                p.sh = b;
              else
                assert(S(2).type=='.', ...
                  '%s access type ''%s'' not supported on shrd', S(2).type);
                assert(~p.is_protected || ...
                  any(strcmp(S(2).subs, fieldnames(ptr.sh))), ...
                  'attempt to add field(s) to protected data', mfilename);
                assert(~any(strcmp(S(2).subs, fieldnames(ptr.pr))), ...
                  '%s: attempt to assign a protected field as shared', mfilename);
                ptr.sh = subsasgn(ptr.sh, S(2:end), b);
              end
              
            otherwise
              % The code below executes faster than the equivalent but
              % commented code that follows.
              method = any(strcmp(S(1).subs, methods(p)));
              assert(~method, '%s: field ''%s'' is a method', ...
                mfilename, S(1).subs);
%               assert(~ismethod(p, S(1).subs), '%s: field ''%s'' is a method', ...
%                 mfilename, S(1).subs);
              if isfield(ptr.sh, S(1).subs)
                ptr.sh = subsasgn(ptr.sh, S, b);
              else
                assert(~p.is_protected, ...
                  '%s: attempt to assign to protected field(s)', mfilename);
                ptr.pr = subsasgn(ptr.pr, S, b);
              end
          end
          
        otherwise
          error('%s: access type ''%s'' not supported', mfilename, S(1).type);
      end
      if ~p.idx
        p.S = ptr;
      end
    end
    
    function disp(p)
      if isscalar(p)
        str = sprintf('id = %d', p.idx);
        if p.is_protected
          str = sprintf('%s [protected]', str);
        else
          if p.share
            str = sprintf('%s [new fields shared]', str);
          else
            str = sprintf('%s [new fields private]', str);
          end
        end
        fprintf('%s\n\n', str);
        if p.idx
          ptr = p.pointers('get', p.idx);
        else
          ptr = p.S;
        end
        empty = true;
        if ~isempty(fieldnames(ptr.pr))
          empty = false;
          disp('pr:')
          disp(ptr.pr);
        end
        if ~isempty(fieldnames(ptr.sh))
          empty = false;
          disp('sh:')
          disp(ptr.sh);
        end
        if empty
          str = sprintf(sprintf('1x1 %s array with no fields\n', mfilename));
          disp(str);
        end
      else
        s    = size(p);
        str  = [ sprintf('%d', s(1)) sprintf('x%d', s(2:end))];
        str  = sprintf(sprintf('%s %s array\n', str, mfilename));
        disp(str);
      end
    end
    
  end
  
  methods
    
    function d = get.data(p)
      if p.idx
        ptr = p.pointers('get', p.idx);
      else
        ptr = p.S;
      end
      d      = ptr.pr;
      fields = fieldnames(ptr.sh);
      for i=1:numel(fields)
        field = fields{i};
        d.(field) = ptr.sh.(field);
      end
    end
    
    function d = get.pr(p)
      if p.idx
        ptr = p.pointers('get', p.idx);
      else
        ptr = p.S;
      end
      d   = ptr.pr;
    end
    
    function d = get.sh(p)
      if p.idx
        ptr = p.pointers('get', p.idx);
      else
        ptr = p.S;
      end
      d   = ptr.sh;
    end
    
    function p = set.data(p,d)
      assert(isstruct(d), '%: data must be a struct', mfilename);
      
      if p.idx
        ptr = p.pointers('get', p.idx);
      else
        ptr = p.S;
      end
      
      all_fields = fieldnames(d);
      
      sdata   = struct();
      sfields = intersect(all_fields, fieldnames(ptr.sh));
      for i=1:numel(sfields)
        field = sfields{i};
        sdata.(field) = d.(field);
      end
      ptr.sh = sdata;
      
      pdata   = struct();
      pfields = setdiff(all_fields, sfields);
      assert(~p.is_protected || isempty(pfields), ...
        '%s: attempt to assign to protected field(s)', mfilename);
      for i=1:numel(pfields)
        field = pfields{i};
        pdata.(field) = d.(field);
      end
      ptr.pr = pdata;
      
      if ~p.idx
        p.S = ptr;
      end
    end
    
    function p = set.pr(p,d)
      if p.idx
        ptr = p.pointers('get', p.idx);
      else
        ptr = p.S;
      end
      
      assert(isstruct(d), '%: prot must be a struct', mfilename);
      assert(~p.is_protected, '%s: attempt to assign to protected field(s)', mfilename);
      assert(isempty(intersect(fieldnames(d), fieldnames(ptr.sh))), ...
        '%s: attempt to assign a shared field as protected', mfilename);
      ptr.pr = d;
      
      if ~p.idx
        p.S = ptr;
      end
    end
    
    function p = set.sh(p,d)
      if p.idx
        ptr = p.pointers('get', p.idx);
      else
        ptr = p.S;
      end
      
      assert(isstruct(d), '%: shrd must be a struct', mfilename);
      assert(~p.is_protected || ...
        isempty(setdiff(fieldnames(d), fieldnames(ptr.sh))), ...
        'attempt to add field(s) to protected data', mfilename);
      assert(isempty(intersect(fieldnames(d), fieldnames(ptr.pr))), ...
        '%s: attempt to assign a protected field as shared', mfilename);
      ptr.sh = d;
      
      if ~p.idx
        p.S = ptr;
      end
    end
    
  end
  
end
