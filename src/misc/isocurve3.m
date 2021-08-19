function h = isocurve3( X,Y,Z,V1,V2,isovalue1,isovalue2,varargin)
% h = isocurve3( X,Y,Z,V1,V2,isovalue1,isovalue2)
% Plot intersection curve between isosurfaces
% X,Y,Z,V1,isovalue1 and X,Y,Z,V2,isovalue2 are input data for isosurface
% Additional arguments are passed to plot3
% h - axes handle
% Leif Persson, Mathematics department, Umeå University, March 2016
diagnostic_mode = false;
%% Compute first isosurface and faces crossing second isosurface
h=isosurface(X,Y,Z,V1,isovalue1);
v=h.vertices; % Vertex coordinates; v(i_v,:) are the x,y,z-coordinates of vertex i_v
f2v=h.faces; % Faces of isosurface 1; f(i_f,:) are the vertex numbers of face i_f
V21 = interp3(X,Y,Z,V2,v(:,1),v(:,2),v(:,3)); % Interpolated V2 values on vertices of isosurface 1
s_f=V21(f2v)-isovalue2; % Second isosurface sign on faces of first surface
is_cf=(s_f(:,1).*s_f(:,2)<=0)|(s_f(:,2).*s_f(:,3)<=0)|(s_f(:,3).*s_f(:,1)<0); % boolean vector for crossing faces
f2v=f2v(is_cf,:); % Restrict to crossing faces
n_v=size(v,1); % Number of vertices
n_f=size(f2v,1); % Number of faces
s_f = s_f(is_cf,:);
%% Restrict to vertices of crossing faces, renumber vertices
old_i_v = unique(f2v);
v = v(old_i_v, :);
V21 = V21(old_i_v);
n_new_v = length(old_i_v);
v2new_v=zeros(size(1,n_v));
v2new_v(old_i_v) = 1:n_new_v; % Translation vector
for i_f=1:n_f,
    f2v(i_f,1)=v2new_v(f2v(i_f,1));
    f2v(i_f,2)=v2new_v(f2v(i_f,2));
    f2v(i_f,3)=v2new_v(f2v(i_f,3));
end
n_v = size(v,1);
%% Alternating edge to vertex mapping
e2v = [];
for i_f=1:n_f,
    if s_f(i_f,1)*s_f(i_f,2)<=0
        e2v = [ e2v; f2v(i_f,[1,2]) ];
    end
    if s_f(i_f,2)*s_f(i_f,3)<=0
        e2v = [ e2v; f2v(i_f,[2,3]) ];
    end
    if s_f(i_f,3)*s_f(i_f,1)<=0
        e2v = [ e2v; f2v(i_f,[3,1]) ];
    end
end
e2v = unique(sort(e2v,2),'rows');
n_e = size(e2v,1);
%% Face to alternating edge mapping
v2e = sparse([e2v(:,1);e2v(:,2)],[e2v(:,2);e2v(:,1)],[(1:n_e)'; (1:n_e)'],n_v,n_v);
f2e = zeros(n_f,2);
for i_f=1:n_f,
    k=1;
    i_e = v2e(f2v(i_f,1),f2v(i_f,2));
    if i_e > 0,
        f2e(i_f, k) = i_e; k=k+1;
    end
    i_e = v2e(f2v(i_f,2),f2v(i_f,3));
    if i_e > 0,
        f2e(i_f, k) = i_e; k=k+1;
    end
    i_e = v2e(f2v(i_f,3),f2v(i_f,1));
    if i_e > 0,
        f2e(i_f, k) = i_e; k=k+1;
    end
end
%% Alternating edge to face mapping
e2f = zeros(n_e,2);
for i_f=1:n_f,
    i_e = f2e(i_f,1);
    if i_e>0,
        if e2f(i_e,1)==0,
            e2f(i_e,1)=i_f;
        elseif e2f(i_e,2)==0
            e2f(i_e,2)=i_f;
        else
            error('Too many faces for edge');
        end
    end
    i_e=f2e(i_f,2);
    if i_e>0,
        if e2f(i_e,1)==0,
            e2f(i_e,1)=i_f;
        elseif e2f(i_e,2)==0
            e2f(i_e,2)=i_f;
        else
            error('Too many faces for edge');
        end
    end
end
%% Build alternating edge-to-edge mapping
e2e = zeros(n_e,2);
for i_e=1:n_e,
    j_f=e2f(i_e,:);
    j_f=j_f(j_f>0); % The faces neighboring edge i_e (may be one or two)
    k_e=f2e(j_f,:); % The edges of those faces
    k_e=reshape(k_e,1,[]);
    k_e=k_e((k_e>0)&(k_e~=i_e)); % The neighboring edges of i_e
    e2e(i_e,1:length(k_e))=k_e;
end
if size(e2e,2)>2,
    error('Incompatible edge-to-edge mapping');
end
%% Compute crossings on alternating edges
tmp=abs(V21(e2v)-isovalue2);
c = zeros(n_e,3);
c(:,1) = (v(e2v(:,1),1).*tmp(:,2)+v(e2v(:,2),1).*tmp(:,1))./(tmp(:,1)+tmp(:,2));
c(:,2) = (v(e2v(:,1),2).*tmp(:,2)+v(e2v(:,2),2).*tmp(:,1))./(tmp(:,1)+tmp(:,2));
c(:,3) = (v(e2v(:,1),3).*tmp(:,2)+v(e2v(:,2),3).*tmp(:,1))./(tmp(:,1)+tmp(:,2));
%% Build adjacency matrix
A = sparse(n_e,n_e);
for i_e=1:n_e,
    tmp = e2e(i_e,:);
    A(i_e, tmp(tmp>0))=1;
end
if any(diag(A)~=0) || any(any(A~=A')),
    error('Incompatible adjacency matrix');
end
%% Compute parts (connected components)
p = []; % Boolean matrix; each row represents nodes of one component
r = ones(1,n_e); % Boolean vector for remaining nodes
while any(r>0)
    x=zeros(1,n_e);
    i_e=find(r,1);
    x(i_e)=1; % Start with first remaining edge
    tmp=x*(A+speye(size(A)))'; % Add nodes
    tmp=(tmp>0);
    while any(tmp~=x)
        x=tmp;
        tmp=x*(A+speye(size(A)))';
        tmp=(tmp>0);
    end
    p=[p;x];
    r = ~any([p; p]>0);
end
p=(p>0);
n_p = size(p,1); % Number of parts (connected components)
%% Part edge numbers (unsorted)
p2e = zeros(size(p));
for i_p=1:n_p,
    tmp=1:n_e;
    tmp = tmp(p(i_p,:));
    p2e(i_p,1:length(tmp))=tmp;
end
n_ep = zeros(1,n_p);
%% Sort edges in each part
sorted_p2e = zeros(size(p,1),size(p,2)+1); % +1 needed for closed curves
for i_p=1:n_p,
    ep=p2e(i_p,:);
    ep=ep(ep>0); % Edge numbers in part
    n_ep(i_p) = length(ep);
    is_open_curve = any(e2f(ep,1)<=0 | e2f(ep,2)<=0);
    if is_open_curve
%        fprintf('Part %d is open\n', i_p);
        i_ep = find(e2f(ep,1)<=0 | e2f(ep,2)<=0);
        ep_start=ep(i_ep(1));
        ep_end  =ep(i_ep(2));
    else                                   % Closed curve
%       fprintf('Part %d is closed\n', i_p);
        n_ep(i_p)=n_ep(i_p)+1;
        ep_start=ep(1);
        ep_end  =ep(1);
    end
    sorted_ep = zeros(1,n_ep(i_p)+1); % +1 needed if closed curve
    is_remaining_ep = true(1,n_e); % Boolean vector indicating remaining edges
    sorted_ep(1) = ep_start;
    is_remaining_ep(ep_start)=false;
    for i_ep=2:n_ep(i_p)-1,
        j_f = e2f(sorted_ep(i_ep-1),:);
        j_f = j_f(j_f>0);
        k_e = f2e(j_f,:);
        k_e = reshape(k_e,1,[]);
        k_e = k_e.*double(is_remaining_ep(k_e));
        k_e = k_e(k_e>0);
        if numel(k_e)==0
            error('Sorting problem...');
        end
        sorted_ep(i_ep)=k_e(1);
        is_remaining_ep(k_e(1))=false;
    end
    i_ep = n_ep(i_p);
    sorted_ep(i_ep) = ep_end;
    sorted_p2e(i_p,1:numel(sorted_ep))=sorted_ep;
end
p2e=sorted_p2e;
%% Plot curves
for i_p=1:n_p,
    plot3(c(p2e(i_p,1:n_ep(i_p)),1),c(p2e(i_p,1:n_ep(i_p)),2),c(p2e(i_p,1:n_ep(i_p)),3),varargin{:});
    hold on;
end
h = gca; % Return handle to current axis
%% Diagnostic plot for edge-face mappings
if diagnostic_mode,
    f_list=containers.Map('KeyType', 'double', 'ValueType', 'double');
    figure;
    hold on;
    view(3);
    i_f=1; f_list(i_f)=0; % Store i_f in face list
    plot_face(i_f, 'LineWidth', 2, 'Color', 'blue','LineStyle', ':');
    i_e=f2e(i_f,:); % Alternating edges of face i_f
    plot_edge(i_e(1), 'LineWidth', 2, 'Color', 'red');
    plot_edge(i_e(2), 'LineWidth', 2, 'Color', 'red');
    j_f=unique(e2f(i_e,:));
    j_f=j_f(j_f>0);
    while any(~isKey(f_list, num2cell(j_f)))
        for j=1:length(j_f),
            if ~isKey(f_list, j_f(j)) break; end
        end
        i_f=j_f(j); f_list(i_f)=0;
        plot_face(i_f, 'LineWidth', 2, 'Color', 'blue','LineStyle', ':');
        i_e=f2e(i_f,:); % Alternating edges of face i_f
        plot_edge(i_e(1), 'LineWidth', 2, 'Color', 'red');
        plot_edge(i_e(2), 'LineWidth', 2, 'Color', 'red');
        j_f=unique(e2f(i_e,:));
        j_f=j_f(j_f>0);
    end
    for i_p=1:n_p,
        plot3(c(p(i_p,:)',1),c(p(i_p,:)',2),c(p(i_p,:)',3),'o','MarkerSize', 10, 'MarkerFaceColor', 'black');
    end
    %plot3(c(p(2,:)',1),c(p(2,:)',2),c(p(2,:)',3),'d','MarkerSize', 10, 'MarkerFaceColor', 'black');
end
%% Auxiliary functions for diagnostic plot
    function plot_face(i_f,varargin)
        ii_e=f2e(i_f,:); % Edges of face i_f
        ii_v=e2v(ii_e,:); % Vertex numbers of edges
        ii_v=unique(ii_v); % Unique vertex numbers
        v_f=v(ii_v,:);
        line(v_f([1,2],1),v_f([1,2],2),v_f([1,2],3),varargin{:});
        line(v_f([2,3],1),v_f([2,3],2),v_f([2,3],3),varargin{:});
        line(v_f([3,1],1),v_f([3,1],2),v_f([3,1],3),varargin{:});
    end
    function plot_edge(i_e,varargin)
        ii_v=e2v(i_e,:); % Vertex numbers of edge
        v_f=v(ii_v,:);   % Vertex coordinates
        line(v_f([1,2],1),v_f([1,2],2),v_f([1,2],3),varargin{:});
        text(0.5*sum(v_f([1,2],1)), 0.5*sum(v_f([1,2],2)), 0.5*sum(v_f([1,2],3)), num2str(i_e));
    end
end
