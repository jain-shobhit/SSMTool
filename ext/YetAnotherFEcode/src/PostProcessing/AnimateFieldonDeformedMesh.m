function AnimateFieldonDeformedMesh(Nodes,Elements,S,varargin)
%% function to animate the displacement snapshots in S (solution cell or matrix) 
% More than one solution signals can be simultaneously animated,
% if S is a cell, then each cell component corresponds to a different
% solution matrices and if S is a matrix then it is a matrix whose columns contain
% the displacement snapshots of a single solution

% Additional Name Value Parameters:
% index: numerical array specifying which DOFs on the node correspond to translational
% displacements in X, Y, Z directions for 3D meshes, X,Y direction for 2D
% meshes and simply a X direction for 1 D meshes
% factor: the factor with which displacements should be scaled  
% filename: for storing animation files (with path)
% framerate: numerical rate of frames to be played per second in the video

[scalefactor,index,filename,framerate] = parse_inputs(varargin{:});

%% video object
nnodes = size(Nodes,1);
myVideo = VideoWriter([filename '.avi']);
myVideo.FrameRate = framerate;

if iscell(S)
    ns = length(S);
    s = zeros(ns,1);
    for j = 1:ns
        s(j) = size(S{j},2);
    end
    nt = min(s); % collect minimal number of snapshots across solutions
    
    for k = 1:ns     % check if all solutions have same number of snapshots and prune if necessary   
        if size(S{k},2) ~= nt
            warning(['Solution cell components do not have same size: truncating at first' num2str(nt) 'snapshots'] )
            S{k} = S{k}(:,1:nt);
        end
    end
else
    nt = size(S,2);
    S = {S};
    ns = 1;
end


nDOFperNode = size(S{1},1)/nnodes;

M(size(S{1},2)) = struct('cdata',[],'colormap',[]);
color = get(groot, 'defaultAxesColorOrder');
color = [0 0 0; color];
for j = 1:nt
    for k = 1:ns
        Solution = S{k};
        meshcolor = color(k,:);
        U = reshape(Solution(:,j),nDOFperNode,[]).';
        disp = U(:,index);
        PlotFieldonDeformedMesh(Nodes,Elements,disp,'factor',scalefactor,'color', meshcolor) ;
        if j*k == 1
            %             set(gca, 'nextplot', 'replacechildren')
            axis manual
            xlim = get(gca,'xlim');
            ylim = get(gca,'ylim');
        end
        set(gca,'xlim',xlim,'ylim',ylim)
    end
    % gif movie
    frame = getframe(gcf);
    M(j) = frame;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if j == 1
        %         set(gca, 'nextplot', 'replacechildren');
        axis manual
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
    end
    cla
end
open(myVideo)
writeVideo(myVideo,M);
close(myVideo)
close(gcf)
end

function [scalefactor,index,filename,framerate] = parse_inputs(varargin)
%% parsing inputs
defaultindex = 1;
defaultFactor = 1;
defaultfilename = 'test'; % plot norm of displacement\
defaultframerate = 100;

p = inputParser;
addParameter(p,'index',defaultindex, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'factor',defaultFactor,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'filename',defaultfilename,@(x)validateattributes(x, ...
                {'char'},{'nonempty'}))
addParameter(p,'framerate',defaultframerate,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );

parse(p,varargin{:});

scalefactor = p.Results.factor;
index = p.Results.index;
filename = p.Results.filename;
framerate = p.Results.framerate;
end