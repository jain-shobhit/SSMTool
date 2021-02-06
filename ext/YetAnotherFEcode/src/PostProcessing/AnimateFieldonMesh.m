function AnimateFieldonMesh(Nodes,Elements,F,filename,varargin)
%% Solution:

myVideo = VideoWriter(['Examples/Results/' filename '.avi']);

if nargin == 7
    myVideo.FrameRate = varargin{1};
else
    myVideo.FrameRate = 100;
end

nt = size(F,2);
M(nt) = struct('cdata',[],'colormap',[]);

for j = 1:nt
        
    PlotFieldonMesh(Nodes,Elements,F(:,j)) ;
    drawnow
    pause(0.1)
    % gif movie
    frame = getframe(gcf);
    M(j) = frame;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if j == 1
        %         set(gca, 'nextplot', 'replacechildren');
        axis manual
        imwrite(imind,cm,['Examples/Results/' filename '.gif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,['Examples/Results/' filename '.gif'],'gif','WriteMode','append');
    end
    cla
end
open(myVideo)
writeVideo(myVideo,M);
close(myVideo)
close(gcf)