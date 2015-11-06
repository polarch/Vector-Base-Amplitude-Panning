function plotTriangulation(mesh, h_ax)
%PLOTTRIANGULATION Plots a triangulated mesh
%
%   INPUTS:
%
%   mesh: mesh structure containing 'vertices' and 'faces' fields in
%       standard Matlab syntax
%   h_ax: axes handle to plot on (e.g. for subplots), gca is the default
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 1/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    h_ax = gca;
end

Nvert = size(mesh.vert,1);

patch('vertices', mesh.vert, 'faces', mesh.faces, 'facecolor','g', 'FaceAlpha',0.9);
axis equal
% number vertices
for i = 1:Nvert
    temps(i,:) = sprintf('%2i', i);
end
text(mesh.vert(:,1), mesh.vert(:,2), mesh.vert(:,3), temps, 'FontSize', 12, 'FontWeight', 'bold');
% indicate cartesian axes
line([0;1.5], [0;0], [0;0],'color','r')
line([0;0], [0;1.5], [0;0],'color','g')
line([0;0], [0;0], [0;1.5],'color','b')
set(gca,'visible','off')
set(findall(h_ax, 'type', 'text'), 'visible', 'on')

end

