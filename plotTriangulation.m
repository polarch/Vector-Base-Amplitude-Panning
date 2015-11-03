function plotTriangulation(mesh)
%PLOTTRIANGULATION Summary of this function goes here
%   Detailed explanation goes here
Nvert = size(mesh.vert,1);

patch('vertices', mesh.vert, 'faces', mesh.faces, 'facecolor','g', 'FaceAlpha',0.9);
axis equal
% number vertices
clear temps
for i = 1:Nvert
    temps(i,:) = sprintf('%2i', i);
end
text(mesh.vert(:,1), mesh.vert(:,2), mesh.vert(:,3), temps, 'FontSize', 12, 'FontWeight', 'bold');
% indicate cartesian axes
line([0;1.5], [0;0], [0;0],'color','r')
line([0;0], [0;1.5], [0;0],'color','g')
line([0;0], [0;0], [0;1.5],'color','b')
set(gca,'visible','off')
set(findall(gca, 'type', 'text'), 'visible', 'on')

end

