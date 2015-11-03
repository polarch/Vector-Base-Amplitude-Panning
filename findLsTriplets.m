function [ls_groups, mesh] = findLsTriplets(ls_dirs, OMIT_LARGE_TRI, aperture_lim)
%TRIANGULATESETUP Summary of this function goes here
%   Detailed explanation goes here

if nargin<2 || isempty(OMIT_LARGE_TRI)
    OMIT_LARGE_TRI = 0;
end

ls_dirs_rad = ls_dirs*pi/180;

% triangulate with delaunay triangulation on the unit sphere (convex hull)
[mesh.vert(:,1), mesh.vert(:,2), mesh.vert(:,3)] = sph2cart(ls_dirs_rad(:,1), ls_dirs_rad(:,2), 1);
mesh.faces = sphDelaunayTriangulation(ls_dirs_rad);

% discard invalid faces, and omit large triangles if asked
mesh = keepValidTriangles(mesh);

% omit large triangles, if asked
if OMIT_LARGE_TRI
    mesh = omitLargeTriangles(mesh, aperture_lim);
end

ls_groups = mesh.faces;

end

function mesh = keepValidTriangles(mesh)
%KEEPVALIDFACES Summary of this function goes here
%   Detailed explanation goes here

valid_faces = [];

for nf=1:size(mesh.faces,1)
    temp = mesh.faces(nf,:);
    vec = mesh.vert(temp,:);
    cvec = cross(vec(2,:)-vec(1,:), vec(3,:)-vec(2,:));
    centroid = mean(vec);
    if acos(dot(centroid, cvec))<pi/2, valid_faces = [valid_faces, nf]; end
end
mesh.faces =  mesh.faces(valid_faces,:);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = omitLargeTriangles(mesh, aperture_lim)
%OMITLARGETRIANGLES Summary of this function goes here
%   Detailed explanation goes here

    valid_faces = [];
    for nf=1:size(mesh.faces,1)
        temp = mesh.faces(nf,:);
        vec = mesh.vert(temp,:);
        a = acos(dot(vec(1,:), vec(2,:)));
        b = acos(dot(vec(2,:), vec(3,:)));
        c = acos(dot(vec(3,:), vec(1,:)));
        abc = [a b c];
        if all(abc<=aperture_lim*pi/180), valid_faces = [valid_faces; nf]; end
    end
    mesh.faces =  mesh.faces(valid_faces,:);    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function faces = sphDelaunayTriangulation(dirs_rad)
%DELAUNAYTRIANGULATION Computes the Delaunay triangulation on the unit sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHDELAUNAYTRIANGULATION.M - 15/7/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to cartesian
N_vert = size(dirs_rad, 1);
[tempx, tempy, tempz] = sph2cart(dirs_rad(:,1), dirs_rad(:,2), ones(N_vert,1));
U_vert = [tempx, tempy, tempz];

% Find the convex hull of the points on the sphere - in this special case
% the result equals the Delaunay triangulation of the points
faces = convhulln(U_vert);

% Invert the triangles
faces = faces(:, 3:-1:1);

% Shift the results to begin each triangle from the smallest entry
for n = 1:size(faces,1)
    tempface = faces(n,:);
    [~, minIdx] = min(tempface);
    faces(n, :) = circshift(tempface, [0 1-minIdx]);
end

% Sort through triangles with smaller entries first
faces = sortrows(faces, 1); % sort through first entry
maxentry = max(faces(:,1)); % sort through second entry
n = 1;
while n <= maxentry
    startIdx = find(faces(:,1) == n, 1, 'first');
    if ~isempty(startIdx)
        endIdx = find(faces(:,1) == n, 1, 'last');
        faces(startIdx:endIdx, :) = sortrows(faces(startIdx:endIdx, :), 2);
        n = n + 1;
    else
        n = n + 1;
    end
end

end