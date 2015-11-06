function U_spread = getSpreadSrcDirs(src_dir, spread, num_src, num_rings_3d)
%GENERATESPREADSRCDIRS Generate virtual sources for VBAP spreading and MDAP
%
%   INPUTS:
%
%   src_dir: panning direction in degrees, [azi] for 2D, [azi elev] for 3D
%   spread: spread angle in degrees defining the extent of the panned source
%   num_src: number of auxiliary sources to use for spreading, default is 8
%   num_rings_3d: number of concentric rings of num_src each to generate
%       inside the spreading surface, default is 1 (valid only for 3D
%       spreading)
%
%   OUTPUTS:
%
%   U_spread: [num_src+1 x 3] matrix of unit vectors pointing to the spread
%       directions, with the last one being the actual panning direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis & Ville Pulkki 1/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% layout dimension
dim = length(src_dir)+1;

% default number of spread sources and concentric rings, if 3D
if nargin<3
    num_src = 8;
    if dim==3
        num_rings_3d = 1;
    end
elseif nargin<4 && dim == 3
    num_rings_3d = 1;
end

% 2D case
if dim == 2
    spread_dirs = (src_dir-spread/2 : spread/num_src : src_dir+spread/2).';
    spread_dirs_rad = spread_dirs*pi/180;
    U_spread = [cos(spread_dirs_rad) sin(spread_dirs_rad)];
    
% 3D case
elseif dim == 3
    % rotation matrix using the axis of rotation-angle definition (around
    % source direction)
    src_dir_rad = src_dir*pi/180;
    [u(1), u(2), u(3)] = sph2cart(src_dir_rad(1), src_dir_rad(2), 1);
    u_x_u = [u(1)^2 u(1)*u(2) u(1)*u(3); u(1)*u(2) u(2)^2 u(2)*u(3); u(1)*u(3) u(2)*u(3) u(3)^2];
    u_x = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
    theta = 2*pi/num_src;
    R_theta = eye(3)*cos(theta) + sin(theta)*u_x + (1-cos(theta))*u_x_u;
    % create a ring of sources on the plane that is purpendicular to the source
    % directions
    spreadbase = zeros(num_src,3);
    elev = src_dir(2);
    if (elev>89) || (elev<-89)
        spreadbase(1,:) = [1 0 0];
    else
        spreadbase(1,:) = cross_prod_scaled(u,[0 0 1]);
    end
    % get ring of directions by rotating the first vector around the source
    for ns = 2:num_src
        spreadbase(ns,:) = R_theta*spreadbase(ns-1,:)';
    end
    % squeeze the purpendicular ring to the desired spread
    spread_rad = (spread/2)*pi/180;
    ring_rad = spread_rad/num_rings_3d;
    U_spread = zeros(num_rings_3d*num_src, 3);
    for nr=1:num_rings_3d
        U_spread((nr-1)*num_src + (1:num_src), :) = repmat(u,[num_src 1]) + spreadbase*tan(ring_rad*nr);
    end
    % normalize vectors to unity
    U_spread_norm = sqrt(sum(U_spread.^2,2));
    U_spread = U_spread./repmat(U_spread_norm, [1 3]);
    
    % append the original source direction at the end
    U_spread(end+1,:) = u;
    
end

end

function x = cross_prod_scaled(a,b)

x = cross(a,b);
x = x/sqrt(sum(x.^2));

end