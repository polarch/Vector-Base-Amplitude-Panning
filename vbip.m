function GainMtx = vbip(src_dirs, ls_groups, ls_invMtx, spread, num_spread_src, num_spread_rings3d)
%VBIP Computes vector-base intensity panning gains for a set of directions
%
%   INPUTS:
%
%   src_dirs: panning direction in degrees, vector for 2D, [Nsrc x 2] matrix 
%       for 3D, in [azi elev] convention
%   ls_groups: valid pairs (for 2D) or triplets (for 3D triplets) returned
%       by findLsPairs() or findLsTriplets()
%   ls_InvMtx: matrix of loudspeaker inversions returned by invertLsMtx()
%   spread: value of spread in degrees of the panning gains for MDAP
%       Additional spreading parameters like number of spread sources and
%       rings can be enabled easily, see this code and getSpreadSrcDirs()
%   num_spread_src: see getSpreadSrcDirs()
%   num_spread_rings3d: see getSpreadSrcDirs()
%
%   OUTPUTS:
%
%   GainMtx: (Nsrc x Nspeaker) VBIP gain matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 1/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    spread = 0;
    num_spread_src = 0;
    num_spread_rings3d = 0;
elseif nargin<5
    num_spread_src = 8;
    num_spread_rings3d = 1;
elseif nargin<6
    num_spread_rings3d = 1;
end

src_num = size(src_dirs,1);
dim = size(ls_groups,2);
ls_num = max(ls_groups(:));

GainMtx = zeros(src_num, ls_num);
% 3D case
if dim == 3
    for ns=1:src_num
        gains = max(0, vbip3(src_dirs(ns,1), src_dirs(ns,2), ls_groups, ls_invMtx, spread, num_spread_src, num_spread_rings3d));
        GainMtx(ns,:) = gains;
    end
    % 2D case
elseif dim == 2
    for ns=1:src_num
        gains = max(0, vbip2(src_dirs(ns), ls_groups, ls_invMtx, spread, num_spread_src));
        GainMtx(ns,:) = gains;
    end
end

end

%%%%%% COMPUTE GAIN FACTORS CORRESPONDING TO INPUT ANGLE, 2D
function gains = vbip2(azi, ls_groups, ls_invMtx, spread, num_spread_src)

if nargin<4, spread = 0; end

ls_num = max(ls_groups(:));
energies = zeros(1,ls_num);

if spread
    
    U_spread = getSpreadSrcDirs(azi, spread, num_spread_src);
    for ns = 1:size(U_spread, 1);
        u_ns = U_spread(ns,:);
        e_ns = zeros(1,ls_num);
        for i=1:size(ls_groups,1);
            e_tmp(1) = ls_invMtx(i,1:2) * u_ns';
            e_tmp(2) = ls_invMtx(i,3:4) * u_ns';
            if min(e_tmp) > -0.001
                e_ns(ls_groups(i,:)) = e_tmp/sum(e_tmp);
            end
        end
        energies = energies + e_ns;
    end 
else
    
    azi_rad = azi*pi/180;
    u = [cos(azi_rad) sin(azi_rad)];
    energies = zeros(1,ls_num);
    
    for i=1:size(ls_groups,1);
        e_tmp(1) = ls_invMtx(i,1:2) * u';
        e_tmp(2) = ls_invMtx(i,3:4) * u';
        if min(e_tmp) > -0.001
            energies(ls_groups(i,:)) = e_tmp/sum(e_tmp);
        end
    end
end
energies=energies./sum(energies);
gains = sqrt(energies);

end



function gains = vbip3(azi, elev, ls_groups, ls_invMtx, spread, num_spread_src, num_spread_rings3d)

if nargin<5, spread = 0; end

ls_num = max(ls_groups(:));
energies = zeros(1,ls_num);

if spread
    
    U_spread = getSpreadSrcDirs([azi elev], spread, num_spread_src, num_spread_rings3d);
    for ns = 1:size(U_spread, 1);
        u_ns = U_spread(ns,:);
        e_ns = zeros(1,ls_num);
        for i=1:size(ls_groups,1);
            e_tmp(1) = ls_invMtx(i,1:3) * u_ns';
            e_tmp(2) = ls_invMtx(i,4:6) * u_ns';
            e_tmp(3) = ls_invMtx(i,7:9) * u_ns';
            if min(e_tmp) > -0.001
                e_ns(ls_groups(i,:)) = e_tmp/sum(e_tmp);
                energies = energies + e_ns;
                break
            end
        end
    end
else   
    azi_rad = azi*pi/180;
    elev_rad = elev*pi/180;
    u = [cos(azi_rad)*cos(elev_rad) sin(azi_rad)*cos(elev_rad) sin(elev_rad)];
    
    energies = zeros(1,ls_num);
    for i=1:size(ls_groups,1);
        e_tmp(1) = ls_invMtx(i,1:3) * u';
        e_tmp(2) = ls_invMtx(i,4:6) * u';
        e_tmp(3) = ls_invMtx(i,7:9) * u';
        if min(e_tmp) > -0.001
            energies(ls_groups(i,:)) = e_tmp/sum(e_tmp);
            break
        end
    end
end
energies = energies/sum(energies);
energies(energies<0) = 0;
gains = sqrt(energies);

end
