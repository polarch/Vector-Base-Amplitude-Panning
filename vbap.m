function GainMtx = vbap(src_dirs, ls_groups, ls_invMtx, spread)
%VBAP Computes vector-base amplitude panning gains for a set of directions
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
%       rings can be enabled easily, see vbip() and getSpreadSrcDirs()
%
%   OUTPUTS:
%
%   GainMtx: (Nsrc x Nspeaker) VBAP/MDAP gain matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Ville Pulkki & Archontis Politis, 1/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    spread = 0;
end

src_num = size(src_dirs,1);
dim = size(ls_groups,2);
ls_num = max(ls_groups(:));

GainMtx = zeros(src_num, ls_num);
% 3D case
if dim == 3
    for ns=1:src_num
        gains = max(0, vbap3(src_dirs(ns,1), src_dirs(ns,2), ls_groups, ls_invMtx, spread));
        GainMtx(ns,:) = gains;
    end
    % 2D case
elseif dim == 2
    for ns=1:src_num
        gains = max(0, vbap2(src_dirs(ns), ls_groups, ls_invMtx, spread));
        GainMtx(ns,:) = gains;
    end
end

end

%%%%%% COMPUTE GAIN FACTORS CORRESPONDING TO INPUT ANGLE, 2D
function gains = vbap2(azi, ls_groups, ls_invMtx, spread)

if nargin<4, spread = 0; end

ls_num = max(ls_groups(:));
gains = zeros(1,ls_num);

if spread
    
    U_spread = getSpreadSrcDirs(azi, spread);
    for ns = 1:size(U_spread, 1);
        u_ns = U_spread(ns,:);
        g_ns = zeros(1,ls_num);
        for i=1:size(ls_groups,1);
            g_tmp(1) = ls_invMtx(i,1:2) * u_ns';
            g_tmp(2) = ls_invMtx(i,3:4) * u_ns';
            if min(g_tmp) > -0.001
                g_ns(ls_groups(i,:)) = g_tmp/sqrt(sum(g_tmp.^2));
            end
        end     
        gains = gains + g_ns;
    end 
else
    
    azi_rad = azi*pi/180;
    u = [cos(azi_rad) sin(azi_rad)];
    gains = zeros(1,ls_num);
    
    for i=1:size(ls_groups,1);
        g_tmp(1) = ls_invMtx(i,1:2) * u';
        g_tmp(2) = ls_invMtx(i,3:4) * u';
        if min(g_tmp) > -0.001
            gains(ls_groups(i,:)) = g_tmp/sqrt(sum(g_tmp.^2));
        end
    end
end
gains=gains./norm(gains);

end



function gains = vbap3(azi, elev, ls_groups, ls_invMtx, spread)

ls_num = max(ls_groups(:));
gains = zeros(1,ls_num);

if spread
    
    U_spread = getSpreadSrcDirs([azi elev], spread);
    for ns = 1:size(U_spread, 1);
        u_ns = U_spread(ns,:);
        g_ns = zeros(1,ls_num);
        for i=1:size(ls_groups,1);
            g_tmp(1) = ls_invMtx(i,1:3) * u_ns';
            g_tmp(2) = ls_invMtx(i,4:6) * u_ns';
            g_tmp(3) = ls_invMtx(i,7:9) * u_ns';
            if min(g_tmp) > -0.001
                g_ns(ls_groups(i,:)) = g_tmp/sqrt(sum(g_tmp.^2));
                gains = gains + g_ns;
                break
            end
        end
    end
else   
    azi_rad = azi*pi/180;
    elev_rad = elev*pi/180;
    u = [cos(azi_rad)*cos(elev_rad) sin(azi_rad)*cos(elev_rad) sin(elev_rad)];
    
    gains = zeros(1,ls_num);
    for i=1:size(ls_groups,1);
        g_tmp(1) = ls_invMtx(i,1:3) * u';
        g_tmp(2) = ls_invMtx(i,4:6) * u';
        g_tmp(3) = ls_invMtx(i,7:9) * u';
        if min(g_tmp) > -0.001
            gains(ls_groups(i,:)) = g_tmp/sqrt(sum(g_tmp.^2));
            break
        end
    end
end
gains = gains/sqrt(sum(gains.^2));

end
