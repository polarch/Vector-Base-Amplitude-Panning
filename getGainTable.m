function [gtable, src_dirs] = getGainTable(ls_dirs, ang_res, spread, method)
%GTABLE Computes a 2D/3D gain table using VBAP/MDAP or VBIP
%
%   INPUTS:
%
%   ls_dirs: vector of loudspeaker directions in degrees, for 2D, or 
%       (Nspeakers x 2) matrix of loudspeaker directions in degrees for 3D
%   ang_res: the angular resolution of the table in degrees, it should be a
%       scalar for 2d VBAP, or a vector [azi_res elev_res] for 3d VBAP.
%   spread: value of spread in degrees of the panning gains for MDAP
%   method: 'vbap' for amplitude panning, or 'vbip' for its energy variant.
%       Default is 'vbap', if not defined.
%
%   OUTPUTS:
%
%   gtable: (Ndirs x Nspeaker) gain matrix. For the indexing and how to 
%       access the gains for a certain direction, check the code and the
%       examples in the included scripts.
%   src_dirs: optional output returning the direction of each entry in the
%       table
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 1/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin<3, spread = 0; method = 'vbap';
    elseif nargin<4, method = 'vbap'; end

    % convert to column vector if not already, to prepare for vbap function
    if size(ls_dirs, 2) > size(ls_dirs, 1)
        ls_dirs = ls_dirs';
    end

    % find dimensionality of array
    if (min(size(ls_dirs)) == 1) || all(ls_dirs(:,2) == 0)
        dim = 2;
    else
        dim = 3;
    end

    % compute VBAP gains for 2D case
    if dim == 2

        % convert to vector if zero azimuth is also defined
        ls_dirs = ls_dirs(:, 1);

        % azimuth resolution in rads
        if nargin < 2 || isempty(ang_res)
            ang_res = 1; % default resolution for the gain table of 1 deg.
        end
        
        az_res = ang_res(1);
        src_dirs = (-180:az_res:180)';

        % find the loudspeaker pairs
        ls_groups = findLsPairs(ls_dirs);
        
    elseif dim == 3

        % Compute directions of the evaluation grid
        if (nargin < 2) || isempty(ang_res)
            ang_res = [2 5]; % default resolution for the gain table of 2 deg. azimuth, 5 deg. elevation
        end        
        
        az_res = ang_res(1); % azimuth resolution in rads
        el_res = ang_res(2); % elevation resolution in rads
        N_azi = round(360/az_res) + 1;
        N_ele = round(180/el_res) + 1;

        azi = (-180:az_res:180)';
        ele = (-90:el_res:90)';

        src_dirs = zeros(N_azi*N_ele, 2);
        src_dirs(:,1) = repmat(azi, N_ele, 1);

        for n = 1:N_ele
            tempIdx = (n-1)*N_azi;
            src_dirs(tempIdx + (1:N_azi), 2) = ele(n);
        end

        % find the loudspeaker triangles
        ls_groups = findLsTriplets(ls_dirs);
    end
    
    % invert matrices
    layoutInvMtx = invertLsMtx(ls_dirs, ls_groups);
    switch method
        case 'vbap'
            % compute vbap gains
            gtable = vbap(src_dirs, ls_groups, layoutInvMtx, spread);
        case 'vbip'
            % compute vbap gains
            gtable = vbip(src_dirs, ls_groups, layoutInvMtx, spread);
    end
    
end

