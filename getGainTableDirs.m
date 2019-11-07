function panning_dirs = getGainTableDirs(ang_res)
%GTABLE Returns the directions of the VBAP gain table for a certain resolution
%
%   INPUTS:
%
%   ang_res: the angular resolution of the table in degrees, it should be a
%       scalar for 2d VBAP, or a vector [azi_res elev_res] for 3d VBAP.
%
%   OUTPUTS:
%
%   panning_dirs: output returning the direction of each entry in the table
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 1/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % dimensionality
    if length(ang_res)==1
        dim = 2;
    elseif length(ang_res)==2
        dim = 3;
    else
        error('wrong size of angular resolution parameters (azi or [azi elev]')
    end

    % compute VBAP gains for 2D case
    if dim == 2

        az_res = ang_res;
        panning_dirs = (-180:az_res:180)';
        
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

        panning_dirs = zeros(N_azi*N_ele, 2);
        panning_dirs(:,1) = repmat(azi, N_ele, 1);

        for n = 1:N_ele
            tempIdx = (n-1)*N_azi;
            panning_dirs(tempIdx + (1:N_azi), 2) = ele(n);
        end
    end
end