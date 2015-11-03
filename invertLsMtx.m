function layoutInvMtx = invertLsMtx(ls_dirs, ls_groups)
%INVERTTRIANGLEMTX Summary of this function goes here
%   Detailed explanation goes here

% Dimensionality, 2d or 3d VBAP
if size(ls_groups, 2) == 2
    dim = 2;
    if isrow(ls_dirs), ls_dirs = ls_dirs'; end
elseif size(ls_groups, 2) == 3
    dim = 3;
end
ls_dirs_rad = ls_dirs*pi/180;

% Number of loudspeakers and loudspeaker groups
N_group = size(ls_groups, 1);

% Convert to cartesian coordinates
if dim == 2
    [U_spkr(:,1), U_spkr(:,2)] = pol2cart(ls_dirs_rad, 1); % size [N_spkr x 2]
elseif dim == 3
    [U_spkr(:,1), U_spkr(:,2), U_spkr(:,3)] = sph2cart(ls_dirs_rad(:,1), ...
        ls_dirs_rad(:,2), 1); % size [N_spkr x 3]
end

% pre-calculate inversions of the speaker groups and store into matrix
layoutInvMtx = zeros(N_group, dim^2);
for n = 1:N_group

    % get the unit vectors for the current group
    tempGroup = U_spkr(ls_groups(n,:), :); % size [3 x 3] for 3d/ [2 x 2] for 2d

    % get inverse of current group
    tempInv = eye(dim) / tempGroup; % size [3 x 3*N_group] for 3d
    tempInv = tempInv(:); % vectorise the inverse matrix by stacking columns
    layoutInvMtx(n, :) = tempInv; % store the vectorized inverse as a row the output
end

end
