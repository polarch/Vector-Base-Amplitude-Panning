function ls_groups = findLsPairs(ls_dirs, OMIT_LARGE_PAIR, aperture_lim)
%FINDLSPAIRS Returns the adjacent loudspeaker pairs for 2D layouts
%
%   INPUTS:
%
%   ls_dirs: vector of loudspeaker directions in degrees
%   OMIT_LARGE_PAIR: {0,1} False to return all pairs, true to keep only
%       ones that have spans less than a limit set in aperture_lim. If not set,
%       all pairs are returned.
%   aperture_lim: limit in degrees above which the respective pairs are
%       omitted
%
%   OUTPUTS:
%
%   ls_groups: (Npairs x 2) matrix. Each row has the speaker indices
%       corresponding to an adjacent pair.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 1/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 || isempty(OMIT_LARGE_PAIR)
    OMIT_LARGE_PAIR = 0;
end
if isrow(ls_dirs), ls_dirs = ls_dirs.'; end

% find the loudspeaker pairs by sorting the angles
[~, idx_sorted] = sort(ls_dirs);
idx_sorted(end + 1) = idx_sorted(1);

ls_groups = zeros(length(idx_sorted)-1, 2);
for n = 1:length(idx_sorted)-1
    ls_groups(n, :) = [idx_sorted(n) idx_sorted(n+1)];
end

% omit pairs with very large angle spans, if asked
if OMIT_LARGE_PAIR
    ls_dirs_sorted = ls_dirs(idx_sorted);
    apertures = ls_dirs_sorted(2:end)-ls_dirs_sorted(1:end-1);
    apertures(end) = 360-abs(apertures(end));
    ls_groups = ls_groups(apertures<aperture_lim,:);
end

end
