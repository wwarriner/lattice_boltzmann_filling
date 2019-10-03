% cube sides are all interval [ -0.5, 0.5 ]
% center is [ 0, 0, 0 ]

function volume = cube_segment_volume( unit_normal, offset_from_center )

unit_normal = sort( abs( unit_normal ) ); % take advantage of cube symmetries

n = dot( unit_normal, [ 0.5 0.5 -0.5 ] ) - offset_from_center;

end

