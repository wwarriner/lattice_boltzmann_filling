% cube sides are all interval [ -0.5, 0.5 ]
% center is [ 0, 0, 0 ]

% normal must be unit normal
% offset is distance from center to plane, negative means plane is on negative
% side of origin, positive means plane on positive side
function volume = cube_segment_volume( normals, offsets )

% setup output
volume = nan( size( offsets ) );

% exploit antisymmetry of offset about origin
flip = 0 < offsets;
offsets( flip ) = -offsets( flip );

% exploit remaining cube symmetries
normals = sort( abs( normals ), 2 );

% case where plane is outside cube
origin_dist = normals * [ -0.5 -0.5 -0.5 ].' - offsets;
plane_outside = 0 <= origin_dist;
consider = plane_outside;
volume( consider ) = 0;
done = consider;

% case where plane does not intercept bottom xy-face inside cube
n = inf( size( offsets ) );
n( ~done ) = normals( ~done, : ) * [ 0.5 0.5 -0.5 ].' - offsets( ~done );
four_verts_inside = n <= 0;
consider = ~done & four_verts_inside;
zi = [ -0.5 -0.5; -0.5 0.5; 0.5 -0.5; 0.5 0.5 ];
z = compute_intercept( normals( consider, [ 1 2 ] ), zi, normals( consider, 3 ), offsets( consider ) );
z = z + 0.5;
volume( consider ) = mean( z, 2 );
done = consider | done;

% case where unit normal has x == 0
flat = normals( :, 1 ) == 0;
consider = ~done & flat;
zi = [ -0.5 -0.5 ];
z = compute_intercept( normals( consider, [ 1 2 ] ), zi, normals( consider, 3 ), offsets( consider ) );
z = max( min( z + 0.5, 1 ), 0 );
yi = [ -0.5 -0.5 ];
y = compute_intercept( normals( consider, [ 1 3 ] ), yi, normals( consider, 2 ), offsets( consider ) );
y = max( min( y + 0.5, 1 ), 0 );
volume( consider ) = y .* z / 2;
done = consider | done;

% remaining case
consider = ~done;
zi = [ -0.5 -0.5; -0.5 0.5; 0.5 -0.5 ];
z = compute_intercept( normals( consider, [ 1 2 ] ), zi, normals( consider, 3 ), offsets( consider ) );
z = z + 0.5;
zc = max( min( z, 1 ), 0 );
yi = [ -0.5 -0.5; -0.5 0.5; 0.5 -0.5 ];
y = compute_intercept( normals( consider, [ 1 3 ] ), yi, normals( consider, 2 ), offsets( consider ) );
y = y + 0.5;
yc = max( min( y, 1 ), 0 );
xi = [ -0.5 -0.5; -0.5 0.5; 0.5 -0.5 ];
x = compute_intercept( normals( consider, [ 2 3 ] ), xi, normals( consider, 1 ), offsets( consider ) );
x = x + 0.5;
xc = max( min( x, 1 ), 0 );
V100 = yc( :, 3 ) .* zc( :, 3 ) .* ( x( :, 1 ) - xc( :, 1 ) );
V010 = xc( :, 3 ) .* zc( :, 2 ) .* ( y( :, 1 ) - yc( :, 1 ) );
V001 = xc( :, 2 ) .* yc( :, 2 ) .* ( z( :, 1 ) - zc( :, 1 ) );
V000 = x( :, 1 ) .* y( :, 1 ) .* z( :, 1 );
volume( consider ) = ( V000 - ( V100 + V010 + V001 ) ) / 6;
done = consider | done;
cases_done = sum( done );

% exploit antisymmetry of offset about origin
volume( flip ) = 1 - volume( flip );

% sanity check
assert( ~any( isnan( volume ) ) );
assert( sum( cases_done ) == numel( offsets ) );
assert( all( 0.0 <= volume ) );
assert( all( volume <= 1.0 ) );

end


function p = compute_intercept( normals, coords, coeff, offset )

if isempty( normals )
    p = [];
    return
end

p = -( normals * coords.' - offset ) ./ coeff;

end

