function indices = sort_points_counterclockwise( points, mean_normals )

% missing points must be at the end of that PAGE in the VEC dimension

% DIMENSION MEANING
VEC = 1;
PAGE = 2;
DIM = 3;

% calculate
a = points( 2 : end, :, : );
c = points( 1, :, : );
p = repmat( mean_normals, [ size( points, 1 ) - 1, 1, 1 ] );
d = sum( cross( a, p, DIM ) .^ 2, DIM );
%d = [ zeros( 1, size( points, PAGE ) ); d ];
[ ~, m ] = max( d, [], VEC );
m = m + ( 0 : size( points, VEC ) : numel( points( :, :, 1 ) ) - 1 );
m = m + ( 0 : numel( points( :, :, 1 ) ) : numel( points ) - 1 ).';
m = m + 1;
p = permute( points( m ), [ 3 2 1 ] ); % indices are DIM PAGE VEC order
p = p - points( 1, :, : ); % by VEC
q = cross( mean_normals, p, DIM );

tr = cross( ( a - c ), repmat( p, [ size( a, 1 ) 1 1 ] ), DIM );
t = mtimesx( permute( mean_normals, [ VEC DIM PAGE ] ), permute( tr, [ DIM VEC PAGE ] ) );
t = permute( t, [ 2 3 1 ] );

ur = cross( ( a - c ), repmat( q, [ size( a, 1 ) 1 1 ] ), DIM );
u = mtimesx( permute( mean_normals, [ VEC DIM PAGE ] ), permute( ur, [ DIM VEC PAGE ] ) );
u = permute( u, [ 2 3 1 ] );

v = atan2( u, t );

[ ~, indices ] = sort( v, VEC );
indices = indices + ( 0 : size( points, VEC ) : numel( points( :, :, 1 ) ) - 1 );
indices = indices + permute( ( 0 : numel( points( :, :, 1 ) ) : numel( points ) - 1 ), [ 1 3 2 ] );
indices = indices + 1;
first = 1 : size( points, VEC ) : numel( points( :, :, 1 ) ) - 1;
first = first + permute( ( 0 : numel( points( :, :, 1 ) ) : numel( points ) - 1 ), [ 1 3 2 ] );
indices = cat( 1, first, indices );

end

