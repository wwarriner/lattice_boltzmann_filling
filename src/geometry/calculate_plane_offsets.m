function offsets = calculate_plane_offsets( normals, fills )

assert( isa( normals, "double" ) );
assert( size( normals, 2 ) == 3 );

assert( isa( fills, "double" ) );
assert( isvector( fills ) );

count = numel( fills );

range = [ -sqrt( 3 ) / 2, sqrt( 3 ) / 2 ];
range = repmat( range, [ count 1 ] );

fn = @(x)cube_segment_volume( normals, x ) - fills;

offsets = bisection_vectorized( fn, range );

end

