function points = determine_surface_points( normals, fills )

assert( isa( normals, "double" ) );
assert( size( normals, 2 ) == 3 );

assert( isa( fills, "double" ) );
assert( isvector( fills ) );

points = calculate_plane_offsets( normals, fills ) .* normals + [ 0.5 0.5 0.5 ];

end

