function calculate_plane_offsets_test()

rng( 217828 );

n = 2 * rand( [ 100 3 ] ) - 1;
n = [ 
    n;
    0 0 1;
    0 1 0;
    1 0 0;
    1 1 0;
    1 0 1;
    0 1 1;
    1 1 1;
    0.0480264745524219 0.684013838616227 0.727886341624542
    ];
n = n ./ vecnorm( n, 2, 2 );
n = sortrows( n );

fills = rand( [ 100 1 ] );

[ n_index, fills ] = meshgrid( 1 : size( n, 1 ), fills );

fprintf( 1, "%i elements calculated in:\n", numel( n_index ) );
t = tic;
offsets = calculate_plane_offsets( n( n_index( : ), : ), fills( : ) );
toc( t );

assert( 0.0 <= all( fills( : ), "all" ) );
assert( all( fills( : ) <= 1.0, "all" ) );

v = cube_segment_volume( n( n_index( : ), : ), offsets );
assert( all( abs( fills( : ) - v ) < 1e-2, "all" ) );

end

