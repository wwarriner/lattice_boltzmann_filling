function determine_surface_points_test()

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

fills = sort( rand( [ 100 1 ] ) );

[ n_index, fills ] = meshgrid( 1 : size( n, 1 ), fills );

fprintf( 1, "%i elements calculated in:\n", numel( n_index ) );
t = tic;
surface_points = determine_surface_points( n( n_index( : ), : ), fills( : ) );
toc( t );

assert( all( 0.0 <= surface_points, "all" ) );
assert( all( surface_points <= 1.0, "all" ) );

end

