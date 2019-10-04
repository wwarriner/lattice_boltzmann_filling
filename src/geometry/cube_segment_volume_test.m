function cube_segment_volume_test()

rng( 314159 );

n = 2 * rand( [ 1000 3 ] ) - 1;
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

off = rescale( rand( [ 1000 1 ] ), -sqrt( 3 ) / 2, 0 );
off = [ -sqrt( 3 ) / 2; -sqrt( 2 ) / 2; off ];
off = [ off; 0.0; -off ];
off = sort( off );

[ n_index, off ] = meshgrid( 1 : size( n, 1 ), off );

fprintf( 1, "%i elements calculated in:\n", numel( n_index ) );
t = tic;
v = cube_segment_volume( n( n_index( : ), : ), off( : ) );
v = reshape( v, size( n_index ) );
toc( t );

assert( all( 0.0 <= v, "all" ) );
assert( all( v <= 1.0, "all" ) );
assert( all( abs( v + flip( v, 1 ) - 1 ) <= 1e-5, "all" ) );

end

