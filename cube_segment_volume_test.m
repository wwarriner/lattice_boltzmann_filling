function cube_segment_volume_test()

n = sort( rand( [ 100 3 ] ), 1 );
n = n ./ vecnorm( n, 2, 2 );

off = sort( rescale( rand( [ 100 1 ] ), -sqrt( 3 ) / 2, 0 ) );
off = [ off; 0.0; flip( -off ) ];

[ n_index, off ] = meshgrid( 1 : size( n, 1 ), off );

v = nan( size( n_index ) );
tic;
for i = 1 : numel( n_index )
    v( i ) = cube_segment_volume( n( n_index( i ), : ), off( i ) );
end
toc;

assert( all( abs( v + flip( v, 1 ) - 1 ) <= 1e-5, "all" ) );

end

