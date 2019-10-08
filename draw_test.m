function canvas = draw_test( interface, int_to_air, int_to_liq, air_to_int, liq_to_int )

canvas = zeros( [ size( squeeze( interface ) ) 3 ] );
canvas = imoverlay( canvas, squeeze( int_to_air ), [ 1 0 0 ] );
canvas = imoverlay( canvas, squeeze( air_to_int ), [ 0 1 0 ] );
canvas = imoverlay( canvas, squeeze( int_to_liq ), [ 0 1 1 ] );
canvas = imoverlay( canvas, squeeze( liq_to_int ), [ 1 0 1 ] );
%canvas = imoverlay( canvas, squeeze( interface ), [ 1 1 1 ] );
imtool( canvas );

end

