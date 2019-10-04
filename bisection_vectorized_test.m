function bisection_vectorized_test()

fn = @sin;
x_range = [ ...
    linspace( -pi/4, -pi/8, 1000 ).' ...
    linspace( pi/8, pi/4, 1000 ).' ...
    ];
x = bisection_vectorized( fn, x_range, 100 );
assert( all( abs( fn( x ) ) <= 1e-10, "all" ) );

fn = @(x)x-0.01;
x_range = [ ...
    linspace( -1, -0.02, 1000 ).' ...
    linspace( 0.02, 1, 1000 ).' ...
    ];
x = bisection_vectorized( fn, x_range, 100 );
assert( all( abs( fn( x ) ) <= 1e-10, "all" ) );

end

