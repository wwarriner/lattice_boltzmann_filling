% cube sides are all interval [ -0.5, 0.5 ]
% center is [ 0, 0, 0 ]

% normal must be unit normal
% offset is distance from center to plane, negative means plane is on negative
% side of origin, positive means plane on positive side
function volume = cube_segment_volume( normal, offset )

% anti-symmetry of offset about origin
flip = false;
if 0 < offset
    offset = -offset;
    flip = true;
end

normal = sort( abs( normal ) ); % take advantage of remaining cube symmetries

origin_dist = dot( normal, [ -0.5 -0.5 -0.5 ] ) - offset;
if 0 <= origin_dist
    volume = 0;
else
    n = dot( normal, [ 0.5 0.5 -0.5 ] ) - offset;
    if n <= 0
        zi = [ -0.5 -0.5; -0.5 0.5; 0.5 -0.5; 0.5 0.5 ];
        z = compute_intercept( normal( [ 1 2 ] ).', zi, normal( 3 ), offset );
        z = z + 0.5;
        volume = mean( z );
        volume = max( min( volume, 1 ), 0 );
    elseif normal( 1 ) == 0
        zi = [ -0.5 -0.5 ];
        z = compute_intercept( normal( [ 1 2 ] ).', zi, normal( 3 ), offset );
        z = max( min( z + 0.5, 1 ), 0 );
        yi = [ -0.5 -0.5 ];
        y = compute_intercept( normal( [ 1 3 ] ).', yi, normal( 2 ), offset );
        y = max( min( y + 0.5, 1 ), 0 );
        volume = y( 1 ) * z( 1 ) / 2;
    else

        zi = [ -0.5 -0.5; -0.5 0.5; 0.5 -0.5 ];
        z = compute_intercept( normal( [ 1 2 ] ).', zi, normal( 3 ), offset );
        z = z + 0.5;
        zc = max( min( z, 1 ), 0 );
        yi = [ -0.5 -0.5; -0.5 0.5; 0.5 -0.5 ];
        y = compute_intercept( normal( [ 1 3 ] ).', yi, normal( 2 ), offset );
        y = y + 0.5;
        yc = max( min( y, 1 ), 0 );
        xi = [ -0.5 -0.5; -0.5 0.5; 0.5 -0.5 ];
        x = compute_intercept( normal( [ 2 3 ] ).', xi, normal( 1 ), offset );
        x = x + 0.5;
        xc = max( min( x, 1 ), 0 );
        V100 = yc( 3 ) * zc( 3 ) * ( x( 1 ) - xc( 1 ) );
        V010 = xc( 3 ) * zc( 2 ) * ( y( 1 ) - yc( 1 ) );
        V001 = xc( 2 ) * yc( 2 ) * ( z( 1 ) - zc( 1 ) );
        V000 = x( 1 ) * y( 1 ) * z( 1 );
        volume = ( V000 - ( V100 + V010 + V001 ) ) / 6;
    end
end

% anti-symmetry of offset about origin
if flip
    volume = 1 - volume;
end

end


function p = compute_intercept( normal, coords, coeff, offset )

p = -( coords * normal - offset ) ./ coeff;

end

