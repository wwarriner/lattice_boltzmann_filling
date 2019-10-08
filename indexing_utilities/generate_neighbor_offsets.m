function offsets = generate_neighbor_offsets( a, include_center )

if nargin < 2
    include_center = false;
end

if isvector( a )
    sz = a;
else
    sz = size( a );
end

strides = [ 1 cumprod( sz( 1 : end - 1 ) ) ];

subs = generate_neighbor_subs( a, include_center );
offsets = sum( subs .* strides, 2 );

end

