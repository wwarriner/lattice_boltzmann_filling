function a = post_rep_fill( a, unary_logical_fn, dim )

if nargin < 3
    dim = find( size( a ) > 1, 1, "first" );
end

p = 1 : ndims( a );
p( dim ) = [];
p = [ dim p ];

a = permute( a, p );
b = unary_logical_fn( a( : ) );
c = a( b );
d = cumsum( b );
a = reshape( c( d ), size( a ) );
a = ipermute( a, p );

end

