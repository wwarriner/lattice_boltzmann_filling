function [ u, offsets ] = revert_sort( s, i, dim )

[ sp, inv ] = rotate_to_dimension( dim, s );
ip = rotate_to_dimension( dim, i );

sz = size( sp );
head = sz( 1 );
tail = sz( 2 : end );
tail_count = prod( tail );

offsets = head * ( 0 : ( tail_count - 1 ) );
offsets = reshape( offsets, [ 1 tail ] );
offsets = repmat( offsets, [ head ones( 1, numel( tail ) ) ] );

[ ~, i_rev ] = sort( ip + offsets, 1 );
offsets = i_rev + offsets;

u = sp( offsets );
u = rotate_from_dimension( u, inv );

end

