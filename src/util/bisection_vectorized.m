function x = bisection_vectorized( ...
    fn, ...
    x_range, ...
    iteration_count ...
    )

if nargin < 3
    iteration_count = 10;
end

assert( isa( fn, "function_handle" ) );

assert( isa( x_range, "double" ) );
assert( size( x_range, 2 ) == 2 );

x_range = sort( x_range, 2 );
x_range = x_range.'; % eases indexing without using sub2ind

ROWS = 2;
cols = size( x_range, 2 );
for i = 1 : iteration_count
    x = mean( x_range, 1 );
    y = fn( x.' ).';
    smaller = y < 0;
    % if smaller, replace index 1 (min)
    % else, replace index 2 (max)
    replace_index = 2 - smaller; % gets column index
    replace_index = replace_index + ROWS .* ( 0 : cols - 1 ); % gets flat index
    x_range( replace_index ) = x;
end
x = x.';

end

