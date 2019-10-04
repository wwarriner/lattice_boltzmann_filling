function mean_normals = calculate_mean_normals( normals )

% missing normals must be at the end of that PAGE in the VEC dimension

% DIMENSION MEANING
VEC = 1;
PAGE = 2;
DIM = 3;

% handle missing
% no influence on mean normal
invalid = any( isnan( normals ), DIM );
count = numel( invalid );
invalid = find( invalid );
invalid = invalid + ( 0 : count : numel( normals ) - 1 );
normals( invalid ) = 0;

% calculate
mean_normals = mean( normals, 1 );
mean_normals = mean_normals ./ vecnorm( mean_normals, 2, DIM );

end

