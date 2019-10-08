function mean_normals = calculate_mean_normals( normals )

% missing normals must be at the end of that PAGE in the VEC dimension

% DIMENSION MEANING
%VEC = 1;
%PAGE = 2;
DIM = 3;

% calculate
mean_normals = nanmean( normals, 1 );
mean_normals = mean_normals ./ vecnorm( mean_normals, 2, DIM );

end

