function k = calculate_surface_curvature( points, normals )

% k = 0 means locally flat
% k < 0 means normal on opposite side of tangent plane as surface
%  (outside local sphere)
% k > 0 means normal on same side of tangent plane as surface
%  (inside local sphere)

% first point in points is central point
% remainder arranged counter-clockwise wrt first normal
% normals associated with points of same index

% missing points must be at the end of that PAGE in the VEC dimension

% DIMENSION MEANING
VEC = 1;
PAGE = 2;
DIM = 3;
DUMMY = 4;

points = post_rep_fill( points, @(x)~isnan(x) );
normals = post_rep_fill( normals, @(x)~isnan(x) );

% determine triangle face areas
p = points( 2 : end, :, : ) - points( 1, :, : ); % by VEC
N = cross( p, -circshift( p, VEC ), DIM );
vn = vecnorm( N, 2, DIM );
areas = 0.5 * vn;
vn = post_rep_fill( vn, @(x)x~=0 );

% compute weighted mean normal
% weights are face areas
%vn( vn == 0 ) = 1; % hack fix vec norms for sliver triangles
N = N ./ vn;
N = sum( areas .* N, VEC );
N = N ./ vecnorm( N, 2, DIM );
% N is one DIM-VEC per PAGE

% compute T_i as
% T_i = ( eye - N*N.' ) * (p_i - p)
% note first part is constant, collapse normals to that value
t_i_pre = eye( 3 ) - mtimesx( permute( N, [ DIM VEC PAGE ] ), permute( N, [ VEC DIM PAGE ] ) );
t_i = mtimesx( t_i_pre, permute( p, [ DIM VEC PAGE ] ) );
t_i = permute( t_i, [ 2 3 1] ); % was in DIM VEC PAGE order
t_i = t_i ./ vecnorm( t_i, 2, DIM );
t_i = mtimesx( permute( t_i, [ DIM DUMMY VEC PAGE ] ), permute( t_i, [ DUMMY DIM VEC PAGE ] ) );
t_i = permute( t_i, [ 3 4 1 2 ] ); % was in DIM DIM VEC PAGE
% t_i is one DIM-DIM matrix per PAGE

% compute k_i as
% k_i = 2*N.' * (p_i - p) ./ vecnorm(p_i - p)^2
% latter part better hand-computed to avoid sqrt
k_i = 2 .* mtimesx( permute( N, [ VEC DIM PAGE ] ), permute( p, [ DIM VEC PAGE ] ) );
k_i = permute( k_i, [ 2 3 1 ] ) ./ sum( p .^ 2, 3 ); % was one face_count page
% k_i is one value per VEC per PAGE

% compute lambda_i as
% lambda_i = [face area of all faces using p_i] / [2 * sum face areas]
% scalar
invalid = areas == 0;
areas = circshift( areas, -1 );
areas = flip( post_rep_fill( flip( areas, VEC ), @(x)x~=0 ), VEC );
areas = circshift( areas, 1 );
lambda_i = areas + circshift( areas, -1 );
lambda_i( invalid ) = 0;
lambda_i = lambda_i ./ sum( lambda_i );
% TODO may be a miscalculation
% lambda_i is one value per face per page

% to vectorize, some points will be "nonsense" points
% in those cases, set product below to 0
% compute M = sum( lambda_i * k_i * T_i * T_i.' );
M = sum( lambda_i .* k_i .* t_i, VEC );
M = squeeze( M );
M = permute( M, [ 2 3 1 ] );
% M is one matrix per page

% find traces
diag_inds = [ 1; 5; 9 ]; % diagonals of 3x3 matrix
diag_inds = diag_inds + ( 0 : 9 : numel( M ) - 1 ); % extended to pages
k = sum( M( diag_inds ) );

end

