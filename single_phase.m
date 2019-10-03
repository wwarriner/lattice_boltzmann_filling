%%
% shanChen.m: Multi-component fluid, using a LB method,
%   based on the Shan-Chen model
% [X.Shan and H.Chen, http://dx.doi.org/10.1103/PhysRevE.47.1815].
%
%
% Lattice Boltzmann sample, written in Matlab
% Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
% Address: EPFL-STI-LIN Station 9
% E-mail: orestis.malaspinas@epfl.ch
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public
% License along with this program; if not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
% Boston, MA  02110-1301, USA.
%

%% D2Q9 LATTICE CONSTANTS
t = [ 16 4 4 4 4 1 1 1 1 ].' ./ 36;
ex = [ 0 1 0 -1 0 1 -1 -1 1 ].';
ey = [ 0 0 1 0 -1 1 1 -1 -1 ].';
opp = [ 1 4 5 2 3 8 9 6 7 ].';
q = numel( t );

%% GEOMETRY
shape = [ 26 101 ];
count = prod( shape );
stride = [ 1 cumprod( shape( 1 : end - 1 ) ) ];
aug_shape = [ 1 shape ];
flat_shape = [ q count ];
full_shape = [ q shape ];
full_count = prod( full_shape );
full_stride = [ 1 cumprod( full_shape( 1 : end - 1 ) ) ];

%% MATERIAL IDS
% AIR = 0;
% SOLID = 1;
% SOLID_BOUNDARY = 2;
% FLUID = 3;
% INTERFACE = 4;

%% BOUNDARIES
SOLID = false( aug_shape );
SOLID( 1, 1, : ) = true;
SOLID( 1, end, : ) = true;
SOLID( 1, :, 1 ) = true;
SOLID( 1, :, end ) = true;
% ub = 1 : stride( 1 ) : stride( 2 );
% lb = ub + count - stride( 2 );
% fb = 1 : stride( 2 ) : count;
% rb = fb + stride( 2 ) - 1;
% SOLID = unique( [ ub lb fb rb ] );

% [x,y] = ind2sub( shape, SOLID );
% x = repmat( x, [ q 1 ] );
% y = repmat( y, [ q 1 ] );
% z = ( 0 : q-1 ).' + ones( 1, numel( SOLID ) );
% FULL_SOLID = sub2ind( full_shape, z, x, y );

%% TIME STEP MANAGEMENT
time_step_count = 1000;

%% UNITS
% tau = nu_lbm / u_s_lbm ^ 2 + 1/2 [ unitless ]
% nu_lbm - viscosity in LBM simulation units [ unitless ]
% u_s_lbm - speed of sound in LBM simulation units, usually 1/3 [ unitless ]
%
% nu_lbm = ( dt / dx^2 ) * ( 1 / Re )
% dt - time step in LBM simulation units (1 / number of iterations) [ unitless ]
% dx - space step in LBM simulation units (1 / number of cells) [ unitless ]
% Re - Reynolds number [ unitless ]
%  - Our case will be very close to zero
%
% Re = u_0_p * ell_0_p / nu_p
% u_0_p - initial velocity in physical units (probably max?) [ m / s ]
%  - For most filling applications this will be 1
%  - Based on recent research we probably want it to be even lower, but
%  that may be impractical.
% ell_0_p - characteristic length in physical units [ m ]
%  - Essentially just the shortest length of the bounding box.
%  - Alternately could use 2 * max EDT value in the casting.
% nu_p - kinematic viscosity in physical units [ m^2 / s ]
%  - Fe will be about 6.944e-7
%
% nu_p = mu_p / rho_p
% mu_p is the dynamic viscosity in physical units [ kg / m*s ]
%  - Fe is about 0.005 in at melt
% rho_p is the density in physical units [ kg / m^3 ]
%  - Fe is about 7200 at melt
%
% g_lbm = ( d_t^2 / d_x ) * g_p
% d_t, d_x as above
% g_p is the acceleration due to gravity [ m / s^2 ]
%  - This is always (approx) 9.8

%% PHYSICAL PARAMETERS
rho_p = 7200; % [ kg / m^3 ]
mu_p = 0.0055; % [ kg / m*s ]
nu_p = rho_p * mu_p; % [ m^2 / s ];
u_0_p = 1.0; % [ m / s ]
ell_0_p = 0.1; % [ m ]
Re = u_0_p * ell_0_p / nu_p; % [ - ]
dt = 1e-7;
dx = 1 / max( shape );
nu_lbm = ( dt / dx^2 ) * ( 1 / Re ); % [ - ]
c_s_lbm = 1.0 / sqrt( 3.0 ); % [ - ]
tau = nu_lbm / c_s_lbm^2 + 1.0 / 2.0; % [ - ]
assert( 0.5 < tau & tau < 5.0 );
% tau = 1.; % relaxation parameter
g_p = 9.8; % [ m / s^2 ]
g_lbm = ( dt^2 / dx ) * g_p;
sg_lbm = sqrt( 2 ) * g_lbm;
dg = -[ 0 0 -g_lbm 0 g_lbm -sg_lbm -sg_lbm sg_lbm sg_lbm ].';
tdg = t .* dg;

ct = ( ( tau - 0.5 ) / 3.0 ) * ( ell_0_p * dx ) ^ 2 / nu_p;
t_p = ct/dt;
u_p_max = ( ell_0_p * dx ) / t_p;

%% INITIAL CONDITIONS
drho = 0.01;
rho_lbm = 1.0;
delta_rho = permute( linspace( -drho, drho, shape( 2 ) ), [ 3 1 2 ] );
delta_rho = repmat( delta_rho, [ q shape( 1 ) 1 ] );
f_in = t .* ( rho_lbm + delta_rho );
f_in( :, :, floor( full_shape(3) *0.25 ) : end ) = 0;
% find interface
LIQUID = 0 < f_in(1,:,:) & ~SOLID;
EMPTY = f_in(1,:,:) == 0 & ~SOLID;
max_conn = conndef( ndims( squeeze( f_in(1,:,:) ) ), "maximal" );
max_conn = reshape( max_conn, [ 1 size( max_conn ) ] );
INTERFACE = LIQUID & imdilate( EMPTY, max_conn );
LIQUID = LIQUID & ~INTERFACE;
finterface = find( INTERFACE );
for i = 1 : q
    f_in( i, finterface ) = t(i) * rho_lbm;
end

%% PLOT SETUP
rho = squeeze( sum( f_in ) );
rho( SOLID ) = 0;
fh = figure();
axh = axes( fh );
ih = imagesc( axh, rho.' );
colorbar( axh );
max_rho = max( rho, [], "all" );
caxis( axh, [ 0 max_rho ] );
title( axh, "Fluid 1 density" );
axis( axh, "equal", "off" );

%% DEBUGGING
total_mass = zeros( time_step_count, 1 );
max_u = zeros( time_step_count, 1 );
max_f = zeros( time_step_count, 1 );

%% MAIN LOOP
mass = zeros( aug_shape );
mass( LIQUID ) = 1;
mass( INTERFACE ) = 0.5;
for cycle = 1 : time_step_count
    %% MACROSCOPIC VARIABLES
    rho = sum( f_in );
    rho( EMPTY ) = 1;
    total_mass( cycle ) = sum( rho, "all" );
    
    ux = reshape( ( ex.' * reshape( f_in, flat_shape ) ), aug_shape );
    uy = reshape( ( ey.' * reshape( f_in, flat_shape ) ), aug_shape );
    u = ux .^ 2 + uy .^ 2;
    max_u( cycle ) = sqrt( max( u, [], "all" ) );
    
    %% COLLISION STEP FLUID 1 AND 2
    cu = 3 * ( ex .* ux + ey .* uy );
    f_eq = t .* ( rho + cu + 0.5 * cu .^ 2 - 1.5 * u );
    %f_out = ( 1 - tau ) * f_in + tau * ( f_eq + rho .* tdg );
    max_f( cycle ) = max( f_out, [], "all" );
    
    %% MASS UPDATE
    dmass = zeros( aug_shape );
    fracs = zeros( aug_shape );
    b = zeros( aug_shape );
    for i = 1 : q
        a = circshift( f_in(i,INTERFACE), [ 0 ex(i) ey(i) ] ) - f_in(i,INTERFACE);
        b( INTERFACE ) = mass( INTERFACE ) ./ rho( INTERFACE );
        b( LIQUID ) = 1;
        b = ( circshift( b, [ 0 ex( i ) ey( i ) ] ) + b ) / 2;
        dmass( INTERFACE ) = dmass( INTERFACE ) ...
            + a.' .* b( INTERFACE );
    end
    mass = mass + dmass;
    
    %% INTERFACE BOUNCE-BACK
    INTERFACE_EMPTY = imdilate( INTERFACE, max_conn ) & EMPTY;
    EMPTY_INTERFACE = imdilate( EMPTY, max_conn ) & INTERFACE;
    for i = 1 : q
        ie = circshift( INTERFACE_EMPTY, [ 0 ex(i) ey(i) ] );
        x = ie & EMPTY & ~SOLID;
        f_eq( i, x ) = t( i ) .* ( 1 + cu( i, x ) + 0.5 * cu( i, x ) .^ 2 - 1.5 * u( 1, x ) );
        f_eq( opp( i ), x ) = f_eq( i, x );
    end
    f_out = ( 1 - tau ) * f_in + tau * ( f_eq + rho .* tdg );
        
    %% OBSTACLE (BOUNCE-BACK)
    for i= 1 : q
        f_out( i, SOLID ) = f_in( opp( i ), SOLID );
    end
    
    %% STREAMING STEP
    for i= 1 : q
        f_in(i,:,:) = circshift( f_out(i,:,:), [ 0 ex(i) ey(i) ] );
    end
    
    %% INTERFACE TRACKING
    % compute mass change for interface (loop over q)
    %  if neighbor n is liquid, f_in_n - f_in_self
    %  if neighbor n is interface, (f_in_n - f_in_self) * (mean m/rho)
    %  this can be summed up over q as we go, no need to track explicitly
    
    % add mass change to mass
    
    % compute fluid fraction as mass / rho
    % compute approx surface normals at interface
    %  weighted sum of negated distance to neighbor vectors which are also interface
    %  weights are 1 for face, sqrt(2) for edge, sqrt(3) for vertex neighbors
    %  then normalize to get unit vector
    
    % reconstruct streaming from empty cells
    %  f_out = f_eq(rho_A,u) + opp f_eq(rho_A,u) - f_in
    % reconstruct streaming to empty cells same way
    % reconstruct streaming to cells for which ci dot normal less than zero
    % FIRST find set of neighbors which are empty and
    
    % determine filled cells
    
    FILLED = INTERFACE & mass > ( 1.001 * rho );
    EMPTIED = INTERFACE & mass < ( -0.001 * rho );
    mass_ex = sum( mass( finterface ) - rho( finterface ) );
    
    LI = LIQUID | INTERFACE;
    rho( ~LI ) = 0;
    rho_avg = zeros( aug_shape );
    ux_avg = zeros( aug_shape );
    uy_avg = zeros( aug_shape );
    n = zeros( aug_shape );
    for i = 2 : q
        rho_avg = rho_avg + circshift( rho, [ 0 ex(i) ey(i) ] );
        ux_avg = ux_avg + circshift( ux, [ 0 ex(i) ey(i) ] );
        uy_avg = uy_avg + circshift( uy, [ 0 ex(i) ey(i) ] );
        n = n + circshift( LI, [ 0 ex(i) ey(i) ] );
    end
    rho_avg = rho_avg ./ n;
    ux_avg = ux_avg ./ n;
    uy_avg = uy_avg ./ n;
    u_avg = ux_avg .^ 2 + uy_avg .^ 2;
    cu_avg = 3 * ( ex .* ux_avg + ey .* uy_avg );
    
    EMPTY = ~( LIQUID | SOLID | INTERFACE );
    LIQUID = ( LIQUID | FILLED ) & ~( EMPTIED | SOLID );
    INTERFACE2 = EMPTY & imdilate( LIQUID, max_conn );
    finterface = find( INTERFACE2 );
    for i = 1 : q
        f_in( i, finterface ) = t( i ) .* ( rho_avg( 1, finterface ) + cu_avg( i, finterface ) + 0.5 * cu_avg( i, finterface ) .^ 2 - 1.5 * u_avg( 1, finterface ) );
    end
    
    %% VISUALIZATION
    rho = squeeze( rho );
    rho( SOLID ) = 0;
    ih.CData = rho';
    max_rho = max( max( rho, [], "all" ), max_rho );
    caxis( axh, [ 0 max_rho ] );
    drawnow( 'limitrate' )
    
    INTERFACE = INTERFACE2;
end
