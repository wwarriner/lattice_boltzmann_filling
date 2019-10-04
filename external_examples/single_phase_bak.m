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
AIR = 0;
SOLID = 1;
SOLID_BOUNDARY = 2;
FLUID = 3;
INTERFACE = 4;

%% BOUNDARIES
ub = 1 : stride( 1 ) : stride( 2 );
lb = ub + count - stride( 2 );
fb = 1 : stride( 2 ) : count;
rb = fb + stride( 2 ) - 1;
SOLID = unique( [ ub lb fb rb ] );

ub = 1 : full_stride( 2 ) : full_stride( 3 );
ub = ub + ( 0 : full_stride( 2 ) - 1 ).';
lb = ub + full_count - full_stride( 3 ) - full_stride( 2 );
fb = 1 : full_stride( 3 ) : full_count;
fb = fb + ( 0 : full_stride( 2 ) - 1 ).';
rb = fb + full_stride( 3 ) - full_stride( 2 );
FULL_SOLID = unique( [ ub lb fb rb ] );

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
dt = 1e-8;
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

ct = ( ( tau - 0.5 ) / 3.0 ) * (ell_0_p*dx)^2/nu_p;
t_p = ct/dt;
u_p_max = (ell_0_p*dx)/t_p;

%% INITIAL CONDITIONS
drho = 0.0001;
rho_lbm = 1.0;
delta_rho = permute( linspace( -drho, drho, shape( 2 ) ), [ 3 1 2 ] );
delta_rho = repmat( delta_rho, [ q shape( 1 ) 1 ] );
f_in = t .* ( rho_lbm + delta_rho );

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
for cycle = 1 : time_step_count
    %% MACROSCOPIC VARIABLES
    rho = sum( f_in );
    total_mass( cycle ) = sum( rho, "all" );
    
    ux = reshape( ( ex.' * reshape( f_in, flat_shape ) ), aug_shape ) ./ rho;
    uy = reshape( ( ey.' * reshape( f_in, flat_shape ) ), aug_shape ) ./ rho;
    max_u( cycle ) = sqrt( max( ux.^2 + uy.^2, [], "all" ) );
    
    %% COLLISION STEP FLUID 1 AND 2
    cu = 3 .* ( ex .* ux + ey .* uy );
    
    % f_eq = t .* ( rho + 3 * cu 
    
    f_eq = rho .* t .* ( 1 + cu + 0.5 .* ( cu .* cu ) - 1.5 .* ( ux .^ 2 + uy .^ 2 ) );
    f_out = f_in - tau .* ( f_in - f_eq ) + tau .* rho .* tdg;
    % NEED to update mass in interface cells
    
    % GAS => cell = 0
    % LIQUID => cell = f_opp_neighbor + f_self
    % SOLID => cell = 0.5 * ( phi_me + phi_neighbor ) * ( 
    
    max_f( cycle ) = max( f_out, [], "all" );
    
    %% OBSTACLE (BOUNCE-BACK)
    for i=1:9
        f_out( i, SOLID ) = f_in( opp( i ), SOLID );
    end
    
    %% STREAMING STEP FLUID 1 AND 2
    for i=1:9
        f_in(i,:,:) = circshift( f_out(i,:,:), [ 0 ex(i) ey(i) ] );
    end
    
    %% VISUALIZATION
    rho = squeeze( rho );
    rho( SOLID ) = 0;
    ih.CData = rho';
    max_rho = max( max( rho, [], "all" ), max_rho );
    caxis( axh, [ 0 max_rho ] );
    drawnow( 'limitrate' )
end
