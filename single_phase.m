%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shanChen.m: Multi-component fluid, using a LB method,
%   based on the Shan-Chen model
% [X.Shan and H.Chen, http://dx.doi.org/10.1103/PhysRevE.47.1815].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample, written in Matlab
% Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
% Address: EPFL-STI-LIN Station 9
% E-mail: orestis.malaspinas@epfl.ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% D2Q9 LATTICE CONSTANTS
t = [ 16 4 4 4 4 1 1 1 1 ].' ./ 36;
cx = [ 0 1 0 -1 0 1 -1 -1 1 ].';
cy = [ 0 0 1 0 -1 1 1 -1 -1 ].';
opp = [ 1 4 5 2 3 8 9 6 7 ].';
q = numel( t );

% GEOMETRY
% SPACE STEP HAS UNITS m
shape = [ 26 101 ];
count = prod( shape );
stride = [ 1 cumprod( shape( 1 : end - 1 ) ) ];
aug_shape = [ 1 shape ];
flat_shape = [ q count ];
full_shape = [ q shape ];
full_count = prod( full_shape );
full_stride = [ 1 cumprod( full_shape( 1 : end - 1 ) ) ];

% BOUNDARIES
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

% TIME STEP MANAGEMENT
time_step_count = 8000;

% UNITS
% tau = nu_lbm/u_s_lbm^2 + 1/2
% nu_lbm - viscosity in LBM simulation units
% u_s_lbm - speed of sound in LBM simulation units, usually 1/3

% PHYSICAL PARAMETERS
tau = 1.; % relaxation parameter
g = 0.001; % gravitational parameter [ s * m/s^2 ]
sg = sqrt( g );
dg = -[ 0 0 -g 0 g -sg -sg sg sg ].';
tdg = t .* dg;

% BOUNDARY CONDITIONS

% RHO has units [ kg / m^3 ]
drho = 1;
delta_rho = permute( linspace( 1 - drho, 1 + drho, shape( 2 ) ), [ 3 1 2 ] );
delta_rho = repmat( delta_rho, [ q shape( 1 ) 1 ] );
f_in = t .* ( 1.0 + delta_rho );

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

% MAIN LOOP (TIME CYCLES)
total_mass = zeros( time_step_count, 1 );
for cycle = 1:time_step_count
    % MACROSCOPIC VARIABLES
    rho = sum( f_in );
    total_mass( cycle ) = sum( rho, "all" );
    ux = reshape( ( cx.' * reshape( f_in, flat_shape ) ), aug_shape ) ./ rho;
    uy = reshape( ( cy.' * reshape( f_in, flat_shape ) ), aug_shape ) ./ rho;
    
    % COLLISION STEP FLUID 1 AND 2
    cu = 3 .* ( cx .* ux + cy .* uy );
    f_eq = rho .* t .* ( 1 + cu + 0.5 .* ( cu .* cu ) - 1.5 .* ( ux .^ 2 + uy .^ 2 ) );
    f_out = f_in - tau .* ( f_in - f_eq ) + tau .* rho .* tdg;
    
    % OBSTACLE (BOUNCE-BACK)
    for i=1:9
        f_out( i, SOLID ) = f_in( opp( i ), SOLID );
    end
    
    % STREAMING STEP FLUID 1 AND 2
    for i=1:9
        f_in(i,:,:) = circshift( f_out(i,:,:), [ 0 cx(i) cy(i) ] );
    end
    
    % VISUALIZATION
    rho = squeeze( rho );
    rho( SOLID ) = 0;
    ih.CData = rho';
    max_rho = max( max( rho, [], "all" ), max_rho );
    caxis( axh, [ 0 max_rho ] );
    drawnow( 'limitrate' )
end
