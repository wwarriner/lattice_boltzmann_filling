function core()
%% LATTICE CONSTANTS
% D2Q9
d = 2;
w = [ 16 4 4 4 4 1 1 1 1 ].' ./ 36;
cx = [ 0 1 0 -1 0 1 -1 -1 1 ].';
cy = [ 0 0 1 0 -1 1 1 -1 -1 ].';
opp = [ 1 4 5 2 3 8 9 6 7 ].';
assert( abs( sum( w ) - 1 ) < 1e-10 );

% D3Q19
% d = 3;
% w = [ 12 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 ].' ./ 36;
% cx = [ 0 1 0 0 -1 0 0 1 -1 -1 1 1 1 -1 -1 0 0 0 0 ];
% cy = [ 0 0 1 0 0 -1 0 1 1 -1 -1 0 0 0 0 1 -1 -1 1 ];
% cz = [ 0 0 0 1 0 0 -1 0 0 0 0 1 -1 -1 1 1 1 -1 -1 ];
% opp = [ 1 5 6 7 2 3 4 10 11 8 9 14 15 12 13 18 19 16 17 ];
% assert( abs( sum( w ) - 1 ) < 1e-10 );

q = numel( w );

%% GEOMETRY
shape = [ 26 101 ];
count = prod( shape );
stride = [ 1 cumprod( shape( 1 : end - 1) ) ];
aug_shape = [ 1 shape ];
flat_shape = [ q count ];
full_shape = [ q shape ];
full_count = prod( full_shape );
full_stride = [ 1 cumprod( full_shape( 1 : end - 1 ) ) ];

n_off = generate_neighbor_offsets( shape );

%% MATERIAL IDS
SOLID = 0;
AIR = 1;
LIQUID = 2;
INTERFACE = AIR + LIQUID;

%% TIME STEP
time_step_count = 1000;

%% PHYSICAL PARAMETERS
rho_p = 7200; % density [ kg / m^3 ]
mu_p = 0.0055; % dynamic viscosity [ kg / m*s ]
nu_p = rho_p * mu_p; % kinematic viscosity [ m^2 / s ]
u_0_p = 1.0; % initial_velocity [ m / s ]
ell_0_p = 0.1; % characteristic length [ m ]

Re = u_0_p * ell_0_p / nu_p; % Reynold's number [ - ]

dt = 1e-8;
dx = 1 / max( shape );

nu_lbm = ( dt / dx ^ 2 ) * ( 1 / Re ); % LBM kinematic viscosity [ - ]
c_s_lbm = 1.0 / sqrt( 3.0 ); % LBM speed of sound, D2Q9, D3Q19 [ - ]
tau = nu_lbm / c_s_lbm^2 + 1.0 / 2.0; % LBM relaxation time [ - ]
assert( 0.5 < tau & tau < 5.0 );

g_p = 9.8; % gravitational acceleration [ m / s^2 ]
g_lbm = ( dt ^ 2 / dx ) * g_p; % LBM gravitational acceleration [ - ]
sg_lbm = g_lbm ./ sqrt( 2 ); % edge diagonal
dg = -g_lbm .* cy; % gravitational force D2Q9
% dg = -[ g_lbm .* cz( 1 : 7 ); sg_lbm .* cz( 8 : end ) ]; % gravitational force D3Q19
wdg = w .* dg;

rho_g = 1.0;

%% VARIABLE SETUP
% macroscopic variables
rho = zeros( aug_shape );
mass = zeros( aug_shape );
fill = zeros( aug_shape );
bubble = zeros( aug_shape );

% microscopic variabls
f_inter = zeros( full_shape );
f = zeros( full_shape );

%% INITIAL CONDITIONS
% GEOMETRY MUST HAVE OUTER SOLID BOUNDARY!

INITIAL_SOLID = false( shape );
INITIAL_SOLID( 1, : ) = true;
INITIAL_SOLID( end, : ) = true;
INITIAL_SOLID( :, 1 ) = true;
INITIAL_SOLID( :, end ) = true;

INITIAL_SPACE = zeros( shape );
INITIAL_SPACE( : ) = AIR;
half = ceil( 0.5 .* shape );
INITIAL_LIQUID = false( shape );
INITIAL_LIQUID( 8 : half( 1 ), 8 : half( 2 ) ) = LIQUID;
INITIAL_INTERFACE = imdilate( INITIAL_LIQUID, conndef( d, "maximal" ) );
INITIAL_INTERFACE = INITIAL_INTERFACE & ~INITIAL_LIQUID;
INITIAL_SPACE( INITIAL_LIQUID ) = LIQUID;
INITIAL_SPACE( INITIAL_INTERFACE ) = INTERFACE;
INITIAL_SPACE( INITIAL_SOLID ) = SOLID;

material = permute( INITIAL_SPACE, [ 3 1 2 ] );

rho( INITIAL_SPACE == SOLID ) = 1;
rho( INITIAL_SPACE == LIQUID ) = 1;
rho( INITIAL_SPACE == INTERFACE ) = 1;
rho( INITIAL_SPACE == AIR ) = 1;

mass( INITIAL_SPACE == LIQUID ) = 1;

fill( INITIAL_SPACE == LIQUID ) = 1;

bubble( INITIAL_SPACE == AIR ) = 1;
bubble( INITIAL_SPACE == INTERFACE ) = 1;

%% PLOT SETUP
fh = figure();
axh = axes( fh );
ih = imagesc( axh, squeeze( fill ).' );
colorbar( axh );
caxis( axh, [ 0 1 ] );
axis( axh, "equal", "off" );
% fh2 = figure();
% axh2 = axes( fh2 );

%% MAIN LOOP
total_mass = zeros( time_step_count + 1, 1 );
total_mass( 1 ) = sum( mass, "all" );
for cycle = 1 : time_step_count
    % calculate surface normals for interface cells
    % stream
    % collide
    % reconstruct missing DFs
    % calculate surface curvature
    % calculate mass transfer
    % determine convert from interface cells
    % determine convert to interface cells
    % determine distribute mass
    % convert and distribute
    
    %% mesh
    liquid = material == LIQUID;
    air = material == AIR;
    interface = material == INTERFACE;
    
    if any( imdilate( material == AIR, conndef( d+1, "maximal" ) ) & material == LIQUID, "all" )
        [];
    end
    
    %% streaming
    for i = 1 : q
        f_inter( i, :, : ) = circshift( f( i, :, : ), [ 0 cx( i ) cy( i ) ] );
    end
    
    %% obstacle-fluid full bounce-back
    for i = 1 : q
        n = circshift( material == SOLID, [ 0 cx( i ) cy( i ) ] );
        POI = ( liquid | interface ) & n;
        %POI_n = circshift( POI, -[ 0 cx( i ) cy( i ) ] );
        f_inter( opp( i ), POI ) = f( i, POI );
    end
    
    %% interface mass flow
    dm = zeros( aug_shape );
    for i = 2 : q
        n = circshift( liquid, [ 0 cx( i ) cy( i ) ] );
        POI = interface & n;
        POI_n = circshift( POI, -[ 0 cx( i ) cy( i ) ] );
        dm( POI ) = dm( POI ) ...
            + ( f( opp( i ), POI_n ) - f( i, POI ) ).';
        n = circshift( interface, [ 0 cx( i ) cy( i ) ] );
        POI = interface & n;
        POI_n = circshift( POI, -[ 0 cx( i ) cy( i ) ] );
        dm( POI ) = dm( POI ) ...
            + 0.5 .* ( fill( POI_n ) + fill( POI ) ) ...
            .* ( f( opp( i ), POI_n ) - f( i, POI ) ).';
    end
    mass = mass + dm;
%     histogram( axh2, mass( interface ) );
%     axh2.XLim = [ -0.1 1.1 ];
    total_mass( cycle + 1 ) = sum( mass, "all" );
    fill( interface ) = mass( interface ) ./ rho( interface );
    
    %% equilibrium
    ux = reshape( ( cx.' * reshape( f_inter, flat_shape ) ), aug_shape );
    uy = reshape( ( cy.' * reshape( f_inter, flat_shape ) ), aug_shape );
    [ f_equil, u2 ] = kernel( rho, ux, uy, w, cx, cy );
    
    %% gas-interface full bounce-back
    int = find( interface );
    f_eq_temp = kernel( rho_g, ux( int ), uy( int ), w, cx, cy );
    for i = 1 : q
        f_inter( opp( i ), int ) = ...
            + f_eq_temp( i, : ) ...
            + f_eq_temp( opp( i ), : ) ...
            - f_inter( i, int );
        % TODO this is more complicated than I thought
        % will take some thinking...
        %
        % need to compute kernel for all interface cells whose i-th neighbor is
        % empty
        % then assign to the empty cell i-th f_inter:
        %  + f_eq i-th
        %  + f_eq opp-th
        %  - f_inter opp-th
    end
    
    %% collision
    f = ( 1 - tau ) .* f_inter + tau .* ( f_equil + rho .* fill .* wdg );
    
    %% update macros
    rho = sum( f );
    u = sqrt( u2 );
    % bubble update
    
    %% convert from interface
    TOL = 1e-3;
    LOWER_FILL_BOUND = -TOL;
    int_to_air = find( fill < LOWER_FILL_BOUND & interface );
    liq_to_int = [];
    for i = 1 : numel( int_to_air )
        index = int_to_air( i );
        neigh = index + n_off;
        inter_neigh = neigh( interface( neigh ) );
        missing = mass( index );
        missing = missing ./ numel( inter_neigh );
        mass( inter_neigh ) = mass( inter_neigh ) + missing;
        mass( index ) = 0.0;
        fill( index ) = 0.0;
        material( index ) = AIR;
        liq_neigh = neigh( liquid( neigh ) );
        liq_to_int = [ liq_to_int; liq_neigh ];
    end
    
    UPPER_FILL_BOUND = 1 + TOL;
    int_to_liq = find( UPPER_FILL_BOUND < fill & interface );
    air_to_int = [];
    for i = 1 : numel( int_to_liq )
        index = int_to_liq( i );
        neigh = index + n_off;
        inter_neigh = neigh( interface( neigh ) );
        excess = mass( index );
        excess = excess ./ numel( inter_neigh );
        mass( inter_neigh ) = mass( inter_neigh ) + excess;
        mass( index ) = 1.0;
        fill( index ) = 1.0;
        material( index ) = LIQUID;
        air_neigh = neigh( air( neigh ) );
        air_to_int = [ air_to_int; air_neigh ];
    end
    
    %% convert to interface
    mean_rho = zeros( size( air_to_int ) );
    mean_ux = zeros( size( air_to_int ) );
    mean_uy = zeros( size( air_to_int ) );
    for i = 1 : numel( air_to_int )
        index = air_to_int( i );
        neigh = index + n_off;
        liq_inter_neigh = [ neigh( interface( neigh ) ); neigh( liquid( neigh ) ) ];
        mean_rho( i ) = mean( rho( liq_inter_neigh ) );
        mean_ux( i ) = mean( ux( liq_inter_neigh ) );
        mean_uy( i ) = mean( uy( liq_inter_neigh ) );
        mass( index ) = 0.0;
        fill( index ) = 0.0;
        material( index ) = INTERFACE;
    end
    f( :, air_to_int ) = kernel( mean_rho, mean_ux, mean_uy, w, cx, cy );
    
    for i = 1 : numel( liq_to_int )
        index = liq_to_int( i );
        mass( index ) = 1.0;
        fill( index ) = 1.0;
        material( index ) = INTERFACE;
    end
    
    
    if any( material == AIR & fill > 0.0, "all" ) || any( material == LIQUID & fill < 1.0, "all" )
        assert( false );
    end
    if any( imdilate( material == AIR, conndef( 3, "maximal" ) ) & material == LIQUID, "all" )
        assert( false );
    end
    
    fprintf( "%i %i %i %i\n", sum( int_to_liq, "all" ), sum( int_to_air, "all" ), sum( liq_to_int, "all" ), sum( liq_to_int, "all" ) );
    
    %% VISUALIZATION
    ih.CData = squeeze( interface ).';
    caxis( axh, [ 0 1 ] );
    drawnow( 'limitrate' )
end

end


function [ f_eq, u2 ] = kernel( rho, ux, uy, w, cx, cy )

if isempty( rho )
    f_eq = [];
    u2 = [];
    return;
end

u2 = ux .^ 2 + uy .^ 2;

if isvector( rho )
    cu = 3.0 .* ( cx .* ux.' + cy .* uy.' );
    f_eq = w .* ( rho.' + cu + 0.5 .* cu .^ 2 - 1.5 .* u2.' );
else
    cu = 3.0 .* ( cx .* ux + cy .* uy );
    f_eq = w .* ( rho + cu + 0.5 .* cu .^ 2 - 1.5 .* u2 );
end

end

