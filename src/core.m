function core()
%% LATTICE CONSTANTS
% D2Q9
d = 2;
w = [ 16 4 4 4 4 1 1 1 1 ].' ./ 36;
cx = [ 0 1 0 -1 0 1 -1 -1 1 ].';
cy = [ 0 0 1 0 -1 1 1 -1 -1 ].';
opp = [ 1 4 5 2 3 8 9 6 7 ].';
same = ( 1 : 9 ).';
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
neighbor_order = [ 5 7 4 2 8 6 1 3 ];
n_off = n_off( neighbor_order );

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

dt = 0.25e-9;
dx = ell_0_p / max( shape );

nu_lbm = ( dt / dx ^ 2 ) * ( 1 / Re ); % LBM kinematic viscosity [ - ]
c_s_lbm = 1.0 / sqrt( 3.0 ); % LBM speed of sound, D2Q9, D3Q19 [ - ]
u_0 = ( dt / dx ) * u_0_p;
Ma = u_0 / c_s_lbm;
tau = nu_lbm / c_s_lbm^2 + 1.0 / 2.0; % LBM relaxation time [ - ]
assert( 0.5 < tau & tau < 5.0 );

g_p = 9.8; % gravitational acceleration [ m / s^2 ]
g_lbm = ( dt ^ 2 / dx ) * g_p; % LBM gravitational acceleration [ - ]
%sg_lbm = g_lbm ./ sqrt( 2 ); % edge diagonal
dg = -g_lbm .* cy; % gravitational force D2Q9
% dg = -[ g_lbm .* cz( 1 : 7 ); sg_lbm .* cz( 8 : end ) ]; % gravitational force D3Q19
wdg = w .* dg;

rho_g = tau;

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

SOLID_NEIGHBOR = bwmorph( ~INITIAL_SOLID, "remove" ); %bwmorph3
SOLID_NEIGHBOR = permute( SOLID_NEIGHBOR, [ 3 1 2 ] );

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
next_material = material;

% rho( INITIAL_SPACE == SOLID ) = 1;
% rho( INITIAL_SPACE == LIQUID ) = 1;
% rho( INITIAL_SPACE == INTERFACE ) = 1;
% rho( INITIAL_SPACE == AIR ) = 1;
rho( : ) = tau;

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
    %% mesh
    liquid = material == LIQUID;
    air = material == AIR;
    interface = material == INTERFACE;
    
    %% streaming
    for i = 1 : q
        f_inter( i, :, : ) = circshift( f( i, :, : ), [ 0 cx( i ) cy( i ) ] );
    end
    
    %% obstacle-fluid full bounce-back
    sn = ( interface | liquid ) & SOLID_NEIGHBOR;
    sn = find( sn );
    for i = 1 : numel( sn )
        index = sn( i );
        neigh = index + n_off;
        solid_neigh = neigh( INITIAL_SOLID( neigh ) );
        dirs = same( INITIAL_SOLID( neigh ) );
        dirs_o = opp( INITIAL_SOLID( neigh ) );
        s = ( solid_neigh - 1 ) .* q + dirs_o;
        f_inter( s ) = f( dirs, index );
    end
    
    %% equilibrium
    ux = reshape( ( cx.' * reshape( f_inter, flat_shape ) ), aug_shape );
    uy = reshape( ( cy.' * reshape( f_inter, flat_shape ) ), aug_shape );
    [ f_equil, u2 ] = kernel( rho, ux, uy, w, cx, cy );
    
    %% gas-interface full bounce-back
    int = find( interface );
    f_eq_temp = kernel( rho_g, ux( int ), uy( int ), w, cx, cy );
    for i = 1 : numel( int )
        index = int( i );
        neigh = index + n_off;
        gas_neigh = neigh( air( neigh ) );
        dirs = same( air( neigh ) );
        dirs_o = opp( air( neigh ) );
        gn = ( gas_neigh - 1 ) .* q + dirs_o;
        f_inter( gn ) = ...
            f_eq_temp( dirs, i ) ...
            + f_eq_temp( dirs_o, i ) ...
            - f_inter( dirs, index );
    end
    
    %% interface mass flow
    int = find( interface );
    dm = zeros( size( int ) );
    for i = 1 : numel( int )
        index = int( i );
        neigh = index + n_off;
        nn = interface( neigh );
        int_neigh = neigh( nn );
        dirs = same( nn );
        dirs_o = opp( nn );
        inds = ( int_neigh - 1 ) .* q + dirs_o;
        dm( i ) = 0.5 .* sum( ...
            ( fill( int_neigh ) + fill( index ) ) ...
            .* ( f( inds ) - f( dirs, index ) ) ...
            );
        nn = liquid( neigh );
        liq_neigh = neigh( nn );
        dirs = same( nn );
        dirs_o = opp( nn );
        inds = ( liq_neigh - 1 ) .* q + dirs_o;
        dm( i ) = dm( i ) + sum( f( inds ) - f( dirs, index ) );
    end
    %assert( ~any( abs( dm ) > 0.1 ), "%.4f", max( abs( dm ), [], "all" ) );
    mass( int ) = mass( int ) + dm;
    total_mass( cycle + 1 ) = sum( mass, "all" );
    fill( interface ) = mass( interface ) ./ rho( interface );
    
    %% collision
    f = ( 1 - tau ) .* f_inter + tau .* ( f_equil + rho .* fill .* wdg );
    
    %% update macros
    rho = sum( f );
    u = sqrt( u2 );
    % bubble update
    
    %% determine interface to air
    TOL = 1e-3;
    LOWER_FILL_BOUND = -TOL;
    int_to_air = find( fill < LOWER_FILL_BOUND & interface );
    next_material( int_to_air ) = AIR;
    
    %% determine interface to liquid
    UPPER_FILL_BOUND = 1 + TOL;
    int_to_liq = find( UPPER_FILL_BOUND < fill & interface );
    next_material( int_to_liq ) = LIQUID;
    
    %% determine liquid to interface
    liq_to_int = [];
    for i = 1 : numel( int_to_air )
        index = int_to_air( i );
        neigh = index + n_off;
        liq_neigh = neigh( next_material( neigh ) == LIQUID );
        liq_to_int = [ liq_to_int; liq_neigh ];
    end
    
    %% determine air to interface
    air_to_int = [];
    for i = 1 : numel( int_to_liq )
        index = int_to_liq( i );
        neigh = index + n_off;
        air_neigh = neigh( next_material( neigh ) == AIR );
        air_to_int = [ air_to_int; air_neigh ];
    end
    
    %% update interface
    fprintf( "%i %i %i %i\n", numel( int_to_liq ), numel( int_to_air ), numel( air_to_int ), numel( liq_to_int ) );
    new_interface = interface;
    new_interface( int_to_air ) = false;
    new_interface( int_to_liq ) = false;
    new_interface( air_to_int ) = true;
    new_interface( liq_to_int ) = true;
    
    %% compute missing mass
    missing = mass( int_to_air );
    
    %% compute excess mass
    excess = mass( int_to_liq ) - 1;
    
    %% convert
    material( int_to_air ) = AIR;
    mass( int_to_air ) = 0.0;
    fill( int_to_air ) = 0.0;
    
    material( int_to_liq ) = LIQUID;
    mass( int_to_liq ) = 1.0;
    fill( int_to_liq ) = 1.0;
    
    material( liq_to_int ) = INTERFACE;
    mass( liq_to_int ) = 1.0;
    fill( liq_to_int ) = 1.0;
    
    mean_rho = zeros( size( air_to_int ) );
    mean_ux = zeros( size( air_to_int ) );
    mean_uy = zeros( size( air_to_int ) );
    for i = 1 : numel( air_to_int )
        index = air_to_int( i );
        neigh = index + n_off;
        nn = interface( neigh ) | liquid( neigh );
        liq_inter_neigh = neigh( nn );
        mean_rho( i ) = mean( rho( liq_inter_neigh ) );
        mean_ux( i ) = mean( ux( liq_inter_neigh ) );
        mean_uy( i ) = mean( uy( liq_inter_neigh ) );
    end
    material( air_to_int ) = INTERFACE;
    f( :, air_to_int ) = kernel( mean_rho, mean_ux, mean_uy, w, cx, cy );
    mass( air_to_int ) = 0.0;
    fill( air_to_int ) = 0.0;
    
    if any( imdilate( material == AIR, conndef( 3, "maximal" ) ) & material == LIQUID, "all" )
        assert( false );
    end
    
    %% redistribute missing mass
    for i = 1 : numel( int_to_air )
        index = int_to_air( i );
        neigh = index + n_off;
        int_neigh = neigh( new_interface( neigh ) );
        mass( int_neigh ) = mass( int_neigh ) + missing( i ) ./ numel( int_neigh );
    end
    
    %% redistribute excess mass
    for i = 1 : numel( int_to_liq )
        index = int_to_liq( i );
        neigh = index + n_off;
        int_neigh = neigh( new_interface( neigh ) );
        mass( int_neigh ) = mass( int_neigh ) + excess( i ) ./ numel( int_neigh );
    end
    
    %% VISUALIZATION
    ih.CData = squeeze( new_interface ).';
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
    f_eq = w .* rho.' .* ( 1 + cu + 1.5 .* cu .^ 2 - 0.5 .* u2.' );
else
    cu = 3.0 .* ( cx .* ux + cy .* uy );
    f_eq = w .* rho .* ( 1 + cu + 1.5 .* cu .^ 2 - 0.5 .* u2 );
end

end

