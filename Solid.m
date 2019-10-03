classdef Solid < Voxels
    
    properties
        value = 1
    end
    
    properties ( Dependent )
        geometric_indices
        phase_indices
    end
    
    methods
        function obj = Solid( phase_shape )
            assert( isa( phase_shape, "double" ) )
            assert( ...
                numel( phase_shape ) == 3 ...
                || numel( phase_shape ) == 4 ...
                );
            
            shape = phase_shape( 2 : end );
            stride = [ 1 prod( shape( 2 : end ) ) ];
            phase_stride = [ 1 prod( shape ) ];
            
            obj.shape = phase_shape( 2 : end );
            obj.stride = stride;
            obj.phase_shape = phase_shape;
            obj.phase_stride = phase_stride;
        end
        
        % function add_voxels( voxels )
        % end
    end
    
    properties ( Access = private )
        shape
        phase_shape
        stride
        phase_stride
    end
    
end

