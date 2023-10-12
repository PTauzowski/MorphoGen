classdef PlaneStressMaterial < Material

    properties
        E, nu, D, dD, M;
    end

    methods
        function obj = PlaneStressMaterial(name)
            obj = obj@Material(name);
        end

        function M = setMassIzoMatrix( m )
            obj.M = [m 0; 0 m];
            M = obj.M;
        end

        function D = setElasticIzo( obj, E, nu )
            D = E/(1-nu*nu)*[ 1  nu 0; ...
                                  nu 1  0;
                      0  0  (1-nu)/2 ];
            obj.D = D;
            obj.nu=nu;
            obj.E=E;

        end

        function dD = setElasticIzoGrad(obj)
            dD(:,:,1) = 1/(1-obj.nu^2)*[ 1  obj.nu 0; ...
                                  obj.nu 1  0;
                      0  0  (1-obj.nu)/2 ];
            dD(:,:,2) = 	[(2*obj.E*obj.nu)/(obj.nu^4-2*obj.nu^2+1) (obj.E*obj.nu^2+obj.E)/(obj.nu^4-2*obj.nu^2+1)	0; ...
	                        (obj.E*obj.nu^2+obj.E)/(obj.nu^4-2*obj.nu^2+1)	(2*obj.E*obj.nu)/(obj.nu^4-2*obj.nu^2+1)	0;...
	                        0	                        0	                        -obj.E/(2*obj.nu^2+4*obj.nu+2)];
            obj.dD = dD;
        end
                
    end
end

