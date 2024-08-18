classdef SolidMaterial < Material
   
    properties
        E, nu, D, dD, M, at, rho;
    end
    
    methods
       function obj = SolidMaterial(name)
            obj = obj@Material(name);
        end
       
        function D = setElasticIzo( obj, E, nu )
            D = E/(1.0+nu)/(1-2*nu)*[ 1-nu nu nu 0 0 0; ...
                                      nu 1-nu nu 0 0 0; ...
                                      nu nu 1-nu 0 0 0;  ...
                                      0 0 0 (1.0 - 2.0 * nu) / 2.0 0 0; ...
                                      0 0 0 0 (1.0 - 2.0 * nu) / 2.0 0; ...
                                      0 0 0 0 0 (1.0 - 2.0 * nu) / 2.0 ];
            obj.E=E;
            obj.nu=nu;
            obj.D=D;
        end

        function D = setElasticIzoThermal( obj, E, nu, at )
            D = setElasticIzo( obj, E, nu );
            obj.at=at;
        end

        function dD = setElasticIzoGrad(obj)
            E=obj.E;
            nu=obj.nu;
            dD(:,:,1) = 1/(1.0+nu)/(1-2*nu)*[ 1-nu nu nu 0 0 0; ...
                                      nu 1-nu nu 0 0 0; ...
                                      nu nu 1-nu 0 0 0;  ...
                                      0 0 0 (1.0 - 2.0 * nu) / 2.0 0 0; ...
                                      0 0 0 0 (1.0 - 2.0 * nu) / 2.0 0; ...
                                      0 0 0 0 0 (1.0 - 2.0 * nu) / 2.0 ];
            
            dD(:,:,2) = [-(2*E*nu^2-4*E*nu)/(4*nu^4+4*nu^3-3*nu^2-2*nu+1)	(2*E*nu^2+E)/(4*nu^4+4*nu^3-3*nu^2-2*nu+1)	(2*E*nu^2+E)/(4*nu^4+4*nu^3-3*nu^2-2*nu+1)	0	0	0; ...
		                (2*E*nu^2+E)/(4*nu^4+4*nu^3-3*nu^2-2*obj.nu+1)	-(2*E*nu^2-4*E*nu)/(4*nu^4+4*nu^3-3*nu^2-2*nu+1)	(2*E*nu^2+E)/(4*nu^4+4*nu^3-3*nu^2-2*nu+1) 	0	0	0; ...
		                (2*E*nu^2+E)/(4*nu^4+4*nu^3-3*nu^2-2*obj.nu+1)	(2*E*nu^2+E)/(4*nu^4+4*nu^3-3*nu^2-2*nu+1)	-(2*E*nu^2-4*E*nu)/(4*nu^4+4*nu^3-3*nu^2-2*nu+1)	0	0	0; ...
		                 0	0	0	-E/(2*nu^2+4*nu+2)	0	0; ...
		                 0	0	0	0	-E/(2*nu^2+4*nu+2)	0; ...
		                 0	0	0	0	0	-E/(2*nu^2+4*nu+2)];
            obj.dD = dD;
        end

        function [ M ] = setMassIzo( m )
            M = [m 0 0; 0 m 0; 0 0 m];
        end

    end
end

