function [ nn ] = closestNode( nodes, pt )
     norms    = sqrt( sum( (nodes-pt).^2, 2 ) );
    [dist,nn] = min( norms ); 
end