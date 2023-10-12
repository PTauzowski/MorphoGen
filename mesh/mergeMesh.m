function [ nodes, elems ] = mergeMesh( nodes1, elems1, nodes2, elems2 )

    [cnodes,i1,i2] = intersect( nodes1, nodes2, 'rows' );
    nidx  = 1:size(nodes2,1);
    nidx( i2 ) = [];
    noi   = 1:size(nidx,2);
    noi = noi + size(nodes1,1);
    nodes = [ nodes1; nodes2(nidx,:) ];
    ninds = zeros( size(nodes2,1), 1);
    ninds(i2)=i1;
    ninds(nidx)=noi;
    el2 = elems2;
    el2(:) = ninds( elems2(:) );
    elems = [ elems1; el2 ];
    
end