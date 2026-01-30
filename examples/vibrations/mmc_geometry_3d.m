function [Phi, PhiDrv, xval, actComp, actDsvb] = calc_Phi_3d(Phi, PhiDrv, xval, i, LSgrid, p, nEhcp, epsilon, actComp, actDsvb, minSz)
%CALC_PHI_3D Compute TDF and derivatives for one MMC component (3D)

    di = xval((i-1)*nEhcp+1 : i*nEhcp);
    x0 = di(1); y0 = di(2); z0 = di(3); l1 = di(4) + eps; l2 = di(5) + eps; l3 = di(6) + eps;
    sa = sin(di(7)); sb = sin(di(8)); sg = sin(di(9));
    ca = cos(di(7)); cb = cos(di(8)); cg = cos(di(9));

    R = [cb*cg, cb*sg, -sb; ...
         sa*sb*cg - ca*sg, sa*sb*sg + ca*cg, sa*cb; ...
         ca*sb*cg + sa*sg, ca*sb*sg - sa*cg, ca*cb];

    xyzLc = [LSgrid.x(:)-x0, LSgrid.y(:)-y0, LSgrid.z(:)-z0];
    xyz   = xyzLc * R';
    x1 = xyz(:,1) + eps; y1 = xyz(:,2) + eps; z1 = xyz(:,3) + eps;

    temp = (abs(x1./l1).^p + abs(y1./l2).^p + abs(z1./l3).^p);
    Phi(:, i) = 1 - temp.^(1/p);

    % deletion tests
    if (l1/minSz < 1.01 && l2/minSz < 1.01) || min(abs(Phi(:, i))) >= epsilon
        Phi(:, i) = -1e3;
        xval((i-1)*nEhcp+[4:6]) = 0;
        actComp = setdiff(actComp, i);
        actDsvb = setdiff(actDsvb, nEhcp*i - nEhcp + 1 : nEhcp*i);
        return;
    end

    Ra = [0 0 0; R(3,1) R(3,2) R(3,3); -R(2,1) -R(2,2) -R(2,3)];
    Rb = [-sb*cg -sb*sg -cb; sa*cb*cg sa*cb*sg -sa*sb; ca*cb*cg ca*cb*sg -ca*sb];
    Rg = [-cb*sg cb*cg 0; -sa*sb*sg - ca*cg sa*sb*cg - sa*sg 0; -ca*sb*sg + sa*cg ca*sb*cg + sa*sg 0];

    dxi = [-R(1,:) zeros(1,6)] + xyzLc * [Ra(1,:); Rb(1,:); Rg(1,:)];
    dyi = [-R(2,:) zeros(1,6)] + xyzLc * [Ra(2,:); Rb(2,:); Rg(2,:)];
    dzi = [-R(3,:) zeros(1,6)] + xyzLc * [Ra(3,:); Rb(3,:); Rg(3,:)];

    temp1 = -temp.^(1/p - 1) .* (x1./l1).^(p-1) / l1;
    temp2 = -temp.^(1/p - 1) .* (y1./l2).^(p-1) / l2;
    temp3 = -temp.^(1/p - 1) .* (z1./l3).^(p-1) / l3;

    dpdl1 = -temp1 .* (x1./l1);
    dpdl2 = -temp2 .* (y1./l2);
    dpdl3 = -temp3 .* (z1./l3);

    dpdL = zeros(size(x1,1), 6); dpdL(:,4:6) = [dpdl1 dpdl2 dpdl3];
    PhiDrv(:, nEhcp*(i-1)+1 : nEhcp*i) = dpdL + temp1 .* dxi + temp2 .* dyi + temp3 .* dzi;
end

