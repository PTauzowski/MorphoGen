function prism = notchedPrism(ip1,ip2,xs2,r,al1,al2,height)
    sfL4= ShapeFunctionL4();
    co=CylinderObject(xs2, r, al1, al2, 0, height);
    pl1 = ShapeObjectRectangular(sfL4,[ip1 0; ip2 0; ip1 height; ip2 height]);
    prism = MorphSpace(pl1,co);
end

