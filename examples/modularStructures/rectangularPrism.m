function prism = rectangularPrism(ip1,ip2,ip3,ip4,height)
sfL4= ShapeFunctionL4();
pl1 = ShapeObjectRectangular(sfL4,[ip1 0; ip2 0; ip1 height; ip2 height]);
pl2 = ShapeObjectRectangular(sfL4,[ip3 0; ip4 0; ip3 height; ip4 height]);
prism = MorphSpace(pl1,pl2);
end

