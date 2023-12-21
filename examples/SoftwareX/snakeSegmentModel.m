function mesh = snakeSegmentModel(width,length, height, r, R,localNodes,resC)
    xs1=[width/2 0];
    xs2=[width/2 length];

    cutAngle=deg2rad(30);

    ip1=[xs2(1)-3*R*cos(cutAngle)/4 length/2];
    ip2=[xs2(1)+3*R*cos(cutAngle)/4 length/2];
    ip3=[xs2(1)-r*cos(cutAngle) xs2(2)-r*sin(cutAngle) ];
    ip4=[xs2(1)+r*cos(cutAngle) xs2(2)-r*sin(cutAngle) ];
    ip5=[xs2(1)-R*cos(cutAngle) length/2 ];
    ip6=[xs2(1)+R*cos(cutAngle) length/2 ];
    ip7=[xs2(1)-R*cos(cutAngle) xs2(2)-R*sin(cutAngle) ];
    ip8=[xs2(1)+R*cos(cutAngle) xs2(2)-R*sin(cutAngle) ];
    ip9=[0 xs1(2)+width/2*tan(cutAngle) ];
    ip10=[width xs1(2)+width/2*tan(cutAngle) ];
    ip11=[0 xs2(2)-width/2*tan(cutAngle) ];
    ip12=[width xs2(2)-width/2*tan(cutAngle) ];
    ip13=[xs1(1)-R*cos(cutAngle) xs1(2)+R*sin(cutAngle) ];
    ip14=[xs1(1)+R*cos(cutAngle) xs1(2)+R*sin(cutAngle) ];

    resY=resC*16;
    resX=resC*3;
    resZ=resC*5;
    resPipe2=resC*30;
    
    mesh=Mesh();
    mp=notchedPrism(ip1,ip2,xs2,r,210,330,height/2);
    mesh.addObjectMesh3D( mp, resY, 6*resX, resZ, localNodes );
    mp=rectangularPrism(ip1,ip3,ip5,ip7,height/2);
    mesh.addObjectMesh3D( mp, resX, resY, resZ, localNodes );
    mp=rectangularPrism(ip2,ip4,ip6,ip8,height/2);
    mesh.addObjectMesh3D( mp, resX, resY, resZ, localNodes );
    mp=notchedPrism(ip6,ip5,xs1,R,30,150,height/2);
    mesh.addObjectMesh3D( mp, resY, 8*resX, resZ, localNodes );
    mp=rectangularPrism(ip8,ip14,ip12,ip10,height/2);
    mesh.addObjectMesh3D( mp,2*resX, 2*resY, resZ, localNodes );
    mp=rectangularPrism(ip13,ip7,ip9,ip11,height/2);
    mesh.addObjectMesh3D( mp, 2*resX, 2*resY, resZ, localNodes )

    mesh1=Mesh();
    mesh1.mergeMesh(mesh);
    mesh1.transformMesh3DDegXY( xs1, 0, [0 0 height/2] );
    mesh.addPipe3D(xs2,r,R,-30,210,height/2,height,resX,resPipe2,resZ,localNodes);
    mesh1.mergeMesh(mesh);
    mesh1.addPipe3D(xs1,r,R, 30,150,0,height/2,10,8*resX,resZ,localNodes);
    %mesh1.addCylinder( [xs1 0], r, height, [18*resX 2*resZ], localNodes );
    mesh=mesh1;
%mesh.addCylinder( [0 0 0] , 3, 5, [20 20 40], sfL8.localNodes );
end

