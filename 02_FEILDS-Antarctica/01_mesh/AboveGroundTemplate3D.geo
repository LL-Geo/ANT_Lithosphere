Mesh.MshFileVersion = 2;
Merge "{MESHFILE}";
Mesh.Algorithm3D = 10;

//dataVolume[]=Extrude {{0,0,DataThickness}}{{Surface{{2}};Layers{{1}};}};
//Physical Volume("DataArea")={{dataVolume[1]}};


Surface Loop(1) = {{87, 105, 91, 86, 98, 82}};
Volume(2) = {{1}};
Physical Volume("Air") = {{2}};

Surface Loop(2) = {{2, 99, 106, 92, 88, 1}};
Volume(4) = {{2}};
Physical Volume("AboveTopo") = {{4}};

Surface Loop(3) = {{1, 100, 107, 93, 89, 3}};
Volume(5) = {{3}};
Physical Volume("Crust") = {{5}};

Surface Loop(4) = {{3, 133, 120, 124, 128, 132}};
Volume(6) = {{4}};
Physical Volume("Mantel1") = {{6}};

Surface Loop(5) = {{133, 155, 142, 146, 150, 154}};
Volume(7) = {{5}};
Physical Volume("Mantel2") = {{7}};

Surface Loop(6) = {{155, 177, 164, 168, 172, 176}};
Volume(8) = {{6}};
Physical Volume("Mantel3") = {{8}};



Surface Loop(7) = {{7, 102, 83, 95, 109, 5, 98, 86, 91, 87, 105, 99, 106, 92, 88, 69, 73, 77, 81}};
Volume(9) = {{7}};
Physical Volume("AirPadding") = {{9}};

Surface Loop(8) = {{96, 110, 103, 84, 7, 8, 93, 107, 100, 89}};
Volume(10) = {{8}};
Physical Volume("CrustPadding") = {{10}};

Surface Loop(9) = {{97, 111, 6, 85, 104, 177, 176, 154, 132, 164, 142, 120, 168, 146, 124, 172, 150, 128, 8}};
Volume(11) = {{9}};
Physical Volume("MantelPadding") = {{11}};

Surface Loop(10) = {{82, 2, 69, 73, 77, 81}};
Volume(12) = {{10}};
Physical Volume("DataArea") = {{12}};

Physical Surface("EarthSurface") = {{1, 7}};
Physical Surface("SurfaceTop") = {{5}};
Physical Surface("SurfaceBottom") = {{6}};




//Physical Surface("EarthSurface") = {{1, 7}};
