Mesh.MshFileVersion = 2;
///Parametrization:
//
// It is assumed that the XY-data arrays are flat and parallel to the surface at a given height and
// corresponding data are not changing vertically across a thin layer. 
//
// .... these are the horizontal coordinates of the lower-left (south-west) end of the data array in the mesh [m]:
//
DataRefX=-3305000.0;
DataRefY=-3305000.0;
SurfaceLevel=-2495.2058438264125;
MohoLevel=-20934.93987128108;
LABLevel=-116793.18875723529;
//
// ... this is the height of the grav and magnetic data above ground [m]
//
DataHeightAboveGround=10000.0;
DataThickness=10000.0;
//
// .... this spacing of the data array [m]:
//
DataSpacingX=10000.0;
DataSpacingY=10000.0;
//
// ... number of data points in east-west (X) and north-south (Y) direction:
//
DataNumX=661;
DataNumY=661;
//
// Note: the resolution specified here should roughly match the resolution of the actual data as input data are interpolated to the resolution in the mesh
//
// ... this is the "thickness" of the data array = the thickness of the vertical layer. 
//
DataMeshSizeVertical=10000.0;
//
// ... this is the thickness of region below the data area. We call this the core region which is the focus of the inversion
//
CoreThickness=400000.0;
//
// ... there is also an air layer and this is its thickness [m] (no updates for density and magnetization here)
//
AirLayerThickness=30000.0;
//
// ... there is padding around the core and air layer. For the subsurface there will be updates in padding region but not for the air layer  
//
PaddingX=100000.0;
PaddingY=100000.0;
PaddingZ=100000.0;
PaddingAir=100000.0;
//
// ... these are factors by which the DataMeshSizeVertical is raised in the air layer and in the core. 
//
MeshSizeAirFactor=1;
MeshSizeCoreFactor=1;
//
// ... these are factors by which the core and air layer mesh size are raised for the padding zone. 
//
MeshSizePaddingFactor=5;
// 
ProgressionType=1;
DataLX=DataSpacingX*(DataNumX-1);
DataLY=DataSpacingY*(DataNumY-1);
Depth=CoreThickness+PaddingZ;
Top=AirLayerThickness+PaddingAir;


MeshSizeAir=1*DataMeshSizeVertical;
MeshSizeCore=1*DataMeshSizeVertical;

MeshSizeBase=5*MeshSizeCore;
MeshSizePaddingAir=5*MeshSizeAir;
MeshSizePaddingSurface=5*MeshSizeCore;

DataMeshGround=DataMeshSizeVertical;




//surface in core region:
Point(100)={DataRefX,       DataRefY, SurfaceLevel, DataMeshGround};
Point(200)={DataRefX+DataLX,DataRefY, SurfaceLevel, DataMeshGround};
Point(300)={DataRefX,DataRefY+DataLY, SurfaceLevel, DataMeshGround};
Point(400)={DataRefX+DataLX,DataRefY+DataLY, SurfaceLevel, DataMeshGround};

//bottom data surface
Point(110)={DataRefX,       DataRefY, DataHeightAboveGround, DataMeshGround};
Point(210)={DataRefX+DataLX,DataRefY, DataHeightAboveGround, DataMeshGround};
Point(310)={DataRefX,DataRefY+DataLY, DataHeightAboveGround, DataMeshGround};
Point(410)={DataRefX+DataLX,DataRefY+DataLY, DataHeightAboveGround, DataMeshGround};
    
//top air layer
Point(120)={DataRefX,       DataRefY, AirLayerThickness, MeshSizeAir};
Point(220)={DataRefX+DataLX,DataRefY, AirLayerThickness, MeshSizeAir};
Point(320)={DataRefX,DataRefY+DataLY, AirLayerThickness, MeshSizeAir};
Point(420)={DataRefX+DataLX,DataRefY+DataLY, AirLayerThickness, MeshSizeAir};
        

//moho
Point(150)={DataRefX,       DataRefY, MohoLevel, MeshSizeCore};
Point(250)={DataRefX+DataLX,DataRefY, MohoLevel, MeshSizeCore};
Point(350)={DataRefX,DataRefY+DataLY, MohoLevel, MeshSizeCore};
Point(450)={DataRefX+DataLX,DataRefY+DataLY, MohoLevel, MeshSizeCore};

//bottom core
//Point(160)={DataRefX,       DataRefY, -CoreThickness, MeshSizeCore};
//Point(260)={DataRefX+DataLX,DataRefY, -CoreThickness, MeshSizeCore};
//Point(360)={DataRefX,DataRefY+DataLY, -CoreThickness, MeshSizeCore};
//Point(460)={DataRefX+DataLX,DataRefY+DataLY, -CoreThickness, MeshSizeCore};

//Padding:
Point(901)={DataRefX-PaddingX,DataRefY-PaddingY, -Depth, MeshSizePaddingSurface};
Point(902)={DataRefX-PaddingX,DataRefY-PaddingY, MohoLevel, MeshSizePaddingSurface};
Point(903)={DataRefX-PaddingX,DataRefY-PaddingY, SurfaceLevel, MeshSizePaddingSurface};
Point(904)={DataRefX-PaddingX,DataRefY-PaddingY, Top, MeshSizePaddingAir};

Point(911)={DataRefX+DataLX+PaddingX,DataRefY-PaddingY, -Depth, MeshSizePaddingSurface};
Point(912)={DataRefX+DataLX+PaddingX,DataRefY-PaddingY, MohoLevel, MeshSizePaddingSurface};
Point(913)={DataRefX+DataLX+PaddingX,DataRefY-PaddingY, SurfaceLevel, MeshSizePaddingSurface};
Point(914)={DataRefX+DataLX+PaddingX,DataRefY-PaddingY, Top, MeshSizePaddingAir};

Point(921)={DataRefX-PaddingX,DataRefY+DataLY+PaddingY, -Depth, MeshSizePaddingSurface};
Point(922)={DataRefX-PaddingX,DataRefY+DataLY+PaddingY, MohoLevel, MeshSizePaddingSurface};
Point(923)={DataRefX-PaddingX,DataRefY+DataLY+PaddingY, SurfaceLevel, MeshSizePaddingSurface};
Point(924)={DataRefX-PaddingX,DataRefY+DataLY+PaddingY, Top, MeshSizePaddingAir};

Point(931)={DataRefX+DataLX+PaddingX,DataRefY+DataLY+PaddingY, -Depth, MeshSizePaddingSurface};
Point(932)={DataRefX+DataLX+PaddingX,DataRefY+DataLY+PaddingY, MohoLevel, MeshSizePaddingSurface};
Point(933)={DataRefX+DataLX+PaddingX,DataRefY+DataLY+PaddingY, SurfaceLevel, MeshSizePaddingSurface};
Point(934)={DataRefX+DataLX+PaddingX,DataRefY+DataLY+PaddingY, Top, MeshSizePaddingAir};

Line(1) = {904, 924};
Line(2) = {924, 934};
Line(3) = {934, 914};
Line(4) = {914, 904};
Line(5) = {120, 320};
Line(6) = {320, 420};
Line(7) = {420, 220};
Line(8) = {220, 120};
Line(9) = {110, 310};
Line(10) = {310, 410};
Line(11) = {410, 210};
Line(12) = {210, 110};
Line(13) = {903, 923};
Line(14) = {923, 933};
Line(15) = {933, 913};
Line(16) = {913, 903};
Line(17) = {100, 300};
Line(18) = {300, 400};
Line(19) = {400, 200};
Line(20) = {200, 100};
Line(21) = {902, 922};
Line(22) = {922, 932};
Line(23) = {932, 912};
Line(24) = {912, 902};
Line(25) = {150, 350};
Line(26) = {350, 450};
Line(27) = {450, 250};
Line(28) = {250, 150};
//Line(29) = {160, 360};
//Line(30) = {360, 460};
//Line(31) = {460, 260};
//Line(32) = {260, 160};
Line(33) = {901, 921};
Line(34) = {921, 931};
Line(35) = {931, 911};
Line(36) = {901, 911};
Line(37) = {300, 310};
Line(38) = {100, 110};
Line(39) = {200, 210};
Line(40) = {400, 410};
Line(41) = {100, 150};
Line(42) = {300, 350};
Line(43) = {200, 250};
Line(44) = {400, 450};
//Line(45) = {150, 160};
//Line(46) = {350, 360};
//Line(47) = {250, 260};
//Line(48) = {450, 460};
Line(49) = {924, 923};
Line(50) = {904, 903};
Line(51) = {934, 933};
Line(52) = {914, 913};
Line(53) = {923, 922};
Line(54) = {903, 902};
Line(55) = {913, 912};
Line(56) = {933, 932};
Line(57) = {922, 921};
Line(58) = {902, 901};
Line(59) = {932, 931};
Line(60) = {912, 911};
Line Loop(1) = {20, 17, 18, 19};
Plane Surface(1) = {1};
Line Loop(2) = {12, 9, 10, 11};
Plane Surface(2) = {2};
Line Loop(3) = {27, 28, 25, 26};
Plane Surface(3) = {3};
//Line Loop(4) = {30, 31, 32, 29};
//Plane Surface(4) = {4};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(5) = {5};
Line Loop(6) = {34, 35, -36, 33};
Plane Surface(6) = {6};
Line Loop(7) = {15, 16, 13, 14};
Plane Surface(7) = {1, 7};
Line Loop(8) = {23, 24, 21, 22};
Plane Surface(8) = {3, 8};

dataVolume[]=Extrude {0,0,DataThickness}{Surface{2};Layers{1};};
//Physical Volume("DataArea")={dataVolume[1]};

//Submantel
Point(950)={DataRefX,       DataRefY, LABLevel, MeshSizeCore};
Point(960)={DataRefX,       DataRefY, LABLevel-(LABLevel+CoreThickness)/2, MeshSizeCore};
Point(970)={DataRefX,       DataRefY, -CoreThickness, MeshSizeCore};

Point(954)={DataRefX,       DataRefY+DataLY, LABLevel, MeshSizeCore};
Point(964)={DataRefX,       DataRefY+DataLY, LABLevel-(LABLevel+CoreThickness)/2, MeshSizeCore};
Point(974)={DataRefX,       DataRefY+DataLY, -CoreThickness, MeshSizeCore};

Point(946)={DataRefX+DataLX,       DataRefY, LABLevel, MeshSizeCore};
Point(956)={DataRefX+DataLX,       DataRefY, LABLevel-(LABLevel+CoreThickness)/2, MeshSizeCore};
Point(966)={DataRefX+DataLX,       DataRefY, -CoreThickness, MeshSizeCore};

Point(945)={DataRefX+DataLX,       DataRefY+DataLY, LABLevel, MeshSizeCore};
Point(955)={DataRefX+DataLX,       DataRefY+DataLY, LABLevel-(LABLevel+CoreThickness)/2, MeshSizeCore};
Point(965)={DataRefX+DataLX,       DataRefY+DataLY, -CoreThickness, MeshSizeCore};



Line(77) = {940, 320};
Line(78) = {936, 120};
Line(79) = {944, 420};
Line(80) = {935, 220};
Line Loop(9) = {49, 14, -51, -2};
Plane Surface(83) = {9};
Line Loop(10) = {53, 22, -56, -14};
Plane Surface(84) = {10};
Line Loop(11) = {22, 59, -34, -57};
Plane Surface(85) = {11};
Line Loop(12) = {77, 6, -79, -64};
Plane Surface(86) = {12};
Line Loop(13) = {8, 5, 6, 7};
Plane Surface(87) = {13};
Line Loop(14) = {37, 10, -40, -18};
Plane Surface(88) = {14};
Line Loop(15) = {18, 44, -26, -42};
Plane Surface(89) = {15};
//Line Loop(16) = {48, -30, -46, 26};
//Plane Surface(90) = {16};
Line Loop(17) = {77, -5, -78, 63};
Plane Surface(91) = {17};
Line Loop(18) = {38, 9, -37, -17};
Plane Surface(92) = {18};
Line Loop(19) = {41, 25, -42, -17};
Plane Surface(93) = {19};
//Line Loop(20) = {45, 29, -46, -25};
//Plane Surface(94) = {20};
Line Loop(21) = {50, 13, -49, -1};
Plane Surface(95) = {21};
Line Loop(22) = {54, 21, -53, -13};
Plane Surface(96) = {22};
Line Loop(23) = {58, 33, -57, -21};
Plane Surface(97) = {23};
Line Loop(24) = {79, 7, -80, -65};
Plane Surface(98) = {24};
Line Loop(25) = {39, -11, -40, 19};
Plane Surface(99) = {25};
Line Loop(26) = {43, -27, -44, 19};
Plane Surface(100) = {26};
//Line Loop(27) = {47, -31, -48, 27};
//Plane Surface(101) = {27};
Line Loop(28) = {51, 15, -52, -3};
Plane Surface(102) = {28};
Line Loop(29) = {15, 55, -23, -56};
Plane Surface(103) = {29};
Line Loop(30) = {59, 35, -60, -23};
Plane Surface(104) = {30};
Line Loop(31) = {8, -78, -62, 80};
Plane Surface(105) = {31};
Line Loop(32) = {38, -12, -39, 20};
Plane Surface(106) = {32};
Line Loop(33) = {20, 41, -28, -43};
Plane Surface(107) = {33};
//Line Loop(34) = {28, 45, -32, -47};
//Plane Surface(108) = {34};
Line Loop(35) = {4, 50, -16, -52};
Plane Surface(109) = {35};
Line Loop(36) = {16, 54, -24, -55};
Plane Surface(110) = {36};
Line Loop(37) = {24, 58, 36, -60};
Plane Surface(111) = {37};

//submantle
Line(113) = {945, 946};
Line(114) = {946, 950};
Line(115) = {950, 954};
Line(116) = {954, 945};

Line(135) = {955, 956};
Line(136) = {956, 960};
Line(137) = {960, 964};
Line(138) = {964, 955};

Line(157) = {965, 966};
Line(158) = {966, 970};
Line(159) = {970, 974};
Line(160) = {974, 965};


Line(119) = {250, 946};
Line(123) = {150, 950};
Line(127) = {350, 954};
Line(118) = {450, 945};

Line(145) = {950, 960};
Line(167) = {960, 970};
Line(149) = {954, 964};
Line(171) = {964, 974};
Line(141) = {946, 956};
Line(163) = {956, 966};
Line(140) = {945, 955};
Line(162) = {955, 965};


Line Loop(38) = {113, 114, 115, 116};
Plane Surface(133) = {38};

Line Loop(39) = {27, 119, -113, -118};
Plane Surface(120) = {39};
Line Loop(40) = {28, 123, -114, -119};
Plane Surface(124) = {40};
Line Loop(41) = {25, 127, -115, -123};
Plane Surface(128) = {41};
Line Loop(42) = {26, 118, -116, -127};
Plane Surface(132) = {42};

Line Loop(43) = {135, 136, 137, 138};
Plane Surface(155) = {43};

Line Loop(44) = {113,141,-135,-140};
Plane Surface(142) = {44};
Line Loop(45) = {114,145,-136,-141};
Plane Surface(146) = {45};
Line Loop(46) = {115,149,-137,-145};
Plane Surface(150) = {46};
Line Loop(47) = {116,140,-138,-149};
Plane Surface(154) = {47};


Line Loop(48) = {157,158,159,160};
Plane Surface(177) = {48};

Line Loop(49) = {135,163,-157,-162};
Plane Surface(164) = {49};
Line Loop(50) = {136,167,-158,-163};
Plane Surface(168) = {50};
Line Loop(51) = {137,171,-159,-167};
Plane Surface(172) = {51};
Line Loop(52) = {138,162,-160,-171};
Plane Surface(176) = {52};



Physical Surface("AllInterfaces") = {1, 2 , 3 , 5,6,7,8 , 83 , 84 , 85 , 86 , 87 , 88 , 89 , 91 , 92 , 93 , 95 , 96 , 97 , 98 , 99 , 100  , 102 , 103 , 104 , 105 , 106 , 107 , 109 , 110 , 111, 82, 2, 77, 69, 73, 81,133,120,124,128,132,155,142,146,150,154,177,164,168,172,176};

//Physical Surface("AllInterfaces") = {1, 2 , 3 , 4 , 5,6,7,8 , 83 , 84 , 85 , 86 , 87 , 88 , 89 , 90 , 91 , 92 , 93 , 94 , 95 , 96 , 97 , 98 , 99 , 100 , 101 , 102 , 103 , 104 , 105 , 106 , 107 , 108 , 109 , 110 , 111, 82, 2, 77, 69, 73, 81,133,120,124,};
// Physical Surface("AllFaces") = {1, 2 , 3 , 4 , 5 , 6 , 7 , 8 , 83 , 84 , 85 , 86 , 87 , 88 , 89 , 90 , 91 , 92 , 93 , 94 , 95 , 96 , 97 , 98 , 99 , 100 , 101 , 102 , 103 , 104 , 105 , 106 , 107 , 108 , 109 , 110 , 111, 82, 2, 77, 69, 73, 81};

