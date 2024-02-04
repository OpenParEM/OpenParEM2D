// Control
//   build=WR90
//   check_limits=true
// EndControl
// RectangularWaveguide
//    name=WR90
//    width=0.02286
//    height=0.01016
//    material=air
//    default_conductor_material=PEC
// EndRectangularWaveguide
Point(1) = {-0.01143,-0.00508,0,0.00508};
Point(2) = {0.01143,-0.00508,0,0.00508};
Point(3) = {0.01143,0.00508,0,0.00508};
Point(4) = {-0.01143,0.00508,0,0.00508};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Surface("air",1) = {1};
