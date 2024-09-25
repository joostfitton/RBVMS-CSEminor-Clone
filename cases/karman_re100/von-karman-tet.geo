//
r = 0.5;
x0 = -8;
x1 = 32.0;

y0 = -8.0;
y1 =  8.0;

// Create outbox
Point(1) = {x0, y0, 0.0, 0.1};
Point(2) = {x1, y0, 0.0, 0.1};
Point(3) = {x1, y1, 0.0, 0.1};
Point(4) = {x0, y1, 0.0, 0.1};

Line(21) = {1,2};
Line(22) = {2,3};
Line(23) = {3,4};
Line(24) = {4,1};

Line Loop(40) = {21,22,23,24};

// Create cylinder
Point(10) = {0.0, 0.0, 0.0, 1.0};

Point(11) = {0.0,  -r, 0.0, 0.05};
Point(12) = {  r, 0.0, 0.0, 0.05};
Point(13) = {0.0,   r, 0.0, 0.05};
Point(14) = { -r, 0.0, 0.0, 0.05};

Circle(31) = {11,10,12};
Circle(32) = {12,10,13};
Circle(33) = {13,10,14};
Circle(34) = {14,10,11};

Line Loop(50) = {31,32,33,34};

// Domain
Surface(100) = {40,50};

Physical Surface(100) = {100};
Physical Line(1) = {24};
Physical Line(2) = {21, 23};
Physical Line(3) = {31,32,33,34};
Physical Line(4) = {22};

Recombine Surface {100};
//Mesh.Algorithm = 8;
