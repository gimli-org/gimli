// Embedded geological body in a 2D crosshole ERT example
// 15.11.2012 fwagner@gfz-potsdam.de

// Mesh sizes
cl1 = 800; // outer box
cl2 = 0.5; // electrodes
cl3 = 1; // body

// Definition of bounding box
Point(1) = {-5000, 0, 0, cl1};
Point(2) = {-5000, 0, -5000, cl1};
Point(3) = {5000, 0, -5000, cl1};
Point(4) = {5000, 0, 0, cl1};

// Definition of inversion region
Point(5) = {0, 0, 0, cl2};
Point(7) = {0, 0, -50, cl2};
Point(8) = {50, 0, -50, cl2};
Point(6) = {50, 0, 0, cl2};

// Connecting points
Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 4};
Line(4) = {4, 3};
Line(5) = {3, 2};
Line(6) = {2, 1};
Line(7) = {5, 7};
Line(8) = {7, 8};
Line(9) = {8, 6};

// Definition of electrodes
For i In {0:9}
Point(newp) = {15,0,-4*(i+2),cl2}; // Borehole 1
Point(newp) = {35,0,-4*(i+2),cl2}; // Borehole 2
EndFor

// Definition of geological body
Point(29) = {40.2, 0, -26.6, cl3};
Point(30) = {42.1, 0, -25.6, cl3};
Point(31) = {41.7, 0, -23.3, cl3};
Point(32) = {35.2, 0, -22.2, cl3};
Point(33) = {27.7, 0, -20.9, cl3};
Point(34) = {20.9, 0, -23.2, cl3};
Point(35) = {17.8, 0, -25.9, cl3};
Point(36) = {15, 0, -26.2, cl3};
Point(37) = {9.4, 0, -27.6, cl3};
Point(38) = {7.2, 0, -30.1, cl3};
Point(39) = {8.5, 0, -33.2, cl3};
Point(40) = {11.9, 0, -34.2, cl3};
Point(41) = {14.9, 0, -34.2, cl3};
Point(42) = {18, 0, -33.9, cl3};
Point(43) = {23.6, 0, -30.8, cl3};
Point(44) = {29.1, 0, -26.7, cl3};
Point(45) = {35.3, 0, -26.3, cl3};

Spline(100) = {31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 29, 30, 31};

// Definition of regions
Line Loop(101) = {3, 4, 5, 6, 1, 7, 8, 9};
Plane Surface(102) = {101}; // Outer box
Line Loop(103) = {100};
Plane Surface(104) = {103}; // Geological body
Line Loop(105) = {2, -9, -8, -7};
Plane Surface(106) = {103, 105}; // Inner box without geological body

// Embedding electrodes in surfaces
Point{10, 12, 14, 16, 20, 22, 24, 26, 28, 9, 11, 13, 15, 17, 23, 25, 27} In Surface{106};
Point{18, 19, 21} In Surface{104}; // Three electrodes lie within the target

// Definition of boundaries
Physical Line(1) = {3, 2, 1}; // Free surface
Physical Line(2) = {4, 5, 6}; // Mixed boundary conditions

// Definition of physical regions
Physical Surface(1) = {102}; // No-inversion region
Physical Surface(2) = {106}; // Inversion region
Physical Surface(3) = {104}; // Geological body

// Definition of electrodes
Physical Point(99) = {9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28}; // Setting electrode marker (99)
