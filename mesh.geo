size = 2; // normal size of mesh elements
refined_size = 0.2; // refined size of mesh elements
r = 3;  // radius of the circle in center

// create points for outer rectangle
Point(1) = { 0, 0, 0, refined_size};
Point(2) = { 40, 0, 0, size};
Point(3) = { 40, 20, 0, size};
Point(4) = { 0, 20, 0, refined_size};

// create poitns for circle
Point(5) = { 20, 10, 0};
Point(6) = { 20+r, 10, 0, refined_size};
Point(7) = { 20-r, 10, 0, refined_size};

// create outer rectangle lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Create circle
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 6};

// Create curves to bound circle and rectangle
Curve Loop(1) = {1, 2, 3, 4}; // rectangle
Curve Loop(2) = {5,6};        // circle

// create the silicone surface
Plane Surface(1) = {1, 2};
Plane Surface(2) = {2};

Physical Line("left border") = {4};
Physical Surface("silicone") = {1};

Mesh 2;

Mesh.SurfaceFaces = 1;
Mesh.Points = 1;

Save "circle_in_rectangle.msh";
