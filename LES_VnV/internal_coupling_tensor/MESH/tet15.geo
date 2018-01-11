Point(1) = {0, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(7) = {0, 1, 1, 1.0};
Point(8) = {0, 0, 1, 1.0};

Line(1) = {3, 7};
Line(5) = {3, 1};
Line(9) = {8, 1};
Line(11) = {8, 7};

Line Loop(23) = {9, -5, 1, -11};
Plane Surface(24) = {23};

surfaceVector[] = Extrude {1, 0, 0} {
  Surface{24};
};

surfaceVector2[] = Extrude {-1, 0, 0} {
  Surface{24};
};

Transfinite Line "*" = 15;
Transfinite Surface "*";
Transfinite Volume "*";
