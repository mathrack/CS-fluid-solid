lx = 6.283;
lz = 3.142;
nx = 32.0;
nz = 32.0;

lc = 1e-7;
dx = lx / nx;
dz = lz / nz;

ypos = 0;
dy0 = 1 / 150;
raison = 1.09;

p1 = newp;
Point(p1) = {0, -1, 0, lc};
autoPoint[] = Translate {dx, 0, 0} { Duplicata{ Point{p1}; } };
p2 = autoPoint[0];
autoPoint[] = Translate {0, 0, dz} { Duplicata{ Point{p1}; } };
p3 = autoPoint[0];
autoPoint[] = Translate {dx, 0, dz} { Duplicata{ Point{p1}; } };
p4 = autoPoint[0];
l1 = newl;
Line(l1) = {p1, p2};
l2 = newl;
Line(l2) = {p1, p3};
l3 = newl;
Line(l3) = {p4, p2};
l4 = newl;
Line(l4) = {p4, p3};
Transfinite Line {l1, l2, l3, l4} = 2 Using Progression 1;
l5 = newll;
Line Loop(l5) = {l1, -l3, l4, -l2} ;
s6 = news;
Plane Surface(s6) = {l5} ;
Transfinite Surface{s6};
Recombine Surface {s6};

dy = dy0;
For i In {0:29}

  autoVol[] = Extrude {0, dy, 0} { Surface{s6}; Layers{1}; Recombine; };

  ypos = ypos + dy;

  s6 = autoVol[5];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  s6 = autoVol[4];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  s6 = autoVol[3];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  s6 = autoVol[2];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  s6 = autoVol[0];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  v6 = autoVol[1];
  Transfinite Volume{v6};

  dy *= raison;

EndFor


dy = (1-ypos);

  autoVol[] = Extrude {0, dy, 0} { Surface{s6}; Layers{1}; Recombine; };

  ypos = ypos + dy;

  s6 = autoVol[5];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  s6 = autoVol[4];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  s6 = autoVol[3];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  s6 = autoVol[2];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  s6 = autoVol[0];
  Transfinite Surface{s6};
  Recombine Surface {s6};
  v6 = autoVol[1];
  Transfinite Volume{v6};