%  Matlab mesh
% MyGeometry_rect, Created by Gmsh
% ASCII
clear msh;
msh.nbNod = 66;
msh.POS = [
0 0 0;
0.5 0 0;
0.5 1 0;
0 1 0;
0.09999999999977846 0 0;
0.1999999999994866 0 0;
0.2999999999994725 0 0;
0.3999999999997362 0 0;
0.5 0.09999999999981414 0;
0.5 0.1999999999995569 0;
0.5 0.299999999999265 0;
0.5 0.3999999999989731 0;
0.5 0.4999999999986921 0;
0.5 0.599999999998945 0;
0.5 0.6999999999992088 0;
0.5 0.7999999999994725 0;
0.5 0.8999999999997362 0;
0.4 1 0;
0.3000000000006937 1 0;
0.2000000000008325 1 0;
0.1000000000004163 1 0;
0 0.8999999999995836 0;
0 0.8 0;
0 0.7000000000006938 0;
0 0.6000000000013874 0;
0 0.5000000000020595 0;
0 0.400000000001665 0;
0 0.3000000000012487 0;
0 0.2000000000008325 0;
0 0.1000000000004162 0;
0.09999999999984224 0.1000000000002958 0;
0.09999999999990604 0.2000000000005774 0;
0.09999999999996981 0.3000000000008519 0;
0.1000000000000336 0.4000000000011267 0;
0.1000000000000973 0.5000000000013861 0;
0.1000000000001612 0.6000000000008989 0;
0.1000000000002249 0.7000000000003967 0;
0.1000000000002887 0.7999999999998946 0;
0.1000000000003525 0.899999999999614 0;
0.1999999999996212 0.1000000000001754 0;
0.1999999999997558 0.2000000000003223 0;
0.1999999999998903 0.3000000000004553 0;
0.200000000000025 0.4000000000005882 0;
0.2000000000001596 0.5000000000007127 0;
0.2000000000002942 0.6000000000004104 0;
0.2000000000004287 0.7000000000000997 0;
0.2000000000005633 0.7999999999997891 0;
0.2000000000006979 0.8999999999996446 0;
0.2999999999995945 0.100000000000055 0;
0.2999999999997168 0.2000000000000672 0;
0.2999999999998389 0.3000000000000585 0;
0.2999999999999611 0.4000000000000499 0;
0.3000000000000832 0.5000000000000392 0;
0.3000000000002053 0.5999999999999218 0;
0.3000000000003274 0.6999999999998028 0;
0.3000000000004494 0.7999999999996834 0;
0.3000000000005715 0.8999999999996752 0;
0.3999999999997627 0.09999999999993456 0;
0.399999999999789 0.1999999999998121 0;
0.3999999999998154 0.2999999999996617 0;
0.3999999999998417 0.3999999999995115 0;
0.3999999999998681 0.4999999999993656 0;
0.3999999999998944 0.5999999999994335 0;
0.399999999999921 0.6999999999995059 0;
0.3999999999999473 0.7999999999995779 0;
0.3999999999999736 0.8999999999997058 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 5 1
 5 6 1
 6 7 1
 7 8 1
 8 2 1
 2 9 3
 9 10 3
 10 11 3
 11 12 3
 12 13 3
 13 14 3
 14 15 3
 15 16 3
 16 17 3
 17 3 3
 3 18 2
 18 19 2
 19 20 2
 20 21 2
 21 4 2
 4 22 3
 22 23 3
 23 24 3
 24 25 3
 25 26 3
 26 27 3
 27 28 3
 28 29 3
 29 30 3
 30 1 3
];
msh.QUADS =[
 1 5 31 30 4
 30 31 32 29 4
 29 32 33 28 4
 28 33 34 27 4
 27 34 35 26 4
 26 35 36 25 4
 25 36 37 24 4
 24 37 38 23 4
 23 38 39 22 4
 22 39 21 4 4
 5 6 40 31 4
 31 40 41 32 4
 32 41 42 33 4
 33 42 43 34 4
 34 43 44 35 4
 35 44 45 36 4
 36 45 46 37 4
 37 46 47 38 4
 38 47 48 39 4
 39 48 20 21 4
 6 7 49 40 4
 40 49 50 41 4
 41 50 51 42 4
 42 51 52 43 4
 43 52 53 44 4
 44 53 54 45 4
 45 54 55 46 4
 46 55 56 47 4
 47 56 57 48 4
 48 57 19 20 4
 7 8 58 49 4
 49 58 59 50 4
 50 59 60 51 4
 51 60 61 52 4
 52 61 62 53 4
 53 62 63 54 4
 54 63 64 55 4
 55 64 65 56 4
 56 65 66 57 4
 57 66 18 19 4
 8 2 9 58 4
 58 9 10 59 4
 59 10 11 60 4
 60 11 12 61 4
 61 12 13 62 4
 62 13 14 63 4
 63 14 15 64 4
 64 15 16 65 4
 65 16 17 66 4
 66 17 3 18 4
];
