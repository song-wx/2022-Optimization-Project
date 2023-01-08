clear 
c1 = [3, 2, 4, 8, 0, 0];
A1 = [-2, 5, 3, -5, 1, 0;
      -1,-2,-5, -6, 0, 1 ];
b1 = [3;-8];
ind = [5,6];
[x1, z1, table1, flag1]=SimplexMin(A1,b1,c1,[],1);%%大M法解题1


c2 = [-3, 1, 3, -1];
A2 = [1, 2, -1, 1;
      1, -1, 2 -1;
      2, -2, 3, 3];
b2 = [0; 6; 9];
[x2, z2, table2, flag2]=SimplexMin(A2,b2,c2,[],2);%%两阶段法解题2


c3 = [2, -1, 0, 0, 0];
A3 = [-2, 1, 1, 1, 0;
      -1, 1, -1, 0, 1];
b3 = [-3; -2];
[x3, w, z3, table3, flag3]=DualSimplex(A3,b3,c3);%%对偶单纯形法解题3


xlswrite('result.xlsx', table1,1,'B2');
xlswrite('result.xlsx', table2,1,'B7');
xlswrite('result.xlsx', table3,1,'B13')
