% standard cube
x = 0.5;
cube_vertices = [-x, -x, x; 
            x -x x; 
            -x x x; 
            x x x; 
            -x x -x; 
            x x -x; 
            -x -x -x; 
            x -x -x];


cube_faces_quad = [1 2 4 3; 
            3 4 6 5; 
            5 6 8 7; 
            7 8 2 1; 
            2 8 6 4; 
            7 1 3 5];
        
cube_faces_tri = [1 2 3;
             3 2 4;
             3 4 5;
             5 4 6;
             5 6 7;
             7 6 8;
             7 8 1;
             1 8 2;
             2 8 4;
             4 8 6;
             7 1 5;
             5 1 3;
             9 10 11;
             11 10 12;
             11 12 13;
             13 12 14;
             13 14 15;
             15 14 16;
             15 16 9;
             9 16 10;
             10 16 12;
             12 16 14;
             15 9 13;
             13 9 11;
            ];
        

m = Mesh(vertices, cube_faces_tri);
write_obj(m, 'test.obj', 'w');

