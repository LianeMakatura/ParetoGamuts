// ========== function to create the lamp geometry ===========
// Create the lamp
module createLampGeometry(
            base, a1, a2, a3, h1, h2, h3,
            radii, stand_h, stand_r,
            shade_end_radius, shade_height, shade_thickness
        ){

    //union(){
        bottom_stand(stand_r, stand_h);

        rod([0,0,0], base, radii[0]); //base
        rod(base, base+a1, radii[1]); //arm 1
        rod(base+a1, base+a1+h1, radii[1]); //hand 1
        rod(base, base+a2, radii[2]); //arm 2
        rod(base+a2, base+a2+h2, radii[2]); //hand 2
        rod(base, base+a3, radii[3]); //arm 3
        rod(base+a3, base+a3+h3, radii[3]); //hand 3

        lampshade_full(base+a1+h1, h1, radii[1], 
                    shade_end_radius, shade_height, shade_thickness);
        lampshade_full(base+a2+h2, h2, radii[2], 
                    shade_end_radius, shade_height, shade_thickness);
        lampshade_full(base+a3+h3, h3, radii[3], 
                    shade_end_radius, shade_height, shade_thickness);
    //}
}





// ============ functions to create lamp components ============
// define a rod between two points
module rod(a, b, r) { 
    translate(a) sphere(r=r); 
    translate(b) sphere(r=r); 

    dir = b-a; 
    h   = norm(dir); 
    if(dir[0] == 0 && dir[1] == 0) { 
        // no transformation necessary 
        cylinder(r=r, h=h); 
    } 
    else { 
        w  = dir / h; 
        u0 = cross(w, [0,0,1]); 
        u  = u0 / norm(u0); 
        v0 = cross(w, u); 
        v  = v0 / norm(v0); 
        multmatrix(m=[[u[0], v[0], w[0], a[0]], 
                      [u[1], v[1], w[1], a[1]], 
                      [u[2], v[2], w[2], a[2]], 
                      [0,    0,    0,    1]]) 
        cylinder(r=r, h=h); 
    } 
} 

module lampshade_full(origin, dir, connection_pt_radius, end_radius, 
                        height, thickness) {
    difference() {
        lampshade(origin, dir, height, connection_pt_radius, end_radius);
        lampshade(origin, dir, height*1.2,      
                    connection_pt_radius*(1-thickness), 
                    end_radius*(1-thickness));
    }
    // add the bulb
    orientation = dir / norm(dir);
    shade_end_local = height*orientation;
//        union() {
            translate(origin + shade_end_local) 
                sphere(r=0.2);
            rod(origin+0.5*shade_end_local, 
                origin + 0.9*shade_end_local,
                0.15);
 //       }
}

module lampshade(origin, dir, height, conn_pt_r, end_r) {
    h = height;
    r1 = conn_pt_r;
    r2 = end_r;

    if(dir[0] == 0 && dir[1] == 0) { 
        // no transformation necessary 
        cylinder(h=h, r1=r1, r2=r2); 
    } 
    else { 
        w  = (dir) / norm(dir); 
        u0 = cross(w, [0,0,1]); 
        u  = u0 / norm(u0); 
        v0 = cross(w, u); 
        v  = v0 / norm(v0); 
        multmatrix(m=[[u[0], v[0], w[0], origin[0]], 
                      [u[1], v[1], w[1], origin[1]], 
                      [u[2], v[2], w[2], origin[2]], 
                      [0,    0,    0,    1]]) 
        cylinder(h=h, r1=r1, r2=r2); 
    } 
}

module bottom_stand(stand_r, stand_h) {
    translate([0,0,-stand_h])
        cylinder(h=stand_h, r=stand_r);
}