use <commonLampModules.scad>

// ====== set the parameter values ======
$fn = 20;

radii = [0.1, 0.1, 0.1, 0.1];
stand_h = 0.25;
stand_r = 1;
shade_end_radius = 0.4;
shade_height = 0.8;
shade_thickness=0.1;

// vectors defined wrt the ancestor endpoint
base = [1, 1, 2];
a1 = [1, 1, 2];
a2 = [0, 1, 2];
a3 = [1, 0, 2];
h1 = [1, 1, -1];
h2 = [-0.4, 1, -1];
h3 = [1, -0.4, -1];

createLampGeometry(base, a1, a2, a3, h1, h2, h3,
            radii, stand_h, stand_r,
            shade_end_radius, shade_height, shade_thickness 
            );

echo("done creating geometry");

