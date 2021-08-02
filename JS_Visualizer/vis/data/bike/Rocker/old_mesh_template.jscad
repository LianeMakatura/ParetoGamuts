/*
Parametric Bike Stay & Rocker 
*/

// ==============================================
//              Define parameters
// ==============================================

var eps = 1e-3;

// ------ Design parameters
//here dp will be the design point array and ap the application point
dp = [1,0,0]; // these are just for testing, they'll be overwritten below
ap = [1];
g_fn = 5; //default resolution, also overwritten 

//INSERT_PARAMETER_VALUES


// -------- shared params & constants
pivotoff_y = lerp(ap[0], -9, 40);

if (ap.length < 2) {
    pivotoff_x = 0; // no  y offset if only 1 application variable
}
else {
    pivotoff_x = lerp(ap[1], 0, -40);
}


pivot = [0,0];
bolt_radius = 3;
bolt_shaft_rad = 7;

// ------- bike rocker params & constants
solid_offset = lerp(dp[0], 25, 45);
thickness = lerp(dp[1], 7, 14);
rocker_depth = lerp(dp[2], 5, 10);

thick_top = thickness;
thick_left = thickness;
thick_right = thickness;
thick_inner = 0.8*thickness;

// these are *always* fixed in same spot; they're correctly relative to the original pivot
rocker_r = [pivot[0] + 100, pivot[1]]; // right bolt hole
rocker_b = [2/3*(rocker_r[0]-pivot[0]), pivot[1] - 50]; // bottom bolt hole
solid_angle = 30; // not used currently


// ------- bike stay params
drop_h = lerp(dp[0], 50, 80);
tube_h = lerp(dp[1], 10, 20);
stay_depth = lerp(dp[2], 6, 12);

// these are *always* fixed in same spot; they're correctly relative to the original pivot
stay_b = [pivot[0] - 180, pivot[1] - 175];
drop_w = 30; // furthest bottom right corner







// ==============================================
//              Define mesh
// ==============================================

function main() {
    const pivotFinal = [pivot[0] + pivotoff_x, pivot[1] + pivotoff_y];

    var rocker = rockerFrame(pivotFinal);
    rocker = rocker.translate([0,0,stay_depth]); // offset to fit with stay
    //const stay = swingArm(pivotFinal);
    
    return [rocker];
}


function swingArm(pivotFinal) {
    var swingArm;
    
    const d = basic_dropout();
    const a = arm(pivotFinal);
    
    var lowercut_corner = [stay_b[0], stay_b[1] + 2*bolt_shaft_rad];
    var upper_corner = [stay_b[0] + 0.7*drop_w, stay_b[1] + drop_h+20];
    var thickness = 50;
    var cutoff = framebar(lowercut_corner, upper_corner, 50);
    
    // offset along the normal
    var n = normalized([-1*(lowercut_corner[1] - upper_corner[1]), lowercut_corner[0] - upper_corner[0]])
    var delta = thickness/2;
    cutoff = cutoff.translate([-delta, 0,0]);
    
    swingArm = difference(
        union(d, a),
        cutoff
        );
    
    swingArm = linear_extrude({height: stay_depth}, swingArm);
    swingArm = add_bolthole(swingArm, pivotFinal, stay_depth + 1); 

    
    return swingArm;
}

function arm(pivotFinal) {
    const ur_drop = [stay_b[0] + drop_w, stay_b[1] + drop_h]; // upper right corner of basic dropout
    
    // push in slightly further so corners are within the dropout rectangle
    var v = [pivotFinal[0] - ur_drop[0], pivotFinal[1] - ur_drop[1]];
    v = normalized(v);
    inner_offset = [ur_drop[0] - v[0]*tube_h*1.3, ur_drop[1] - v[1]*tube_h*1.3]
    inner_offset_up =  [pivotFinal[0] - v[0]*bolt_shaft_rad*1.5, pivotFinal[1] - v[1]*bolt_shaft_rad*1.5];
    
    const a = framebar(inner_offset_up, inner_offset, tube_h);
    const amorph = frametrapezoid(inner_offset_up, pivotFinal, tube_h, 2*bolt_shaft_rad); 

    return union(a, amorph);
}

function basic_dropout() {
    var base = square({size: [drop_w, drop_h], center: false}); 
    base = base.translate([stay_b[0], stay_b[1], 0]);

    // add the hole for fixed position
    const holeloc = [stay_b[0] + drop_w *0.35, stay_b[1]];
    const holewidth = bolt_radius;    
    const inset_y = bolt_shaft_rad;
    var holesq = square({size: [holewidth*2, inset_y*2], center: true});
    holesq = holesq.translate([holeloc[0], holeloc[1], 0])
    var hole = circle({r: holewidth, center: true});
    hole = hole.translate([holeloc[0], holeloc[1] + inset_y, 0]);
    
    var drop = difference(
        base,
        union(
            hole,  holesq
        )
    );
    
    return drop;
}



function rockerFrame(pivotFinal) {
    var rocker;
    
    // construct outer frame
    const f1 = framebar(pivotFinal, rocker_r, thick_top);
    const f2 = framebar(pivotFinal, rocker_b, thick_left);
    const f3 = framebar(rocker_r, rocker_b, thick_right);
    const f =  union(f1, f2, f3);
    rocker = linear_extrude({height: rocker_depth}, f);

    // construct inner frame
    const uppervec = [rocker_r[0] - pivotFinal[0], rocker_r[1] - pivotFinal[1]];
    const frac = 0.55;
    const barpt = [pivotFinal[0] + uppervec[0]*frac, pivotFinal[1] + uppervec[1]*frac];
    const f_inner = framebar(rocker_b, barpt, thick_inner);
    
    const tri = polygon([pivotFinal, rocker_r, rocker_b]);
    var solid = square({size: [solid_offset*2, solid_offset*3], center: true}).translate([pivotFinal[0], pivotFinal[1], 0]);
    // solid = solid.rotateZ(-solid_angle);
    const solid_inset = intersection(tri, solid);
    
    var inner_frame = union(f_inner, solid_inset);
    const inner_depth = rocker_depth;// *2/3;
    inner_frame = linear_extrude({height: inner_depth}, inner_frame);
    
    rocker = union(
        rocker,
        inner_frame.translate([0, 0, rocker_depth/2 - inner_depth / 2])
    );
    
    // add the boltholes
    rocker = add_bolthole(rocker, pivotFinal, rocker_depth + 1); 
    rocker = add_bolthole(rocker, rocker_r, rocker_depth + 1);
    rocker = add_bolthole(rocker, rocker_b, rocker_depth + 1);

    return rocker;
}


function add_bolthole(base, location2D, depth) {
    const circle_outer_r = bolt_shaft_rad;
    const circle_inner_r = bolt_radius;
    
    const outer = circle({r: circle_outer_r, center: true});
    const outercyl = linear_extrude({height: depth}, outer);
    const inner = circle({r: circle_inner_r, center: true});
    const innercyl = linear_extrude({height: depth}, inner);
    
    const result = difference(
        union(
            base, 
            outercyl.translate([location2D[0], location2D[1], 0])
        ), 
        innercyl.translate([location2D[0], location2D[1], 0])
    );
    return result;
}

function framebar(p1, p2, width) {
    const delta = width / 2;
    const dx = p2[0] - p1[0];
    const dy = p2[1] - p1[1]; 
    var n = [-dy, dx];
    const nlength = Math.sqrt(n[0] * n[0] + n[1] * n[1]);
    n = [n[0] / nlength, n[1] / nlength];
    
    const offset = [n[0]*delta, n[1]*delta];

    const v1 = [p1[0] + offset[0], p1[1] + offset[1]];
    const v2 = [p1[0] - offset[0], p1[1] - offset[1]];
    const v3 = [p2[0] - offset[0], p2[1] - offset[1]];
    const v4 = [p2[0] + offset[0], p2[1] + offset[1]];

    return polygon([v1, v2, v3, v4]);
}

function frametrapezoid(p1, p2, width1, width2) {
    const delta1 = width1 / 2;
    const delta2 = width2 / 2;
    const dx = p2[0] - p1[0];
    const dy = p2[1] - p1[1]; 
    var n = [-dy, dx];
    const nlength = Math.sqrt(n[0] * n[0] + n[1] * n[1]);
    n = [n[0] / nlength, n[1] / nlength];
    
    const offset1 = [n[0]*delta1, n[1]*delta1];
    const offset2 = [n[0]*delta2, n[1]*delta2];

    const v1 = [p1[0] + offset1[0], p1[1] + offset1[1]];
    const v2 = [p1[0] - offset1[0], p1[1] - offset1[1]];
    const v3 = [p2[0] - offset2[0], p2[1] - offset2[1]];
    const v4 = [p2[0] + offset2[0], p2[1] + offset2[1]];

    return polygon([v1, v2, v3, v4]);
}


// ========================================================
//      Generic helpers: rescaling inputs & math functions
// ========================================================

// functions to scale the meshes
function lerp(t, minval, maxval) {
    return minval + t* (maxval - minval);
}

function scaleParam(tvec, minvals, maxvals, scaleFactor) {
    var res = [];
    for (i=0; i<3; i++) {
        res.push(lerp(tvec[i]) * scaleFactor, minvals[i], maxvals[i]);
    }
    return new CSG.Vector3D(res);
}

function normalized(n) {
    const nlength = Math.sqrt(n[0] * n[0] + n[1] * n[1]);
    nnormed = [n[0] / nlength, n[1] / nlength];
    return nnormed;
}

