/*
Parametric Bike Stay & Rocker 
*/

// ==============================================
//              Define parameters
// ==============================================

var eps = 1e-3;

// ------ Design parameters
//here dp will be the design point array and ap the application point
dp = [0,1,0.5]; // these are just for testing, they'll be overwritten below
ap = [1];
g_fn = 5; //default resolution, also overwritten 

//INSERT_PARAMETER_VALUES


// --------- bike rocker
const bulginess = dp[0];
const rocker_tube_h = lerp(dp[1], 6, 9);
const holiness = (dp[1] / 4) + 0.5; // revisit
const rocker_depth = lerp(dp[2], 5, 8);


// --------- bike stay
const drop_material = dp[0];
const stay_tube_h = lerp(dp[1], 6, 9);
const lower_height = lerp(dp[1], 6, 9);
const stay_depth = lerp(dp[2], 5, 8);


// --------- shared param
const ap_normalized = ap[0];
const pivotoff_y = lerp(ap[0], -9, 40);
if (ap.length < 2) {
    pivotoff_x = 0; // no  y offset if only 1 application variable
}
else {
    pivotoff_x = lerp(ap[1], 0, -40);
}



// ==============================================
//              Define mesh
// ==============================================

function main(params) {
    let parts = [];
    parts.push(left_part(pivotoff_y, stay_tube_h, drop_material, stay_depth, lower_height));
    // parts.push(top_part(pivotoff_y, holiness, 20*(ap_normalized/2 + 0.5), bulginess, rocker_depth, lower_height));

    parts.push(bottom_part(stay_depth, lower_height));
    // parts.push(right_part(stay_depth));

    return parts;
}

function right_part(stay_depth) {
    const p5 = unproject([560.5, 610.6]);
    const p8 = unproject([480.8, 390.1]);
    
    const total_height = length(v_sub(p5, p8));
    const top_cylinder_rad = stay_depth / 2 * 0.5;
    const top_cylinder = cylinder({start: lerp_point(0.1, p5, p8), end: lerp_point(0.5, p5, p8), r1: top_cylinder_rad, r2: top_cylinder_rad});
    
    const bottom_cylinder_rad = stay_depth / 10;
    const bottom_cylinder = cylinder({start: lerp_point(0.1, p5, p8), end: p8, r1: bottom_cylinder_rad, r2: bottom_cylinder_rad});
    
    const spiral_length_per = 0.8;
    let spiral = linear_extrude({ height: total_height * spiral_length_per, twist: 360*9, slices: 150}, circle().translate([stay_depth / 4,0,0]) );
    spiral = spiral.rotateX(90);
    const theta_rad = Math.acos(dot(normalize(v_sub(p5, p8)), [0,-1]));
    const theta = 180 - theta_rad * 180/Math.PI;
    spiral = spiral.rotateZ(-theta).translate(lerp_point(1 - spiral_length_per, p5 ,p8));
    
    let shock = union(top_cylinder, bottom_cylinder, spiral);
    return shock.translate([0, 0, rocker_depth/2]);
}

function bottom_part(stay_depth, width) {
    const p3 = unproject([109.39, 323.95]);
    const p6 = unproject([353.6, 340.9]);
    const p7 = unproject([440.1, 374.8]);
    
    let a1, a2, a3, a4;
    {
        let d;
        [a1, a4, a2, a3, d] = cage_at_line_end(p3, p6, width / 2);
    }
    
    let b1, b2;
    {
        let d,x,y;
        [b1, b2, x, y, d] = cage_at_line_end(p6, p7, width / 2);
    }
    
    let c1, c2, c3, c4;
    {
        let d;
        [c1, c4, c2, c3, d] = cage_at_line_end(p7, p6, width / 2);
    }
    
    let path = new CSG.Path2D([a1], false);
    path = path.appendBezier([a2, a3, a4]);
    path = path.appendBezier([b2]);
    path = path.appendBezier([c1]);
    path = path.appendBezier([c2, c3, c4]);
    path = path.appendBezier([b1]);
    path = path.appendBezier([a1]);
    
    path = path.close();
    
    let part = linear_extrude({height: stay_depth}, path.innerToCAG());
    part = add_bolthole(part, p3, stay_depth + 1);
    part = add_bolthole(part, p7, stay_depth + 1);
    
    return part;
}

function top_part(pivotoff_y, holiness, bulge_ext, bulginess, stay_depth, width) {
    const p1 = unproject([353.6, 566.5 + pivotoff_y]);
    const p4 = unproject([482.5, 568.2]);
    const p5 = unproject([560.5, 610.6]);
    
    let a1, a2, a3, a4;
    {
        let d;
        [a1, a4, a2, a3, d] = cage_at_line_end(p1, p5, width / 2);
    }
    
    let b1, b2, b3, b4;
    {
        let d;
        [b1, b4, b2, b3, d] = cage_at_line_end(p5, p1, width / 2);
        b2 = lerp_point(0.3, b2, b3);
    }
    
    let c1, c2, c3;
    {
        c2 = v_add(p4, [bulge_ext/3, -bulge_ext]);
        const c = project_onto_line(a1, b4, p4);
        c1 = lerp_point(lerp(bulginess, -0.3, 0.7), c, b4);
        c3 = lerp_point(lerp(bulginess, -0.3, 0.7), c, a1);
    }
    
    const q = lerp_point(0.8, a4, b1);
    const a5 = v_add(q, [0, bulge_ext/5]);
    
    let path = new CSG.Path2D([a1], false);
    path = path.appendBezier([a2, a3, a4]);
    path = path.appendBezier([a5, b1]);
    path = path.appendBezier([b2, b3, b4]);
    path = path.appendBezier([c1, c2, c3, a1]);
    
    path = path.close();
    
    let part = linear_extrude({height: stay_depth}, path.innerToCAG());
    
    part = add_bolthole(part, p1, stay_depth + 1);
    part = add_bolthole(part, p5, stay_depth + 1);
    const p4_hole_radius = 2;
    const inner = circle({r: p4_hole_radius*3/4, center: true});
    const innercyl = linear_extrude({height: stay_depth}, inner);
    part = difference(part, innercyl.translate([p4[0], p4[1], 0]));
    
    // the hole cannot be so large that it swollows the p4 hole.
    // compute the projection of p4 onto line
    const pp = project_onto_line(p1, p5, p4);
    // compute the distance between p4 and pp and substract the radius of the p4 hole
    const max_hole_ext = (length(v_sub(p4, pp)) - p4_hole_radius) * 0.6;
    
    const side_interpolant = lerp(holiness, 0.75, 0.5);
    const wide_interpolant = lerp(holiness, width/8, Math.min(width, max_hole_ext * 1.5));
    const pc = lerp_point(0.5, p1, p5);
    let hole_path = new CSG.Path2D([lerp_point(side_interpolant, p1, pc), lerp_point(side_interpolant, p5, pc)], false);
    const hole_part = hole_path.rectangularExtrude(wide_interpolant,stay_depth,64, true);
    
    return difference(part, hole_part);
}

function left_part(pivotoff_y, tube_h, drop_material, stay_depth, lower_height) {
    const p1 = unproject([353.6, 566.5 + pivotoff_y]);
    const p2 = unproject([0, 340.9]);
    const p3 = unproject([109.39, 323.95]);
    
    let a1, a2, a3, a4;
    let off1;
    {
        let d;
        [a1, a2, a3, a4, d] = cage_at_line_end(p1, p2, tube_h / 2);
        const offset = v_sub(a2, p1);
        
        off1 = [p2[0] + offset[0], p2[1] + offset[1]];
    }
    
    let b1, b2, b3;
    {
        const d2 = v_sub(p2, p1);
        b1 = v_add(a1, d2);
        b2 = v_add(p2, [-lower_height/2, -lower_height/2]);
        b3 = v_add(p2, [0, -lower_height/2]);
    }
    
    let c1, c2, c3, c4;
    let off2;
    {
        let d;
        [c4, c1, c3, c2, d] = cage_at_line_end(p3, p2, lower_height / 2);
        const offset = v_sub(c4, p3);
        
        off2 = [p2[0] + offset[0], p2[1] + offset[1]];
    }
    
    let e1, e2, e3;
    {
        e2 = intersectLines(off1, a2, c3, off2);
        const v1 = lerp(drop_material, 0.2, 0.9);
        e1 = lerp_point(v1, e2, c4);
        const v2 = lerp(drop_material, 0.2, 0.5);
        e3 = lerp_point(v2, e2, a2);
    }
    
    let path = new CSG.Path2D([a1, b1], false);
    path = path.appendBezier([b2, b3]);
    path = path.appendBezier([null, c1]);
    path = path.appendBezier([c2, c3, c4]);
    path = path.appendBezier([e1, e2, e3]);
    path = path.appendBezier([a2]);
    path = path.appendBezier([a4, a3, a1]);
    
    path = path.close();
    
    let part = linear_extrude({height: stay_depth}, path.innerToCAG());
    
    part = add_bolthole(part, p1, stay_depth + 1);
    part = add_bolthole(part, p2, stay_depth + 1);
    part = add_bolthole(part, p3, stay_depth + 1);
 
    return part;
}

function add_bolthole(base, location2D, depth) {
    const bolt_radius = 2;
    const bolt_shaft_rad = 5;
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

// ========================================================
//      Generic helpers: rescaling inputs & math functions
// ========================================================

function unproject(p) {
    return [p[0] / 5 - 60, p[1] / 5 - 100];
}


function dot(a, b) {
    return a[0] * b[0] + a[1] * b[1];
}

function dot3(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function length(n) {
    return Math.sqrt(dot(n, n));
}

function length3(n) {
    return Math.sqrt(dot3(n, n));
}

function normalize(n) {
    return v_mul(n, [1.0 / length(n), 1.0 / length(n)]);
}

function normaliz3(n) {
    return v_mul3(n, 1.0 / length3(n));
}

function v_add(lhs, rhs) {
    return [lhs[0] + rhs[0], lhs[1] + rhs[1]];
}

function v_add3(lhs, rhs) {
    return [lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]];
}

function v_sub(lhs, rhs) {
    return [lhs[0] - rhs[0], lhs[1] - rhs[1]];
}

function v_sub3(lhs, rhs) {
    return [lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]];
}

function v_mul(lhs, rhs) {
    return [lhs[0] * rhs[0], lhs[1] * rhs[1]];
}

function v_mul3(lhs, rhs) {
    return [lhs[0] * rhs[0], lhs[1] * rhs[1], lhs[2] * rhs[2]];
}

// functions to scale the meshes
function lerp(t, minval, maxval) {
    return minval + t* (maxval - minval);
}

function lerp_point(t, p1, p2) {
    return [lerp(t, p1[0], p2[0]), lerp(t, p1[1], p2[1])];
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

// line intercept math by Paul Bourke http://paulbourke.net/geometry/pointlineplane/
// Determine the intersection point of two line segments
// Return FALSE if the lines don't intersect
function intersect(x1, y1, x2, y2, x3, y3, x4, y4) {
    // Check if none of the lines are of length 0
      if ((x1 === x2 && y1 === y2) || (x3 === x4 && y3 === y4)) {
          return false
      }

      denominator = ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1))

    // Lines are parallel
      if (denominator === 0) {
          return false
      }

      let ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator
      let ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator

    // is the intersection along the segments
    //   if (ua < 0 || ua > 1 || ub < 0 || ub > 1) {
    //       return false
    //   }

    // Return a object with the x and y coordinates of the intersection
      let x = x1 + ua * (x2 - x1)
      let y = y1 + ua * (y2 - y1)

      return [x, y]
}

function intersectLines(p1_l1, p2_l1, p1_l2, p2_l2) {
    return intersect(...p1_l1, ...p2_l1, ...p1_l2, ...p2_l2);
}

function project_onto_line(A, B, p) {
    const x1=A[0], y1=A[1], x2=B[0], y2=B[1], x3=p[0], y3=p[1];
    const px = x2-x1, py = y2-y1, dAB = px*px + py*py;
    const u = ((x3 - x1) * px + (y3 - y1) * py) / dAB;
    const x = x1 + u * px, y = y1 + u * py;
    return [x,y]; //this is D
}

// creates offset points on the sides of p1, distance delta away, along normal n (normal to p1 - p2)
// also creates cage points for a bezier
// also returns vector from p2 to p1
function cage_at_line_end(p1, p2, delta) {
    const dx = p1[0] - p2[0];
    const dy = p1[1] - p2[1]; 
    let n = [-dy, dx];
    const nlength = Math.sqrt(n[0] * n[0] + n[1] * n[1]);
    n = [n[0] / nlength, n[1] / nlength];
    
    const offset = [n[0]*delta, n[1]*delta];
    
    const a = [p1[0] + offset[0], p1[1] + offset[1]];
    const b = [p1[0] - offset[0], p1[1] - offset[1]];
    const c = v_add(a, [dx / nlength * delta, dy / nlength * delta]);
    const d= v_add(b, [dx / nlength * delta, dy / nlength * delta]);
    
    return [a, b, c, d, [dx, dy]];
}

// takes a 2d vector, returns basis in R3
function ONB(n) {
    const a = 1.0 / (1.0 + 0.0);
    const b = -n[0] * n[1] + a;
    const b1 = [1 - n[0] * n[0] * a, b, -n[0]];
    const b2 = [b, 1 - n[1] * n[1] * a, -n[1]];
    return [b1, b2, n];
}

function helix(t, r, p) {
    return [r * Math.sin(t), r * Math.cos(t), p * t];
}

function to_world(basis, p, x) {
    const [t,b,n] = basis;
    return [
        t[0] * x[0] + b[0] * x[1] + n[0] * x[2] + p[0],
        t[1] * x[0] + b[1] * x[1] + n[1] * x[2] + p[1],
        t[2] * x[0] + b[2] * x[1] + n[2] * x[2] + p[2]
        ];
}