/*
Parametric LBracket design
Author: Anonymous, February 2020
*/

// ==============================================
//              Define parameters
// ==============================================

var eps = 1e-3;

// ------ Design parameters
//here dp will be the design point array and ap the application point
dp = [0,0,0]; // these are just for testing, they'll be overwritten below
ap = [1];
g_fn = 5; //default resolution, also overwritten 

//INSERT_PARAMETER_VALUES


var scaleF = 10;

var l = scaleF * lerp(dp[0], 0.4, 1);
var w = scaleF * lerp(dp[1], 0.2, 1);
if (dp.length < 3) {
    var t = scaleF * 0.1;
}
else {
    var t = scaleF * lerp(dp[2], 0.1, 0.3);
}

var angle = lerp(ap[0], 45 + eps, 90 - eps);



// ==============================================
//              Define mesh
// ==============================================

function main() {
    return Lbracket();
}

function Lbracket() {
    base = base_sketch();
    
    bracket = linear_extrude({height:w, center: false, twist: 0, slices: 1}, base);
    return bracket.rotateX(90);
}

function base_sketch () {
    // define invariant side
    v3 = new CSG.Vector2D([0,0]);
    v5 = v3.plus(new CSG.Vector2D([l, 0]));
    v6 = v5.plus(new CSG.Vector2D([0, t]));

    // define pts to be rotated
    v1temp = new CSG.Vector2D(v5);
    v2temp = v1temp.minus(new CSG.Vector2D([0, t]));

    // rotate the "vertical" beam as desired
    var rot = CSG.Matrix4x4.rotationX(-angle);
    rotv1 = rot.rightMultiply1x3Vector(new CSG.Vector3D([0, v1temp.x, v1temp.y]));
    rotv2 = rot.rightMultiply1x3Vector(new CSG.Vector3D([0, v2temp.x, v2temp.y]));

    v1 = new CSG.Vector2D([rotv1.y, rotv1.z]);
    v2 = new CSG.Vector2D([rotv2.y, rotv2.z]);

    m_rot = slope2D(v1, v3);
    m_horz = slope2D(v3, v5); //should always be 0
    v4 = intersectLines2D(v2, m_rot, v6, m_horz);

    p = polygon([v1, v2, v4, v6, v5, v3]);
    return p;
}

function slope2D(p1, p2) {
    return (p2.y - p1.y) / (p2.x - p1.x);
}

function intersectLines2D(p1, m1, p2, m2) {
    //Find intersection of two lines, given the slope and a point (not
    //necessarily the intercept). Derived from the two-point form.
    p1x = p1.x; p1y = p1.y;
    p2x = p2.x; p2y = p2.y;
    
    x = (p2y - p1y + m1*p1x - m2*p2x) / (m1 - m2);
    y = m2*(x - p2x) + p2y;
    return new CSG.Vector2D([x,y]);
}

// ==============================================
//      Generic helpers: rescaling inputs
// ==============================================

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

