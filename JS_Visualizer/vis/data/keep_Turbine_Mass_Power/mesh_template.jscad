/*
Parametric turbine design, based on the OpenSCAD model by Denise Lee (12March2017)
*/

// ==============================================
//              Define parameters
// ==============================================


//here dp will be the design point array and ap the application point
dp = [1,0.5,1]; // these are just for testing, they'll be overwritten below
ap = [1];
g_fn = 5; //default resolution, also overwritten - not used right now

//INSERT_PARAMETER_VALUES

// === Design parameters
var blade_radius = lerp(dp[0], 40, 100);
var turbine_height = lerp(dp[1], 6, 15);
var pitch_angle = lerp(dp[2], 5, 20); //Converting pitch to twist angle for linear extrude

// === contextual parameters (none related to design)


// === Other parameters, fixed for now
// continuous
var blade_thickness = 4; //thickness of blade
var shaft_fit_r = 2.6; //Your connecting motor pin to turbine 
var stem_top_r = 4+shaft_fit_r;
var stem_bot_r = 5+shaft_fit_r;

var blade_arch = blade_radius/11; // Curvature of blades, selected straight edge or use ARC module 
var percent_offset = 40; //Percent of turbine height with blades at blade radius


// === Discrete choices
var apply_curvature =1; //1 Curve blades, 0 for straight edge
var rotation_dir = -1; //CCW = -1, CW = 1
var num_blades = 3;
var curve_profile = 1; // convex if 1, none if 0, concave if -1

// === Fabrication choices
var printer_layer_height = 1; //for preview (F5) x20, set this to 0.1 for rendering for export STL(F6)!!!
var slicing = turbine_height / printer_layer_height; //equal printing slicing height
var offset_slicing = Math.round(slicing*(1/((100-percent_offset)/100)));  //offset_slicing must be greater or equal slicing




// ==============================================
//              Define mesh
// ==============================================

function main(){
    // should set the parameters here as well
    return turbine();
}

function turbine() {
    var turbine;
    // create the blade profile sketches
    for (i=0; i < num_blades; i++) {
        b = blade(i);
        
        if (i == 0) {
            turbine = b;
        } else {
            turbine = union(turbine, b);
        }
    }
    
    // twisting extrude
    var blade_cirf = 2*Math.PI*blade_radius;
    var twist_angle = 360*turbine_height/(blade_cirf*Math.tan(pitch_angle* Math.PI/180));
    turbine = linear_extrude({
        height: turbine_height,
        center: false,
        twist: twist_angle*rotation_dir,
        slices: slicing}, turbine);
        
    // adjust curve profile -- intersection is time intensive
    if (curve_profile != 0) {
        if (curve_profile == 1) {
            curve = edge_curve(1); // convex
        } else if (curve_profile == -1) {
            curve = edge_curve(0); // concave
        }
        turbine = intersection(turbine, curve);
    }
    
    // include center and pin hole
    var stem = cylinder({
        h: turbine_height, 
        r1: stem_bot_r, 
        r2: stem_top_r,
        center: false, 
        fn: 20});
    var pin = cylinder({
        h: turbine_height + 0.1, 
        r1: shaft_fit_r, 
        r2: shaft_fit_r,
        center: false, 
        fn: 20});
    turbine = union(turbine, stem);
    turbine = difference(turbine, pin);
    
    return turbine;
}

function blade(id) {
    var blade;
    if(apply_curvature != 1) { 
        blade = square([blade_radius, blade_thickness]); 
    } else {
        // used curved blade
        blade = arc(blade_radius, blade_thickness, blade_arch); 
        if(rotation_dir == -1) {
            mirror([0,1,0], blade);
        }
    }
    blade = translate([0,-blade_thickness/2], blade);
    blade = rotate([0, 0, id*360/num_blades], blade);
    
    return blade;
}


//length and breadth of inner arc
function arc(length, width, arch_height){
    //r = (l^2 + 4*b^2)/8*b 
    radius = (pow(length,2) + 4*pow(arch_height,2))/(8*arch_height);

    var res1 = difference(
        // ring from concentric circles 
        difference(
            circle({r:radius+width, center:true, fn:100}),
            circle({r:radius, center:true, fn:100})
        ).translate([0,-radius+arch_height,0]),
    
        // bounding square to truncate arc 
        square(
            [(radius+width)*2,(radius+width)*2]
        ).translate([-(radius+width),-(radius+width)*2,0])
    )

    // two rectangles to trim the edges of the arc
    var res2 = union(
        square(
            [length/2,arch_height*2]
        ).translate([-length,-arch_height]),
        square(
            [length/2,arch_height*2]
        ).translate([length/2,-arch_height])
    );
    
    resf = difference(res1, res2).translate([length/2,0,0]);

    return resf;
}


function edge_curve(convex) {
    var layer_h = turbine_height/offset_slicing;
    const epsilon = 0.005;
    
    if (convex == 1) {
        //Curve convex: y = 1 - 1/exp(x)
        exp_pow = 3;
        end_point = 1 - 1/Math.exp(exp_pow*1); 
    } else {
        //Curve concave: y = base^x - 1
        base = 3;
        end_point = pow(base,1)-1;
    }

    var res;
    var c;

    var bot_r = blade_radius/6; //Bottom blade radius
    var top_r= blade_radius;
    var delta_r = top_r - bot_r;
    for (i=0; i < offset_slicing; i++) {
        // create the slice
        if (i < slicing) { 
            if (convex == 1) { 
                offset_r = delta_r*((1-1/Math.exp(exp_pow*(i/slicing)))/end_point); //Normalised N
                offset_r_increment = delta_r*((1-1/Math.exp(exp_pow*((i+1)/slicing)))/end_point); //Normalised N+1  
            } else {
                offset_r = delta_r*((pow(base,i/slicing)-1)/end_point); //Normalised N 
                offset_r_increment = delta_r*((pow(base,(i+1)/slicing)-1)/end_point);  //Normalised N+1    
            }
            c = cylinder({
                h: layer_h + epsilon, 
                r1: bot_r + offset_r, 
                r2: bot_r + offset_r_increment, 
                center: false, 
                fn: 20});
            c = translate([0,0,i*layer_h], c);
        } else {
            c = cylinder({
                h: layer_h + epsilon, 
                r1: top_r, 
                r2: top_r, 
                center: false, 
                fn:20});    
            c = translate([0,0,i*layer_h], c);
        }  
        
        // union it with others
        if (i==0) {
            res = c;
        } else {
            res = union(res, c);
        }
    }
    
    return res;
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

