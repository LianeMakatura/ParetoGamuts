var g_fn = 5;

function createLampGeometry(
    base, a1, a2, a3, h1, h2, h3,
    radii, stand_h, stand_r,
    shade_end_radius, shade_height, shade_thickness) {
    var res = [];

    res = res.concat(bottom_stand(stand_r, stand_h));

    res = res.concat(rod(new CSG.Vector3D([0,0,0]), base, radii[0])); //base
    res = res.concat(rod(base, base.plus(a1), radii[1])); //arm 1
    res = res.concat(rod(base.plus(a1), base.plus(a1.plus(h1)), radii[1])); //hand 1
    res = res.concat(rod(base, base.plus(a2), radii[2])); //arm 2
    res = res.concat(rod(base.plus(a2), base.plus(a2.plus(h2)), radii[2])); //hand 2
    res = res.concat(rod(base, base.plus(a3), radii[3])); //arm 3
    res = res.concat(rod(base.plus(a3), base.plus(a3.plus(h3)), radii[3])); //hand 3

    res = res.concat(lampshade_full(base.plus(a1.plus(h1)), h1, radii[1],
                shade_end_radius, shade_height, shade_thickness));
    res = res.concat(lampshade_full(base.plus(a2.plus(h2)), h2, radii[2],
                shade_end_radius, shade_height, shade_thickness));
    res = res.concat(lampshade_full(base.plus(a3.plus(h3)), h3, radii[3],
                shade_end_radius, shade_height, shade_thickness));

    return union(res).translate([0,0,stand_h]);
}

function cylinder_correct(start, end, r1, r2, fn) {
    var dir = end.minus(start);
    var c_length = dir.length();
    var c = cylinder({start: [0,0,0], end: [0,0,c_length], r1: r1, r2: r2, fn: fn});
    var t = new CSG.Matrix4x4.translation(start);
    if(!(dir.x == 0 && dir.y == 0)) {
        var w  = dir.dividedBy(c_length);
        var u0 = w.cross(new CSG.Vector3D([0,0,1]));
        var u  = u0.dividedBy(u0.length());
        var v0 = w.cross(u);
        var v  = v0.dividedBy(v0.length());
        var r = new CSG.Matrix4x4([u.x, u.y, u.z, 0,
                                   v.x, v.y, v.z, 0,
                                   w.x, w.y, w.z, 0,
                                   0,   0,   0,   1]);
        return c.transform(r.multiply(t));
    }
    else return c.transform(t);
}

function rod(a, b, r) {
    var res = [];
    res.push(cylinder_correct(a,b,r,r,g_fn));
    res.push(sphere({r: r, fn: g_fn}).translate(a));
    res.push(sphere({r: r, fn: g_fn}).translate(b));
    return res;
}

function lampshade_full(origin, dir, connection_pt_radius, end_radius,
                        height, thickness) {
    var q1 = lampshade(origin, dir, height, connection_pt_radius, end_radius);
    var q2 = lampshade(origin, dir, height*1.2,
                    connection_pt_radius*(1-thickness),
                    end_radius*(1-thickness));
    var part1 = difference(q1[0], q2[0]);

    // add the bulb
    var orientation = dir.dividedBy(dir.length());
    var shade_end_local = orientation.times(height);
    var q1 = sphere({r: 0.2, fn: g_fn}).translate(origin.plus(shade_end_local));
    var q2 = rod(origin.plus(shade_end_local.times(0.5)),
        origin.plus(shade_end_local.times(0.9)),
                0.15);
    var part2 = union(q1, q2);

    return [part1, part2];
}

function lampshade(origin, dir, height, conn_pt_r, end_r) {
    //return [cylinder({start: origin, end: origin.plus(dir.times(height)), r1: conn_pt_r, r2: end_r, fn: g_fn})];
    return [cylinder_correct(origin, origin.plus(dir.times(height)), conn_pt_r, end_r, g_fn)]
}

function bottom_stand(stand_r, stand_h) {
    return [cylinder({h:stand_h, r:stand_r, fn: g_fn}).translate([0,0,-stand_h])];
}


function lerp(minval, maxval, t) {
    return minval + t* (maxval - minval);
}

function scaleParam(tvec, minvals, maxvals, scaleFactor) {
    var res = [];
    for (i=0; i<3; i++) {
        res.push(lerp(minvals[i], maxvals[i], tvec[i]) * scaleFactor);
    }
    return new CSG.Vector3D(res);
}

function main() {
    //here dp will be the design point array and ap the application point
    //INSERT_PARAMETER_VALUES
    
    var scaleFactor = 4;

    var r = 0.1 * scaleFactor;
    var radii = [r, r, r, r];
    var stand_h = 0.25* scaleFactor;
    var stand_r = 1* scaleFactor;
    var shade_end_radius = 0.4* scaleFactor;
    var shade_height = 0.8;
    var shade_thickness=0.1;


    // vectors defined wrt the ancestor endpoint
    var base = scaleParam(dp.slice(0, 3), [0,0,2], [1,1,6], scaleFactor);
    var a1 = scaleParam(dp.slice(3, 6), [0,0,1], [1,1,2], scaleFactor);
    var a2 = scaleParam(dp.slice(6, 9), [-1, 0, 1], [0,1,2], scaleFactor);
    var a3 = scaleParam(dp.slice(9, 12), [0,-1,1], [1, 0,2], scaleFactor);
    var h1 = scaleParam(dp.slice(12, 15), [0.4, 0.4, -2], [1, 1, -1], scaleFactor);
    var h2 = scaleParam(dp.slice(15, 18), [-1,0.4,-2], [-0.4, 1, -1], scaleFactor);
    var h3 = scaleParam(dp.slice(18, 21), [0.4, -1, -2], [1, -0.4, -1], scaleFactor);

    return createLampGeometry(base, a1, a2, a3, h1, h2, h3,
                radii, stand_h, stand_r,
                shade_end_radius, shade_height, shade_thickness
                );
}