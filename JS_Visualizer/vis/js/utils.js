color_values = ['#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'];

function get_distinct_color(i) {
    return color_values[i % color_values.length];
}

function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
}

function rgbToHex(r, g, b) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

/* accepts parameters
 * h  Object = {h:x, s:y, v:z}
 * OR 
 * h, s, v
*/
//https://stackoverflow.com/a/17243070
function HSVtoRGB(h, s, v) {
    var r, g, b, i, f, p, q, t;
    if (arguments.length === 1) {
        s = h.s, v = h.v, h = h.h;
    }
    i = Math.floor(h * 6);
    f = h * 6 - i;
    p = v * (1 - s);
    q = v * (1 - f * s);
    t = v * (1 - (1 - f) * s);
    switch (i % 6) {
        case 0: r = v, g = t, b = p; break;
        case 1: r = q, g = v, b = p; break;
        case 2: r = p, g = v, b = t; break;
        case 3: r = p, g = q, b = v; break;
        case 4: r = t, g = p, b = v; break;
        case 5: r = v, g = p, b = q; break;
    }
    return [
        Math.round(r * 255),
        Math.round(g * 255),
        Math.round(b * 255)
    ];
}

function isArray(x) {
    return Array.isArray(x);
}

function distance_points_sqr(lhs, rhs) {
    var e = 0;
    for(var i = 0; i < lhs.length; i++)
        e += (lhs[i] - rhs[i]) * (lhs[i] - rhs[i]);
    return e;
}

function file_extension(filename) {
    return filename.split('.').pop();
}

//https://stackoverflow.com/a/40475362
function makeArr(startValue, stopValue, cardinality) {
    var arr = [];
    var step = (stopValue - startValue) / (cardinality - 1);
    for (var i = 0; i < cardinality; i++) {
      arr.push(startValue + (step * i));
    }
    return arr;
}

function range(n) {
    return [...Array(n).keys()];
}

function lerp(t, min, max) {
    return min + (max - min) * t;
}

function random_id() {
    return Math.random().toString(36).replace(/[^a-z]+/g, '').substr(2, 10);
}

String.prototype.replaceAll = function(search, replacement) {
    var target = this;
    return target.replace(new RegExp(search, 'g'), replacement);
};

function remove_children(parent) {
    while (parent.firstChild) {
        parent.firstChild.remove();
    }
}

function scale_axis_to_unity(value, axis_info) {
    return (value - axis_info.min) / (axis_info.max - axis_info.min);
}

function parse_mapping_entry(m) {
    const shortcuts = {
        "d": "design",
        "p": "performance",
        "a": "application"
    };
    if(m[0] in shortcuts)
        return [shortcuts[m[0]], m[1]]
    else return m;
}

function parse_mapping(options) {
    return options.map(m => parse_mapping_entry(m));
}

function evaluate_mapping_entry(entry, point) {
    const part1 = point[entry[0] + "_point"];
    if(entry[1] != null)
        return part1[entry[1]];
    else return part1;
}

function evaluate_mapping(mapping, point) {
    return mapping.map(entry => evaluate_mapping_entry(entry, point));
}

//mesh = [triangle], triangle = [v1, v2, v3], vi = [x, y ,z]
function compute_mesh_center_of_mass(mesh) {
    let totalVolume  = 0;
    let center = [0, 0, 0];

    for (let i = 0; i < mesh.length; i++) {
        let tri = mesh[i];
        let currentVolume = (tri[0].x*tri[1].y*tri[2].z - tri[0].x*tri[2].y*tri[1].z - tri[1].x*tri[0].y*tri[2].z + tri[1].x*tri[2].y*tri[2].z + tri[2].x*tri[0].y*tri[1].z - tri[2].x*tri[1].y*tri[0].z) / 6;
        center[0] += ((tri[0].x + tri[1].x + tri[2].x) / 4) * currentVolume;
        center[1] += ((tri[0].y + tri[1].y + tri[2].y) / 4) * currentVolume;
        center[2] += ((tri[0].z + tri[1].z + tri[2].z) / 4) * currentVolume;
        totalVolume += currentVolume;
    }

    center = [center[0] / totalVolume, center[1] / totalVolume, center[2] / totalVolume];
    return {center, totalVolume};
}

async function load_js_file(path) {
    return new Promise(function(resolve, reject) {
        var script = document.createElement('script');
        script.type = "text/javascript";
        script.src = path;
        script.onload = () => resolve();
        script.onerror = () => reject();
        document.head.appendChild(script);
    });
}

function geometry_from_data(vertices, triangles) {
    var geom = new THREE.Geometry();

    for(let vertex_index = 0; vertex_index < vertices.length; vertex_index++) {
        let v = vertices[vertex_index];
        geom.vertices.push(new THREE.Vector3(v[0],v[1],v[2]));
    }

    for(let triangle_index = 0; triangle_index < triangles.length; triangle_index++) {
        let t = triangles[triangle_index];
        geom.faces.push( new THREE.Face3( t[0], t[1], t[2] ) );
    }

    return geom;
}

// extracts a list of triangles from a mesh [[[x,y,z], [x,y,z], [x,y,z]]]
async function data_from_mesh(mesh) {
    let geometry = mesh.geometry;

    if(geometry.faces == null || geometry.vertices == null) {
        geometry = new THREE.Geometry().fromBufferGeometry(geometry);
    }

    let faces = geometry.faces;
    let vertices = geometry.vertices;

    let triangles = [];
    for ( let face_idx = 0; face_idx < faces.length; face_idx++ ) {
        let face = geometry.faces[face_idx];
        triangles.push([
            vertices[face.a],
            vertices[face.b],
            vertices[face.c]
        ]);
    }

    return triangles;
}

function download_blob(filename, blob) {
    var downloadAnchorNode = document.createElement('a');
    downloadAnchorNode.setAttribute("href",     window.URL.createObjectURL(blob)     );
    downloadAnchorNode.setAttribute("download", filename);
    document.body.appendChild(downloadAnchorNode); // required for firefox
    downloadAnchorNode.click();
    downloadAnchorNode.remove();
}