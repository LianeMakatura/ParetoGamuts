function color_map_to_plotjs(points, point_map, group_map, showLine) {
    const dataset_map = {};
    const datasets = [];
    for(var color_grp_idx in group_map) {
        dataset_map[color_grp_idx] = datasets.length;
        const color = group_map[color_grp_idx];
        datasets.push({
            data: [],
            showLine: showLine,
            fill: false,
            pointBackgroundColor: color //"rgba(238,36,36,1)" //--fixed color for all points
        });
    }

    for(var idx = 0; idx < points.length; idx++) {
        const color_grp_idx = point_map[idx];
        const dataset_idx = dataset_map[color_grp_idx];
        const p = points[idx].point;
        datasets[dataset_idx].data.push({x: p[0], y: p[1]});
    }

    return datasets;
}

//computes a map from point_idx to color_grp_idx and a map color_grp_idx to color string
class GroupIdColorMap {
    compute_map(points) {
        const point_map = [];
        const group_map = [];
        const patch_map = {};

        for(var point_idx = 0; point_idx < points.length; point_idx++) {
            const point_grp_id = points[point_idx].attributes.patch_id;

            if(!(point_grp_id in patch_map)) {
                patch_map[point_grp_id] = group_map.length;
                group_map.push(get_distinct_color(point_grp_id));
            }

            const group_id = patch_map[point_grp_id];
            point_map.push(group_id);
        }

        return [point_map, group_map];
    }
}

class AttributeColorMap {
    constructor(options) {
        this.mapping = parse_mapping(options.mapping);
        this.color_mode = options.color_mode;
    }

    compute_map(points) {
        const point_map = [];
        const group_map = [];

        for(var point_idx = 0; point_idx < points.length; point_idx++) {
            //TODO p needs to be normalized
            const vals = evaluate_mapping(this.mapping, points[point_idx].attributes);
            let rgb = [0,0,0];

            if(this.color_mode == null || this.color_mode == "rgb") {
                for(let c = 0; c < Math.min(vals.length, 3); c++)
                    rgb[c] = Math.floor(vals[c] * 255);
            }
            else if(this.color_mode == "hsv") {
                const hsv = [0,0,0];
                for(let c = 0; c < Math.min(vals.length, 3); c++)
                    hsv[c] = vals[c];
                rgb = HSVtoRGB(...hsv);
            }
            else throw new Error(`Invalid color mode: ${this.color_mode}`);

            const color = rgbToHex(...rgb);
            point_map.push(point_idx);
            group_map.push(color);
        }

        return [point_map, group_map];
    }
}

class RainbowColorMap {
    constructor(options) {
        this.mapping = parse_mapping(options.mapping);
        this.s_value = options.s_value != null ? options.s_value : 1;
        this.v_value = options.v_value != null ? options.v_value : 1;
    }

    compute_map(points) {
        const point_map = [];
        const group_map = [];

        for(var point_idx = 0; point_idx < points.length; point_idx++) {
            //TODO p needs to be normalized
            const vals = evaluate_mapping(this.mapping, points[point_idx].attributes);

            const rgb = HSVtoRGB(vals[0], this.s_value, this.v_value);

            const color = rgbToHex(...rgb);
            point_map.push(point_idx);
            group_map.push(color);
        }

        return [point_map, group_map];
    }
}

function create_grayscale_color_map(options) {
    const opt2 = Object.assign({}, options);
    opt2.mapping = [options.mapping, options.mapping, options.mapping];
    return new AttributeColorMap(opt2);
}

class DiscreteColorMap {
    compute_map(points) {
        const point_map = [];
        const group_map = [];

        for(var point_idx = 0; point_idx < points.length; point_idx++) {
            point_map.push(point_idx);
            var color = get_distinct_color(point_idx);
            group_map.push(color);
        }

        return [point_map, group_map];
    }
}

function create_color_map(options) {
    if(options == null || options.type == null)
        options = {type: "GroupIdColorMap"};

    if(options.type == "GroupIdColorMap")
        return new GroupIdColorMap();
    if(options.type == "AttributeColorMap")
        return new AttributeColorMap(options);
    if(options.type == "RainbowColorMap")
        return new RainbowColorMap(options);
    if(options.type == "GrayscaleColorMap")
        return create_grayscale_color_map(options);
    if(options.type == "DiscreteColorMap")
        return new DiscreteColorMap();
    throw new Error("Unknown color map type: " + options.type);
}

//points is a list of {point, color} where point is {point, attributes}
//color can either be null which then applies the color_map to this point
//or a color string directly
//this function returns the same as color_map.compute_map but also a reordered list of points
function apply_color_map_partial(points, color_map) {
    const mappable_points = [], fixed_points = [];
    for(let idx = 0; idx < points.length; idx++) {
        if(points[idx].color == null)
            mappable_points.push(points[idx].point);
        else fixed_points.push(points[idx]);
    }

    let [point_map, group_map] = color_map.compute_map(mappable_points);

    for(let q of fixed_points) {
        mappable_points.push(q.point);

        let group_idx = group_map.findIndex(x => x == q.color);
        if(group_idx == -1) {
            group_idx = group_map.length;
            group_map.push(q.color);
        }

        point_map.push(group_idx);
    }

    return [mappable_points, point_map, group_map];
}