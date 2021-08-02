function is_selection_point_index(item) {
    return Number.isInteger(item);
}

function resolve_selection_point(item, points) {
    if(is_selection_point_index(item))
        return points[item];
    else return item;
}

function make_highlight_point(point, color, use_color) {
    return {
        point,
        color,
        use_color
    };
}

class HighlightGroup {
    //list of {point: index or {design, performance, application, ...}, color: null or string, use_color: bool}
    points = null;

    constructor(p) {
        this.points = p;
    }
}

class SelectionItem {
    //either index or {design, performance, application, ...}
    selected_point = null;
    selected_point_color = null;
    highlighted_groups = [];

    constructor(selected_point, selected_point_color, highlighted_groups = []) {
        this.selected_point = selected_point;
        this.selected_point_color = selected_point_color;
        this.highlighted_groups = highlighted_groups;
    }

    iterate_points(include_selected, include_groups, points, clb) {
        if(include_selected)
            clb(is_selection_point_index(this.selected_point) ? resolve_selection_point(this.selected_point, points) : this.selected_point, this.selected_point_color);

        if(include_groups)
            this.highlighted_groups.forEach(g => {
                g.points.forEach(q => clb(is_selection_point_index(q.point) ? resolve_selection_point(q.point) : q.point, q.use_color ? q.color : this.selected_point_color));
            });
    }

    map_points(include_selected, include_groups, points, clb) {
        let res = [];
        this.iterate_points(include_selected, include_groups, points, (point, color) => {
            res.push(clb(point, color));
        })
        return res;
    }

    map_selection(map_idx, map_point) {
        const apply = (p_item) => is_selection_point_index(p_item) ? (map_idx != null ? map_idx(p_item) : p_item) : (map_point != null ? map_point(p_item) : p_item);
        return new SelectionItem(apply(this.selected_point), this.selected_point_color, this.highlighted_groups.map(g => new HighlightGroup(g.points.map(q => ({point: apply(q.point), color: q.color, use_color: q.use_color})))));
    }
}

function iterate_points(selection_info, include_selected, include_groups, points, clb) {
    selection_info.forEach(selection_item => {
        selection_item.iterate_points(include_selected, include_groups, points, (point, color) => {
            clb(point, color, selection_item);
        });
    });
}

function map_points(selection_info, include_selected, include_groups, points, clb) {
    let res = [];
    iterate_points(selection_info, include_selected, include_groups, points, (point, color, selection_item) => {
        res.push(clb(point, color, selection_item));
    })
    return res;
}

function map_selection_info(selection_info, map_idx, map_point) {
    return selection_info.map(selection_item => selection_item.map_selection(map_idx, map_point));
}

class DataSourceControl {
    html_template_chart = `
        <div id="CHART-ID" style="flex: 1;overflow: hidden">
        </div>
    `;
    html_template_selectors = `
        <div id="SEL-ID" style="flex-grow:0; margin:3px;">
        </div>
    `;
    html_template_selector = `
        <label for="SEL-ID" style="margin:0;">AXIS-LABEL</label> <select id="SEL-ID"></select>
    `;

    constructor(div, options, axes, select_clb) {
        this.axes = axes;
        this.options = options;
        this.select_clb = select_clb;
        this.input = {points: null, filtered_point_indices: null, point_index_map: null, selection_info: null, application_point: null};
        const n_dims_chart = this.get_chart_class().NUM_DIMS;

        if(this.options.data_source.show_filtered_points == null)
            this.options.data_source.show_filtered_points = true;

        const iterate_axes = (type) => range(axes[type].length).map(i => [type, i]);
        const d_axes = iterate_axes("design");
        const p_axes = iterate_axes("performance");
        const a_axes = iterate_axes("application");

        this.available_axes_ids = d_axes.concat(p_axes).concat(a_axes);
        if(options.data_source.type == "design")
            this.available_axes_ids = d_axes;
        if(options.data_source.type == "performance")
            this.available_axes_ids = p_axes;
        if(options.data_source.type == "application")
            this.available_axes_ids = a_axes;

        if(options.data_source.initial != null) {
            this.used_axes_ids = parse_mapping(options.data_source.initial);
            if(n_dims_chart != null && this.used_axes_ids.length != n_dims_chart)
                throw new Error(`Invalid intital mapping, expected ${n_dims_chart} axis but got ${this.used_axes_ids.length}!`);
        }
        else {
            this.used_axes_ids = n_dims_chart == null ? this.available_axes_ids : this.available_axes_ids.slice(0, n_dims_chart);
        }

        if(div != null) {
            this.chart_id = random_id();
            div.insertAdjacentHTML("beforeend", this.html_template_chart.replace("CHART-ID", this.chart_id));
        }

        if(n_dims_chart != null && this.available_axes_ids.length != n_dims_chart) {
            const sel_group_id = random_id();
            div.insertAdjacentHTML("beforeend", this.html_template_selectors.replace("SEL-ID", sel_group_id));

            const sel_group = document.getElementById(sel_group_id);
            for(let sel_idx = 0; sel_idx < n_dims_chart; sel_idx++) {
                const sel_id = random_id();
                sel_group.insertAdjacentHTML("beforeend", this.html_template_selector.replaceAll("SEL-ID", sel_id).replace("AXIS-LABEL", `Axis-${sel_idx}`));
                const sel = document.getElementById(sel_id);

                for(let a_idx = 0; a_idx < this.available_axes_ids.length; a_idx++) {
                    const [type, idx] = this.available_axes_ids[a_idx];
                    const axes_name = this.axes[type][idx].name;

                    var opt = document.createElement('option');
                    opt.value = a_idx;
                    opt.innerHTML = axes_name;
                    opt.selected = this.used_axes_ids[sel_idx][0] == type && this.used_axes_ids[sel_idx][1] == idx;
                    sel.appendChild(opt);
                }

                sel.onchange = () => {
                    const id = this.available_axes_ids[sel.value];
                    this.used_axes_ids[sel_idx] = id;

                    this.create_chart();
                    this.set_points(this.input.points);
                    this.set_filtered_points(this.input.filtered_point_indices, this.input.point_index_map, this.input.application_point);
                    this.set_selected_points(this.input.selection_info);
                };
            }
        }

        this.create_chart();
    }

    get_chart_class() {
        const plot_options = this.options.chart;
        var chart_class = null;
        if(plot_options.type == "chart2d")
            chart_class = ChartJS2DViewer;
        else if(plot_options.type == "radar")
            chart_class = RadarViewer;
        else if(plot_options.type == "pca")
            chart_class = PCAViewer;
        else if(plot_options.type == "triangle")
            chart_class = TriangleViewer;
        else if(plot_options.type == "chart3d")
            chart_class = Scatter3DViewer;
        else if(plot_options.type == "error")
            chart_class = ErrorChartViewer;
        else if(plot_options.type == "mesh")
            chart_class = MeshViewer;
        else throw new Error("Error! Unknown plot type: " + plot_options.type);
        return chart_class;
    }

    create_chart() {
        const chart_class = this.get_chart_class();

        const axes = this.used_axes_ids.map(id => this.axes[id[0]][id[1]]);
        this.chart = new chart_class(document.getElementById(this.chart_id), this.options.title, random_id(), axes, (idx) => this.onPointSelect(idx), this.options.chart);
        this.chart.resize();
    }

    map_local_idx_to_global(local_point_idx) {
        return this.options.data_source.show_filtered_points ? this.input.filtered_point_indices[local_point_idx] : local_point_idx;
    }

    map_global_idx_to_local(global_point_idx) {
        return this.options.data_source.show_filtered_points ? this.input.point_index_map[global_point_idx] : global_point_idx;
    }

    onPointSelect(local_point_idx) {
        const global_point_idx = this.map_local_idx_to_global(local_point_idx);
        this.select_clb(global_point_idx);
    }

    //takes a {design, performance, application, patch_id}
    map_point(point) {
        return {point: evaluate_mapping(this.used_axes_ids, point), attributes: point};
    }

    map_points(points) {
        return points.map(x => this.map_point(x));
    }

    set_points(points) {
        this.input.points = points;
        if(!this.options.data_source.show_filtered_points)
            this.chart.set_points(this.map_points(points));
    }

    set_filtered_points(filtered_point_indices, point_index_map, application_point) {
        this.input.filtered_point_indices = filtered_point_indices;
        this.input.point_index_map = point_index_map;
        this.input.application_point = application_point;
        if(this.options.data_source.show_filtered_points) {
            this.chart.set_points(this.map_points(filtered_point_indices.map(idx => this.input.points[idx])));
        }

        if(this.chart.set_fixed_axis_info != null) {
            const fixed_values = this.used_axes_ids.map(id => id[0] == "application" ? application_point[id[1]] : null);
            this.chart.set_fixed_axis_info(fixed_values);
        }
    }

    set_selected_points(selection_info) {
        this.input.selection_info = selection_info;

        const mapped_selection_info = map_selection_info(selection_info, idx => this.map_global_idx_to_local(idx), p => this.map_point(p));

        this.chart.set_selected_points(mapped_selection_info);
    }

    get_additional_selections(selection_info) {
        if(this.chart.get_additional_selections == null)
            return selection_info;

        const local_info = map_selection_info(selection_info, idx => this.map_global_idx_to_local(idx), p => this.map_point(p));
        this.chart.get_additional_selections(local_info);

        return map_selection_info(local_info, idx => this.map_local_idx_to_global(idx), p => p);
    }

    get_chart() {
        return this.chart;
    }
}