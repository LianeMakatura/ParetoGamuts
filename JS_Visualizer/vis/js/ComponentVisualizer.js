class VisiblePointSet {
    filtered_point_indices = null;
    inv_filtered_indices = null;
    selected_point_indices = [];

    //both points can contain null
    //null means automatic acceptance in a dimension
    compare_application_points(point_a, point_b, threshold) {
        for(var i = 0; i < point_a.length; i++)
            if(point_a[i] != null && point_b[i] != null && Math.abs(point_a[i] - point_b[i]) > threshold)
                return false;
        return true;
    }

    distance_design_performance(p1, p2) {
        let a = 0;//distance_points_sqr(p1.design_point, p2.design_point);
        let b = distance_points_sqr(p1.performance_point, p2.performance_point);
        return a + b;
    }

    //sets selected point to a random new one
    //application point has to have the correct length but can contain null values
    filter_points(all_points, application_point, threshold) {
        this.filtered_point_indices = [];
        this.inv_filtered_indices = [];
        for(var point_idx = 0; point_idx < all_points.length; point_idx++) {
            if(this.compare_application_points(all_points[point_idx].application_point, application_point, threshold)) {
                this.inv_filtered_indices.push(this.filtered_point_indices.length);
                this.filtered_point_indices.push(point_idx);
            }
            else this.inv_filtered_indices.push(null);
        }

        let old_selected_point_indices = this.selected_point_indices;
        this.selected_point_indices = [];
        for(let selected_point_idx = 0; selected_point_idx < old_selected_point_indices.length; selected_point_idx++) {
            let old_point = all_points[old_selected_point_indices[selected_point_idx]];
            let best_new_point_idx = -1;

            for(let new_point_index = 0; new_point_index < this.filtered_point_indices.length; new_point_index++) {
                let best_dist = best_new_point_idx == -1 ? 10e10 : this.distance_design_performance(old_point, all_points[best_new_point_idx]);
                let new_dist = this.distance_design_performance(old_point, all_points[this.filtered_point_indices[new_point_index]]);

                if(new_dist < best_dist) {
                    best_new_point_idx = this.filtered_point_indices[new_point_index];
                }
            }

            if(best_new_point_idx != -1) {
                this.selected_point_indices.push(best_new_point_idx);
            }
        }

        // make sure that points are not selected twice
        if (this.selected_point_indices.length > 0)
            this.selected_point_indices = Array.from(new Set(this.selected_point_indices));

        // if we weren't able to pick a single point just chose the first one
        if(this.selected_point_indices.length == 0)
            this.selected_point_indices = this.filtered_point_indices.length > 0 ? [this.filtered_point_indices[0]] : [];
    }

    get_filtered_point_indices() {
        return this.filtered_point_indices;
    }

    get_filtered_point_inv_map() {
        return this.inv_filtered_indices;
    }

    get_selected_point_indices() {
        return this.selected_point_indices;
    }

    set_selected_point_indices(q) {
        this.selected_point_indices = q;
    }
}

class ComponentVisualizer {
    point_set = null;
    component_info = null;
    plots = [];

    html_string = `
        <div id="XYZ-top" class="split content" style="position: relative; display: flex; flex-direction: column;">

        </div>
        <div id="XYZ-bottom" class="split content" style="position: relative; display: flex; flex-direction: column;">

        </div>
    `;

    html_string2 = `
        <div id="XYZ" style="width: W_PER%; height: 100%; flex: 1; display: flex; flex-flow: column;">
        </div>
    `;

    constructor(parent_div, additional_display_div, component_info, select_clb) {
        this.select_clb = select_clb;
        this.parent_div = parent_div;
        this.additional_display_div = additional_display_div;
        this.point_set = new VisiblePointSet();
        this.component_info = component_info;
        const id_base = random_id();
        parent_div.innerHTML = this.html_string.replaceAll("XYZ", id_base);

        const add_options = (options) => {
            options.chart.mesh_provider = this.component_info.mesh_provider;
            options.chart.problem_path = visualizationDataFolder + this.component_info.name + "/";
            options.chart.three_window = this.component_info.three_window;
            return options;
        };

        const axes = {
            design: component_info.problem_data.design_variables,
            performance: component_info.problem_data.performance_metrics,
            application: component_info.problem_data.application_variables
        };
        const create_plot = (options, local_parent_div, n_plots_in_row) => {
            const id_plot = random_id();
            local_parent_div.lastChild.insertAdjacentHTML("beforeend", this.html_string2.replace("XYZ", id_plot).replace("W_PER", Math.floor(100 / n_plots_in_row)));
            const plot = new DataSourceControl(document.getElementById(id_plot), add_options(options), axes, (point_idx) => this.selected_point(point_idx));
            plot.set_points(this.component_info.pareto_points);
            return plot;
        };

        const get_plot_info = (variable_type, div_id_part, chart_options) => {
            const opt_ext = {
                data_source: {
                    type: variable_type
                },
                chart: chart_options,
                title: variable_type + " Space"
            };

            const div = document.getElementById(id_base + '-' + div_id_part);

            return [opt_ext, div];
        };

        //add standard plots
        const top_plot_info = component_info.problem_data.design_chart_options.map(opt => get_plot_info("design", "top", opt));
        const bottom_plot_info = component_info.problem_data.performance_chart_options.map(opt => get_plot_info("performance", "bottom", opt));

        //add additional plots
        if(component_info.problem_data.top_charts != null)
            component_info.problem_data.top_charts.forEach(opt => {
                top_plot_info.push([opt, document.getElementById(id_base + '-top')]);
            });
        if(component_info.problem_data.bottom_charts != null)
            component_info.problem_data.bottom_charts.forEach(opt => {
                bottom_plot_info.push([opt, document.getElementById(id_base + '-bottom')]);
            });

        //actually create the plots
        document.getElementById(id_base + '-top').insertAdjacentHTML("beforeend", `<div style="display: flex; height: 100%; align-items: center"></div>`);
        const top_plots = top_plot_info.map(x => create_plot(x[0], x[1], top_plot_info.length));
        document.getElementById(id_base + '-bottom').insertAdjacentHTML("beforeend", `<div style="display: flex; height: 100%; align-items: center"></div>`);
        const bottom_plots = bottom_plot_info.map(x => create_plot(x[0], x[1], bottom_plot_info.length));

        //create the extra plots if the user asked for them
        const extra_plots = [];
        if(component_info.problem_data.extra_charts != null) {
            if(additional_display_div == null)
                throw new Error("User requested extra charts but there was no space fr them");
            additional_display_div.innerHTML = component_info.problem_data.extra_charts.html_layout;

            for(const [plot_name, opt] of Object.entries(component_info.problem_data.extra_charts.charts)) {
                const plot_div = document.getElementById(plot_name);
                const parent_div = plot_div.parentElement;
                if($(parent_div).find("[name='qqq']").length == 0) {
                    parent_div.insertAdjacentHTML("beforeend", `<div name="qqq" style="display: flex; height: 100%; align-items: center"></div>`);
                }
                const num_plots_in_row = parent_div.children.length - 1;
                plot_div.remove();
                const plot = create_plot(opt, parent_div, num_plots_in_row);
                extra_plots.push(plot);
            }
        }

        //add the mesh viewer
        const mesh_plot = new DataSourceControl(null, add_options(get_plot_info("design", null, {type: "mesh"})[0]), axes, (point_idx) => this.selected_point(point_idx));
        mesh_plot.set_points(this.component_info.pareto_points);

        //save all plots
        this.plots = top_plots.concat(bottom_plots).concat([mesh_plot]).concat(extra_plots);

        Split(['#' + id_base + '-top', '#' + id_base + '-bottom'], {
            direction: 'vertical',
            sizes: [50, 50],
            gutterSize: 8,
            cursor: 'row-resize',
            onDragEnd: () => this.resize()
        });
    }

    // returns [[mesh]], that is for each selected point a list of meshes
    async get_meshes(use_high_quality = false) {
        return await Promise.all(this.point_set.get_selected_point_indices().map(async idx => {
            const point_design = this.component_info.pareto_points[idx].design_point;
            const point_application = this.component_info.pareto_points[idx].application_point;

            return await this.component_info.mesh_provider.load_meshes(point_design, point_application, {"high-quality": use_high_quality});
        }));
    }

    //returns [{design_point, application_point, performance_point}]
    get_selected_points() {
        return this.point_set.get_selected_point_indices().map(idx => {
            const design_point = this.component_info.pareto_points[idx].design_point;
            const application_point = this.component_info.pareto_points[idx].application_point;
            const performance_point = this.component_info.pareto_points[idx].performance_point;
            return {design_point, application_point, performance_point};
        });
    }

    show_selected() {
        let selection_info = this.point_set.get_selected_point_indices().map((idx, i) => new SelectionItem(idx, i == 0 ? "red" : get_distinct_color(i)));

        for(const plot of this.plots) {
            selection_info = plot.get_additional_selections(selection_info);
        }

        for(const plot of this.plots)
            plot.set_selected_points(selection_info);

        // var patch_id_str = "Patch-IDs: " + JSON.stringify(this.point_set.get_selected_point_indices().map(i => this.component_info.pareto_points[i].patch_id));
        // document.getElementById("patch-label").innerHTML = patch_id_str;
    }

    selected_point(point_idx) {
        if(point_idx != null){
            if(cntrlIsPressed)
                this.point_set.get_selected_point_indices().push(point_idx);
            else this.point_set.set_selected_point_indices([point_idx]);
        }
        else this.point_set.set_selected_point_indices([]);

        this.show_selected();

        if(this.select_clb != null)
            this.select_clb();
    }

    set_application_point(application_point, threshold) {
        this.point_set.filter_points(this.component_info.pareto_points, application_point, threshold);

        for(const plot of this.plots)
            plot.set_filtered_points(this.point_set.get_filtered_point_indices(), this.point_set.get_filtered_point_inv_map(), application_point);

        this.show_selected();
    }

    delete_visualizer() {
        for(const plot of this.plots)
            plot.dispose();
        remove_children(this.parent_div);
        if(this.additional_display_div != null)
            remove_children(this.additional_display_div);
    }

    async take_screenshots() {
        return await Promise.all(this.plots.map(async plot => {
            const data = plot.get_chart().get_screenshot_data();

            if(data.canvas != null)
                //return data.canvas; //this results in black backgrounds for chartjs
                return html2canvas(data.canvas, {scale: 4});
            else if(data.div != null)
                return html2canvas(data.div, {scale: 4});
            else throw new Error("Plot did not return screenshot data: " + data);
        }));
    }

    resize() {
        for(const plot of this.plots)
            plot.get_chart().resize();
    }
}