class ErrorChartViewer {
    html_string = `
        <canvas id="left-top-canvas">
        </canvas>
    `;

    static NUM_DIMS = null;

    constructor(div, chart_name, chart_id, axis_info, selection_callback, opts) {
        this.chart_id = chart_id;
        this.color_map = create_color_map(opts.color_map);
        this.error_mode = opts.error_mode != null ? opts.error_mode : "absolute";
        //absolute for |x-x*|, relative for |(x-x*)/(1/2*(x+x*))|, nadier for ...
        this.autoscale_values = opts.autoscale_values != null ? opts.autoscale_values : false;

        var html_chart = this.html_string.replaceAll("left-top-canvas", chart_id + "-canvas");
        div.innerHTML = html_chart;

        var ctx = document.getElementById(chart_id + "-canvas").getContext('2d');
        this.chart = Chart.Scatter(ctx, {
            type: 'scatter',
            data: {
                datasets: []
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                legend: {
                    display: false
                },
                title: {
                    display: true,
                    text: "Percentage of Optimal Performance"
                },
                animation: chartjs_animation
            }
        });

        this.chart.options.scales.xAxes[0].ticks.min = 0;
        this.chart.options.scales.xAxes[0].ticks.max = 1;
        this.chart.options.scales.xAxes[0].scaleLabel = { display: true, labelString: "Context" };
        this.chart.options.scales.yAxes[0].ticks.min = 0;
        this.chart.options.scales.yAxes[0].ticks.max = 100
        this.chart.options.scales.yAxes[0].scaleLabel = { display: true, labelString: "Performance (as % of optimal)" };

        //load error data
        this.error_data = null;
        fetch(opts.problem_path + opts.data_file).then(async res => {
            this.error_data = await res.json();
        });
    }

    dispose() {

    }

    get_screenshot_data() {
        return {canvas: document.getElementById(this.chart_id + "-canvas")};
    }

    resize() {
        this.chart.resize();
    }

    set_points(points) {
        this.points = points;
    }

    get_closest_point(design_point) {
        var best_idx = null;
        var best_distance = 1e5;
        for(var idx = 0; idx < this.error_data.length; idx++) {
            const d = distance_points_sqr(design_point, this.error_data[idx].designPt);
            if(d < best_distance) {
                best_idx = idx;
                best_distance = d;
            }
        }

        return best_idx;
    }

    compute_value(x, x_star) {
        if(this.error_mode == "absolute") {
            return Math.sqrt(distance_points_sqr(x, x_star));
        }
        else if(this.error_mode == "relative") {
            const nominator = Math.sqrt(distance_points_sqr(x, x_star));
            const denominator = 0.5 * Math.sqrt(x.reduce((total, _, i) => total + (x[i] + x_star[i]) * (x[i] + x_star[i]), 0));
            return nominator / denominator;
        }
        else if(this.error_mode == "nadier") {
            const volume_f = (p) => p.reduce((total, _, i) => total * (1 - p[i]), 1);
            const volume_optimal = volume_f(x_star);
            const volume_pint = volume_f(x);
            return volume_pint / volume_optimal * 100;
        }
        else throw new Error(`Invalid error mode: ${this.error_mode}`);
    }

    compute_values(entry) {
        //perfEvalsPerContext, closestPtPerContext
        return entry.closestPtPerContext.map((_, i) => this.compute_value(entry.perfEvalsPerContext[i], entry.closestPtPerContext[i]));
    }

    set_selected_points(selection_info) {
        this.chart.data.datasets = [];

        if(this.error_data != null) {
            let scale_factor = 1;
            if(this.autoscale_values) {
                scale_factor = 0;
                iterate_points(selection_info, true, false, this.points, (point, _) => {
                    const best_idx = this.get_closest_point(point.attributes.design_point);
                    if(best_idx != null) {
                        const y_values = this.compute_values(this.error_data[best_idx]);
                        scale_factor = Math.max(scale_factor, Math.max(...y_values));
                    }
                });
                scale_factor = 1 / scale_factor * 0.8;
            }

            this.chart.data.datasets = map_points(selection_info, true, false, this.points, (point, color) => {
                const best_idx = this.get_closest_point(point.attributes.design_point);

                if(best_idx != null) {
                    const y_values = this.compute_values(this.error_data[best_idx]);
                    const x_values = makeArr(0, 1, y_values.length);
                    const points = range(y_values.length).map(i => ({x: x_values[i], y: y_values[i] * scale_factor}));

                    let [_, a, b] = apply_color_map_partial([{point, color}], this.color_map);
                    const chart_color = b[a[0]];
                    return {
                        label: "Scaled Error",
                        data: points,
                        backgroundColor: Chart.helpers.color(chart_color).alpha(0.2).rgbString(),
                        borderColor: chart_color,
                        pointBackgroundColor: chart_color,
                        showLine: true,
                        fill: false,
                    };
                }
                else return null;
            });
        }

        this.chart.update();
    }

    get_additional_selections(selection_info) {
        if(this.error_data == null)
            return [];

        iterate_points(selection_info, true, false, this.points, (point, color, selection_item) => {
            if(point.attributes.design_point != null) {
                const best_idx = this.get_closest_point(point.attributes.design_point);

                if(best_idx != null) {
                    const group_points = [];
                    const performance_points = this.error_data[best_idx].perfEvalsPerContext;
                    const app_values = makeArr(0, 1, performance_points.length);
                    for(var idx = 0; idx < performance_points.length; idx++) {
                        group_points.push(make_highlight_point({
                            design_point: this.error_data[best_idx].designPt,
                            performance_point: performance_points[idx],
                            application_point: [app_values[idx]]
                        }, null, false));
                    }
                    selection_item.highlighted_groups.push(new HighlightGroup(group_points));
                }
            }
        });
    }
}