var gfopt = {
	zoom: {
		// Container for pan options
		pan: {
			// Boolean to enable panning
			enabled: true,

			// Panning directions. Remove the appropriate direction to disable
			// Eg. 'y' would only allow panning in the y direction
			mode: 'xy',

			rangeMin: {
				// Format of min pan range depends on scale type
				x: null,
				y: null
			},
			rangeMax: {
				// Format of max pan range depends on scale type
				x: null,
				y: null
			}
		},

		// Container for zoom options
		zoom: {
			// Boolean to enable zooming
			enabled: true,

			// Enable drag-to-zoom behavior
			drag: false,

			// Drag-to-zoom rectangle style can be customized
			// drag: {
			// 	 borderColor: 'rgba(225,225,225,0.3)'
			// 	 borderWidth: 5,
			// 	 backgroundColor: 'rgb(225,225,225)'
			// },

			// Zooming directions. Remove the appropriate direction to disable
			// Eg. 'y' would only allow zooming in the y direction
			mode: 'xy',

			rangeMin: {
				// Format of min zoom range depends on scale type
				x: null,
				y: null
			},
			rangeMax: {
				// Format of max zoom range depends on scale type
				x: null,
				y: null
			},

			// Speed of zoom via mouse wheel
			// (percentage of zoom on a wheel event)
			speed: 0.1
		}
	}
};

var chartjs_animation = {
    duration: 0
};

class ChartJS2DViewer {
    static NUM_DIMS = 2;

    html_string = `
        <canvas id="left-top-canvas">
        </canvas>
    `;

    constructor(div, chart_name, chart_id, axis_info, selection_callback, opts) {
        this.chart_id = chart_id;
        this.color_map = create_color_map(opts.color_map);

        this.selection_dataset = [];
        this.point_dataset = [];

        this.selection_callback = selection_callback;
        this.axis_info = axis_info;
        var html_chart = this.html_string.replaceAll("left-top-canvas", chart_id + "-canvas");
        div.innerHTML = html_chart;

        var ctx = document.getElementById(chart_id + "-canvas").getContext('2d');
        this.chart = Chart.Scatter(ctx, {
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
                    text: chart_name
                },
              onClick: (a,b) => this.plot_click(a,b),//have to use the lambda, otherwise this will be different
              plugins: gfopt,
              animation: chartjs_animation
            }
        });

        var var_x = this.axis_info[0];
        var var_y = this.axis_info[1];

        this.chart.options.scales.xAxes[0].ticks.min = var_x.min;
        this.chart.options.scales.xAxes[0].ticks.max = var_x.max;
        this.chart.options.scales.xAxes[0].scaleLabel = { display: true, labelString: var_x.name };
        this.chart.options.scales.yAxes[0].ticks.min = var_y.min;
        this.chart.options.scales.yAxes[0].ticks.max = var_y.max;
        this.chart.options.scales.yAxes[0].scaleLabel = { display: true, labelString: var_y.name };

        this.app_vis_dataset = [];
    }

    set_fixed_axis_info(axis_values) {
        const fixed_axis = axis_values.filter(x => x != null).length;
        if(fixed_axis == 2) {
            this.app_vis_dataset = {
                data: [{x: axis_values[0], y: axis_values[1]}],
                showLine: false,
                fill: false,
                pointBackgroundColor: "black"
            };
        }
        else if(fixed_axis == 1) {
            const fixed_idx = axis_values[0] != null ? 0 : 1;
            const p1 = [0,0], p2 = [1,1];
            p1[fixed_idx] = axis_values[fixed_idx];
            p2[fixed_idx] = axis_values[fixed_idx];
            this.app_vis_dataset = {
                data: [{x: p1[0], y: p1[1]}, {x: p2[0], y: p2[1]}],
                showLine: true,
                fill: false,
                pointBackgroundColor: "black"
            };
        }
    }

    dispose() {

    }

    get_screenshot_data() {
        return {canvas: document.getElementById(this.chart_id + "-canvas")};
    }

    resize() {
        this.chart.resize();
    }

    plot_click(a, b) {
        var valueX = null, valueY = null;
        for (var scaleName in this.chart.scales) {
            var scale = this.chart.scales[scaleName];
            if (scale.isHorizontal()) {
                valueX = scale.getValueForPixel(event.offsetX);
            } else {
                valueY = scale.getValueForPixel(event.offsetY);
            }
        }

        //build complete point out of click coordinates
        var best_point_idx = null;
        var best_distance = 1e10;
        for(var idx = 0; idx < this.points.length; idx++) {
            var point = this.points[idx].point;
            var px = point[0];
            var py = point[1];
            var distance = (valueX - px) * (valueX - px) + (valueY - py) * (valueY - py);
            if(distance < best_distance) {
                best_distance = distance;
                best_point_idx = idx;
            }
        }

        if(best_point_idx != null)
            this.selection_callback(best_point_idx);
    }

    set_points(points) {
        this.points = points;

        this.point_dataset = [];
        if(this.points != null) {
            let [point_map, group_map] = this.color_map.compute_map(this.points);

            this.point_dataset = color_map_to_plotjs(this.points, point_map, group_map, true);
        }

        this.update_plot();
    }

    set_selected_points(selection_info) {
        const point_color = map_points(selection_info, true, true, this.points, (point, color) => ({point, color}));
        let [new_points, point_map, group_map] = apply_color_map_partial(point_color, this.color_map);
        this.selection_dataset = color_map_to_plotjs(new_points, point_map, group_map, false);

        this.update_plot();
    }

    update_plot() {
        this.chart.data.datasets = this.selection_dataset.concat(this.point_dataset).concat([this.app_vis_dataset]);
        this.chart.update();
    }
}