class TriangleViewer {
    static NUM_DIMS = 3;
    TrianglePoints = [[0, 0], [1,0], [0.5,1]];

    html_string = `
        <canvas id="left-top-canvas">
        </canvas>
    `;

    constructor(div, chart_name, chart_id, axis_info, selection_callback, opts) {
        this.color_map = create_color_map(opts.color_map);
        this.chart_id = chart_id;

        this.selection_dataset = [];
        this.point_dataset = [];

        this.selection_callback = selection_callback;
        this.x_name = "x-drop-" + chart_id;
        this.y_name = "y-drop-" + chart_id;
        this.z_name = "z-drop-" + chart_id;
        this.axis_info = axis_info;
        var html_chart = this.html_string.replaceAll("left-top-canvas", chart_id + "-canvas");
        div.innerHTML = html_chart;

        var ctx = document.getElementById(chart_id + "-canvas").getContext('2d');
        this.chart = Chart.Scatter(ctx, {
            data: {
                datasets: [{
                    pointBackgroundColor: 'red',
                    data: [{
                        x: -1,
                        y: -1
                    }]
                }]
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

        this.triangle_dataset = [{
            pointBackgroundColor: 'red',
            showLine: true,
            fill: false,
            lineTension: 0,
            borderColor: "black",
            borderWidth: 3,
            data: [{
                x: this.TrianglePoints[0][0],
                y: this.TrianglePoints[0][1]
            },
            {
                x: this.TrianglePoints[1][0],
                y: this.TrianglePoints[1][1]
            },
            {
                x: this.TrianglePoints[2][0],
                y: this.TrianglePoints[2][1]
            },
            {
                x: this.TrianglePoints[0][0],
                y: this.TrianglePoints[0][1]
            }]
        }];
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

        var best_point_idx = null;
        var best_distance = 1e10;
        for(var idx = 0; idx < this.projected_points.length; idx++) {
            var point = this.projected_points[idx].point;
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
            this.projected_points = this.points.map(x => ({point: this.map_to_triangle(x.point)}));

            let [point_map, group_map] = this.color_map.compute_map(this.points);

            this.point_dataset = color_map_to_plotjs(this.projected_points, point_map, group_map, true);
        }

        this.update_plot();
    }

    set_selected_points(selection_info) {
        const point_color = map_points(selection_info, true, true, this.points, (point, color) => {
            return {
                point: {
                    point: this.map_to_triangle(point.point),
                    attributes: point.attributes
                },
                color
            };
        });
        let [new_points, point_map, group_map] = apply_color_map_partial(point_color, this.color_map);
        this.selection_dataset = color_map_to_plotjs(new_points, point_map, group_map, false);

        this.update_plot();
    }

    update_plot() {
        this.chart.data.datasets = this.selection_dataset.concat(this.point_dataset).concat(this.triangle_dataset);
        this.chart.update();
    }

    map_to_triangle(point) {
        //get the selected normalized coordinates
        var p = [scale_axis_to_unity(point[0], this.axis_info[0]),
                 scale_axis_to_unity(point[1], this.axis_info[1]),
                 scale_axis_to_unity(point[2], this.axis_info[2])]

        //project p onto the (1,1,1) plane
        //by intersecting the line p * lambda with the (1,1,1) plane getting
        var lambda = 1 / (p[0] + p[1] + p[2]);
        //the intersecting point p' is also its own barycentric coordinate for the unit triangle in the (1,1,1) plane
        var p_prime = [lambda * p[0], lambda * p[1], lambda * p[2]];

        //convert to 2d point in the plot by using the barycentric coordinates to find the same point in a different triangle
        return [
            this.TrianglePoints[0][0] * p_prime[0] + this.TrianglePoints[1][0] * p_prime[1] + this.TrianglePoints[2][0] * p_prime[2],
            this.TrianglePoints[0][1] * p_prime[0] + this.TrianglePoints[1][1] * p_prime[1] + this.TrianglePoints[2][1] * p_prime[2]
        ];
    }
}