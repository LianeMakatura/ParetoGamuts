class Scatter3DViewer {
    static NUM_DIMS = 3;

    html_string = `
        <div id="left-top-canvas" >

        </div>
    `;

    constructor(div, chart_name, chart_id, axis_info, selection_callback, opts) {
        this.color_map = create_color_map(opts.color_map);

        var html_chart = this.html_string.replaceAll("left-top-canvas", chart_id + "-canvas");
        div.innerHTML = html_chart;
        this.chart_div = document.getElementById(chart_id + "-canvas");
        this.axis_info = axis_info;
        this.chart_name = chart_name;
        this.selected_point = null;

        this.layout = {
          title: this.chart_name,
          autosize: false,
          margin: {
            t: 0, //top margin
            l: 0, //left margin
            r: 0, //right margin
            b: 0, //bottom margin
            autosize: true,
            pad:4,
          },
          showlegend: false,
          scene: {
            aspectratio: {x: 1, y: 1, z: 1},
            xaxis: {
              title: this.axis_info[0].name,
              range: [this.axis_info[0].max, this.axis_info[0].min],
              autorange: false
            },
            yaxis: {
              title: this.axis_info[1].name,
              range: [this.axis_info[1].max, this.axis_info[1].min],
              autorange: false
            },
            zaxis: {
              title: this.axis_info[2].name,
              range: [this.axis_info[2].min, this.axis_info[2].max],
              autorange: false
            }
          }
        };

        var data = [
            //points
            {
                type: 'scatter3d',
                mode: 'markers',
                x: [0],
                y: [0],
                z: [0],
                marker: {
                    size: 3.5,
                    color: []
                }
            },
            //selected points
            {
                type: 'scatter3d',
                mode: 'markers',
                marker: {
                    size: 8,
                    color: []
                }
            },
            //highlight points
            {
                type: 'scatter3d',
                mode: 'markers',
                marker: {
                    size: 5,
                    color: []
                }
            },
            //fixed axis points
            {

            }
        ];

        var config = {
            displayModeBar: false,
        };

        Plotly.newPlot(this.chart_div, data, this.layout, config);
    }

    //shows either a point, line or rectangle based on how many application variables are shown
    set_fixed_axis_info(axis_values) {
        let update = {};

        const fixed_axis = axis_values.filter(x => x != null).length;
        if(fixed_axis == 3) {
            update = {
                type: 'scatter3d',
                mode: 'markers',
                x: [[axis_values[0]]],
                y: [[axis_values[1]]],
                z: [[axis_values[2]]],
                marker: {
                    size: 5,
                    color: "black"
                }
            }
        }
        else if(fixed_axis == 2) {
            const free_idx = axis_values[0] == null ? 0 : (axis_values[1] == null ? 1 : 2);
            const fixed_idx_1 = (free_idx + 1) % 3;
            const fixed_idx_2 = (free_idx + 2) % 3;

            const create_point = v => {
                const p = [null, null, null];
                p[free_idx] = v;
                p[fixed_idx_1] = axis_values[fixed_idx_1];
                p[fixed_idx_2] = axis_values[fixed_idx_2];
                return p;
            };

            const p1 = create_point(0), p2 = create_point(1);

            update = {
                type: 'scatter3d',
                mode: 'lines',
                x: [[p1[0], p2[0]]],
                y: [[p1[1], p2[1]]],
                z: [[p1[2], p2[2]]],
                marker: {
                    size: 5,
                    color: "black"
                }
            }
        }
        else if(fixed_axis == 1) {
            const fixed_idx = axis_values[0] != null ? 0 : (axis_values[1] != null ? 1 : 2);
            const free_idx_1 = (fixed_idx + 1) % 3;
            const free_idx_2 = (fixed_idx + 2) % 3;
            const create_point = (x,y) => {
                const p = [null,null,null];
                p[fixed_idx] = axis_values[fixed_idx];
                p[free_idx_1] = x;
                p[free_idx_2] = y;
                return p;
            };
            const points = [
                create_point(0, 0),
                create_point(1, 0),
                create_point(1, 1),
                create_point(0, 1),
            ];
            update = {
                type: 'scatter3d',
                mode: 'lines',
                x: [[]],
                y: [[]],
                z: [[]],
                marker: {
                    size: 5,
                    color: "black"
                }
            }
            const qq = ["x", "y", "z"];
            for(let idx = 0; idx < 4; idx++) {
                for(let c = 0; c < 3; c++) {
                    update[qq[c]][0].push(points[idx][c]);
                    update[qq[c]][0].push(points[(idx + 1) % 4][c]);
                }
            }
        }

        Plotly.restyle(this.chart_div, update, 3);
    }

    dispose() {

    }

    get_screenshot_data() {
        return {div: this.chart_div};
    }

    resize() {
        this.layout.width = this.chart_div.parentNode.getBoundingClientRect().width;
        this.layout.height = this.chart_div.parentNode.getBoundingClientRect().height;
        Plotly.relayout(this.chart_div, this.layout);
    }

    set_points(points) {
        this.points = points;

        var x = points != null ? points.map(x => x.point[0]) : [];
        var y = points != null ? points.map(x => x.point[1]) : [];
        var z = points != null ? points.map(x => x.point[2]) : [];

        let [point_map, group_map] = points != null ? this.color_map.compute_map(points) : [[], []];
        var c = point_map.map((v, i) => group_map[v]);

        var update = {
            x: [x],
            y: [y],
            z: [z],
            'marker.color' : [c],
        };

        Plotly.restyle(this.chart_div, update, 0);
    }

    set_selected_points(selection_info) {
        const f = (include_selected, include_highlight) => {
            var x = [], y = [], z = [], c = [];
            iterate_points(selection_info, include_selected, include_highlight, this.points, (point, color) => {
                let [_, a, b] = apply_color_map_partial([{point, color}], this.color_map);
                const chart_color = b[a[0]];
                x.push(point.point[0]);
                y.push(point.point[1]);
                z.push(point.point[2]);
                c.push(chart_color);
            });
            return {
                x: [x],
                y: [y],
                z: [z],
                'marker.color' : [c],
            };;
        };

        const selected_update = f(true, false);
        const highlight_update = f(false, true);

        Plotly.restyle(this.chart_div, selected_update, 1);
        Plotly.restyle(this.chart_div, highlight_update, 2);
    }
}