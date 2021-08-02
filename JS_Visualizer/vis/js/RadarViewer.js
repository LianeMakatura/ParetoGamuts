class RadarViewer {
    html_string = `
        <canvas id="XYZ">
        </canvas>
    `;

    static NUM_DIMS = null;

    constructor(div, chart_name, chart_id, axis_info, selection_callback, opts) {
        this.chart_id = chart_id;
        this.color_map = create_color_map(opts.color_map);

        div.insertAdjacentHTML("beforeend", this.html_string.replace("XYZ", chart_id+"-canvas"));
        var ctx = document.getElementById(chart_id + "-canvas").getContext('2d');
        this.chart = new Chart(ctx, {
            type: 'radar',
            data: {
                labels: axis_info.map(x => x.name),
                datasets: []
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                legend: {
                    display: false
                },
                scale: {
                    ticks: {
                      beginAtZero: true
                    }
                },
                title: {
                    display: true,
                    text: chart_name + " Point"
                },
                animation: chartjs_animation
            }
        });

        this.chart.update();
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

    set_selected_points(selection_info) {
        this.chart.data.datasets = map_points(selection_info, true, false, this.points, (point, color) => {
            let [_, a, b] = apply_color_map_partial([{point, color}], this.color_map);
            const chart_color = b[a[0]];
            return {
                label: "Selected Point",
                data: point.point,
                backgroundColor: Chart.helpers.color(chart_color).alpha(0.2).rgbString(),
                borderColor: chart_color,
                pointBackgroundColor: chart_color,
            }
        });

        this.chart.update();
    }
}