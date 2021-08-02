class PCAViewer {
    static NUM_DIMS = 2;

    constructor(div, chart_name, chart_id, axis_info, selection_callback, opts) {
        var pca_axis_info = [...Array(PCAViewer.NUM_DIMS).keys()].map(i => ({name: "Dim " + i, min: -10, max: 10}));
        this.chart = new ChartJS2DViewer(div, chart_name, chart_id, pca_axis_info, selection_callback, opts);
    }

    dispose() {

    }

    get_screenshot_data() {
        return this.chart.get_screenshot_data();
    }

    resize() {
        this.chart.resize();
    }

    set_points(points) {
        if(points.length > 0) {
            var point_values = points.map(x => x.point);
            this.pca_vectors = PCA.getEigenVectors(point_values).slice(0, PCAViewer.NUM_DIMS);
            var adData = PCA.computeAdjustedData(point_values, ...this.pca_vectors);
            this.pca_avg_data = adData.avgData[0];
        }

        this.chart.set_points(points.map(x => ({point: this.map_to_pca(x.point), attributes: x.attributes})));
    }

    set_selected_points(selection_info) {
        this.chart.set_selected_points(map_selection_info(selection_info, null, p => ({point: this.map_to_pca(p.point), attributes: p.attributes})));
    }

    map_to_pca(std_point) {
        var vectors = this.pca_vectors.map(function(v){return v.vector});
        var matrixMinusMean = PCA.subtract([std_point], [this.pca_avg_data]);
        var adjustedData = PCA.multiply(vectors, PCA.transpose(matrixMinusMean));
        return PCA.transpose(adjustedData)[0];
    }
}