var application_variables = null;

var component_data = null;
var component_visualizers = null;
var custom_center_div = null;
var custom_scene_generator = null;
var three_window = null;

var cntrlIsPressed = false;

function resize_renderer_split() {
    for(var c_idx in component_visualizers) {
        component_visualizers[c_idx].resize();
    }
}

function get_application_point() {
    var application_point = Array(application_variables.length).fill(null);
    for(var idx = 0; idx < application_point.length; idx++) {
        application_point[idx] = parseFloat($("#z-slider-" + idx)[0].value);
    }
    return application_point;
}

function change_value_application() {
    var application_point = get_application_point();

    for(var c_idx in component_visualizers) {
        component_visualizers[c_idx].set_application_point(application_point, parseFloat($("#z-threshold-input")[0].value));
    }

    for(var idx = 0; idx < application_point.length; idx++) {
        $("#z-slider-label-" + idx)[0].innerHTML = application_variables[idx].name + ": " + Math.round(application_point[idx] * 100) / 100;
    }

    selected_point_callback();

    if(custom_scene_generator != null) {
        custom_scene_generator.update(application_point);
    }
}

async function load_component_data(folder_path) {
    var problem_res = await fetch(visualizationDataFolder + folder_path + "problem.json");
    var problem_data = await problem_res.json();

    var fronts_res = await fetch(visualizationDataFolder + folder_path + "fronts.json");
    var fronts = await fronts_res.json();
    var nd = problem_data.design_variables.length, na = problem_data.application_variables.length, np = problem_data.performance_metrics.length;
    var pareto_points = fronts.map(x => ({design_point: x.slice(0, nd), application_point: x.slice(nd, nd + na), performance_point: x.slice(nd + na, nd + na + np), patch_id: x[nd + na + np]}))

    var mesh_provider = await create_mesh_provider(visualizationDataFolder + folder_path, problem_data.mesh_provider);

    return {name: folder_path.replace(".json", ""), problem_data, pareto_points, mesh_provider, three_window: three_window};
}

async function selected_point_callback() {
    if (custom_center_div == null)
        return;

    let selection_info = [];
    for (let i = 0; i < component_visualizers.length; i++) {
        let selected_points = component_visualizers[i].get_selected_points();
        let nested_meshes = null;
        if (custom_center_div.needs_meshes != undefined && custom_center_div.needs_meshes()) {
            nested_meshes = await component_visualizers[i].get_meshes();
        }
        let selected_points_info = selected_points.map((x,i) => {
            x.meshes = nested_meshes == null ? null : nested_meshes[i];
            return x;
        });
        selection_info.push({component_name: component_data[i].name, selected_points_info});
    }

    custom_center_div.set_points(selection_info);
}

function component_selector(selector_div, display_div, additional_display_div, hash_name, idx, create_header) {
    var s_idx = hash_name + "-component-selector";

    //gets the index currently selected (even if no header was created)
    var f_get_component_idx = () => {
        if(create_header)
            return parseInt($("#" + s_idx)[0].value);
        else return 0;
    }

    var f_create_visualizer = () => {
        component_visualizers[idx].delete_visualizer();
        component_visualizers[idx] = new ComponentVisualizer(display_div, additional_display_div, component_data[f_get_component_idx()], selected_point_callback);
        component_visualizers[idx].set_application_point(get_application_point(), parseFloat($("#z-threshold-input")[0].value));
    };

    if(create_header) {
        html_string = `
        <label for="XYZ" style="margin:0;">Component: </label> <select id="XYZ"></select>
        `;

        selector_div.innerHTML = html_string.replaceAll("XYZ", s_idx);
        for(var c_idx = 0; c_idx < component_data.length; c_idx++) {
            var opt = document.createElement('option');
            opt.value = c_idx;
            opt.innerHTML = component_data[c_idx].name;
            $("#" + s_idx)[0].appendChild(opt);
        }

        $("#" + s_idx)[0].onchange = f_create_visualizer;
        //select a component
        $("#" + s_idx)[0].value = idx;
    }

    component_visualizers[idx] = new ComponentVisualizer(display_div, additional_display_div, component_data[idx], selected_point_callback);
}

function create_split_UI(num_columns) {
    var cols = ['#left-column', '#middle-column'];
    var sizes = [60, 40];

    if(num_columns == 2) {
        cols.push('#right-column');
        sizes = [33, 34, 33];
    }

    Split(cols, {
        gutterSize: 8,
        cursor: 'col-resize',
        sizes: sizes,
        onDragEnd: () => resize_renderer_split()
    });

    let space_middle_top = custom_center_div == null ? 12 : 22;

    Split(['#middle-top', '#middle-bottom'], {
        direction: 'vertical',
        sizes: [space_middle_top, 100 - space_middle_top],
        gutterSize: 8,
        cursor: 'row-resize',
        onDragEnd: () => resize_renderer_split()
    });
}

async function download_meshes(high_quality = true) {
    const replace_invalid_chars = filename => filename.replace(/[/\\?%*:|"<>]/g, '-');

    // we need to do two flattens. One over the outer array of components.
    // one over the inner array of selected points
    var stl_name_data = await Promise.all(component_visualizers.map(async (c,i) => {
        const nested_meshes = await c.get_meshes(high_quality);
        return nested_meshes.map((meshes, selected_point_idx) => {
            return meshes.map((mesh, mesh_idx) => {
                const suffix = `-${selected_point_idx}-${mesh_idx}.stl`;
                const name = c.component_info.name ? replace_invalid_chars(c.component_info.name) + suffix
                                                   : i + suffix;
                var exporter = new THREE.STLExporter();
                var stl_data = exporter.parse( mesh, { binary: true } );
                return [name, stl_data];
            });
        });
    }));
    stl_name_data = stl_name_data.flat(2);

    var download_data = null;
    if(stl_name_data.length == 1) {
        download_data = [stl_name_data[0][0], new Blob([stl_name_data[0][1].buffer], {type: 'application/octet-stream'})];
    }
    else {
        var zip = new JSZip();
        for(var idx = 0; idx < stl_name_data.length; idx++) {
            zip.file(stl_name_data[idx][0], stl_name_data[idx][1].buffer);
        }
        const file = await zip.generateAsync({type: "blob"});
        download_data = ["meshes.zip", file];
    }

    download_blob(download_data[0], download_data[1]);
}

async function take_screenshots() {
    const files = [];
    for(var component_idx = 0; component_idx < component_visualizers.length; component_idx++) {
        const screenshots = await component_visualizers[component_idx].take_screenshots();
        screenshots.forEach((canvas, chart_idx) => {
            const blob = canvas.toDataURL('image/jpeg', 1.0);
            const filename = `${component_idx}-${chart_idx}.jpg`;
            files.push([filename, blob]);
        });
    }

    var zip = new JSZip();
    for(var idx = 0; idx < files.length; idx++) {
        zip.file(files[idx][0], files[idx][1].split('base64,')[1], {base64: true});
    }
    const file = await zip.generateAsync({type: "blob"});

    download_blob("screnshots.zip", file);
}

async function initialize_custom_center_div(parent_div) {
    try {
        await load_js_file(visualizationDataFolder + "/ui.js");
    }
    catch (e) {
        return;
    }

    custom_center_div = new CustomCenterDiv(parent_div);
}

async function initialize_custom_scene_generator() {
    try {
        await load_js_file(visualizationDataFolder + "/scene_generator.js");
    }
    catch (e) {
        return;
    }

    custom_scene_generator = new CustomSceneGenerator();
    custom_scene_generator.init(visualizationDataFolder, three_window.scene);
}

async function initialize_UI() {
    //load the root component
    var problem_res = await fetch(visualizationDataFolder + "problem.json");
    var global_problem_data = await problem_res.json();

    //create the extra div content if the problem specifies it
    //wait for this because code further down checks the global variables set by it
    await initialize_custom_center_div(document.getElementById('middle-middle'));

    //create the three js window
    three_window = new ThreeWindow(document.getElementById('middle-bottom'));

    //this is either a proper component or only a list of filenames
    component_data = [];
    if(Array.isArray(global_problem_data)) {
        create_split_UI(2);
        component_data = global_problem_data.map(async x => load_component_data(x + "/"));
        component_data = await Promise.all(component_data);
    }
    else {
        let num_columns = 1;
        if(typeof  OVERWRITE_NUM_COMPONENTS_UI !== 'undefined')
            num_columns = OVERWRITE_NUM_COMPONENTS_UI;
        create_split_UI(num_columns);
        component_data.push(await load_component_data(""));
    }

    //the application variables for each component have to be the same!
    application_variables = component_data[0].problem_data.application_variables;

    //this needs the three js scene to be available
    initialize_custom_scene_generator().then(() => {
        // check if we should add the standard stuff
        const add_standard_objects = custom_scene_generator == null ||
                                     custom_scene_generator.add_standard_objects == null ||
                                     custom_scene_generator.add_standard_objects();
        if(add_standard_objects)
            three_window.add_standard_objects();
    });

    //only when we have multiple components create both ui elements
    if(component_data.length == 1) {
        component_visualizers = [null];
        component_selector($("#left-header")[0], $("#left-content")[0], $("#right-content")[0], "left", 0, false);
    }
    else {
        component_visualizers = [null, null];
        const sheight = "30px";
        $("#left-header")[0].setAttribute("style","height:"+sheight);
        $("#right-header")[0].setAttribute("style","height:"+sheight);
        $("#left-content")[0].setAttribute("style",`height:calc(100% - ${sheight})`);
        $("#right-content")[0].setAttribute("style",`height:calc(100% - ${sheight})`);
        component_selector($("#left-header")[0], $("#left-content")[0], null, "left", 0, true);
        component_selector($("#right-header")[0], $("#right-content")[0], null, "right", 1, true);
    }

    for(let i = 0; i < application_variables.length; i++) {
        var id = "z-slider-" + i;
        var slider = document.createElement('input');
        slider.id = id;
        slider.type = 'range';
        slider.min = application_variables[i].min;
        slider.max = application_variables[i].max;
        slider.value = 0.0;
        slider.step = (slider.max - slider.min) / (parseInt($("#z-num-slider-stops")[0].value) - 1);
        slider.onchange = () => change_value_application();
        var newDiv = document.createElement("div");
        newDiv.innerHTML = '<label id="z-slider-label-' + i + '" for="' + id + '" style="margin:0;">' + application_variables[i].name + '</label>';
        $("#z-slider-group")[0].appendChild(newDiv.childNodes[0]);
        $("#z-slider-group")[0].appendChild(slider);
    }

    var z_range = Math.max(...application_variables.map(x => x.max)) - Math.min(...application_variables.map(x => x.min));
    $("#z-threshold-input")[0].onchange = () => change_value_application();
    $("#z-threshold-input")[0].step = z_range / 1000;

    $(document).on("keypress", function (e) {
        // use e.which
        if(e.key.toLowerCase() === 'h') {
            three_window.set_mesh_quality(e.shiftKey);
        }
        if(e.key.toLowerCase() === 'd') {
            download_meshes();
        }
        if(e.key.toLowerCase() === 'i') {
            take_screenshots();
        }
    });

    $(document).keydown(function(event) {
        if (event.which == "17" || event.which == "91")
            cntrlIsPressed = true;
    });

    $(document).keyup(function(event) {
        if (event.which == "17" || event.which == "91")
            cntrlIsPressed = false;
    });

    var is_checked_str = localStorage.getItem('automatic-mesh-reload');
    if(is_checked_str != null && is_checked_str != undefined)
    {
        let v = JSON.parse(is_checked_str).value;
        $("#ck-automatic-mesh-load")[0].checked = v;
        three_window.update_meshes_continuously = v;
    }
    $("#ck-automatic-mesh-load")[0].onclick = () => {
        let v = $("#ck-automatic-mesh-load")[0].checked;
        three_window.update_meshes_continuously = v;
        localStorage.setItem('automatic-mesh-reload', JSON.stringify({value: v}));
    };

    change_value_application();
    resize_renderer_split();
    window.onresize = resize_renderer_split;
}

window.onload = () => {
    initialize_UI();
};