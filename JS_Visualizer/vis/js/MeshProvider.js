class MeshFileProvider {
    constructor(data_folder, options) {
        this.data_folder = data_folder;
        this.list_file_path = options.list_file == null ? "meshes-index.json" : options.list_file;
    }

    async load() {
        let res = await fetch(this.data_folder + this.list_file_path);
        this.meshes_list = await res.json();
    }

    async load_meshes(point_design, point_application, options) {
        //this searches through the list of meshes, not the front points!
        let best_distance = 1e10;
        let best_idx = -1;
        for(var i = 0; i < this.meshes_list.length; i++) {
            const d = distance_points_sqr(this.meshes_list[i].design, point_design) + distance_points_sqr(this.meshes_list[i].application, point_application);
            if(d < best_distance) {
                best_distance = d;
                best_idx = i;
            }
        }

        if(best_idx == -1)
            return null;

        const filename = this.meshes_list[best_idx].filename;
        let mesh_res = await fetch(data_folder + filename);
        let mesh_data = await mesh_res.text();
        let loader = file_extension(filename).toLowerCase() == "stl" ? new THREE.LegacyJSONLoader() : new THREE.OBJLoader();
        let object = loader.parse( mesh_data );
        // extract meshes from the object
        let meshes = [];
        object.traverse(child => {
            if(child instanceof THREE.Mesh) {
                meshes.push(child);
            }
        });
        return meshes;
    }
}

class JscadMeshProvider {
    constructor(data_folder, options) {
        this.data_folder = data_folder;
        this.template_file_path = options.template_file == null ? "mesh_template.jscad" : options.list_file;
    }

    async load() {
        var jscad_res = await fetch(this.data_folder + this.template_file_path);
        this.jscad_template = await jscad_res.text();
    }

    async load_meshes(point_design, point_application, options) {
        var json_text = "var dp = " + JSON.stringify(point_design) + ";\nvar ap = " + JSON.stringify(point_application) + ";";
        if("high-quality" in options && options["high-quality"])
            json_text = json_text + "\ng_fn = 10;";
        var jscad_evaluated_code = this.jscad_template.replace("//INSERT_PARAMETER_VALUES", json_text);
        var xhr = new XMLHttpRequest();
        xhr.overrideMimeType('text/plain; charset=x-user-defined');
        xhr.responseType = "arraybuffer";
        xhr.open("POST", JSCAD_SERVER_URL, true);

        let promise = new Promise((resolve, reject) => {
            xhr.onload = function(event) {
                let mesh_data = event.target.response || event.target.responseText;
                let loader = new THREE.STLLoader();
                let buffer_geometry = loader.parse( mesh_data ); //specular: 0x111111, shininess: 200
                let geometry = new THREE.Geometry().fromBufferGeometry(buffer_geometry);
                geometry.computeFlatVertexNormals();
                let material = new THREE.MeshPhongMaterial({
                    color: new THREE.Color('#54668e'),
                    specular: new THREE.Color('#fff'),
                    opacity: 0.95,
                    reflectivity: 0.6,
                  });
                let mesh = new THREE.Mesh( geometry, material );
                resolve([mesh]);
            };
            xhr.send(jscad_evaluated_code);
        });
        return promise;
    }
}

class JSMeshGeneratorProvider {
    mesh_generator = null;

    constructor(data_folder, options) {
        this.data_folder = data_folder;
        this.js_file_path = options.js_file == null ? "js_mesh_generator.js" : options.list_file;
    }

    async load() {
        await load_js_file(this.data_folder + this.js_file_path);
        this.mesh_generator = new JSMeshGenerator();
    }

    async load_meshes(point_design, point_application, options) {
        return this.mesh_generator.generate_mesh(point_design, point_application);
    }
}

async function create_mesh_provider(data_folder, options) {
    let provider = null;
    if(options.type == "files")
        provider = new MeshFileProvider(data_folder, options);
    else if(options.type == "jscad")
        provider =  new JscadMeshProvider(data_folder, options);
    else if(options.type == "js_mesh_generator")
        provider = new JSMeshGeneratorProvider(data_folder, options);
    else throw new Error("Invalid mesh provider type: " + options);
    await provider.load();
    return provider;
}