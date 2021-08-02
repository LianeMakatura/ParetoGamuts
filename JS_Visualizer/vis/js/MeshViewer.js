//creates the ui component
class ThreeWindow {
    use_high_quality_meshes = false;
    update_meshes_continuously = true;
    max_radius = 10;

    constructor(container) {
        this.container = container;
        THREE.Object3D.DefaultUp.set(0.0, 0.0, 1.0);
        this.renderer = new THREE.WebGLRenderer({antialias:true});
        var w = container.offsetWidth;
        var h = container.offsetHeight;
        this.renderer.setSize(w, h-10);
        container.appendChild(this.renderer.domElement);
        //Add Camera
        this.camera = new THREE.PerspectiveCamera( 75, w / h, 2, 1000 );
        this.camera.position.x = 100;
        this.camera.position.y = 100;
        this.camera.position.z = 100;
        this.scene = new THREE.Scene();
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.addEventListener( 'change', () => this.render() );

        // Configure renderer clear color
        this.renderer.setClearColor("#6495ED");

        const max_radius = this.max_radius;
        this.camera.position.set(max_radius, max_radius, max_radius);
        this.camera.near = max_radius / 1000;
        this.camera.far = max_radius * 1000;
        this.camera.updateProjectionMatrix();

        this.viewers = [];
    }

    add_standard_objects() {
        var light = new THREE.AmbientLight( 0x404040 ); // soft white light
        this.scene.add( light );

        var pointLight  = new THREE.PointLight(0xffffff, 1, 10000);
        pointLight.position.set( 1000, 1000, 1000 );
        this.scene.add( pointLight );

        var material = new THREE.MeshPhongMaterial( {
            color: 0xaaaaaa,
            polygonOffset: true,
            polygonOffsetFactor: 1, // positive value pushes polygon further away
            polygonOffsetUnits: 1
        } );
        const max_radius = this.max_radius;
        var geometry = new THREE.PlaneGeometry( max_radius*2, max_radius*2, 32, 32);
        var plane = new THREE.Mesh( geometry, material );
        plane.translateZ(-0.005);
        // this.scene.add( plane );
        // wireframe
        var geo = new THREE.WireframeGeometry( geometry ); // or WireframeGeometry
        var mat = new THREE.LineBasicMaterial( { color: 0x333333, linewidth: 1 } );
        var wireframe = new THREE.LineSegments( geo, mat );
        plane.add( wireframe );
        //coordinate system dummy
        var axesHelper = new THREE.AxesHelper( max_radius * 1.5 );
        this.scene.add( axesHelper );
    }

    render() {
        this.renderer.render(this.scene, this.camera);
    }

    resize() {
        this.camera.aspect = this.container.offsetWidth / this.container.offsetHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize( this.container.offsetWidth, this.container.offsetHeight-10 );
        this.render();
    }

    register_viewer(viewer) {
        this.viewers.push(viewer);
    }

    set_mesh_quality(use_high_quality) {
        this.viewers.forEach(x => x.load_models(use_high_quality));
    }
}

class MeshViewer {
    meshes_shown = [];

    constructor(div, chart_name, chart_id, axis_info, selection_callback, opts) {
        this.mesh_provider = opts.mesh_provider;
        this.three_window = opts.three_window;
        this.three_window.register_viewer(this);
    }

    dispose() {
        this.remove_last_mesh();
    }

    get_screenshot_data() {
        return {div: this.three_window.container};
    }

    resize() {
        this.three_window.resize();
    }

    set_points(points) {
        this.points = points;
    }

    set_selected_points(selection_info) {
        if(!this.three_window.update_meshes_continuously)
            return;
        this.last_selection_info = selection_info;
        this.load_models(false);
    }

    async load_models(use_high_quality) {
        let sprite = new THREE.TextSprite({
            fillStyle: 'red',
            fontFamily: '"Times New Roman", Times, serif',
            fontSize: 5,
            fontStyle: 'normal',
            text: "Loading Mesh",
        });
        sprite.material.depthTest = false;
        sprite.material.depthWrite = false;
        this.three_window.scene.add(sprite);
        this.three_window.render();

        // this is [[mesh]]
        const meshes_nested = await Promise.all(map_points(this.last_selection_info, true, false, this.points, async point => {
            if(point.attributes.design_point == null || point.attributes.application_point == null)
                return null;
            return await this.mesh_provider.load_meshes(point.attributes.design_point, point.attributes.application_point, {"high-quality": use_high_quality});
        }));

        this.remove_last_mesh();
        this.three_window.scene.remove(sprite);

        this.meshes_shown = meshes_nested.flat().map(mesh => {
            var bbox = new THREE.Box3().setFromObject(mesh);
            const center = new THREE.Vector3();
            bbox.getCenter(center);
            this.three_window.scene.add(mesh);

            return mesh;
        });

        this.three_window.render();
    }

    remove_last_mesh() {
        const q = this.meshes_shown;
        this.meshes_shown = [];
        q.forEach(x => {
            this.three_window.scene.remove(x);
        });
        this.three_window.render();
    }
}