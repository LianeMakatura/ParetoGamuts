class CustomSceneGenerator {
    Ni = 10;
    numFaces = this.Ni - 1;
    scalefactor = 1;
    roof_vis_scale_factor = 12
    houseInset = 1;

    Init() {
    }

    async init(data_folder, scene) {
        var mtlLoader = new THREE.MTLLoader();
        mtlLoader.setPath( data_folder + 'meshes/' );
        var url = "modern_house_base_v5_lowpoly.mtl";
        mtlLoader.load( url, materials => {
            materials.preload();
            var objLoader = new THREE.OBJLoader();
            for (const [key, value] of Object.entries(materials.materials)) {
                value.side = THREE.DoubleSide;
            }
            objLoader.setMaterials( materials );
            objLoader.load(
                // resource URL
                data_folder + "meshes/modern_house_base_v5_lowpoly.obj",
                (o) => {
                    scene.add( o );
                    this.house = o;
                },
                // called when loading is in progresses
                function ( xhr ) {
                    console.log( ( xhr.loaded / xhr.total * 100 ) + '% loaded' );
                },
                // called when loading has errors
                function ( error ) {
                    console.log( 'An error happened' );
                }
            );
        });
        // var mtlLoader = new THREE.MTLLoader();
        // mtlLoader.setPath( data_folder + 'meshes/' );
        var url = "morning_arrows.mtl";
        mtlLoader.load( url, materials => {
            materials.preload();
            var objLoader = new THREE.OBJLoader();
            for (const [key, value] of Object.entries(materials.materials)) {
                value.side = THREE.DoubleSide;
            }
            objLoader.setMaterials( materials );
            objLoader.load(
                // resource URL
                data_folder + "meshes/sun_arrows.obj",
                (o) => {
                    var vector = new THREE.Vector3(-0.5, -0.5, -0.6);
                    var axis = new THREE.Vector3(0, 1, 0);
                    o.quaternion.setFromUnitVectors(axis, vector.clone().normalize());
                    o.position.copy(vector.clone().multiplyScalar(-200));
                    scene.add( o );
                    this.morning_arrows = o;
                },
                // called when loading is in progresses
                function ( xhr ) {
                    console.log( ( xhr.loaded / xhr.total * 100 ) + '% loaded' );
                },
                // called when loading has errors
                function ( error ) {
                    console.log( 'An error happened' );
                }
            );
        });
        var url = "evening_arrows.mtl";
        mtlLoader.load( url, materials => {
            materials.preload();
            var objLoader = new THREE.OBJLoader();
            for (const [key, value] of Object.entries(materials.materials)) {
                value.side = THREE.DoubleSide;
            }
            objLoader.setMaterials( materials );
            objLoader.load(
                // resource URL
                data_folder + "meshes/sun_arrows.obj",
                (o) => {
                    var vector = new THREE.Vector3(-0.5, 0.5, -0.6);
                    var axis = new THREE.Vector3(0, 1, 0);
                    o.quaternion.setFromUnitVectors(axis, vector.clone().normalize());
                    o.position.copy(vector.clone().multiplyScalar(-200));
                    scene.add( o );
                    this.evening_arrows = o;
                },
                // called when loading is in progresses
                function ( xhr ) {
                    console.log( ( xhr.loaded / xhr.total * 100 ) + '% loaded' );
                },
                // called when loading has errors
                function ( error ) {
                    console.log( 'An error happened' );
                }
            );
        });
        var url = "stationary_landscape_v5_lowpoly.mtl";
        mtlLoader.load( url, materials => {
            materials.preload();
            var objLoader = new THREE.OBJLoader();
            for (const [key, value] of Object.entries(materials.materials)) {
                value.side = THREE.DoubleSide;
            }
            objLoader.setMaterials( materials );
            objLoader.load(
                // resource URL
                data_folder + "meshes/stationary_landscape_v5_lowpoly.obj",
                (o) => {
                    console.log('Landscape loaded')
                    scene.add( o );
                    this.landscape = o;
                },
                // called when loading is in progresses
                function ( xhr ) {
                    console.log( ( xhr.loaded / xhr.total * 100 ) + '% loaded' );
                },
                // called when loading has errors
                function ( error ) {
                    console.log( 'An error happened' );
                }
            );
        });

        // track last rotation to undo before new one on each application var change
        this.prev_app_value = 0;
    }

    update(application_point) {
        console.log("App point" + application_point)
        let offset = (this.numFaces / 2 - this.houseInset) * this.roof_vis_scale_factor;

        //this is not setting the translation but applying one -- need to do incremental update
        this.house.translateX(offset );
        // this.house.translateY(offset );

        this.house.rotateZ(-1* this.prev_app_value  * Math.PI*2);
        this.house.rotateZ(application_point[0] * Math.PI*2);

        this.house.translateX(- offset );
        // this.house.translateY(- offset );

        this.prev_app_value = application_point[0];
    }
}