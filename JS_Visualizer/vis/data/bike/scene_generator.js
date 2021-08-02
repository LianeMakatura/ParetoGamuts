class CustomSceneGenerator {
    Init() {
    }

    async init(data_folder, scene) {
        var mtlLoader = new THREE.MTLLoader();
        mtlLoader.setPath( data_folder + 'meshes/' );
        var url = "shock.mtl";
        mtlLoader.load( url, materials => {
            materials.preload();
            var objLoader = new THREE.OBJLoader();
            for (const [key, value] of Object.entries(materials.materials)) {
                value.side = THREE.DoubleSide;
            }
            objLoader.setMaterials( materials );
            objLoader.load(
                // resource URL
                data_folder + "meshes/shock.obj",
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
    }

    update(application_point) {

    }
}