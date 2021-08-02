class JSMeshGenerator {
	house_z = -1.5;
	houseInset = 1; // roof overhangs house connection by this many squares

	dp_scalefactor = 2.75; // in general would need a lerp, but [0,5] is our range
	vis_scale_factor = 12

	generate_mesh(dp, ap) {
		let Ni = Math.sqrt(dp.length);// +2; // assume a perfect square, +2 for implicit BC (zeros on edges)
		let Nj = Ni;

		// create the roof vertices
		let [roof_verts, roof_tris] = this.createRoof(dp, Ni, Nj);

		// create the box underneath for attaching to house		
		let [house_verts, house_tris] = this.createHouseConnection(roof_verts, Ni, Nj);

		roof_verts = this.rotate_and_translate(ap[0], roof_verts, Ni-1, Nj-1); // -1 because counting faces not vertices
		house_verts = this.rotate_and_translate(ap[0], house_verts, Ni-1, Nj-1); // -1 because counting faces not vertices

		let m1 = this.meshFromData(roof_verts, roof_tris, 0x3c3c3c);
		let m2 = this.meshFromData(house_verts, house_tris, 0xffffffff);
		return [m1, m2];
	}

	meshFromData(verts, tris, color) {
		let geom = geometry_from_data(verts, tris);
		// var material = new THREE.MeshPhongMaterial({color: color, specular: 0x111111, shininess: 60, flatShading: true, side: THREE.DoubleSide});
		var material = new THREE.MeshPhongMaterial({color: color, specular: 0x111111, shininess: 60, flatShading: true});
		geom.computeFlatVertexNormals();
		return new THREE.Mesh( geom, material );
	}

	createHouseConnection(roof_pts, Ni, Nj) {
		let Bi = [this.houseInset, Ni - this.houseInset-1]; // extents of boundary edges
		let Bj = [this.houseInset, Nj - this.houseInset-1];

		// construct the walls of the building
		let house_pts = [];
		let house_tris = [];
		let hID = 0;
		for (let j=Bj[0]; j<Bj[1]; j++) {
			hID = this.addWallSegment(Bi[0], j, Ni, house_pts, house_tris, hID, roof_pts, true, false);
			hID = this.addWallSegment(Bi[1], j, Ni, house_pts, house_tris, hID, roof_pts, true, true);
		}
		for (let i=Bi[0]; i<Bi[1]; i++) {
			hID = this.addWallSegment(i, Bj[0], Ni, house_pts, house_tris, hID, roof_pts, false, true);
			hID = this.addWallSegment(i, Bj[1], Ni, house_pts, house_tris, hID, roof_pts, false, false);
		}

		return [house_pts, house_tris];
	}

	addWallSegment(i, j, Ni, house_pts, house_tris, hID, roof_pts, increment_j, reverse_normals) {
		// ids of roof indices + duplicate points for separate mesh
		let r_ij = this.get_linearID(i, j, Ni);
		house_pts.push(roof_pts[r_ij]);
		r_ij = hID;
		hID++;

		let r_ij1 = this.get_linearID(i+1, j, Ni);
		if (increment_j) {
			r_ij1 = this.get_linearID(i, j+1, Ni);
		}
		house_pts.push(roof_pts[r_ij1]);
		r_ij1 = hID;
		hID++;


		// wall vertices
		house_pts.push([i, j, this.house_z]);
		let w_i_j = hID;
		hID++;

		if (increment_j) {
			house_pts.push([i, j+1, this.house_z]);
		} else {
			house_pts.push([i+1, j, this.house_z]);
		}
		let w_ij1 = hID;
		hID++;

		// wall triangles	
		if (reverse_normals) {			
			this.add_square_reverse(w_i_j, w_ij1, r_ij, r_ij1, house_tris);
		} else {
			this.add_square(w_i_j, w_ij1, r_ij, r_ij1, house_tris);
		}

		return hID;
	}

	rotate_and_translate(thetaIn, pts, Ni, Nj) {
		let theta = thetaIn;//-0.0015;
		console.log("Theta" + theta);

		theta = 2*Math.PI * theta;
		let cosT = Math.cos(theta);
		let sinT = Math.sin(theta);

		let rotatedPts = [];
		for (let i=0; i<pts.length; i++) {
			let pt = pts[i];
			pt = [pt[0] - Ni/2, pt[1] - this.houseInset, pt[2]];

			let x = cosT * pt[0] - sinT * pt[1];
			let y = sinT * pt[0] + cosT * pt[1];
			let newPt = [x + Ni/2 - this.houseInset, 
				y, 
				pt[2] - this.house_z]; // align corner with (0,0) and move up to move on top of house base
			rotatedPts.push([newPt[0]*this.vis_scale_factor, newPt[1]*this.vis_scale_factor, newPt[2]*this.vis_scale_factor]); 
		}

		return rotatedPts;
	}

	// create the roof
	createRoof(dp, Ni, Nj) {
		let [roof_verts, roof_verts_upper] = this.roofVerts(dp, Ni, Nj);
		let numVerts = roof_verts.length;

		let roof_tris = [];
		let upper_tris = [];
		for (let i=0; i < Ni-1; i++) {
			for (let j=0; j < Nj-1; j++) {
				let i_j = this.get_linearID(i, j, Ni);
				let i_j1 = this.get_linearID(i, j+1, Ni);
				let i1_j = this.get_linearID(i+1, j, Ni);
				let i1_j1 = this.get_linearID(i+1, j+1, Ni);
				
				this.add_square_reverse(i_j, i_j1, i1_j, i1_j1, roof_tris);
				this.add_square(numVerts+i_j, numVerts+i_j1, numVerts+i1_j, numVerts+i1_j1, upper_tris);
			}
		}

		// connect the roof layers
		let Bi = [0, Ni -1]; // extents of boundary edges
		let Bj = [0, Nj -1];
		let i, j, ij_low, ij1_low, ij_high, ij1_high;
		for (let j=Bj[0]; j<Bj[1]; j++) {
			// small i edge
			i = Bi[0];
			ij_low = this.get_linearID(i, j, Ni);
			ij1_low = this.get_linearID(i, j+1, Ni);
			ij_high = numVerts + ij_low;
			ij1_high = numVerts + ij1_low;
			this.add_square(ij_low, ij1_low, ij_high, ij1_high, upper_tris);

			// large i edge
			i = Bi[1];
			ij_low = this.get_linearID(i, j, Ni);
			ij1_low = this.get_linearID(i, j+1, Ni);
			ij_high = numVerts + ij_low;
			ij1_high = numVerts + ij1_low;
			this.add_square_reverse(ij_low, ij1_low, ij_high, ij1_high, upper_tris);
		}
		for (let i=Bi[0]; i<Bi[1]; i++) {
			// small j edge
			j = Bj[0];
			ij_low = this.get_linearID(i, j, Ni);
			ij1_low = this.get_linearID(i+1, j, Ni);
			ij_high = numVerts + ij_low;
			ij1_high = numVerts + ij1_low;
			this.add_square_reverse(ij_low, ij1_low, ij_high, ij1_high, upper_tris);

			// large j edge
			j = Bj[1];
			ij_low = this.get_linearID(i, j, Ni);
			ij1_low = this.get_linearID(i+1, j, Ni);
			ij_high = numVerts + ij_low;
			ij1_high = numVerts + ij1_low;
			this.add_square(ij_low, ij1_low, ij_high, ij1_high, upper_tris);
		}

		roof_verts = roof_verts.concat(roof_verts_upper);
		roof_tris = roof_tris.concat(upper_tris);

		return [roof_verts, roof_tris];
	}


	// create the zero-padded roof vertices from the zheights in dp
	roofVerts(dp, Ni, Nj) {
		let verts = [];
		let verts_upper = [];
		for (let j=0; j < Nj; j++) {
			for (let i=0; i<Ni; i++) {
				// let z = 0;
				// if (i==0 || i==Ni-1 || j==0 || j==Nj-1) {
				// 	z = 0;
				// }
				// else {
				// 	z = dp[this.get_linearID(i-1, j-1, Ni-2)] * dp_scalefactor;
				// }
				let z = dp[this.get_linearID(i, j, Ni)] * this.dp_scalefactor;
				verts.push([i, j, z]);
				verts_upper.push([i, j, z+0.25]);
			}
		}

		return [verts, verts_upper];
	}

	// takes linear IDs of the four courners
	add_square(i_j, i_j1, i1_j, i1_j1, tris) {
		tris.push([i_j, i1_j, i_j1]);
		tris.push([i_j1, i1_j, i1_j1]);
	}

	// takes linear IDs of the four courners -- normals reversed relative to add_square
	add_square_reverse(i_j, i_j1, i1_j, i1_j1, tris) {
		tris.push([i_j, i_j1, i1_j]);
		tris.push([i_j1, i1_j1, i1_j]);
	}

	get_linearID(i, j, Ni) {
		return (j * Ni) + i;
	}

}