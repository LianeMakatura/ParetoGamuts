class JSMeshGenerator {
	// returns [[[x, y, z]], [[i1, i2, i3]]]
	generate_mesh(dp, ap) {
		Ni = sqrt(dp.length) +2; // assume a perfect square, +2 for implicit BC (zeros on edges)
		Nj = Ni;

		dp_scalefactor = 5; // in general would need a lerp, but [0,5] is our range

		const verts = roofVerts(dp, Ni, Nj);

		var tris = [];
		for (i=0; i<Ni-1; i++) {
			for (j=0; j<Nj-1; j++) {
				i_j = get_linearID(i, j);
				i_j1 = get_linearID(i, j+1);
				i1_j = get_linearID(i+1, j);
				i1_j1 = get_linearID(i+1, j+1);
				
				add_square(i_j, i_j1, i1_j, i1_j1, tris);
			}
		}

		return [verts, tris];
	}

	// create the zero-padded roof vertices from the zheights in dp
	roofVerts(dp, Ni, N, dp_scalefactor) {
		verts = [];
		for (i=0; i < Ni-1; i++) {
			for (j=0; j<Nj-1; j++) {
				if (i==0 || i==Ni || j==0 || j==Nj) {
					z = 0;
				}
				else {
					z = dp[get_linearID(i,j)] * dp_scalefactor;
				}
				verts.push([i, j, z]);
			}
		}

		return verts;
	}


	// takes linear IDs of the four courners
	add_square(i_j, i_j1, i1_j, i1_j1, tris) {
		tris.push([i_j, i_j1, i1_j]);
		tris.push([i_j1, i1_j1, i1_j]);
	}

	get_linearID(i, j, Ni) {
		return (j * Ni) + i;
	}

}