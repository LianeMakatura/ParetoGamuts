import autograd.numpy as anp
import numpy as np
from numpy import linalg as LA

from pymoo.util.misc import stack
from pymoo.model.problem import Problem


class Lamp(Problem):

	def __init__(self, focal_z):
		super().__init__(n_var=21,
					n_obj=3,
					n_constr=0,
					xl=anp.array(21*[0]),
					xu=anp.array(21*[1]))

		# other constants
		self.bMax = 4.0;
		self.aMax = 2.0;
		self.hMax = 1.0;
		self.wMax = 0.05;
		self.endVolume = 0.3**3;
		self.radii = np.ones([4])*self.wMax*2;
		scale = self.bMax/2+self.aMax; 			## TODO: we can't achieve anything taller

		# set the context variable
		self.focalPoint_z = self.lerp(focal_z, 1, 5) * scale;

		self.focalPoint = np.array([2*scale, 2*scale, self.focalPoint_z]);


	def _evaluate(self, x, out, *args, **kwargs):
		scaledx = self.rescaleParams(x)
		base = scaledx[:, 0:3]
		a1 = scaledx[:, 3:6]
		a2 = scaledx[:, 6:9]
		a3 = scaledx[:, 9:12]
		h1 = scaledx[:, 12:15]
		h2 = scaledx[:, 15:18]
		h3 = scaledx[:, 18:21]

		f1 = self.f1(base, a1, a2, a3, h1, h2, h3)
		f2 = self.f2(base, a1, a2, a3, h1, h2, h3)
		f3 = self.f3(base, a1, a2, a3, h1, h2, h3)
		out["F"] = anp.column_stack([f1, f2, f3])


	# instability
	def f1(self, base, a1, a2, a3, h1, h2, h3):
		COM = self.get_centerOfMass(base, a1, a2, a3, h1, h2, h3);
		f1 = COM[:, 0]**2 + COM[:, 1]**2;   

		# normalize so all possible outputs are between 0 and 1
		f1 = f1 / 15.0; #observed via jscad
		return f1

	# mass
	def f2(self, base, a1, a2, a3, h1, h2, h3):
		minMass = (LA.norm(np.array([0.0,0.0,1.0]))*(self.aMax*3 + self.bMax) + LA.norm(np.array([0.4,0.4,1.0]) * self.hMax*3) )*2*(self.wMax);
		maxMass  = LA.norm(np.array([1.0,1.0,2.0]))*(self.hMax*3 + self.aMax*3 + self.bMax)*(2*self.wMax);

		mass = self.get_mass(base, a1, a2, a3, h1, h2, h3);
		# normalize so all possible outputs are between 0 and 1
		f2 = (mass -minMass)/(maxMass - minMass);
		return f2

	# distance to focal point
	def f3(self, base, a1, a2, a3, h1, h2, h3):
		# average distance to the focal point
		numpts = base.shape[0];
		focalRep = self.repRow(self.focalPoint, numpts);
		endPosCost = self.rownorm(base+a1+h1 -focalRep) \
					+ self.rownorm(base+a2+h2 -focalRep) \
					+ self.rownorm(base+a3+h3 -focalRep);
		endPosCost = endPosCost / 3.0;

		#normalize so outputs between 0 and 1; estimated in jscad viewer 
		f3 = endPosCost/20.0;#(endPosCost/3- 10)/ 8;
		return f3


	def get_mass(self, base, a1, a2, a3, h1, h2, h3):
		barMasses = self.get_elementMasses(base, a1, a2, a3, h1, h2, h3);
		# undo the square to match Schulz et al implementation
		barMasses[:, 0] = barMasses[:, 0] / self.radii[0];
		barMasses[:, 1] = barMasses[:, 1] / self.radii[1];
		barMasses[:, 2] = barMasses[:, 2] / self.radii[2];
		barMasses[:, 3] = barMasses[:, 3] / self.radii[3];
		barMasses[:, 4] = barMasses[:, 4] / self.radii[1];
		barMasses[:, 5] = barMasses[:, 5] / self.radii[2];
		barMasses[:, 6] = barMasses[:, 6] / self.radii[3];
		# totalMass = self.rowsum(barMasses) + 3*self.endVolume;
		return self.rowsum(barMasses);


	def get_elementMasses(self, base, a1, a2, a3, h1, h2, h3):
		numpts = base.shape[0];
		masses = np.zeros([numpts, 7]); # one column for each element
		masses[:, 0] = self.rownorm(base) * self.radii[0]**2;
		masses[:, 1] = self.rownorm(a1) * self.radii[1]**2;
		masses[:, 2] = self.rownorm(a2) * self.radii[2]**2;
		masses[:, 3] = self.rownorm(a3) * self.radii[3]**2;
		masses[:, 4] = self.rownorm(h1) * self.radii[1]**2;
		masses[:, 5] = self.rownorm(h2) * self.radii[2]**2;
		masses[:, 6] = self.rownorm(h3) * self.radii[3]**2;
		return masses

	def get_centerOfMass(self, base, a1, a2, a3, h1, h2, h3):
		masses = self.get_elementMasses(base, a1, a2, a3, h1, h2, h3);
		totalMass = self.rowsum(masses) + 3*self.endVolume;

		centerOfMass = base * 0.5 * self.repCol(masses[:, 0], 3) \
					+ (base + a1*0.5) * self.repCol(masses[:, 1], 3) \
					+ (base + a2*0.5) * self.repCol(masses[:, 2], 3) \
					+ (base + a3*0.5) * self.repCol(masses[:, 3], 3) \
					+ (base + a1 + h1*0.5) * self.repCol(masses[:, 4], 3) \
					+ (base + a2 + h2*0.5) * self.repCol(masses[:, 5], 3) \
					+ (base + a3 + h3*0.5) * self.repCol(masses[:, 6], 3) \
					+ (base + a1 + h1) * self.endVolume \
					+ (base + a2 + h2) * self.endVolume \
					+ (base + a3 + h3) * self.endVolume;
		return centerOfMass / self.repCol(totalMass, 3); 

	def repCol(self, vec, numcols): # creates an nxCols matrix for elementwise mult
		numrows = vec.shape[0];
		a = vec;
		for i in range(0, numcols - 1):
			a = np.concatenate((a, vec));
		a = np.reshape(a, [numrows, numcols], 'F'); # 'F' treats it as a column-major ordering
		return a

	def repRow(self, vec, numrows): # creates an nx3 matrix for elementwise mult
		numcols = vec.shape[0];
		a = vec;
		for i in range(0, numrows - 1):
			a = np.concatenate((a, vec));
		a = np.reshape(a, [numrows, numcols], 'C'); # 'C' treats it as a row-major ordering
		return a

	def rownorm(self, mat):
		return LA.norm(mat, axis=1);

	def rowsum(self, mat):
		return np.sum(mat, axis=1);



	# ====== parameter setup

	def rescaleParams(self, x):
		scaledx = np.empty_like(x)
		#base
		scaledx[:, 0] = self.lerp(x[:, 0], 0, 1) * self.bMax;
		scaledx[:, 1] = self.lerp(x[:, 1], 0, 1) * self.bMax;
		scaledx[:, 2] = self.lerp(x[:, 2], 1, self.bMax) * self.bMax;
		#arm 1
		scaledx[:, 3] = self.lerp(x[:, 3], 0, 1) * self.aMax;
		scaledx[:, 4] = self.lerp(x[:, 4], 0, 1) * self.aMax;
		scaledx[:, 5] = self.lerp(x[:, 5], 1, 2) * self.aMax;
		#arm 2
		scaledx[:, 6] = self.lerp(x[:, 6], -1, 0) * self.aMax;
		scaledx[:, 7] = self.lerp(x[:, 7], 0, 1) * self.aMax;
		scaledx[:, 8] = self.lerp(x[:, 8], 1, 2) * self.aMax;
		#arm 3
		scaledx[:, 9] = self.lerp(x[:, 9], 0, 1) * self.aMax;
		scaledx[:, 10] = self.lerp(x[:, 10], -1, 0) * self.aMax;
		scaledx[:, 11] = self.lerp(x[:, 11], 1, 2) * self.aMax;
		#hand 1
		scaledx[:, 12] = self.lerp(x[:, 12], 0.4, 1) * self.hMax;
		scaledx[:, 13] = self.lerp(x[:, 13], 0.4, 1) * self.hMax;
		scaledx[:, 14] = self.lerp(x[:, 14], -2, -1) * self.hMax;
		#hand 2
		scaledx[:, 15] = self.lerp(x[:, 15], -1, -0.4) * self.hMax;
		scaledx[:, 16] = self.lerp(x[:, 16], 0.4, 1) * self.hMax;
		scaledx[:, 17] = self.lerp(x[:, 17], -2, -1) * self.hMax;
		#hand 3
		scaledx[:, 18] = self.lerp(x[:, 18], 0.4, 1) * self.hMax;
		scaledx[:, 19] = self.lerp(x[:, 19], -1, -0.4) * self.hMax;
		scaledx[:, 20] = self.lerp(x[:, 20], -2, -1) * self.hMax;

		return scaledx

	def lerp(self, t, min, max):
		return min + t*(max-min)


