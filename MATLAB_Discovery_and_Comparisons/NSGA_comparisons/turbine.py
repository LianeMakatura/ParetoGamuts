import autograd.numpy as anp
import numpy as np
from pymoo.util.misc import stack
from pymoo.model.problem import Problem


class Turbine(Problem):

	def __init__(self, windv):
		super().__init__(n_var=3,
							n_obj=2,
							n_constr=0,
							xl=anp.array([0,0,0]),
							xu=anp.array([1,1,1]))

		scaledv = self.lerp(windv, 4, 6)
		self.v = scaledv # context parameter

	def _evaluate(self, x, out, *args, **kwargs):
		scaledx = self.rescaleParams(x)
		radius = scaledx[:, 0]
		height = scaledx[:, 1]
		pitch = scaledx[:, 2]

		f1 = self.f1(radius, height, pitch)
		f2 = self.f2(radius, height, pitch)
		out["F"] = anp.column_stack([f1, f2])

	## mass 
	def f1(self, radius, height, pitch): 
		f1 = radius * height * 1/pitch 

		# normalize so all possible outputs are between 0 and 1	
		min_mass = 229.1831; # 0,0,1
		max_mass = 2.1486e+04; #1,1,0
		f1 = (f1 - min_mass) / (max_mass - min_mass);
		return f1

	## power output
	def f2(self, radius, height, pitch):
		betz = False;

		fluidDensity = 1; # can really disregard, assumed constant
		sweptArea = np.pi * radius**2;

		if betz:
			Cp = 16/27; # BetzCoeff; theoretical limit, for ideal only
			min_power = -2.0106e+06; # for experiment with wind in 4-6m/s
			max_power = -9.5318e+04;
		else:
			Cp = self.performanceCoefficient(radius);
			min_power = -1096672.537507; # for experiment with wind in 4-6m/s
			max_power = 28929.459366;

		P_out = Cp * 1/2 * fluidDensity * sweptArea * self.v**3; 
		f2 = -P_out; # want to maximize power, so we minimize the negative

		# normalize so all possible outputs are between 0 and 1
		f2 = (f2 - min_power) / (max_power - min_power);
		return f2 

	def performanceCoefficient(self, radius):
		TSR = self.tipSpeedRatio(radius)
		# copied from the results of matlabs polyfit in Turbine
		degree = 6;
		p = [-1.05348512795590e-07,	7.08209729939267e-06,	-0.000140277525378244,	0.000307565692335347,	0.0118771972725566,	-0.0352202780490948,	0.0160943595028349]

		evald = 0;
		for i in range(0, degree+1):
			evald = evald + p[i]*TSR**(degree - i);
		return evald;

	def tipSpeedRatio(self, radius):
		#function of tip-speed ratio (TSR, aka lambda)
		rpm = 8; #chosen to preserve tip ratio range
		omegam = 2 * np.pi * rpm / 60;
		tsr = radius * omegam / self.v;
		return tsr


	def rescaleParams(self, x):
		scaledx = np.empty_like(x)
		scaledx[:, 0] = self.lerp(x[:, 0], 40, 100) #radius [m]
		scaledx[:, 1] = self.lerp(x[:, 1], 6, 15) 	# height [m]
		scaledx[:, 2] = self.lerp(x[:, 2], 5, 20) * (np.pi / 180) # pitch angle [rad]
		return scaledx

	def lerp(self, t, min, max):
		return min + t*(max-min)

