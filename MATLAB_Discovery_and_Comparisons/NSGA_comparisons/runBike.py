import autograd.numpy as anp
import numpy as np
from pymoo.util.misc import stack
from pymoo.model.problem import Problem
from pymoo.factory import get_problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
from pymoo.util.termination.default import MultiObjectiveDefaultTermination
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter

import matplotlib.pyplot as plt
from pymoo.performance_indicator.hv import Hypervolume
from ticTocGenerator import *
import time


class BikeRocker(Problem):

    def __init__(self, pivot_offset):
        super().__init__(n_var=3,
                         n_obj=2,
                         n_constr=0,
                         xl=np.array(3*[0]),
                         xu=np.array(3*[1]))

        scaled_pivot_offset = self.lerp(pivot_offset, 0, -40)
        self.pivot_offset = scaled_pivot_offset
        self.pivotbase = [0,0]
        self.rocker_r = anp.array([self.pivotbase[0] + 100, self.pivotbase[1]])
        self.rocker_b = anp.array([2/3*(self.rocker_r[0] - self.pivotbase[0]), self.pivotbase[1] - 50])
        self.solid_angle = 30
        self.finalPivot = anp.array([self.pivotbase[0] + self.pivot_offset, self.pivotbase[1]])

    def _evaluate(self, x, out, *args, **kwargs):
        scaled_x = self.rescaleParams(x)
        solid_offset = scaled_x[:, 0]
        thickness = scaled_x[:, 1]
        depth = scaled_x[:, 2]
        f1 = self.f1(solid_offset, thickness, depth)
        f2 = self.f2(solid_offset, thickness, depth)
        out["F"] = np.column_stack([f1, f2])

    def f1(self, solid_offset, thickness, depth):
        # import pdb; pdb.set_trace()
        top = thickness * np.sqrt(np.sum((self.finalPivot - self.rocker_r) ** 2, axis=-1))
        left = thickness * np.sqrt(np.sum((self.finalPivot - self.rocker_b) ** 2, axis=-1))
        right = thickness * np.sqrt(np.sum((self.rocker_r - self.rocker_b) ** 2, axis=-1))
        front_area = top + left + right + solid_offset ** 2
        f1 = front_area * depth
        f1 = (f1 - 1.1645e+04) / (6.4756e+04 - 1.1645e+04)
        return f1

    def f2(self, solid_offset, thickness, depth):
        top = thickness * np.sqrt(np.sum((self.finalPivot - self.rocker_r) ** 2, axis=-1))
        left = thickness * np.sqrt(np.sum((self.finalPivot - self.rocker_b) ** 2, axis=-1))
        right = thickness * np.sqrt(np.sum((self.rocker_r - self.rocker_b) ** 2, axis=-1))
        front_area = top + left + right + solid_offset ** 2
        mass = front_area * depth
        axisLen = np.sqrt(np.sum((self.finalPivot - self.rocker_r) ** 2, axis=-1))
        
        f2 = axisLen / mass; 
        f2 = (f2 - 0.0015) / (0.01 - 0.0015)
        return f2

    def lerp(self, t, min_, max_):
        return min_ + t*(max_-min_)

    def rescaleParams(self, x):
        scaledx = np.empty_like(x)
        scaledx[:, 0] =  self.lerp(x[:, 0], 25, 45)
        scaledx[:, 1] =  self.lerp(x[:, 1], 7, 14)
        scaledx[:, 2] =  self.lerp(x[:, 2], 5, 10)
        return scaledx


def runBike(run, z_range):
	# ## test the output in comparison to matlab code
	# problem = Turbine(0.1)
	# points = np.array([[0.0,0.0,0.0],
	# 					[1.0,1.0,1.0]])
	# out = {}
	# print(problem._evaluate(points, out))
	# print(out["F"])

	tic(); 
	algorithm = NSGA2(
		pop_size=200,
		sampling=get_sampling("real_random"),
		crossover=get_crossover("real_sbx", prob=0.9, eta=15),
		mutation=get_mutation("real_pm", eta=20),
		eliminate_duplicates=True
	)
	# termination = get_termination("n_gen", 50)

	funcEvals = 0;
	time = 0;
	numSteps = len(np.arange(z_range[0], z_range[1], z_range[2]))
	allRes = []
	hypervolume = []
	for i in np.arange(z_range[0], z_range[1], z_range[2]):
		print('run:', run, 'context:', i)
		pivot_offset = i # need [0,1] inclusive; for spacing of 0.1, need 11 fronts
		problem = BikeRocker(pivot_offset)
		res = minimize(problem,
					algorithm,
					# seed=1,
					pf=problem.pareto_front(use_cache=False),
					save_history=True,
					verbose=True)
		allRes.append(res)
		funcEvals += res.algorithm.evaluator.n_eval
		time += res.exec_time
		import pdb; pdb.set_trace()

		# plot the hypervolume over time
		# create the performance indicator object with reference point (1,1)
		metric = Hypervolume(ref_point=np.array([1.0, 1.0]))

		# collect the population in each generation
		pop_each_gen = [a.pop for a in res.history]

		# receive the population in each generation
		obj_and_feasible_each_gen = [pop[pop.get("feasible")[:,0]].get("F") for pop in pop_each_gen]
		# calculate for each generation the HV metric
		hv = [metric.calc(f) for f in obj_and_feasible_each_gen]
		hypervolume.append(hv[-1])


	toc()
	print("Total number of function evals over", numSteps, "Pareto fronts:", funcEvals)
	print("Total execution time from solver:", time)
	return time, funcEvals, hypervolume

	# scipy.io.savemat('./{}.mat'.format(os.path.basename(__file__).split('.')[0]), stats)

	# Visualize
	# Design Space
	plot = Scatter(title = "Design Space", axis_labels="x")
	plot.add(res.X)
	plot.do()
	# plot.apply(lambda ax: ax.set_xlim(-0.5, 1.5))
	# plot.apply(lambda ax: ax.set_ylim(-2, 2))
	plot.show()

	# Augmented Objective Space
	plot = Scatter(title = "Objective Space", labels=["Mass", "Power", "Wind Speed"])

	for i in range(0, numSteps):
		vwind = i / (numSteps-1);

		front = allRes[i].F;
		numPts = front.shape[0];
		z = vwind * np.ones([numPts, 1]);
		FCfront = np.concatenate((front, z), 1);
		plot.add(FCfront);

	plot.do()
	plot.apply(lambda ax: ax.set_xlim(0, 1))
	plot.apply(lambda ax: ax.set_ylim(0, 1))
	plot.apply(lambda ax: ax.set_zlim(0, 1))
	plot.show()

	import pdb; pdb.set_trace()

	# function evaluations at each snapshot
	n_evals = np.array([a.evaluator.n_eval for a in res.history])

	# visualze the convergence curve
	plt.plot(n_evals, hv, '-o')
	plt.title("Convergence")
	plt.xlabel("Function Evaluations")
	plt.ylabel("Hypervolume")
	plt.show()
    
if __name__ == '__main__':
	time, funcEvals, hypervolume = runBike('BikeRocker', [0.1, 1.0, 0.2])
	print(time, funcEvals, hypervolume)