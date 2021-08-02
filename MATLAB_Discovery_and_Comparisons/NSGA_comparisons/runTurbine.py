import numpy as np
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter

import matplotlib.pyplot as plt
from pymoo.performance_indicator.hv import Hypervolume

from ticTocGenerator import *
from turbine import Turbine 
import scipy.io
import os

def runTurbine():
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
	numSteps = 6;
	allRes = []
	hypervolume = []

	for i in range(0, numSteps):
		vwind = i / (numSteps-1); # need [0,1] inclusive; for spacing of 0.1, need 11 fronts
		problem = Turbine(vwind)
		res = minimize(problem,
					algorithm,
					seed=1,
					pf=problem.pareto_front(use_cache=False),
					save_history=True,
					verbose=True)
		allRes.append(res)
		funcEvals += res.algorithm.evaluator.n_eval
		time += res.exec_time

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
