import numpy as np
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
from pymoo.util.termination.default import MultiObjectiveDefaultTermination
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter

import matplotlib.pyplot as plt
from pymoo.performance_indicator.hv import Hypervolume

from ticTocGenerator import *
from lamp import Lamp 

def runLamp():
	# ## test the output in comparison to matlab code
	# problem = Lamp(0.1)
	# points = np.array([21*[0.0],
	# 					21*[1.0]])
	# # points = np.array([21*[0.0]])
	# out = {}
	# problem._evaluate(points, out)
	# print(out["F"])


	# set up the solver
	tic(); 
	algorithm = NSGA2(
		pop_size=200,
		sampling=get_sampling("real_random"),
		crossover=get_crossover("real_sbx", prob=0.9, eta=15),
		mutation=get_mutation("real_pm", eta=20),
		eliminate_duplicates=True
	)

	termination = MultiObjectiveDefaultTermination(
		x_tol=1e-8,
		cv_tol=1e-6,
		f_tol=0.0025,
		nth_gen=5,
		n_last=30,
		n_max_gen=10000,
		n_max_evals=1000000
	)


	## run the problem 
	funcEvals = 0;
	time = 0;
	numSteps = 11;
	allRes = [];
	hypervolume = []

	for i in range(0, numSteps):
		focalz = i / (numSteps-1); # need [0,1] inclusive; for spacing of 0.1, need 11 fronts
		problem = Lamp(focalz)
		res = minimize(problem,
					algorithm,
					termination,
					seed=1,
					pf=problem.pareto_front(use_cache=False),
					save_history=True,
					verbose=True)
		allRes.append(res);

		funcEvals += res.algorithm.evaluator.n_eval
		time += res.exec_time;

		# plot the hypervolume over time
		# create the performance indicator object with reference point (1,1)
		metric = Hypervolume(ref_point=np.array([1.0, 1.0, 1.0]))

		# collect the population in each generation
		# pop_each_gen = [a.pop for a in res.history]
		pop_each_gen = res.history[-1].pop

		# receive the population in each generation
		# obj_and_feasible_each_gen = [pop[pop.get("feasible")[:,0]].get("F") for pop in pop_each_gen]
		obj_and_feasible_each_gen = pop_each_gen[pop_each_gen.get("feasible")[:,0]].get("F")

		# calculate for each generation the HV metric
		# hv = [metric.calc(f) for f in obj_and_feasible_each_gen]
		hv = metric.calc(obj_and_feasible_each_gen)
		hypervolume.append(hv)

	toc()

	print("Total number of function evals over", numSteps, "Pareto fronts:", funcEvals)
	print("Total execution time from solver:", time)
	# return time, funcEvals, hypervolume

	# Visualize
	# Design Space
	plot = Scatter(title = "Design Space", axis_labels="x")
	plot.add(res.X[:, 0:2])
	plot.do()
	# plot.apply(lambda ax: ax.set_xlim(-0.5, 1.5))
	# plot.apply(lambda ax: ax.set_ylim(-2, 2))
	plot.show()

	# Augmented Objective Space
	plot = Scatter(title = "Objective Space", labels=["Instability", "Mass", "Focal_Z"])

	for i in range(0, numSteps):
		focalz = i / (numSteps-1);

		front = allRes[i].F[:, 0:2];
		import pdb; pdb.set_trace()
		numPts = front.shape[0];
		z = focalz * np.ones([numPts, 1]);
		FCfront = np.concatenate((front, z), 1);
		plot.add(FCfront);

	plot.do()
	plot.apply(lambda ax: ax.set_xlim(0, 1))
	plot.apply(lambda ax: ax.set_ylim(0, 1))
	plot.apply(lambda ax: ax.set_zlim(0, 1))
	plot.show()

	# plot the hypervolume over generations
	# create the performance indicator object with reference point (1,1)
	metric = Hypervolume(ref_point=np.array([1.0, 1.0, 1.0]))

	# collect the population in each generation
	pop_each_gen = [a.pop for a in res.history]

	# receive the population in each generation
	obj_and_feasible_each_gen = [pop[pop.get("feasible")[:,0]].get("F") for pop in pop_each_gen]

	# calculate for each generation the HV metric
	hv = [metric.calc(f) for f in obj_and_feasible_each_gen]

	# function evaluations at each snapshot
	n_evals = np.array([a.evaluator.n_eval for a in res.history])

	# visualze the convergence curve
	plt.plot(n_evals, hv, '-o')
	plt.title("Convergence")
	plt.xlabel("Function Evaluations")
	plt.ylabel("Hypervolume")
	plt.show()


if __name__ == '__main__':
	runLamp()