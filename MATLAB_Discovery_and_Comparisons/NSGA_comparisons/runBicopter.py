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


class Bicopter(Problem):

    def __init__(self, r):
        super().__init__(n_var=32,
                         n_obj=2,
                         n_constr=0,
                         xl=np.array(32*[0]),
                         xu=np.array(32*[1]))
        self.anchor = anchor
        # anchored_r = self.anchor[-1] + 0.001 * r
        anchored_r = r
        self.r =  self.lerp(anchored_r, 0.5, 1)

        mass_density = 1.0
        self.m = self.r * mass_density
        self.q_init = np.zeros(6)
        self.dt = 0.25
        self.position_goal = np.array([1, 0, 0, 0, 0, 0])

    def _evaluate(self, x, out, *args, **kwargs):
        scaledu = self.rescaleParams(x)
        f1, f2 = self.bicopter_dynamics(scaledu)
        out["F"] = np.column_stack([f1, f2])

    def bicopter_dynamics(self, u):
        g = -9.81
        last_q = np.tile(self.q_init, [u.shape[0], 1])
        position_goal = np.tile(self.position_goal, [u.shape[0], 1])
        I = self.inertia_func(self.m, self.r)
        n = u.shape[1] // 2
        f1 = np.zeros(u.shape[0])
        f2 = np.zeros(u.shape[0])
        shaped_u = np.reshape(u, [u.shape[0], n, 2])

        for i in range(n):
            x = last_q[:, 0]
            y = last_q[:, 1]
            theta = last_q[:, 2]
            old_vel = last_q[:, 3:]
            delat_x = old_vel * self.dt
            run_val = 0.1 * 0.5 * np.linalg.norm(u[:, i*2:i*2+2], axis=1) ** 2 * 0.7
            f2 = f2 + run_val

            vx = - (shaped_u[:, i, 0] + shaped_u[:, i, 1]) * np.sin(theta) / self.m * self.dt
            vy = (shaped_u[:, i, 0] + shaped_u[:, i, 1]) * np.cos(theta) - self.m * g * self.dt
            vtheta = self.r * (shaped_u[:, i, 0] - shaped_u[:, i, 1]) / I * self.dt

            last_q = last_q + np.concatenate([delat_x, np.expand_dims(vx, -1), np.expand_dims(vy, -1), np.expand_dims(vtheta, -1)], axis=1)

        final_val = 0.5 * np.linalg.norm(last_q-position_goal, axis=1) ** 2 #/ 300
        f1 = f1 + final_val
        return f1/1e3, f2/5

    def inertia_func(self, m, r):
        return 1/12 * m * ((2*r) ** 2)

    def lerp(self, t, min_, max_):
        return min_ + t*(max_-min_)

    def rescaleParams(self, x):
        scaledx = np.empty_like(x)
        anchored_x = x #0.001 * x + np.tile(self.anchor[:32], [x.shape[0],1])
        for i in range(x.shape[1]):
            scaledx[:, i] = self.lerp(anchored_x[:,i], -5, 5)
        return scaledx

tic()

# / 300

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
funcEvals = 0
time = 0
numSteps = 6
allRes = []
for i in range(0, numSteps-1):
    r = 0.1 + i / (numSteps-1) # need [0,1] inclusive; for spacing of 0.1, need 11 fronts
    problem = Bicopter(r)
    algorithm = NSGA2(
        pop_size=200,
        n_offsprings=10,
        sampling=np.random.random((10, 32)),
        crossover=get_crossover("real_sbx", prob=0.9, eta=15),
        mutation=get_mutation("real_pm", eta=20),
        eliminate_duplicates=True
    )
    res = minimize(problem,
                algorithm, 
                termination,
                seed=1,
                pf=problem.pareto_front(use_cache=False),
                save_history=True,
                verbose=True)
    allRes.append(res)
    funcEvals += res.algorithm.evaluator.n_eval
    time += res.exec_time
toc()

print("Total number of function evals over", numSteps, "Pareto fronts:", funcEvals)
print("Total execution time from solver:", time)

# Design Space
plot = Scatter(title = "Design Space", axis_labels="x")
plot.add(res.X[:, 0:2])
plot.do()
plot.apply(lambda ax: ax.set_xlim(-0.5, 1.5))
plot.apply(lambda ax: ax.set_ylim(-2, 2))
plot.show()

hv_all = []
# Objective Space
plot = Scatter(title = "Objective Space", labels=["Distance to goal", "Energy", "Length"])
for i in range(0, numSteps-1):
    # r = 0.0025 + i / (numSteps-1) * 0.6
    r = 0.1 + i / (numSteps-1) 

    front = allRes[i].F
    numPts = front.shape[0]
    z = r * np.ones([numPts, 1])
    FCfront = np.concatenate((front, z), 1)
    plot.add(FCfront[np.where(FCfront[:,0] < 1)])
    metric = Hypervolume(ref_point=np.array([1.0, 1.0]))
    hv_all.append(metric.calc(front[np.where(front[:,0] < 1)]))

plot.do()
plot.apply(lambda ax: ax.set_xlim(0, 1))
plot.apply(lambda ax: ax.set_ylim(0.3, 1))
plot.apply(lambda ax: ax.set_zlim(0, 1))
plot.show()

print(hv_all)

