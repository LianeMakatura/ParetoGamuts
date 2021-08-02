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


class Bike(Problem):

    def __init__(self):
        super().__init__(n_var=7,
                         n_obj=4,
                         n_constr=0,
                         xl=np.array(7*[0]),
                         xu=np.array(7*[1])) 
        # n_var: 0: solid_offset; 1: thickness; 2: depth; 3: drop_h; 4: tube_h; 5: stay_depth; 6: pivot_offset x; 7 (optional): pivot_offset y
        # scaled_pivot_offset = self.lerp(pivot_offset, 0, -40)
        # self.pivot_offset = scaled_pivot_offset
        self.pivotbase = [0,0]
        self.rocker_r = anp.array([self.pivotbase[0] + 100, self.pivotbase[1]])
        self.rocker_b = anp.array([2/3*(self.rocker_r[0] - self.pivotbase[0]), self.pivotbase[1] - 50])
        self.solid_angle = 30
        self.drop_w = 30
        self.stay_b = [self.pivotbase[0] - 180, self.pivotbase[1] - 175]
        # self.finalPivot = anp.array([self.pivotbase[0] + self.pivot_offset, self.pivotbase[1]])

    def _evaluate(self, x, out, *args, **kwargs):
        scaled_x = self.rescaleParams(x)
        # scaledx[:, 0] =  self.lerp(x[:, 0], 25, 45)
        # scaledx[:, 1] =  self.lerp(x[:, 1], 7, 14)
        # scaledx[:, 2] =  self.lerp(x[:, 2], 5, 10)
        # scaledx[:, 3] =  self.lerp(x[:, 3], 50, 80)
        # scaledx[:, 4] =  self.lerp(x[:, 4], 10, 20)
        # scaledx[:, 5] =  self.lerp(x[:, 5], 6, 16)
        # scaledx[:, 6] =  self.lerp(x[:, 6], 0, -40)
        solid_offset = scaled_x[:, 0] # (25, 45)
        thickness = scaled_x[:, 1] # (7, 14)
        depth = scaled_x[:, 2] # (5, 10)
        drop_h = scaled_x[:, 3] # (50, 80)
        tube_h = scaled_x[:, 4] # (10, 20)
        stay_depth = scaled_x[:, 5] # (6, 16)
        pivot_offset_x = scaled_x[:, 6] # (0, -40)
        finalPivot = anp.array([self.pivotbase[0] + pivot_offset_x, self.pivotbase[1]]) # [(0, -40), 0]
        # Rocker
        f1 = self.f1(finalPivot, solid_offset, thickness, depth)
        f2 = self.f2(finalPivot, solid_offset, thickness, depth)
        # Stay
        ur_drop = anp.array([self.stay_b[0] + self.drop_w, self.stay_b[1] + drop_h]) # [-150, (-125, -95)]
        tube_length = anp.array(anp.sqrt(anp.sum((finalPivot - ur_drop) ** 2, axis=-1)))
        f3 = self.f3(drop_h, stay_depth, tube_h, tube_length)
        f4 = self.f4(tube_length, tube_h, stay_depth)
        out["F"] = np.column_stack([f1, f2, f3, f4])

    def f1(self, finalPivot, solid_offset, thickness, depth):
        # solid_offset: (25, 45)
        # thickness: (7, 14)
        # depth: (5, 10)
        # finalPivot: [(0, -40), 0]
        # rocker_r = [100, 0]
        # rocker_b = [2/3 * 100, -50]
        top = thickness * np.sqrt(np.sum((finalPivot - self.rocker_r) ** 2, axis=-1)) # (0, 140)
        left = thickness * np.sqrt(np.sum((finalPivot - self.rocker_b) ** 2, axis=-1)) # (large positive scalar, 4/9 * 10000)
        right = thickness * np.sqrt(np.sum((self.rocker_r - self.rocker_b) ** 2, axis=-1)) # scalar
        front_area = top + left + right + solid_offset ** 2
        f1 = front_area * depth # []
        f1 = (f1 - 1.1645e+04) / (6.4756e+04 - 1.1645e+04)
        return f1 

    def f2(self, finalPivot, solid_offset, thickness, depth):
        top = thickness * np.sqrt(np.sum((finalPivot - self.rocker_r) ** 2, axis=-1))
        left = thickness * np.sqrt(np.sum((finalPivot - self.rocker_b) ** 2, axis=-1))
        right = thickness * np.sqrt(np.sum((self.rocker_r - self.rocker_b) ** 2, axis=-1))
        front_area = top + left + right + solid_offset ** 2
        mass = front_area * depth
        axisLen = np.sqrt(np.sum((finalPivot - self.rocker_r) ** 2, axis=-1))
        
        f2 = axisLen / mass; 
        f2 = (f2 - 0.0015) / (0.01 - 0.0015)
        return f2

    def f3(self, drop_h, stay_depth, tube_h, tube_length):
        # drop_h = scaled_x[:, 3] # (50, 80)
        # tube_h = scaled_x[:, 4] # (10, 20)
        # stay_depth = scaled_x[:, 5] # (6, 16)
        drop_mass = drop_h * drop_h * stay_depth
        tube_mass = tube_h * tube_length * stay_depth
        f3 = drop_mass + tube_mass
        f3 = (f3 - 2.35e+04) / (1.16e+05 - 2.35e+04)
        return f3
    
    def f4(self, tube_length, tube_h, stay_depth):
        f4 = tube_length / (tube_h * stay_depth)
        f4 = (f4 - 41/64) / (55/12 - 41/64)
        return f4

    def lerp(self, t, min_, max_):
        return min_ + t*(max_-min_)

    def rescaleParams(self, x):
        # n_var: 0: solid_offset; 1: thickness; 2: depth; 3: drop_h; 4: tube_h; 5: stay_depth; 6: pivot_offset x; 7 (optional): pivot_offset y
        scaledx = np.empty_like(x)
        scaledx[:, 0] =  self.lerp(x[:, 0], 25, 45)
        scaledx[:, 1] =  self.lerp(x[:, 1], 7, 14)
        scaledx[:, 2] =  self.lerp(x[:, 2], 5, 10)
        scaledx[:, 3] =  self.lerp(x[:, 3], 50, 80)
        scaledx[:, 4] =  self.lerp(x[:, 4], 10, 20)
        scaledx[:, 5] =  self.lerp(x[:, 5], 6, 16)
        scaledx[:, 6] =  self.lerp(x[:, 6], 0, -40)
        return scaledx


def runBike():
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
    allRes = []
    hypervolume = []
    problem = Bike()
    res = minimize(problem,
                algorithm,
                # seed=1,
                pf=problem.pareto_front(use_cache=False),
                save_history=True,
                verbose=True)
    allRes.append(res)
    funcEvals += res.algorithm.evaluator.n_eval
    time += res.exec_time

    # plot the hypervolume over time
    # create the performance indicator object with reference point (1,1)
    metric = Hypervolume(ref_point=np.array([1.0, 1.0, 1.0, 1.0]))

    # collect the population in each generation
    pop_each_gen = [a.pop for a in res.history]

    # receive the population in each generation
    obj_and_feasible_each_gen = [pop[pop.get("feasible")[:,0]].get("F") for pop in pop_each_gen]
    # calculate for each generation the HV metric
    hv = [metric.calc(f) for f in obj_and_feasible_each_gen]
    hypervolume.append(hv[-1])


    toc()
    # print("Total number of function evals over", numSteps, "Pareto fronts:", funcEvals)
    print("Total execution time from solver:", time)
    print("time, funcEvals:", time, funcEvals)

    # scipy.io.savemat('./{}.mat'.format(os.path.basename(__file__).split('.')[0]), stats)

    # Visualize
    # Design Space
    # plot = Scatter(title = "Design Space", axis_labels="x")
    # plot.add(res.X)
    # plot.do()
    # # plot.apply(lambda ax: ax.set_xlim(-0.5, 1.5))
    # # plot.apply(lambda ax: ax.set_ylim(-2, 2))
    # plot.show()

    # Augmented Objective Space
    plot = Scatter(title = "Objective Space", labels=["Mass", "Power", "Wind Speed"])

    import pdb; pdb.set_trace()
    front = allRes[0].F;
    numPts = front.shape[0];
    # z = vwind * np.ones([numPts, 1]);
    FCfront = front# np.concatenate((front, z), 1);
    plot.add(FCfront);

    plot.do()
    plot.apply(lambda ax: ax.set_xlim(0, 1))
    plot.apply(lambda ax: ax.set_ylim(0, 1))
    plot.apply(lambda ax: ax.set_zlim(0, 1))
    plot.show()



    # function evaluations at each snapshot
    n_evals = np.array([a.evaluator.n_eval for a in res.history])

if __name__ == '__main__':
    runBike()
