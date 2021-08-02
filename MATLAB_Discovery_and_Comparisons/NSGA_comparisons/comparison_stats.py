import numpy as np
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.util.termination.default import MultiObjectiveDefaultTermination

import matplotlib.pyplot as plt
from pymoo.performance_indicator.hv import Hypervolume

from ticTocGenerator import *
from turbine import Turbine 
from lamp import Lamp
import scipy.io
import os


'''
Compute a single fixed-context pareto front
'''
def FC_NSGA(problem):
    tic(); 
    algorithm = NSGA2(pop_size=200,
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

    res = minimize(problem,
        algorithm,
        termination,
        seed=1,
        pf=problem.pareto_front(use_cache=False),
        save_history=True,
        verbose=True)

    fc_feval = res.algorithm.evaluator.n_eval
    fc_time = res.exec_time
    
    # plot the hypervolume over time (REVISIT: we only need HV at last time point)
    # create the performance indicator object with reference point (1,1)
    metric = Hypervolume(ref_point=np.array([1.0, 1.0])) # for turbine
    # metric = Hypervolume(ref_point=np.array([1.0, 1.0, 1.0])) #for lamp
    pop_each_gen = [a.pop for a in res.history] # collect the population in each generation
    obj_and_feasible_each_gen = [pop[pop.get("feasible")[:,0]].get("F") for pop in pop_each_gen]  # receive the population in each generation

    # calculate for each generation the HV metric
    all_hv = [metric.calc(f) for f in obj_and_feasible_each_gen]
    fc_hv = all_hv[-1] # only save final population's hv

    return fc_time, fc_feval, fc_hv

'''
Compute statistics for a single approximation of the gamut 
(repeated calls to the fixed context method, at the desired contextual sampling frequency)
'''
def gamutApproximation(problemClass, contextRange, numFCFronts):
    zmin = contextRange[0]
    zmax = contextRange[1]
    contextStepSize = (zmax-zmin)/(numFCFronts);

    zmin = zmin + 0.5*contextStepSize;
    zmax = zmax - 0.5*contextStepSize; # center the samples in the range

    stats = dict()
    stats['time'] = []
    stats['num_evals'] = []
    stats['HV'] = []
    
    totalFuncEvals = 0
    totalTime=0
    for i in range(numFCFronts):
        z = zmin + contextStepSize*i #current context
        
        # run a single fixed-context instance of the problem with NSGA
        problemFunc = problemClass(z)
        time, funcEvals, hypervolume = FC_NSGA(problemFunc)
        
        # collect stats
        stats['time'].append(time)
        stats['num_evals'].append(funcEvals)
        stats['HV'].append(hypervolume)
        totalFuncEvals += funcEvals
        totalTime += time

    stats['totalTime'] = totalTime
    stats['total_num_evals'] = totalFuncEvals

    return stats


'''
Collect full statistics for this problem over N repeated experiments
'''
def collectNSGAStats(problemClass, experimentRepetitions, numFCFronts, contextRange):
    time = np.empty([experimentRepetitions, numFCFronts])
    nevals = np.empty([experimentRepetitions, numFCFronts])
    hv = np.empty([experimentRepetitions, numFCFronts])
    totalTime = np.empty([experimentRepetitions, 1])
    totalEvals = np.empty([experimentRepetitions, 1])

    for i in range(experimentRepetitions):
        stats = gamutApproximation(problemClass, contextRange, numFCFronts)
        
        time[i, :] = stats['time']
        nevals[i, :] = stats['num_evals']
        hv[i, :] = stats['HV']
        totalTime[i] = stats['totalTime']
        totalEvals[i] = stats['total_num_evals']
        
    finalStats = dict()
    finalStats['time'] = time
    finalStats['num_evals'] = nevals
    finalStats['HV'] = hv
    finalStats['totalTime'] = totalTime
    finalStats['total_num_evals'] = totalEvals

    filename = problemClass.__name__ + "_numFronts" + str(numFCFronts) + "_NSGA"
    scipy.io.savemat('./{}.mat'.format(filename), finalStats)


'''
SET PARAMETERS
'''

# parser = argparse.ArgumentParser(description='Collect statistics for NSGA')
# parser.add_argument('--problemName', dest='problemClass', type=str, required=True,
#                     help='name of problem to evaluate (Turbine, Lamp, Bicopter)')
# parser.add_argument('--numTrials', dest='experimentRepetitions', type=int,
#                     default=1, help='number of times to repeat gamut identification')
# parser.add_argument('--numFronts', dest='numFCFronts', type=int,
#                     default=10, help='number of fixed context fronts to compute')
# parser.add_argument('--contextMin', dest='zmin', type=float,
#                     default=0.0, help='minimum context value')
# parser.add_argument('--contextMax', dest='zmax', type=float,
#                     default=1.0, help='maximum of context value')

# args = parser.parse_args()

# problemDict = {
#     "Turbine" : Turbine, 
#     "Lamp" : Lamp
# }
# if args.problemClass in problemDict:
#     probClass = problemDict[args.problemClass]
#     contextRange = [args.zmin, args.zmax]
#     collectNSGAStats(probClass, args.experimentRepetitions, args.numFCFronts, contextRange)

# else:
#     print('Invalid problem class. Valid options are Turbine, Lamp')

experimentRepetitions = 1
contextRange = [0, 1.0]
problemClass = Turbine
# problemClass = Lamp

# for i in range(5, 201, 5): %% FOR FULL COMPARISON, TAKES A LONG TIME
for i in range(5, 6, 5):
    numFCFronts = i
    collectNSGAStats(problemClass, experimentRepetitions, numFCFronts, contextRange)
