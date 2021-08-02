Pareto Gamuts: Exploring Optimal Designs Across Varying Contexts
SIGGRAPH 2021

# Running an existing example
-----------------------------------

All examples can be run through `examplesMeshBased.m` or `examplesAnalytic.m`

The input/output of each problem call is detailed below:
    @param problemID
        the problemIDs correspond to various example files, written in analytic/ and meshBased/   
        AppMappingFunction details the problemID needed to run each example (but they're already pre-filled in the exampleX.m files)

    @return stats
    	statistics about the run including number of function evals, timings in seconds, hypervolume measurements over time, etc.

    @return buffer
        the discovery buffer, containing all patches and corresponding information found along the way. If interested, see `AppBuffer.m` for more details; else, feel free to ignore this.

    @return des_figh   
        figure handle for the final design space plot, so you can change it after the fact (eg, add comparison to other algorithms' results)  

    @return perf_figh   
        figure handle for the final performance space plot, to alter after the fact   

    @return pts  
        n x (D+A+d) array containing all points on the extracted pareto front   
        Each row corresponds to a single pareto optimal solution, specified as [design, application, performance]   

    @return labels   
        n x 1 array containing the family/class/patch label for each point in pts   
        labels(i) corresponds to the solution given by pts(i, :)   

    @ filledPtsArray
    	points that were filled in using the sampled continuous patches after extracting the pointwise gamut from the buffer

	@idxToVis
		dictionary of variable indices to visualize along each axis of the pareto gamut

Each problem call finds the pareto front, then generates 2 plots showing the extracted pareto gamut + its pre-image



### Before Running
- you will need to change the path `basewd` near the top of `AppDiscover.m`  
- ensure that `prepJSvis = false` in `AppDiscover.m` 
   - the code for the JS visualizer is provided with pre-populated data to explore the results from the video. If you wish to move computed data from MATLAB to the JS viewer automatically, you will have to change the JS vis directory path variable in meshBased/prepForInteractiveVis.m . We do not provide any additional documentation on this process, but it should be discoverable using the samples provided.
- you may want to change the userParams (e.g., stepSize or maxRuns)

  



# Creating a new problem description
-----------------------------------

1. create a folder with your new example (probably in either `analytic/` or `meshBased/`, although you can make a new one too)
   - inside this folder, create a file/class called ProblemFolderName.m   
   -copy the contents of `analytic/ProblemTemplate.m`, and follow the instructions to define your problem
2. create a new functionID/entry for your example in `AppMappingFunction.m`
    
### Hints

- You will need to define the design + application variables, performance metrics and (if applicable) constraints 
- You may also need to change the createProblem function (member of the problem class)

- the variables (eg xSymb, zSymb) must be between [0,1] in the optimization
   - create the aliased variables + use the lerp function to map this to your desired range  
- performance metrics (or at least the points along the front) must evaluate between [0,1] as well
   - this is so scalarization makes sense, and it's fully contained in the buffer
   - to get a sense of a metric's spread/normalization factor, you can run the Discover function 
     and look at the plot of all discovered patches to get a sense of the spread/normalization factor. 
     The extracted pareto front and statistics may be incorrect if the points are outside [0,1])  
- your performance metrics should be able to take in and evaluate *multiple* points on a single call    
    - if you use analytic functions (without matrix equations) this is done for you; don't worry about it.  
    - if you have matrix equations or otherwise can't evaluate multiple at once, 
            you need to make sure the points are evaluated one at a time

