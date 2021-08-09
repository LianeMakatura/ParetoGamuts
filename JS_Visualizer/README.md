# Pareto Gamuts: Exploring Optimal Designs Across Varying Contexts

Liane Makatura, Minghao Guo, Adriana Schulz, Justin Solomon, and Wojciech Matusik

SIGGRAPH 2021

This repository contains the code to accompany our SIGGRAPH 2021 paper, 
[Pareto Gamuts: Exploring Optimal Designs Across Varying Contexts](http://paretogamuts.csail.mit.edu). 
For more details, please refer to the paper.

## Running the visualizer

The JS visualizer (shown in the video) comes pre-loaded with all necessary data to recreate the examples in our paper, so you do not need to run our MATLAB code.
The visualizer is best used in Chrome. It does not work in Safari.

The visualizer is interactive or nearly interactive for all problems once it's loaded; loading may take a few seconds.
The turbine is not quite interactive because of the large amount of data and semi-intensive jscad model; it takes a few seconds to update between point selections. 

## Prepare the data files

To conserve space, we have compressed all the data files (e.g. `fronts.json(.zip)`, `designExp.json(.zip)`) inside each problem class. To uncompress them, cd into this directory (`JS_Visualizer`) and run:

`python prepare.py decompress`

This will automatically loop over the compressed files and extract the full version inside the proper directory.

## Set up docker

Since our gamuts contain thousands of individual designs, it is intractable to generate and/or share the corresponding meshes a priori. Instead, we use JSCAD to generate the corresponding meshes on-the-fly. We provide simple instructions to set up the docker image in `./jscad-server`

## Selecting a problem

You can toggle between the problems by (un)commenting lines in `vis/js/settings.js`

## Launching the visualizer

We recommend using VSCode with the LiveServer extension; simply hit "Go Live" in the bottom right corner to start the viewer (and copy the link to chrome if that's not your default browser)

## Using the interface

You can select points to examine from any 2D design or performance plot (not 3D).
To view multiple designs simultaneously, press Control. To go back to single design mode, press Control.
To download a mesh, use `d`

## JS Visualizer: Contributors 

This code was written by Hannes Hergeth and Liane Makatura.