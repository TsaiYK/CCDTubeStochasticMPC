# Control Co-Design (CCD) using Tube-based Stochastic MPC

These codes provide the stochastic model predictive control (SMPC) solution with the generated images and animations.
![Stochastic tube-based MPC on a state plane.](results_mpc_test/GIF_plane.gif)

1. 'main_sim.m': The script implements a Stochastic Model Predictive Control (SMPC) simulation for a linear system with disturbances, including state and control constraint tightening, dynamic iteration, and visualization of results through plots and animated GIFs.
2. 'ConstraintTigtening.m': This MATLAB class handles constraint tightening for Model Predictive Control (MPC) with stochastic disturbances. It computes nominal constraint parameters, time-varying constraint sets, and terminal constraint sets for state and input variables. The class includes methods for finding nominal constraints, constructing time-varying constraint sets, and calculating the nominal terminal set.
3. 'ellipsedata.m': This function generates data points for ellipses representing contour curves of Gaussian distributions with specified covariance and mean. It's useful for visualizing confidence intervals in 2D space. The function takes inputs for covariance matrix, center point, number of points, and sigma rule, returning a matrix of ellipse data points.
4. 'find_multidim_contour.m': This function calculates confidence intervals (CIs) for 90%, 95%, and 99% in 2D or 4D space, given a mean vector and covariance matrix. For 2D, it uses the ellipsedata function to generate ellipses. For 4D, it implements a custom contour4D function to create hyperspheres. The output provides coordinates for plotting these confidence intervals.
5. 'PlotGenerator.m': This function generates various plots for visualizing MPC results, including state and control trajectories, confidence intervals, and constraint violations. It creates figures showing nominal trajectories, actual trajectories (if specified), and confidence intervals for both state and control variables. The function also calculates and returns probabilities of constraint violations for states and controls.
6. 'animated_gif_creator.m': This function creates an animated GIF from a series of image files, allowing customization of loop count, delay times between frames, and providing a progress bar during the creation process

You can install the MPT3 package and get the folder 'tbxmanager' via https://www.mpt3.org/.

For more details about the method and case study, please read the paper: Tsai, Y.-K. and Malak, R., 2025, “Control Co-Design with Performance-Robustness Trade-Off using Tube-Based Stochastic Model Predictive Control,” *ASME Letters in Dynamic Systems and Control*.

If you have any further question, please contact the developer by yktsai0121@tamu.edu.
