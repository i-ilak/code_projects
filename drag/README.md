# Drag

This project was motivated by [this picture](https://static.boredpanda.com/blog/wp-content/uploads/2018/02/dwayne-the-rock-johnson-skyscraper-jump-funny-reactions-1-5a7ab25b4d416__700.jpg) of "The Rock" jumping between skyscrapers. We'd like to know if this jump is possible and what are limiting factors. Note that we will be reusing the `integrate.hpp` from the [n_body](https://github.com/i-ilak/code_projects/tree/main/n_body/src) project, since we will need an ODE solver.

### Current status of the project:
Runs and produces reasonable results, but there seems to be a numerical instability: depending on the initial conditions we get only roughly half the possible solutions (all permissable angles bigger or smaller than zero, but not the whole range). Some more investigation is necessary here. Also, need to add testing. 
