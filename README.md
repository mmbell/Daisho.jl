# Daisho.jl

Daisho is a data analysis and assimilation package written in Julia for atmospheric science applications. The code is being developed to provide functionality similar to the FRACTL and SAMURAI C++ software with additional polarimetric radar processing capabilities. The name comes from the Japanese term meaning "large and small" and is the name of matched pair of traditionally made Japanese swords (nihonto) worn by the samurai in ancient Japan. The software is in active development and is still in pre-release. 

### Installation

After cloning this repository, start Julia using Daisho.jl as the project directory. This can be done on the command line using `julia --project` or set the JULIA_PROJECT environmental variable:

`export JULIA_PROJECT=/path/to/Daisho.jl`

To install Daisho, in the REPL, go into Package mode by pressing `]`. You will see the REPL change color and indicate `pkg` mode. 

If you are actively developing or modifying Daisho then you can install the module using `dev /path/to/Daisho.jl` in `pkg` mode. This will update the module as changes are made to the code. You should see the dependencies being installed, and then the Daisho package will be precompiled. Exit Package mode with ctrl-C.

If you wish to just install a static version of the latest release, run `activate` to activate the package environment. Then, run `instantiate` to install the necessary dependencies. Exit Package mode with ctrl-C.

Test to make sure the precompilation was successful by running `using Daisho` in the REPL. If everything is successful then you should get no errors and it will just move to a new line.

### Running Daisho

The primary current functionality is a novel gridding approach for weather radars on moving platforms that takes into account beam characteristics. Documentation is limited at the current time as the software is still in its early stages.

### Future plans
Multi-Doppler and variational analysis are in development, as is integration with the Springsteel.jl framework. Interested users are welcome to contribute to improve the code. Stay tuned for more functionality!
