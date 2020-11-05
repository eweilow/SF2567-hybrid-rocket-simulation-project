# SF2567-hybrid-rocket-simulation-project

SF2567 - Project Course in Scientific Computing, simulating a hybrid rocket engine

## Dependencies

To manage dependencies in a nice way, this project is based around a single Dockerfile. One must thus run all the code inside Docker, but that's a useful tradeoff to manage the native dependencies. The ones that come with the Dockerfile are:

- Native dependencies:
  - Build tools:
    - wget
    - cmake + make
    - gcc + gfortran
    - libopenblas-dev + liblapack-dev
    - [sundials](https://computing.llnl.gov/projects/sundials)
  - [coolprop](http://www.coolprop.org/) (installed with pip)
  - [NASA CEA](https://www.grc.nasa.gov/WWW/CEAWeb/ceaHome.htm) (installed with pip)
- pip dependencies:
  - [numpy](https://pypi.org/project/numpy/)
  - [scipy](https://pypi.org/project/scipy/)
  - [matplotlib](https://pypi.org/project/matplotlib/)
  - [scikits.odes](https://pypi.org/project/scikits.odes/)
  - [coolprop](https://pypi.org/project/coolprop/)
  - [rocketcea](https://pypi.org/project/rocketcea/)
  - [Cython](https://pypi.org/project/Cython/)

## Get started

Start by building the main project Dockerfile. This has all the dependencies necessary already installed. The tag (`-t sf2567-python`) is important here!

```bash
docker build -t sf2567-python .
```

The build should now be saved as a Docker image on your computer, to be used by the different parts of this project.

To run each part of the simulation, one can use the command

```bash
bash run.sh && python plot.py
```

This should be mostly cross-compatible between platforms (and works in Windows if one has WSL enabled).
