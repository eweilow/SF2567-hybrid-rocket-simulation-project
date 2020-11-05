FROM python:3.9.0-buster

WORKDIR /

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install wget -y
RUN apt-get install cmake make -y
RUN apt-get install gcc gfortran -y
RUN apt-get install libopenblas-dev liblapack-dev -y

## Build sundials
RUN wget -c https://computing.llnl.gov/projects/sundials/download/sundials-5.1.0.tar.gz -O - | tar -xz
# Jump into the sundials build directory
WORKDIR /build-sundials-5.1.0
RUN cmake -DLAPACK_ENABLE=ON -DSUNDIALS_INDEX_SIZE=64 -DCMAKE_INSTALL_PREFIX=/sundials-install /sundials-5.1.0
RUN make install


## Jump back to root
WORKDIR /


## Common dependencies
RUN pip install numpy
RUN pip install scipy
RUN pip install matplotlib


## Install scikits.odes (which is used mainly for DAEs)
ENV SUNDIALS_INST=/sundials-install
RUN pip install scikits.odes
# Install nose, which is used to test scikits.ode
RUN pip install nose
# Set the right library path for sundials
ENV LD_LIBRARY_PATH=$SUNDIALS_INST/lib:$LD_LIBRARY_PATH
# Print the version of scikits.odes
RUN python -c 'import pkg_resources; print(pkg_resources.get_distribution("scikits.odes").version)'
# Run scikits.odes tests
RUN python -c 'import scikits.odes as od; od.test()'


## rocketcea is used for chemical combustion modelling
RUN pip install rocketcea


## CoolProp is used for propellant properties
# Cython is necessary to install CoolProp
RUN pip install Cython
RUN git clone https://github.com/CoolProp/CoolProp --recursive
WORKDIR /CoolProp/wrappers/Python
RUN apt-get install git g++ p7zip libpython-dev -y
RUN python setup.py install


# Set the main working directory for the Dockerfile
WORKDIR /project