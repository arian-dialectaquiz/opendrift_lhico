# opendrift_lhico installing
LHiCo's implementation of OpenDrift to Cananeia and Santos  (ECOM) and Amazon (ROMS) estuaries and South Brazil Bight (ECOM and ROMS).
Opendrift version 1.11


1. install Mambaforge:

    https://github.com/conda-forge/miniforge#mambaforge

2. Clone this github repository:

    git clone https://github.com/arian-dialectaquiz/opendrift_lhico.git

4. Create environment with the required dependencies and install OpenDrift

    cd opendrift_lhico/
    mamba env create -f environment.yml
    conda activate opendrift
    pip install --no-deps -e .

5. Choose the eulerian model to be used:
    ecom
    roms

 4.1. Check if the necessary grid data has been downloaded. This is needed to mask land data in ECOM and to make angle rotation in ROMS.

    cd ecom/opendrift/readers/ECOM_dependencies

    cd roms/opendrift/readers/roms_grid

  4.2. Understand the readers and change as will, but keep the stable version unchanged.

    cd ecom/opendrift/readers/
      reader_ECOM_S2Z_many_grids_stable
    cd roms/opendrift/readers/
      reader_roms_lhico
      
5. See the run cases to further examples, or the 'How to make a run' file to the first steps.

6. This framework was based on OpenDrift model (https://opendrift.github.io/index.html). 
    
    Check the original github to theoretical background, and subscribe to the forum, where the developers are constantly interacting and solving issues.
