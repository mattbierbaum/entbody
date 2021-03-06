This is the full version of entbody with all features included.

Options in the Makefile so that it works easily
with other systems (edit as necessary):
- CUDA   - compile and run on CUDA device
- DOPLOT - display with opengl
- FPS    - calculate frames per second (not portable)
- POINTS - make the opengl display points, not fancy circles
- OPENMP - use threading through openmp
- TEMPS  - since we depend on preprocessors, output the source before gcc

How to build:
    * need ENTBODY env variable set to the root directory
        export ENTBODY=`pwd`

    * make and run the examples (base program doesn't compile) 
        make

Other addons:
    * if you would like to use plotting, install freeglut, OpenGL 
        sudo apt-get install freeglut3-dev

    * image save and initial conditions come from OpenIL (DevIL)
        sudo apt-get install libdevil-dev

    * sound output comes from OpenAL
        sudo apt-get install libopenal-dev

    * multithreading is available through OpenMP
        > already installed with gcc

    * advanced computing on GPU through OpenCL
        > install OpenCL as appropriate for your system

