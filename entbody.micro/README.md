Introduction
------------
This is the smallest implementation of entbody that I could come up with.
It is a single file at the smallest, still including parts of interaction
that do not always need to be there.

Features
--------
- 2D only with periodic boundary conditions determined by an array pbc = [0,0]
- Fast plotting library through OpenGL (freeglut) and OpenIL (saving images)
- Keyboard interactions easy to add
- Neighborlist creation. Forces are not assumed to be Newtonian
- Fancy coloring based on local stress

How to use
----------
To compile use make.  The options for plotting and FPS are provided in the Makefile.
