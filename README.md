Gravitational lensing simulated using discrete time newtonian model.

Build with: gcc -fopenmp -g -lm -o lensing lensing.c linalg.c tracing.c

Render resolution preferences can be found in tracing.h, black hole masses and schwarzchild radii must be manually altered in lensing.c. Be sure to alter N_BODIES and N_OBJECTS when adding black holes and render objects, found in lensing.h and tracing.h respectively.
