# VegasTables
C++ code to perform bicubic interpolation (using GSL) on tabular data to accurately interface micropysics calculations with hydro codes.  It was specifically developed to utilize the tables [here](http://www.physics.unlv.edu/astro/xstartables.html).

The code is self-contained in a .h file so that it can easily be called within a hydro code.  For example, to interface with Athena++, create folder named 'user' in /src and place vegas_tables.h there, and then call it from any pgen using
```
#include "../user/vegas_tables.h"  
```

The example driver bicubic_speed_test.cpp compares the speed of 1000 random evalutions using both bilinear and bicubic interpolation on the sample table Blondin_100by200.dat.  Compile and run it using
```
g++ bicubic_speed_test.cpp -o place_bets -lgsl -lm
./place_bets
```

For further documentation, refer to the header information in vegas_tables.h

There is a C-wrapper for this code - contact me to get a hold of that.
