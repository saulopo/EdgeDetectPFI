# EdgeDetectPFI
A Fortran 90 implementation of the paper *EdgeDetectPFI: an algorithm for  automatic edge detection in potential field anomaly images - application to  dike-like magnetic structures*, published in Computers &amp; Geosciences

[http://dx.doi.org/]

**Authors:**  
* Saulo P. Oliveira [saulopo@ufpr.br]
* Francisco J. F. Ferreira [francisco.ferreira@ufpr.br]
* Jeferson de Souza [jdesouza@ufpr.br]
         
## Instructions:

1. Compile the file EdgeDetect.f90. Please use an optimization compiler  option, if available. In a linux distribution with gfortran installed, use

 gfortran -O3 EdgeDetectPFI.f90

2. Make sure the file input.dat provided in this package is in the same folder as the fortran file.
   If you intend to use your own input data, please save it as a CSV (Comma Separated Values) file with four columns, namely:
   - Column 1: x coordinates of the i-th grid point 
   - Column 2: y coordinates of the i-th grid point
   - Column 3: Signum transform of Mz at the i-th grid point
   - Column 4: Signum transform of Mz-|Mh| at the i-th grid point

3. execute a.out:

 ./a.out

4. The executable file will generate the CSV file output.csv. View this file on matlab/octave with the auxiliary script ViewCSV.m. If you intend to use your own graphic viewer, the output file is organized as follows:
   - Column 1: x coordinates of the i-th grid point 
   - Column 2: y coordinates of the i-th grid point
   - Column 3: source depth at the i-th grid point
   - Column 4: source width at the i-th grid point
   
   Please note that the grid points are given at the sources only.
