# sr
Fooling around with SR encodings. Create a conda environment as
```
conda create -n sr cmake fmt
```
then, after activation and in the repo main folder:
```
mkdir build
cd build
cmake ../
make
```

Executables:
## koza
this presents the basic idea of using the new genotype phenotype mapping to get the koza quintic.

Usage:
```
./koza n_trials verbosity
```
Example: 
```
./main 200  0
```
Will produce the ERT over 200 runs

## autodiff
this contains also an idea to implement automated differentiation up to second order in the most efficient (albeit tedious) way.

Usage:
```
./autodiff
```
