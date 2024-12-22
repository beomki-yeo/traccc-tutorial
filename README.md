# traccc-tutorial

### Description

This is a repository for traccc tutorial. 
In case one having a problem in compilation, please let the developers know by filing an issue or sending an email.

### Prerequisites
- gcc with C++ 20 support
- CMake >= 3.22
- (Optional) CUDA >= 12.4

### CMake Build Options

| Option | Description | Default |
| --- | --- | --- |
| BUILD_CUDA  | Build the CUDA tutorials | OFF |

### Setup in Perlmutter

[Perlmutter](https://docs.nersc.gov/systems/perlmutter/architecture/) users can compile the tutorial project by running the following commands in the login node.
          
```              
module load gcc/12.2.0
module load cmake/3.30.2
module load cudatoolkit/12.4
```

To use A100, configure the project with a right architecture number:

```
cmake <project_directory> -DBUILD_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=80
```
