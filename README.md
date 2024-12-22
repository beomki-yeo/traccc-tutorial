# traccc-tutorial

### Description

This is a repository for traccc tutorial. 
In case one having a problem in compilation, please let the devs know by filing an issue or sending an email.

### Prerequisites
- gcc with C++ 20 support
- CMake >= 3.22
- (Optional) CUDA >= 12.4

[Perlmutter](https://docs.nersc.gov/systems/perlmutter/architecture/) users can compile the tutorial project by running the following commands (CPU-only)
          
```              
module load gcc/12.2.0
module load cmake/3.30.2
```

### CMake Build Options

| Option | Description |
| --- | --- |
| BUILD_CUDA  | Build the CUDA tutorials |
