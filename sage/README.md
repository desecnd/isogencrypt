# Isogencrypt SageMath Package

SageMath auxiliary package and collection of tools for isogencrypt - isogeny-based cryptographic research library.

## SageMath

This package requires `SageMath` version `10.4` installed on the system to run properly. For specific instructions how to obtain SageMath on your system, refer to [this guide](https://doc-gitlab.sagemath.org/html/en/installation/index.html).


## Structure

```bash
sage
├── benches
├── isogencrypt_sage
├── poc
├── scripts
├── tests
├── verifiers
└── pyproject.toml
```

1. `benches` - Benchmarks that measure the performance of isogeny protocols. Each entry corresponds with `.c` benchmark implementation for comparison
2. `isogencrypt_sage` - Flat python package with `sage` as dependency, that contains core isogeny functions and protocol implementations, used by other submodules.
3. `poc` - "Proof of Concept" SageMath-only scripts (with no other dependencies, should run with plain `sage` call) for each of the implemented isogeny-based protocol
4. `scripts` - Helper functionalities like generating benchmark testcases, converting the output formats, plotting the results, etc.
5. `tests` - Unit Tests for `isogencrypt_sage` 
6. `verifiers` - Servers both as unit tests (executed by pytest) and as a template for comparison with low-level `.c` implementation
7. `pyproject.toml` contains the configuration properties of the python package. 

## Run

In order to correctly run all scripts that depend on the `isogencrypt_sage` module with `sage -python` it is recommended to add directory with the package to `PYTHONPATH`. 

```bash
# Add current directory to PYTHONPATH
$ export PYTHONPATH=$PYTHONPATH:`pwd` 

# Check if Sage python executable can import the package.
$ cd scripts
$ sage -python -c "import isogencrypt_sage"
$ sage -c "import isogencrypt_sage"
```

Additional solution is to create virtualenv and install the package in editable mode, but it may lead to conflicts with `sage -python` version and their dependencies.

## Testing

Both unit tests in `tests` and `verifiers` can be executed with `pytest` module invocation.

```bash
# Run the pytest module with sagemath's python executable
$ sage -python -m pytest
```