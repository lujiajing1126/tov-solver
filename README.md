# TOV Solver

> This is a high performance Python C extension to solve TOV equation

## Requirements

  - `Python 3`
  - `numpy`
  - `gsl`

## Install

```
$ pip install git+https://github.com/lujiajing1126/tov-solver.git@master
```


## Use
```python
from tov import solve
# rho(numpy.ndarray): baryon density
# pressure(numpy.ndarray): pressure of the beta-equilibrium matter at specfic baryon density
# eps(numpy.ndarray): total energy density at specfic baryon density
# Every array should be in reverse order, and with the same dimension
# Pressure should be monotonically decreasing
solve(rho, pressure, eps)
```

## Refs

  - (2008) S. L. Shapiro                   978-3-527-61767-8
  - (2006) O. E. Nicotra et al.            [10.1051/0004-6361:20053670](https://doi.org/10.1051/0004-6361:20053670)
  - (2008) Z. H. Li and H.-J. Schulze      [10.1103/PhysRevC.78.028801](https://doi.org/10.1103/PhysRevC.78.028801)
  - (2010) G. F. Burgio and H.-J. Schulze  [10.1051/0004-6361/201014308](https://doi.org/0.1051/0004-6361/201014308)
  - (2012) Z. H. Li et al.                 [10.1088/0256-307X/29/1/012101](https://doi.org/10.1088/0256-307X/29/1/012101)
  - (2015) B. K. Sharma et al.             [10.1051/0004-6361/201526642](https://doi.org/10.1088/0256-307X/29/1/012101)

## LICENSE

MIT