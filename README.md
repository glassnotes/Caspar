# Caspar 
A Python implementation of the SU(_n_) factorization scheme of [citation forthcoming].

If you run into any trouble, or find a bug, feel free to drop me a line at
odimatte@uwaterloo.ca or create an issue in the repo.
                                                                                   
## Dependencies                                                                    
                                                                                   
You will need the numpy package. Caspar was developed using Python 3.5.2 but
also was tested and runs smoothly in 2.7.10. 

## Installation

In the main directory of the folder, type `python setup.py install`.
                                                                                   
## Usage (basic)

There are two important files: `factorization_script.py`, and `user_matrix.py`. 
Enter the SU(_n_) matrix you wish to factorize in the variable `SUn_mat` in 
`user_matrix.py`. Then, factorize it by running 
```
python factorization_script.py
```

The output of this script will be a series of lines in the following format, for example:
```
4,5    [-2.8209, 2.5309, 2.3985]
3,4    [-1.7534, 1.4869, -1.753]
...
```
This is a sequence of SU(2) transformations. The first two integers indicate 
the modes on which the transformation acts. The set of three floats are the 
parameters of the transformation (see parametrization below). 
The original matrix `SUn_mat` is obtained by embedding each SU(2) transformation 
into the indicated modes of an SU(_n_) transformation, and multiplying them 
together from top to bottom of the list (with each transformation added to 
the product on the right, e.g. _U_ = _U_<sub>45</sub> _U_<sub>34</sub>...).  

_Important note_: at the current time, Caspar works well for up to about _n_ = 10,
if the matrix is a Haar-random, dense special unitary. After this point, it 
begins to suffer from issues due to precision. However, for sparser matrices,
e.g. the m-qubit Paulis, it has seen success for up to 6 qubits (n = 2^6 = 64,
and this is just as high as I tested).


## Usage (detailed)
                                                         
An arbitrary element of SU(_n_) can be fully expressed using at most 
_n_<sup>2</sup> - 1 parameters. We put forth a factorization scheme that 
decomposes elements of SU(_n_) as a sequence of SU(2) transformations. 
SU(2) transformations require in general 3 parameters, [_a_, _b_, _g_], 
written in matrix form as  [[_e_<sup>_i_(_a_+_g_)/2</sup> cos(_b_/2), 
-_e_<sup>_i_(_a_-_g_)/2</sup> sin(_b_/2)], 
[_e_<sup>-_i_(_a_-_g_)/2</sup> sin(_b_/2), _e_<sup>_i_(_a_+_g_)/2</sup> cos(_b_/2)]].
                                                                                   
There are two main functions: `sun_factorization` and `sun_reconstruction`, 
each contained in the appropriately named files.
                                                                                   
The function `sun_factorization` takes an SU(_n_) matrix (as a numpy matrix) 
and decomposes it into a sequence of _n_(_n_-1)/2 such SU(2) transformations. 
The full set of _n_<sup>2</sup> - 1 parameters is returned as a list of tuples 
of the form ("_i_,_i_+1", [_a_<sub>_k_</sub>, _b_<sub>_k_</sub>, _g_<sub>_k_</sub>]) 
where _i_ and _i_+1 indicate the modes on which the transformation acts (our
factorization uses transformations only on adjacent modes).

The following code snippet can be used to factorize the SU(3) matrix below.

```python
import numpy as np
from caspar import sun_factorization 

n = 3

SUn_mat = np.matrix([[0., 0., 1.],                                               
                    [np.exp(2 * 1j * np.pi/ 3), 0., 0.],                        
                    [0., np.exp(-2 * 1j * np.pi / 3), 0.]])    

# Perform the decomposition
parameters = sun_factorization(SUn_mat)

# The output produced is 
#
# Factorization parameters: 
#   2,3    [2.0943951023931953, 0.0, 2.0943951023931953]
#   1,2    [0.0, 3.1415926535897931, 0.0]
#   2,3    [0.0, 3.1415926535897931, 0.0]
```

It is also possible to reconstruct an SU(_n_) transformation based on a list 
of parameters for SU(2) transformations given in the form 
("_i_,_i_+1", [_a_<sub>_k_</sub>, _b_<sub>_k_</sub>, _g_<sub>_k_</sub>]). 
The matrix is computed by multiplication on the right. At the moment only 
adjacent mode transformations are supported.
```python
from caspar import sun_reconstruction

parameters = [("1,2", [1.23, 2.34, 0.999]), 
              ("4,5", [0.3228328, 0.23324, -0.2228])]

new_SUn_mat = sun_reconstruction(6, parameters)
```
