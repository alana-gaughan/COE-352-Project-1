# COE-352-Project-1

## Results
When the code for SVD decomposition is run, it compares my find_SVD() function with functions from various other Python libraries. It tests the results on 3 matrices of 3 different sizes. Below is the reformatted printed output.

Comparing my Function vs SciPy Routines:

```
Matrix:
[[ 3  3  2]
 [ 2  3 -2]]
 
My SVD:
U=[[ 0.7815437  0.6238505] S=[[5.54801894 0.         0.        ] V=[[ 0.64749817  0.10759258  0.75443354]
   [ 0.6238505 -0.7815437]]   [0.         2.86696457 0.        ]]   [ 0.7599438  -0.16501062 -0.62869461]
                                                                    [ 0.05684667  0.9804057  -0.18860838]]


SciPy SVD:
U=[[-0.7815437 -0.6238505]  S=[5.54801894 2.86696457] V=[[-0.64749817 -0.7599438  -0.05684667]
   [-0.6238505  0.7815437]]                              [-0.10759258  0.16501062 -0.9804057 ]
                                                         [-0.75443354  0.62869461  0.18860838]]

My Inverse:
[[ 0.11462451  0.04347826]
 [ 0.07114625  0.13043478]
 [ 0.22134387 -0.26086957]]

Scipy Inverse:
[[ 0.11462451  0.04347826]
 [ 0.07114625  0.13043478]
 [ 0.22134387 -0.26086957]]

My Condition Number: 1.93515434599516
Numpy Condition Number: 1.9351543459951601

Matrix:
[[ 3  3]
 [ 2  3]
 [-2  2]]

My SVD:
U=[[-0.7599438  -0.16501062  0.62869461] S=[[5.54801894 0.        ] V=[[-0.6238505 -0.7815437]
   [-0.64749817  0.10759258 -0.75443354]    [0.         2.86696457]    [-0.7815437  0.6238505]]
   [-0.05684667  0.9804057   0.18860838]]   [0.         0.        ]]


SciPy SVD:
U=[[-0.7599438   0.16501062  0.62869461] S=[5.54801894 2.86696457] V=[[-0.6238505 -0.7815437]
   [-0.64749817 -0.10759258 -0.75443354]                              [ 0.7815437 -0.6238505]]
   [-0.05684667 -0.9804057   0.18860838]]  

My Inverse:
[[ 0.13043478  0.04347826 -0.26086957]
 [ 0.07114625  0.11462451  0.22134387]]

Scipy Inverse:
[[ 0.13043478  0.04347826 -0.26086957]
 [ 0.07114625  0.11462451  0.22134387]]

My Condition Number: 1.93515434599516
Numpy Condition Number: 1.9351543459951595

Matrix:
[[3 4]
 [0 5]]
 
My SVD:
U=[[-0.70710678 -0.70710678] S=[[6.70820393 0.        ] V=[[-0.31622777 -0.9486833 ]
   [-0.70710678  0.70710678]]   [0.         2.23606798]]   [-0.9486833   0.31622777]]

SciPy SVD:
U=[[ 0.70710678 -0.70710678]  S=[6.70820393 2.23606798] V=[[ 0.31622777  0.9486833 ]
   [ 0.70710678  0.70710678]]                              [-0.9486833   0.31622777]]

My Inverse:
[[ 3.33333333e-01 -2.66666667e-01]
 [ 8.82771154e-17  2.00000000e-01]]

Scipy Inverse:
[[ 3.33333333e-01 -2.66666667e-01]
 [-1.98965378e-17  2.00000000e-01]]

My Condition Number: 2.9999999999999987
```
Numpy Condition Number: 3.000000000000001



