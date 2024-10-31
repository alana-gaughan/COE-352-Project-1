# COE-352-Project-1

## Results
When the code for SVD decomposition is run, it compares my find_SVD() function with functions from various other Python libraries. It tests the results on 3 matrices of 3 different sizes. Below is the reformatted printed output.

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
```
```
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
```
```
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
Numpy Condition Number: 3.000000000000001
```
For the next problem, we solved the equilibrium spring mass system with different boundary conditions. All cases below have all three massed = 1, and all spring constants = 1 as well (though different cases have a different number of springs so that the problem makes sense). Below are the properly formatted results:
```
Case 1: Two fixed ends:
Displacement of masses = [14.715 19.62  14.715]
Elongation of springs = [ 14.715   4.905  -4.905 -14.715]
Internal Stresses = [ 14.715   4.905  -4.905 -14.715]
Condition Number = 5.828427124746198

Case 2: Top end fixed, bottom end free:
Displacement of masses = [29.43 49.05 58.86]
Elongation of springs = [29.43 19.62  9.81]
Internal Stresses = [29.43 19.62  9.81]
Condition Number = 16.393731622284225

Case 3: Top end free, bottom end fixed:
Displacement of masses = [58.86 49.05 29.43]
Elongation of springs = [ -9.81 -19.62 -29.43]
Internal Stresses = [ -9.81 -19.62 -29.43]
Condition Number = 16.393731622284296

Case 4: Two free ends:
Displacement of masses = [ 2.21109691e-16 -7.84937407e-16  9.82722686e-16]
Elongation of springs = [-1.00604710e-15  1.76766009e-15]
Internal Stresses = [-1.00604710e-15  1.76766009e-15]
Condition Number = 187965912.2651711
```
From these cases we can see some interesting results. First of all, the result for case 1 matches the result from homework 2 problem 6. Secondly, the displacement result from case 2 is the opposite of the result from case 3, which is expected because these are the same systems algebraically except flipped. This doesn't make sense physically for case 3, but I suppose this means that the spring stiffness in case 3 is just not stiff enough. Additionally, in the case where there are two free ends, the displacement, elongation, and internal stresses are all close to 0. That is because the system is free falling. You can also see that the condition number is quite high for this matrix.
