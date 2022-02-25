This is the first RPI trial, it contains the following:
1) An UBPF algorithm based on newton-raphson that reads the Ybus,P load, Q load and busnames from text files that were extracted offline from openDSS.
2) A python script that calls openDSS to solve the same problem, solved by the newton raphson algorithm.
3) Both python scripts have the capability to view the amount of memory used when solving the problem.
