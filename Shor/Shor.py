#####################################################
#
# Simulate the execution of Shor's algorithm 
#
# Copyright (c) 2018 christianb93
# Permission is hereby granted, free of charge, to 
# any person obtaining a copy of this software and 
# associated documentation files (the "Software"), 
# to deal in the Software without restriction, 
# including without limitation the rights to use, 
# copy, modify, merge, publish, distribute, 
# sublicense, and/or sell copies of the Software, 
# and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice 
# shall be included in all copies or substantial 
# portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY 
# OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# If you want to run this on a machine without X-server,
# do a
# export MPLBACKEND="AGG"
#####################################################


import numpy as np
import matplotlib.pyplot as plt
import tempfile
import argparse



########################################################
# Given a rational number p / q, this function will
# calculate the (finite) continuous fraction expansion
########################################################


def cf(p,q):
    X = []
    while (q != 0):
        a = p // q
        r = p % q
        X.append(a)
        p = q
        q = r
    return X
    
    
    
########################################################
# This function directly returns the first convergent of
# the continous fraction expansion of a ratio p / q
# which has a numerator exceeding b or the m-th convergent,
# whichever condition occurs first
########################################################

def convergent(p, q, b = None, m = None):
    X = []
    n = 0
    p_n_minus2 = 0
    p_n_minus1 = 1
    q_n_minus1 = 0
    q_n_minus2 = 1
    while (q != 0):
        #
        # Compute a_n
        #
        a_n = p // q
        r = p % q
        p = q
        q = r
        n = n + 1
        #
        # Compute p_n and q_n
        #
        p_n = a_n * p_n_minus1 + p_n_minus2
        q_n = a_n * q_n_minus1 + q_n_minus2
        #
        # If we have reached the precision or
        # the desired number of steps return result
        #
        if (b != None):
            if q_n > b:
                return p_n_minus1, q_n_minus1
        if m != None: 
            if n > m:
                return p_n, q_n
        #
        # Shift variables
        #
        q_n_minus2 = q_n_minus1
        q_n_minus1 = q_n
        p_n_minus2 = p_n_minus1
        p_n_minus1 = p_n
    return p_n, q_n
    
    
########################################################
# Use the Euclidian algorithm to calculate the greatest
# common divisor of two numbers a and b
########################################################

def gcd(a,b):
    while b != 0:
        t = b
        b = a % b
        a = t
    return a

##############################################
# Use exponentation by squaring to calculate
# a power x**k mod M
##############################################

def power(x,k,M):
    if k == 0:
        return 1
    p = 1
    while k > 0:
        if 0 == (k % 2):
            x = x**2 % M
        else:
            p = p * x % M
            x = x**2 % M
        k = k // 2
    return p


########################################################
# Prepare the initial state for Shor's algorithm
# The state is a linear combination of states
# |k> |x mod M>
# which we can therefore describe as a matrix with
# N rows and M columns. 
########################################################

def init(N,M, x):
    state = np.zeros((N,M),dtype=np.complex)
    for k in range(N):
        y = power(x,k,M)
        state[k,y] = np.complex(1, 0)
    return state / np.sqrt(N)
        

        

########################################################
# Given a state
# S = \sum_a_{st} |s> |t mod M>
# represented by a matrix S with N rows and M columns,
# this function will calculate the quantum Fourier 
# transform of this state, acting on the first
# register
########################################################
        
def QF(S):
    N = S.shape[0]
    #
    # The t-th column of our new state is the
    # discrete Fourier transform of the t-th 
    # column of the old state. As the numpy.fft
    # module uses slightly different conventions, we have
    # to use the numpy.fft.iffy routine and multiply by
    # the square root of N
    #
    return np.sqrt(N) * np.fft.ifft(S, axis = 0)
        

##############################################
# Draw from a finite distribution defined
# by a vector p with elements adding up
# to one. We return a number between
# 0 and the number of elements of p minus 1
# and i is returned with probability p[i]
##############################################
def draw(p):
    u = np.random.uniform()
    x = 0
    n = 0
    for i in p:
        x = x + i
        if x >= u:
            return n
        n += 1
    return i-i        
        

##############################################
# Use Shor's algorithm to find the period of
# a number x mod M. We assume that x and M
# are coprime
##############################################
def findPeriod(x,M):
    #
    # We need to find the number of bits that M**2 has
    #
    r = 0
    n = int(np.log2(M**2)) + 1
    N = 2**n
    print("Using ",n,"qubits")
    success = 0
    while (0 == success):
        #
        # Initialize the state
        #
        print("Setting up initial state")
        state = init(N,M,x)
        #
        # Apply QFT
        #
        print("Running quantum Fourier transform")
        state = QF(state)
        #
        # Now we compute, for each value of s between 0 and N, the probability
        # to measure s
        #
        P = np.zeros(N)
        for s in range(N):
            for i in range(M):
                P[s] = P[s] + abs(state[s,i])**2
        #
        # Now we measure, i.e. we draw a value from the probability distribution P. 
        #
        s = draw(P)
        print("Measured s = ", s)
        #
        # For debugging, print the convergents
        #
        X = cf(s, N)
        for i in range(len(X)):
            print("Convergent ", i, ":", convergent(s, N, b = None, m = i))
        #
        # Get the approximation to r by continued fractions
        #
        c,r = convergent(s, N, b = M, m = None)
        print("Guess for period: r = ", r, ", c = ", c)
        if (1 == power(x,r,M)):
            print("Quantum part successful")
            success = 1
        else:
            print("Algorithm failed, repeating")
    
        
    return r, P
    
    
####################################################
# Parse arguments
####################################################
        
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--M", 
                    type=int,
                    default=21,
                    help="Number that we want to factor")
    args=parser.parse_args()
    return args
    
####################################################
#
# Main
#
####################################################

#
# Get the argument - this is the number M that we want to factor
#
M = get_args().M
print("Trying to factor M = ",M)
#
# Next we need to choose a number x which is coprime to M
# We pick x randomly and start over if x is not coprime to M
# (in a real-world application, we could stop if this happens
# as we then have found a factor)
#
x = M
while (gcd(x,M) != 1):  
    x = np.random.randint(low = 2, high = M - 1)
print("Using x = ", x)    

#
# Run quantum part
#
r, P = findPeriod(x, M)

if (1 == (r % 2)):
    print("Period r is not even, please run script again")
    exit

#
# Now use this to find a divisor. We calculate x**(r/2)  - 1
# 
a = (power(x, r // 2, M) - 1) % M
print("Calculating gcd(",a, ",",M,")")
d = gcd(a % M,M)
if (0 != (M % d)):
    print("Ups, ",d,"should be a factor but is not - something went wrong")
else:
    print("Found factor d = ", d)


#
# If we get to this point, the algorithm was sucessful. Plot P
#
fig = plt.figure(figsize=(12,5))
axis = fig.add_subplot(1,1,1)
axis.plot(P)
outfile = tempfile.gettempdir() + "/ShorSampleOutput.png"
print("Saving plot of probability distribution to ", outfile)
fig.savefig(outfile)
plt.show()
