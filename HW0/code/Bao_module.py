# Bao's module contains Bao's functions

import numpy as np


# computing n!
def factorial(n):
    if type(n)==int:
        if n==0:
            return 1
        else:
            fact=1
            for i in np.arange(1,n+1):
                fact=fact*i
            return fact
    else:
        print('Wrong type of argument!!')

# check if an integer is prime
def isPrime(n):
    if type(n)==int:
        s=0
        for i in np.arange(2,n):   
            if n%i==0:
                s=s+1
            else:
                s=s
        if s!=0:
            print('not prime')
        else:
            print('prime')
    else:
        print('Wrong type of argument!!')