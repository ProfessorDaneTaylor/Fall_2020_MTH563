# Bao's module contains Bao's functions

import numpy as np


# computing n!
def factorial(n):
    """
    Bao's first function.

    Args:
        n (int): Some postive integer.

    Returns:
        int: Returns :math:`n!`
    """
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
    """
    Bao's second function.

    Args:
        n (int): Some integer.

    Returns:
        string: Returns "prime" if n is prime, "not prime" if n is not prime
    """
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