#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 11:17:29 2023

@author: khd2
"""
import numpy as np
import os
import sys

import time

import gc
gc.collect()

regionname = "7115"#sys.argv[1]
startfl=0#int(sys.argv[2])
endfl=1#int(sys.argv[3])
ly=7#int(sys.argv[4])
lx=3#int(sys.argv[5])



a = 1 
x = np.linspace(0,1,lx)
y = np.linspace(0,1,ly)



def G_fn_potential(M, lx,ly):
    
    B = np.zeros((lx,ly), dtype=float)
    B = B +np.reshape(M,(lx,ly))

    Bxp=np.zeros((lx,ly), dtype=float)
    Byp=np.zeros((lx,ly), dtype=float)
    
    for i in range(lx):
#        print(i)
#        print_(i)
        for j in range(ly):            
            xi1 = x[0]
            c=0
            for i1 in range(lx):
                yi1 = y[0]
                for j1 in range(ly):
                    
#                    dx = abs(x[i1]) - abs(xi1)   # x seperation
#                    dy = abs(y[j1]) - abs(yi1)   # y seperation

                    Bxp[i,j] += ( B[i1,j1]*  (x[i] - x[i1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
                    Byp[i,j] += (  B[i1,j1]* (y[j] - y[j1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)


                    yi1 = y[j1]
                    
                xi1 = x[i1]

    return Bxp, Byp



def GP_newV(M,x,y,lx,ly):

    B = np.zeros((lx,ly), dtype=float)
    B = B +np.reshape(M,(lx,ly))
    
    Bxp=np.zeros((lx,ly), dtype=float)
    Byp=np.zeros((lx,ly), dtype=float)
    
    for i in range(lx//2):
        for j in range(ly//2):
            for i1 in range(lx):
                for j1 in range(ly):
                    X = x[i] - x[i1]
                    X_1 = x[i+1] - x[i1]
                    
                    Y = y[j] - y[j1]
                    Y_1 = y[j+1] - y[j1]
                    
                    r = (X**2 + Y**2 + a**2)**(3/2)
                    r_X_1 = (X_1**2 + Y**2 + a**2)**(3/2)
                    r_Y_1 = (X**2 + Y_1**2 + a**2)**(3/2)

# fill top and bottom of x-axis (rows) in the matric                    
                    Bxp[i,j] += ( B[i1,j1] * X ) / r
                    Bxp[-(i+1),-(j+1)] += ( B[-(i1+1),-(j1+1)]*-1*(X) ) / r

# 
                    if lx%2 == 0:             # if lx is even number
#                   ''' ##  fill rows from right and left sids of the matrix '''
                        if ((i+1) != lx//2):         ## i is not at the middle of lx
                            Bxp[i,-(j+1)] += ( B[i1,-(j1+1)]* X ) / r      # fill rows from right
                            Bxp[-(i+1),j] += ( B[-(i1+1),j1]*-1*X ) / r    # fill rows from left
                        else:
                            Bxp[i,-(j+1)] += ( B[i1,-(j1+1)]* X ) / r
                            
                            Bxp[i+1,j] +=  (B[i1,j1]*X_1 ) / r_X_1
                    else:                     # if lx is odd number
                        Bxp[i,-(j+1)] += ( B[i1,-(j1+1)]* X ) / r
                        Bxp[-(i+1),j] += ( B[-(i1+1),j1]*-1*X ) / r
                        
 #                       Bxp[i+1,j] +=  (B[i1,j1]*X_1 ) / r_X_1
 #                       Bxp[i+1,-(j+1)] +=  (B[i1,-(j1+1)]*X_1 ) / r_X_1
                        
# fill the middle of x-axis in the matrix if lx is odd number'''
                    if (((i+1)==lx//2) and (lx%2 != 0) ):
#                        X_1 = x[i+1] - x[i1]
                        Bxp[i+1,j] +=  (B[i1,j1]*X_1 ) / r_X_1
                        Bxp[i+1,-(j+1)] +=  (B[i1,-(j1+1)]*X_1 ) / r_X_1


# fill the middle of y-axis in the matrix if ly is odd number'''
                    if ly%2 !=0:
                        if ((j+1) == ly//2):
                            Bxp[i,j+1] +=  (B[i1,j1]*X ) / r_Y_1
                            Bxp[-(i+1),j+1] +=  (B[-(i1+1),j1]*-1*X ) / r_Y_1

## fill the center of the matrix if ly and lx are odd numbers                    
                    if ((ly%2 !=0) and (lx%2 !=0) and ((i+1)==lx//2)  and ((j+1)==ly//2)):
                            Bxp[lx//2,ly//2] +=  (B[i1,j1]*(x[lx//2] - x[i1]) ) / ((x[lx//2] - x[i1])**2 + (y[ly//2] - y[j1])**2 + a**2)**(3/2)

##################### Byp term ################################
                    Byp[i,j] += ( B[i1,j1] * Y ) / r
                    Byp[-(i+1),-(j+1)] += ( B[-(i1+1),-(j1+1)]*-1*Y ) / r

                    if lx%2 == 0:             # if lx is even number
#                   ''' ##  fill rows from right and left sids of the matrix '''
                        if ((i+1) != lx//2):         ## i is not at the middle of lx
                            Byp[i,-(j+1)] += ( B[i1,-(j1+1)]* -Y ) / r      # fill rows from right
                            Byp[-(i+1),j] += ( B[-(i1+1),j1]* Y ) / r    # fill rows from left
                        else:
                            Byp[i,-(j+1)] += ( B[i1,-(j1+1)]* -Y ) / r
                            Byp[-(i+1),j] += ( B[-(i1+1),j1]* Y ) / r

                    else:                     # if lx is odd number
                        Byp[i,-(j+1)] += ( B[i1,-(j1+1)]* -Y ) / r
                        Byp[-(i+1),j] += ( B[-(i1+1),j1]* Y ) / r

# fill the middle of x-axis in the matrix if lx is odd number'''
                    if (((i+1)==lx//2) and (lx%2 != 0)):
#                        X_1 = x[i+1] - x[i1]
                        Byp[i+1,j] +=  (B[i1,j1]*Y ) / r_X_1
                        Byp[i+1,-(j+1)] +=  (B[i1,-(j1+1)]* -Y ) / r_X_1

# fill the middle of y-axis in the matrix if ly is odd number'''
                    if ly%2 !=0:
                        if ((j+1) == ly//2):
                            Byp[i,j+1] +=  (B[i1,j1]*Y_1 ) / r_Y_1
                            Byp[-(i+1),j+1] +=  (B[-(i1+1),j1]*Y_1 ) / r_Y_1
                            
#    ## fill the center of the matrix if ly and lx are odd numbers                    
                    if ((ly%2 !=0) and (lx%2 !=0) and ((i+1)==lx//2)  and ((j+1)==ly//2)):
                        Byp[lx//2,ly//2] +=  (B[i1,j1]* (y[ly//2] - y[j1]) ) / ((x[lx//2] - x[i1])**2 + (y[ly//2] - y[j1])**2 + a**2)**(3/2)


    return Bxp, Byp



M = np.arange(0,63765,1)
lx = 585
ly = 109

x = np.linspace(0,1,lx)
y = np.linspace(0,1,ly)

import gc
gc.collect()

import time

start_time = time.time()

Bxp_new, Byp_new = GP_newV(M,x,y,lx,ly)

end = time.time()
print('Duration1:', ((end - start_time))/60, "minutes")

#end_time = datetime.now()
#print('Duration1: {}'.format(end_time - start_time))

#end = time.time()
#print('time 1 = ',((end - start)*688**3)/3600)


start_time = time.time()

Bxp, Byp = G_fn_potential(M, lx,ly)
end = time.time()
print('Duration2:', ((end - start_time))/60, "minutes")
