#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 11:17:29 2023

@author: khd2
"""
import numpy as np
import os
import sys

# def print_(i):
#     if i%20 == 0:
#         print(i)
import time

import gc
gc.collect()

regionname = "7115"#sys.argv[1]
startfl=0#int(sys.argv[2])
endfl=1#int(sys.argv[3])
ly=7#int(sys.argv[4])
lx=3#int(sys.argv[5])

# regionname = sys.argv[1]
# startfl=int(sys.argv[2])
# endfl=int(sys.argv[3])
# ly=int(sys.argv[4])
# lx=int(sys.argv[5])
sample_x, sample_y = 1,1  # x and y steps for fast calculations
# print(lx,ly)

a = 1 

#ly=4#int(sys.argv[4])
#lx=6#int(sys.argv[5])

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
#                    Byp[i,j] += (   (y[j] - y[j1]) ) #/ ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#                    Bxp[i,j] +=  (B[i1,j1]*(x[i] - x[i1])) /((x[i] - x[i1])**2 + (y[j] - y[j1])**2 +1 )**(3/2)
                    Byp[i,j] += (  B[i1,j1]* (y[j] - y[j1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)


                    yi1 = y[j1]
                    
                xi1 = x[i1]

    return Bxp, Byp

### rotate the matrix for grid
#rotated = list(zip(*rotated))[::-1]
#rotated = list(zip(*rotated))[::-1]
#z = np.linspace(0,0,ly).reshape(1,-1).tolist()

# def GP_new(x,y,lx,ly):

#     Bxp=np.zeros((lx,ly), dtype=float)
    
#     if lx%2 != 0:        ## if lx is not even number, put 0 in the middel
#         for j in range(ly):
#             Bxp[(lx//2),j] = 0
#     for i in range(lx//2):
#         for j in range(ly):
#             for i1 in range(lx):
#                 for j1 in range(ly):
#                     Bxp[i,j] += ( (x[i] - x[i1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#                     Bxp[-(i+1),-(j+1)] += ( -1*(x[i] - x[i1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#     return Bxp

# def cal_r(list_,indx,length):
#     r = []
#     for i in range(length):
#         r.append((list_[indx]-list_[i]) )
#     return np.array(r)

# def GP_newV(M,x,y,lx,ly):

#     B = np.zeros((lx,ly), dtype=float)
#     B = B +np.reshape(M,(lx,ly))
    
#     Bxp=np.zeros((lx,ly), dtype=float)
#     Byp=np.zeros((lx,ly), dtype=float)
    
#     for i in range(lx//2):
#         for j in range(ly//2):
#             for i1 in range(lx):
#                 for j1 in range(ly):
# #                    Bxp[i,j] += ( B[i1,j1]* (x[i] - x[i1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
# #                    Bxp[-(i+1),-(j+1)] += ( B[-(i1+1),-(j1+1)]*-1*(x[i] - x[i1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)

# #                    if ((i+1)==lx//2):
# #                        Bxp[i+1,j] +=  (B[i1,j1]*(x[i+1] - x[i1]) ) / ((x[i+1] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)

# #                    Byp[i,j] += (B[i1,j1]* (y[j] - y[j1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
# #                    Byp[-(i+1),-(j+1)] += ( B[-(i1+1),-(j1+1)]*-1*(y[j] - y[j1]) ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
                    
# #                   ''' ##  fill left of the matrix '''
#                     X = x[i] - x[i1]
#                     X_1 = x[i+1] - x[i1]
                    
#                     Bxp[i,j] += ( B[i1,j1]* (X) ) / (X**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#                     Bxp[-(i+1),-(j+1)] += ( B[-(i1+1),-(j1+1)]*-1*(X) ) / (X**2 + (y[j] - y[j1])**2 + a**2)**(3/2)

#                     if lx%2 == 0:
# #                   ''' ##  fill rigth of the matrix '''
#                         if ((i+1) != lx//2):       ## i is not at the middel of lx
#                             Bxp[i,-(j+1)] += ( B[i1,-(j1+1)]* X ) / (X**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#                             Bxp[-(i+1),j] += ( B[-(i1+1),j1]*-1*X ) / (X**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#                         else:
#                             Bxp[i,-(j+1)] += ( B[i1,-(j1+1)]* X ) / (X**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#                     else:
#                         Bxp[i,-(j+1)] += ( B[i1,-(j1+1)]* X ) / ((x[i] - x[i1])**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#                         Bxp[-(i+1),j] += ( B[-(i1+1),j1]*-1*X ) / (X**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
                     
# ########
# #                   ''' ##  fill the meddile of the matrix if lx is odd number'''
#                     if ((i+1)==lx//2):
# #                        X_1 = x[i+1] - x[i1]
#                         Bxp[i+1,j] +=  (B[i1,j1]*X_1 ) / (X_1**2 + (y[j] - y[j1])**2 + a**2)**(3/2)
#                         Bxp[i+1,-(j+1)] +=  (B[i1,-(j1+1)]*X_1 ) / (X_1**2 + (y[j] - y[j1])**2 + a**2)**(3/2)


# #                    ''' ##  fill the meddile of y in the matrix if ly is odd number'''
#                     if ly%2 !=0:
#                         if ((j+1) == ly//2):
#                             Bxp[i,j+1] +=  (B[i1,j1]*X ) / (X**2 + (y[j+1] - y[j1])**2 + a**2)**(3/2)
#                             Bxp[-(i+1),j+1] +=  (B[-(i1+1),j1]*-1*X ) / (X**2 + (y[j+1] - y[j1])**2 + a**2)**(3/2)
                    
#                     if ((ly%2 !=0) and (lx%2 !=0) and ((i+1)==lx//2)):
#                             Bxp[lx//2,ly//2] +=  (B[i1,j1]*(x[lx//2] - x[i1]) ) / ((x[lx//2] - x[i1])**2 + (y[ly//2] - y[j1])**2 + a**2)**(3/2)
# ########
#    return Bxp, Byp

def GP_newV(M,x,y,lx,ly):

    B = np.zeros((lx,ly), dtype=float)
    B = B +np.reshape(M,(lx,ly))
    
    Bxp=np.zeros((lx,ly), dtype=float)
    Byp=np.zeros((lx,ly), dtype=float)
    
    for i in range(lx//2):
        for j in range(ly//2):
            for i1 in range(lx):
                X = x[i] - x[i1]
                X_1 = x[i+1] - x[i1]
                for j1 in range(ly):
                    
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
#original                            Bxp[lx//2,ly//2] +=  (B[i1,j1]*(x[lx//2] - x[i1]) ) / ((x[lx//2] - x[i1])**2 + (y[ly//2] - y[j1])**2 + a**2)**(3/2)
#test
                            Bxp[lx//2,ly//2] +=  (B[i1,j1]*X_1 ) / (X_1**2 + Y_1**2 + a**2)**(3/2)

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
#original                        Byp[lx//2,ly//2] +=  (B[i1,j1]* (y[ly//2] - y[j1]) ) / ((x[lx//2] - x[i1])**2 + (y[ly//2] - y[j1])**2 + a**2)**(3/2)
#test
                        Byp[lx//2,ly//2] +=  (B[i1,j1]* Y_1 ) / (X_1**2 + Y_1**2 + a**2)**(3/2)


    return Bxp, Byp



M = np.arange(0,63765,1)
lx = 585
ly = 109#int(M.size/lx)

x = np.linspace(0,1,lx)
y = np.linspace(0,1,ly)
M = np.array([5,4,6,3,2,5,
            5,4,6,3,7,5,
            5,8,6,10,2,5,
            5,4,6])

import gc
gc.collect()
#X = grid(M,x,y, lx,ly)
#Bxp_new = GP_new(x,y,lx,ly)

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

# s1=sys.argv[6]+"/bz_"+regionname+"_"
# #p = "/media/khd2/Spare-Data-Disk/ARTop_winding/testNew_output/Data/AR_7115/output"
# #s1=p+"/bz_"+regionname+"_"

# ##==========================
# #s1=sys.argv[6]+"/bz_"+regionname+"_"
# #print(s1)
# s3=".txt"
# L=lx*ly
# ##==========================

# for i in range(startfl,endfl):
#     print("pot field ",i,"\n")
#     s2=str(i)#, base = 10, pad = 1)
#     path=s1+s2+s3
# #    print(path)
#     if os.path.isfile(path)==True:
#         pol=np.loadtxt(path, delimiter=" ", dtype=np.float64)
#         M=np.array(pol)
# #        print(M.shape, 'done')
        
#         Bxp, Byp = G_fn_potential(M, lx,ly)

#         Bxpf=np.reshape(Bxp,(L,1))
#         Bypf=np.reshape(Byp,(L,1))


        
#         Nx="Bxp_"+regionname+"_"
#         Ny="Byp_"+regionname+"_"

#         z1=sys.argv[6]+"/"
# #        print(z1)
# #        z1=p+"/";

#         z3="Gf.txt"          
        
#         a_file = (z1+Nx+s2+z3)
#         np.savetxt(a_file, Bxpf[:])

#         a_file = (z1+Ny+s2+z3)
#         np.savetxt(a_file, Bypf[:])        

