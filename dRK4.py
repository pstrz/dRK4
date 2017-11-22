# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 13:11:43 2017

@author: philstrzelecki
"""

import numpy as np

def RK4(x, t, h):
    f = L96
    for i in range(1, t.size+1):
        k1 = f(x[:,i-1])
        k2 = f(x[:,i-1] + 1/2*h*k1)
        k3 = f(x[:,i-1] + 1/2*h*k2)
        k4 = f(x[:,i-1] + h*k3)
        y = x[:,i-1] + 1/3*h*(1/2*k1 + k2 + k3 + 1/2*k4)
        y = np.reshape(y,(10,1))
        x = np.append(x, y, axis=1)
    return x

def dRK4(x, t, h):
    f = L96
    df = dL96
    dx = np.eye(x.size)
    
    x = RK4(x, t, h)
    for i in range(1,t.size+1):
        k1 = f(x[:,i-1])
        dk1 = df(x[:,i-1]) @ dx
        k2 = f(x[:,i-1] + 1/2*h*k1)
        dk2 = df(x[:,i-1] + 1/2*h*k1) @ (dx + 1/2*h*dk1)
        k3 = f(x[:,i-1] + 1/2*h*k2)
        dk3 = df(x[:, i-1] + 1/2*h*k2) @ (dx + 1/2*h*dk2)
        k4 = f(x[:,i-1] + h*k3)
        dk4 = df(x[:,i-1] + h*k3) @ (dx + h*dk3)
        dx = dx + 1/3*h*(1/2*dk1 + dk2 + dk3 + 1/2*dk4)
    return dx
    
def L96(x):
    F = 8
    G = F * np.ones(x.size)
    
    f = (np.roll(x,(-1,0)) - np.roll(x,(2,0))) * np.roll(x,(1,0)) - x + G
    return f
    
def dL96(x):
    X = np.diag(x)
    df = np.roll(X,-1, axis=0) - np.roll(X,(2,3), axis=(0,1)) + np.roll(X, (2,1), axis=(0,1)) - np.roll(X, (-1,1), axis=(0,1)) - np.eye(x.size)
    return df
    
h = 0.01
t = np.arange(0, 0.2, 0.01)

x = np.array([[1], [2], [3], [4], [5], [6], [7], [8], [9], [10]])

dM = dRK4(x,t,h)
    