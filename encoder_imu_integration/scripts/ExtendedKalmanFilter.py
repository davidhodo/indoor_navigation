# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 18:03:52 2012

@author: hododav
"""

from numpy import matrix, vstack, hstack
# import matlib which gives matrix versions of eye, zeros, etc.
from numpy.matlib import eye, zeros, arange
# import linalg for matrix exponential function (expm)
from scipy.linalg import expm

# Helper functions

def c2dm(A, Bu, Bw, Q, dt):
    """KFC2DM Converts a continuous system with process noise to discrete 
    
       [Ad,Bd,Qd]=KFc2dm(A,Bu,Bw,Q,dt)
           Params:
               - A = continuous time A matrix
               - Bu - continuous time input matrix
               - Bw - continuous time process noise input matrix
               - Q - continuous time process noise covariance matrix
               - ts = sample time
           Outputs:
               - Ad = state transition matrix
               - Bd = state estimate prediction after time update
    
       Assuming a continuous time system of the form:
           dx/dt=A x(t) + B u(t) + Bw w(t)
       where w is a Guassian, white noise with statistics ~ N(0,Q)
    
       The system is discretized to give:
           x[k+1]=Ad x[k] + Bd u[k] + w[k]
       where w ~ N(0,Qd) 
    """

    #example    
    #A=eye(3)
    #Bu=arange(6).reshape([3,2])
    #Bw=Bu
    #Q=1.2*eye(2)
    #dt=0.05
    
    n=A.shape[0]
    
    # use van Loan method
    # form matrix:
    #     AA=[-A Bw*Q*Bw';
    #         zeros(n)   A'];
    AA=vstack((hstack((-A,Bw*Q*Bw.T)),hstack((zeros((n,n)),A))))

    # calculate matrix exponential
    BB=expm(AA*dt);
    # extract state transition matrix and discrete process covariance
    Ad=BB[n:2*n][:,n:2*n].T
    Qd=Ad*BB[0:n][:,n:2*n]
    
    return (Ad, Qd)

class ExtendedKalmanFilter:
    def __init__(self, time_step, num_states, num_measurements, Qk, Rk,
                 evolution_function, evolution_jacobian, 
                 observation_function, observation_jacobian):
        self._num_states=num_states;  # number of filter states
        self._num_measurements=num_measurements
        self._dt=time_step
        self._Xk=zeros(self._num_states) # current state estimate (+)
        self._Xkm=zeros(self._num_states) # current state estimate (-)
        self._Pk=zeros(self._num_states,self._num_states) # current estimate error covariance (+)
        self._Pkm=zeros(self._num_states,self._num_states) # current estimate error covariance (-)
        self._Qk=Qk  # process noise covariance matrix
        self._Rk=Rk  # measurement noise covariance matrix
        
        self.evolution_function = evolution_function
        self.evolution_jacobian = evolution_jacobian
        self.observation_function = observation_function
        self.observation_jacobian = observation_jacobian        
            
    def initialize(self, X_0, P_0 ):
        self._Xk=X_0
        self._Pk=P_0
    
    def prediction_step(self,Uk):
        self._Uk=Uk
        self._Xkm = self.evolution_function(self._dt, self._Xk, self._Uk)
        self._Phik = self.evolution_jacobian(self._dt, self._Xk, self._Uk)
        self._Pkm = self._Phik * self._Pk * self._Phik.T + self._Qk
        #print "After predictio nstep : ", self.xkm

    def correction_step(self,Yk):
        self._Yk=Yk
        self._Hk = self.observation_jacobian(self._dt, self._Xkm, self._Ukm)
        self._Kk = self._Pkm * self._Hk.T * (self._Hk * self._Pkm * self._Hk.T+self._Rk).I
        self._Pk = (matrix(identity(self.num_states)) - self.Kk * self.Hk) * self.Pkm
        self._Xk = self._Xkm + self._Kk * (self.yk - self.observation_function(self.dt,self.xkm))

    def linear_propagation(self,x,u,params):
        self._Phik = self.evolution_jacobian(self._dt, self._Xk, self._Uk)
        return x + self._Phik*x
    
    def linear_measurement(self,x,u):
        self._Hk = self.observation_jacobian(self._dt, self._Xkm)
