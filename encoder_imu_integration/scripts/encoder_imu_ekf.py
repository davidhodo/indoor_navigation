#!/usr/bin/env python

"""
Provides an EncoderImuEkf class that implements an
Extended Kalman Filter (EKF) for blending measurements from
a set of wheel encoders and a 

"""

import numpy

class EncoderImuEkf:
    def __init__(self):
        self.num_states=5;
        self.x=numpy.zeros(self.num_states)
        self.P=numpy.zeros(self.num_states,self.nnum_states)
        
        pass
    
    def TimeUpdate(self,t,u):
        pass
    
    def MeasurementUpdate(self,t,y):
        pass
    
    