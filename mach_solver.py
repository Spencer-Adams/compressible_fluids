import json 
import numpy as np 
import math
import matplotlib.pyplot as plt 

class nozzle_solver:
    """This class includes all of the the functions associated with finding and plotting things about a nozzle"""
    def __init__(self, gamma, expansion_ratio):
        self.gamma = gamma 
        self.expansion_ratio = expansion_ratio


    def calc_f_m_j(self, m_j):
        """This function returns a value for F(M_j)"""
        m_j_sq = m_j**2
        # first, second, third, and exponent are just portions of the equation for f(M_j). 
        # Splitting it up this way allows for a code that's easier to write and read.
        first = 1/(m_j)
        second = 2/(self.gamma + 1)
        third = 1+(((self.gamma - 1)/2)*m_j_sq)
        exponent = (self.gamma + 1) / (2*(self.gamma - 1))
        F_M_j = (first)*(second * third)**(exponent) - self.expansion_ratio
        return F_M_j


    def calc_del_f_del_M(self, m_j):
        """This function returns a value for the derivative of F with respect to M, del_F/del_M"""
        m_j_sq = m_j**2
        # first, second, third, and exponent are just portions of the equation for del_F/del_M. 
        # Splitting it up this way allows for a code that's easier to write and read.
        first = 2**((1-3*self.gamma)/(2-2*self.gamma))
        second = (m_j_sq - 1)/((m_j_sq)*(2+m_j_sq*(self.gamma-1)))
        third = (1+((self.gamma - 1)/2)*m_j_sq)/(self.gamma+1)
        exponent = (self.gamma+1)/(2*(self.gamma-1))
        del_F_DEL_M = first*second*(third**exponent)
        return del_F_DEL_M
    
    
    def calc_epsilon(self, m_j):
        """This function retrieves the value that can trigger the stopping condition in newton solver"""
        m_j_sq = m_j**2
        # first, second, third, and exponent are just portions of the equation for del_F/del_M. 
        # Splitting it up this way allows for a code that's easier to write and read.
        first = 1/m_j
        second = 2/(self.gamma + 1)
        third = 1 + ((self.gamma - 1)/2)*m_j_sq
        exponent = (self.gamma + 1)/(2*(self.gamma-1))
        epsilon = abs(((first * (second * third)**exponent)) - self.expansion_ratio)/self.expansion_ratio
        return epsilon


    def newton_solver(self, start_mach, end_epsilon):
        """This function numerically solves for the mach number given the"""
        m_j = start_mach
        epsilon = 10 #arbitrary placeholder
        while epsilon >= end_epsilon:
            print("mach_number: ",m_j)
            f_m_j = self.calc_f_m_j(m_j)
            del_f_del_M = self.calc_del_f_del_M(m_j)
            M_j_plus_one = m_j - (f_m_j/del_f_del_M)
            epsilon = self.calc_epsilon(m_j)
            m_j = M_j_plus_one

        
        return m_j


if __name__ == "__main__":
    """This is where the code actually runs"""
    gamma = 1.25 # input gamma
    expansion_ratio = 77.5 # input expansion ratio 
    start_mach = 3.5 # if anticipating Mach<1, initialize slightly greater than 0, else initialize at greater than 3 
    end_epsilon = 0.001 # stopping criteria 
    nozzle = nozzle_solver(gamma, expansion_ratio) # this line initializes the nozzle class 

    mach_number = nozzle.newton_solver(start_mach, end_epsilon) # this line runs the newton solver for finding the mach number. 
    # print(mach_number)

