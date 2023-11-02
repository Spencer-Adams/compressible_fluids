import json 
import numpy as np 
import math
import matplotlib.pyplot as plt 

class ramjet_solver:
    """This class contains important nozzle attributes that can assist in nozzle design"""
    def __init__(self, json_file):
        self.json_file = json_file 

    # First, here are some functions that allow us to take the information 
    # from the text file so that we can plot it. 
    def load_project_one_json(self):
        """This function loads the json file that references the text file"""
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.gamma = input_vals["gamma"]
            self.P_inf = input_vals["P_inf"]*10**3  # in Pa
            self.T_inf = input_vals["T_inf"] # in K
            self.m_weight = input_vals["m_weight"] # in g/mol=kg/kmol
            self.R_g = 8314.41/self.m_weight 
            self.expansion_ratio_diffuser = input_vals["diffuser_exp_ratio"] # this value is calculated in the Mathcad file
            self.expansion_ratio_nozzle = input_vals["nozzle_exp_ratio"] # this value is calculated in the Mathcad file
            self.end_epsilon = input_vals["end_epsilon"]


    # Now, here are the functions associated with the newton solver (HW 4, problem 1)
    def calc_f_m_j(self, m_j, expansion_ratio):
        """This function returns a value for F(M_j)"""
        m_j_sq = m_j**2
        # first, second, third, and exponent are just portions of the equation for f(M_j). 
        # Splitting it up this way allows for a code that's easier to write and read.
        first = 1/(m_j)
        second = 2/(self.gamma + 1)
        third = 1+(((self.gamma - 1)/2)*m_j_sq)
        exponent = (self.gamma + 1) / (2*(self.gamma - 1))
        F_M_j = (first)*(second * third)**(exponent) - expansion_ratio
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


    def calc_epsilon(self, m_j, expansion_ratio): #### could be a non class-specific function. 
        """This function retrieves the value that can trigger the stopping condition in newton solver"""
        epsilon = abs(self.calc_f_m_j(m_j, expansion_ratio))/expansion_ratio
        return epsilon
    

    def isentropic_newton_solver(self, start_mach, expansion_ratio):
        """This function numerically solves for the mach number given the"""
        m_j = start_mach
        epsilon = 10 #arbitrary placeholder
        while epsilon >= self.end_epsilon:
            # print("mach_number: ",m_j)
            f_m_j = self.calc_f_m_j(m_j, expansion_ratio)
            del_f_del_M = self.calc_del_f_del_M(m_j)
            M_j_plus_one = m_j - (f_m_j/del_f_del_M)
            epsilon = self.calc_epsilon(m_j, expansion_ratio)
            m_j = M_j_plus_one
        self.m_e = m_j
        return m_j
    

    def run_isentropic_solver(self):
        """This function runs all of the necessary functions to run the isentropic solver"""
        print("Ramjet Stuff:")
        self.load_project_one_json()
        m_3 = self.isentropic_newton_solver(0.1, self.expansion_ratio_diffuser) ## 0.5 is a subsonic starting guess because we're on the right side of diffuser throat
        print("m_3: ", m_3)
        m_5 = self.isentropic_newton_solver(3, self.expansion_ratio_nozzle) ## 3 is a supersonic starting guess because we're on the right side of the nozzle throat
        print("m_5", m_5)
   

if __name__ == "__main__":
    """This is where the code actually runs"""
    nozzle_isentropic = ramjet_solver("ramjet.json")
    run_isentropic = nozzle_isentropic.run_isentropic_solver()
    

