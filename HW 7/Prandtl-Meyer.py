import json 
import numpy as np 
import math
import matplotlib.pyplot as plt 

class Prandtl_Meyer:
    """This class contains important nozzle attributes that can assist in nozzle design"""
    def __init__(self, json_file):
        self.json_file = json_file
        self.load_json() 

    # First, here are some functions that allow us to take the information 
    # from the text file so that we can plot it. 
    def load_json(self):
        """This function loads the json file that references the text file"""
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.M_1 = input_vals["M_1"]
            self.gamma = input_vals["gamma"]
            self.one_ATM = input_vals["one_ATM"] 
            self.P_1 = input_vals["P_1"]*self.one_ATM # in Pa
            self.T_1 = input_vals["T_1"] # in K
            self.theta_1 = np.deg2rad(input_vals["theta_1"]) # value in input should be given in degrees 
            self.theta_2 = np.deg2rad(input_vals["theta_2"]) # value in input should be given in degrees 
            
            # self.R_g = 8314.41/self.m_weight 
            self.end_epsilon = input_vals["end_epsilon"]

    def calc_mu_1(self):
        """This function calculates the forward mach angle"""
        self.mu_1 = np.arcsin(1/self.M_1)


    def calc_v_M1(self): # this is an angle value
        """This function calculates the angle as a function of the entrance mach"""
        self.v_M1 = ((self.gamma+1)/(self.gamma-1))**0.5*np.arctan((((self.gamma-1)/(self.gamma+1))*(self.M_1**2-1))**0.5)-np.arctan((self.M_1**2-1)**0.5)

    
    def calc_v_M2(self, M_2guess): # this is an angle value
        """This function calculates the angle as a function of the mach after the turn"""
        v_M2 = ((self.gamma+1)/(self.gamma-1))**0.5*np.arctan((((self.gamma-1)/(self.gamma+1))*(M_2guess**2-1))**0.5)-np.arctan((M_2guess**2-1)**0.5) 
        return v_M2


    def calc_del_v_M2(self, M_2guess):
        """This function is the derivative of the Prandtl-meyer equation"""
        del_v_M2 = (1/M_2guess)*(((M_2guess**2-1)**0.5)/(1+((self.gamma-1)/2)*(M_2guess**2)))
        return del_v_M2


    def M_2_newton_solver(self, m_j):
        """This function numerically solves for the mach number given the"""
        epsilon = 10 #arbitrary placeholder
        while epsilon >= self.end_epsilon:
            v_m_j = self.calc_v_M2(m_j)
            del_v_del_Mj = self.calc_del_v_M2(m_j)
            delta = ((self.theta_1 + self.v_M1 - v_m_j) /del_v_del_Mj) 
            m_j = m_j + delta
            # print("m_j",m_j)
            epsilon = abs(delta)
            # print("epsilon", epsilon)
        self.M_2 = m_j
        return m_j
    
    
    def backward_newton_solver(self, m_j):
        """This function numerically solves for the mach number given the"""
        epsilon = 10 #arbitrary placeholder
        while epsilon >= self.end_epsilon:
            v_m_j = self.calc_v_M2(m_j)
            del_v_del_Mj = self.calc_del_v_M2(m_j)
            delta = ((np.radians(20) + v_m_j - 1.8) /del_v_del_Mj) 
            m_j = m_j + delta
            # print("m_j",m_j)
            epsilon = abs(delta)
            # print("epsilon", epsilon)
        back_m2 = m_j
        return back_m2
    
    
    def calc_mu_2(self):
        """This function solves for mu_2 given M_2"""
        self.mu_2 = np.arcsin(1/self.M_2)

    
    def calc_exit_mach_wave_angle(self):
        """This function calculates the exit mach wave angle"""
        self.exit_mach_wave_angle = self.mu_2 - self.theta_1


    def calc_P_o1(self):
        """Solves the stagnation pressure given the static pressure and mach number at the entrance"""
        self.P_o1 = self.P_1*((1 + ((self.gamma-1)/2)*self.M_1**2))**(self.gamma/(self.gamma-1))


    def calc_T_o1(self):
        """Solves the stagnation pressure given the static pressure and mach number at the entrance"""
        self.T_o1 = self.T_1*((1 + ((self.gamma-1)/2)*self.M_1**2))


    def calc_P_o2(self):
        """ISENTROPIC"""
        self.P_o2 = self.P_o1


    def calc_T_o2(self):
        """ISENTROPIC"""
        self.T_o2 = self.T_o1

    
    def calc_P_2(self):
        """Calculates the static pressure based on the stagnation pressure"""
        self.P_2 = self.P_o2/(((1 + ((self.gamma-1)/2)*self.M_2**2))**(self.gamma/(self.gamma-1)))

    
    def calc_T_2(self):
        """Calculates the static temperature based on stag temp and mach"""
        self.T_2 = self.T_o2/((1 + ((self.gamma-1)/2)*self.M_2**2))

    def calc_theta_two_total(self):
        """calculates sum of two thetas"""
        self.theta_two_total = self.theta_1 + self.theta_2

    ### now, calculate the beta angle to do the shock part of the question (section 2 to 3)
    def calc_a_b_c(self):
        """the a part of the cubic function"""
        self.a = (1 + ((self.gamma-1)/2)*self.M_2**2)*np.tan(self.theta_two_total) #### theta only works again because it's the same for the second station
        self.b = (self.M_2**2-1)
        self.c = (1 + ((self.gamma+1)/2)*self.M_2**2)*np.tan(self.theta_two_total) #### theta only works again because it's the same for the second station

    
    ### calculate beta explicitly
    def calc_lambda(self):
        """the lambda coefficient used for solving beta"""
        self.lamb = (self.b**2- 3*self.a*self.c)**0.5

    
    def calc_chi(self):
        """the chi coefficient used for solving beta"""
        self.chi = (self.b**3 - 9*self.a*(self.a+((self.gamma+1)/4)*self.M_2**4*np.tan(self.theta_two_total)))/(self.lamb**3)
    

    def calc_beta(self):
        """beta value used for calculating M3"""
        delta = 1
        numerator = self.b + 2*self.lamb*np.cos((4*np.pi*delta + np.arccos(self.chi))/(3))
        denominator = 3*self.a
        self.beta = np.arctan(numerator/denominator)

    ## now, calculate Mn2 
    def calc_Mn2(self):
        """This function is used to calculate properties as they cross normal to the oblique shock"""
        self.Mn2 = self.M_2*np.sin(self.beta)

    
    def calc_Mn3(self):
        """This function is used to calculate properties as they cross normal to the oblique shock"""
        numerator = 1 + ((self.gamma-1)/2)*(self.M_2*np.sin(self.beta))**2
        denominator = self.gamma*(self.M_2*np.sin(self.beta))**2 - ((self.gamma-1)/2)
        self.Mn3 = np.sqrt((numerator)/(denominator))
        
    
    def calc_M_3(self):
        """this function calculates M_3 using Mn3"""
        self.M_3 = (self.Mn3)/(np.sin(self.beta-self.theta_two_total))

    
    def calc_P_3(self):
        """this function calculates P_3 using Mn1"""
        self.P_3 = self.P_2*(1+(2*self.gamma/(self.gamma+1))*(self.Mn2**2-1))

    
    def calc_T_3(self):
        """this function calculates T_3 using Mn1"""
        self.T_3 = self.T_2*(1+(2*self.gamma/(self.gamma+1))*(self.Mn2**2-1))*((2+(self.gamma-1)*self.Mn2**2)/((self.gamma+1)*self.Mn2**2))


    def calc_P_o3(self):
        """This function calculates P_o3"""
        self.P_o3 = self.P_3*((1 + ((self.gamma-1)/2)*self.M_3**2))**(self.gamma/(self.gamma-1))

    def calc_T_o3(self):
        """This function calcs T_o3"""
        self.T_o3 = self.T_3*((1 + ((self.gamma-1)/2)*self.M_3**2))

    
    



    def run_prandtl_meyer(self):
        """This function runs all of the necessary functions to run the isentropic solver"""
        self.calc_mu_1()
        print("\nEntry Mach Wave Angle:\n", np.degrees(self.mu_1), "degrees\n")
        self.calc_v_M1()
        print("\nv(M1):\n", np.degrees(self.v_M1), " degrees\n")
        if self.M_1 > 1:
            M_2 = self.M_2_newton_solver(self.M_1)
        else:
            print("Can't find M_2 for a non supersonic M_1 because expansion fans only exist in supersonic flows")
        print("\nM2:\n", M_2, "\n")

        self.backward_newton_solver(1.8)
        print("backwards_m_2")

        self.calc_mu_2()
        print("\nmu_2:\n", np.degrees(self.mu_2), "degrees\n")
        self.calc_exit_mach_wave_angle()
        print("\nExit Mach Wave Angle:\n", np.degrees(self.exit_mach_wave_angle), "degrees\n")
        self.calc_P_o1()
        print("\nP_o1:\n", self.P_o1, "Pa\n")
        self.calc_T_o1()
        print("\nT_o1:\n", self.T_o1, "K\n")
        self.calc_P_o2()
        print("\nP_o2:\n", self.P_o2, "Pa\n")
        self.calc_T_o2()
        print("\nT_o2:\n", self.T_o2, "K\n")
        self.calc_P_2()
        print("\nP_2:\n", self.P_2, "Pa\n")
        self.calc_T_2()
        print("\nT_2:\n", self.T_2, "K\n")
        self.calc_theta_two_total()
        print("\ntheta 2 total:\n", np.degrees(self.theta_two_total), "degrees\n")
        self.calc_a_b_c()
        print("\na:\n", self.a, "\n")
        print("\nb:\n", self.b, "\n")
        print("\nc:\n", self.c, "\n")
        self.calc_lambda()
        print("\nLambda:\n", self.lamb, "\n")
        self.calc_chi()
        print("\nChi:\n", self.chi, "\n")
        self.calc_beta()
        print("\nShock angle (Beta):\n", np.degrees(self.beta), "degrees\n")
        self.calc_Mn2()
        print("\nMn2:\n", self.Mn2, "\n")
        self.calc_Mn3()
        print("\nMn3:\n", self.Mn3, "\n")
        self.calc_M_3()
        print("\nM_3:\n", self.M_3, "\n")
        self.calc_P_3()
        print("\nP_3:\n", self.P_3, "Pa\n")
        self.calc_T_3()
        print("\nT_3:\n", self.T_3, "K\n")
        self.calc_P_o3()
        print("\nP_o3:\n", self.P_o3, "Pa\n")
        self.calc_T_o3()
        print("\nT_o3:\n", self.T_o3, "K\n")


if __name__ == "__main__":
    """This is where the code actually runs"""

    expander = Prandtl_Meyer("Prandtl-Meyer.json") # this line initializes the nozzle class 
    expander.run_prandtl_meyer()

    

