import json 
import numpy as np 
import math
import matplotlib.pyplot as plt 

class normal_oblique:
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
            self.M_inf = input_vals["M_inf"]
            self.gamma = input_vals["gamma"]
            self.one_ATM = input_vals["one_ATM"] 
            self.P_inf = input_vals["P_inf"]*self.one_ATM # in Pa
            self.T_inf = input_vals["T_inf"] # in K
            self.theta= np.deg2rad(input_vals["theta"]) # value in input should be given in degrees
            self.r = input_vals["r"] # in m  
            self.end_epsilon = input_vals["end_epsilon"]


    def calc_P_oinf(self):
        """Solves the stagnation pressure given the static pressure and mach number at the entrance"""
        self.P_oinf = self.P_inf*((1 + ((self.gamma-1)/2)*self.M_inf**2))**(self.gamma/(self.gamma-1))


    def calc_T_oinf(self):
        """Solves the stagnation pressure given the static pressure and mach number at the entrance"""
        self.T_oinf = self.T_inf*((1 + ((self.gamma-1)/2)*self.M_inf**2))


    ### now, calculate the beta angle to do the shock part of the question (section 2 to 3)
    def calc_a_b_c(self):
        """the a part of the cubic function"""
        self.a = (1 + ((self.gamma-1)/2)*self.M_inf**2)*np.tan(self.theta) 
        self.b = (self.M_inf**2-1)
        self.c = (1 + ((self.gamma+1)/2)*self.M_inf**2)*np.tan(self.theta)


        ### calculate beta explicitly
    def calc_lambda(self):
        """the lambda coefficient used for solving beta"""
        self.lamb = (self.b**2- 3*self.a*self.c)**0.5

    
    def calc_chi(self):
        """the chi coefficient used for solving beta"""
        self.chi = (self.b**3 - 9*self.a*(self.a+((self.gamma+1)/4)*self.M_inf**4*np.tan(self.theta)))/(self.lamb**3)
    

    def calc_beta(self):
        """beta value used for calculating M3"""
        delta = 1
        numerator = self.b + 2*self.lamb*np.cos((4*np.pi*delta + np.arccos(self.chi))/(3))
        denominator = 3*self.a
        self.beta = np.arctan(numerator/denominator) 

    
    ## now, calculate Mninf 
    def calc_Mninf(self):
        """This function is used to calculate properties as they cross normal to the oblique shock"""
        self.Mninf = self.M_inf*np.sin(self.beta)

    
    def calc_Mn1(self):
        """This function is used to calculate properties as they cross normal to the oblique shock"""
        numerator = 1 + ((self.gamma-1)/2)*(self.M_inf*np.sin(self.beta))**2
        denominator = self.gamma*(self.M_inf*np.sin(self.beta))**2 - ((self.gamma-1)/2)
        self.Mn1 = np.sqrt((numerator)/(denominator))
        
    
    def calc_M_1(self):
        """this function calculates M_3 using Mn3"""
        self.M_1 = (self.Mn1)/(np.sin(self.beta-self.theta))

    
    def calc_P_1(self):
        """this function calculates P_1 using Mninf"""
        self.P_1 = self.P_inf*(1+(2*self.gamma/(self.gamma+1))*(self.Mninf**2-1))

    
    def calc_T_1(self):
        """this function calculates T_3 using Mn1"""
        self.T_1 = self.T_inf*(1+(2*self.gamma/(self.gamma+1))*(self.Mninf**2-1))*((2+(self.gamma-1)*self.Mninf**2)/((self.gamma+1)*self.Mninf**2))


    def calc_P_o1(self):
        """This function calculates P_o1"""
        self.P_o1 = self.P_1*((1 + ((self.gamma-1)/2)*self.M_1**2))**(self.gamma/(self.gamma-1))


    def calc_T_o1(self):
        """This function calcs T_o1"""
        self.T_o1 = self.T_1*((1 + ((self.gamma-1)/2)*self.M_1**2))

    
    def calc_Mn2(self):
        """This finds Mn2"""
        numerator = 1 + ((self.gamma-1)/2)*(self.M_1*np.sin(np.pi/2))**2 # 
        denominator = self.gamma*(self.M_1*np.sin(np.pi/2))**2 - ((self.gamma-1)/2)
        self.Mn2 = np.sqrt((numerator)/(denominator))


    def calc_M_2(self):
        """M_2 is equal to Mn2 because you're dividing by one (beta is 90 and theta is zero)"""
        self.M_2 = self.Mn2

    
    def calc_P_2(self):
        """this function calculates P_2 using Mn1"""
        self.P_2 = self.P_1*(1+(2*self.gamma/(self.gamma+1))*(self.M_1**2-1))

    
    def calc_T_2(self):
        """this function calculates T_2 using Mn1"""
        self.T_2 = self.T_1*(1+(2*self.gamma/(self.gamma+1))*(self.M_1**2-1))*((2+(self.gamma-1)*self.M_1**2)/((self.gamma+1)*self.M_1**2))

    
    def calc_P_o2(self):
        """This function calculates P_o2"""
        self.P_o2 = self.P_2*((1 + ((self.gamma-1)/2)*self.M_2**2))**(self.gamma/(self.gamma-1))

    
    def calc_T_o2(self):
        """This function calcs T_o1"""
        self.T_o2 = self.T_2*((1+((self.gamma-1)/2)*self.M_2**2))

    
    def calc_stag_recovery(self):
        """This function calculates the ratio of Po2 over Poinfty"""
        self.stag_recovery = self.P_o2/self.P_oinf

    
    def run_normal_oblique(self):
        """This function runs all of the necessary functions to run the isentropic solver"""
        print("\nM_inf:\n", self.M_inf, "\n")
        self.calc_P_oinf()
        print("\nP_oinf:\n", self.P_oinf, "Pa\n")
        self.calc_T_oinf()
        print("\nT_oinf:\n", self.T_oinf, "K\n")
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
        print("\ntheta:\n", np.degrees(self.theta), "degrees\n")
        self.calc_Mninf()
        print("\nMninf:\n", self.Mninf, "\n")
        self.calc_Mn1()
        print("\nMn1:\n", self.Mn1, "\n")
        self.calc_M_1()
        print("\nM_1:\n", self.M_1, "\n")
        self.calc_P_1()
        print("\nP_1:\n", self.P_1, "Pa\n")
        self.calc_T_1()
        print("\nT_1:\n", self.T_1, "K\n")
        self.calc_P_o1()
        print("\nP_o1:\n", self.P_o1, "Pa\n")
        self.calc_T_o1()
        print("\nT_o1:\n", self.T_o1, "K\n")
        self.calc_Mn2()
        print("\nMn2:\n", self.Mn2, "\n")
        self.calc_M_2()
        print("\nM_2:\n", self.M_2, "\n")
        self.calc_P_2()
        print("\nP_2:\n", self.P_2, "Pa\n")
        self.calc_T_2()
        print("\nT_2:\n", self.T_2, "K\n")
        self.calc_P_o2()
        print("\nP_o2:\n", self.P_o2, "Pa\n")
        self.calc_T_o2()
        print("\nT_o2:\n", self.T_o2, "K\n")
        self.calc_stag_recovery()
        print("\nStag_recovery:\n", self.stag_recovery, "\n")


    def plotter(self):
        for i in range(1, 90):
            self.theta = np.deg2rad(i)
            self.run_normal_oblique()
            plt.scatter(np.degrees(self.theta), self.stag_recovery)
        plt.title("Stagnation recovery vs turning angle")
        plt.xlabel("turning angle [degrees]")
        plt.ylabel("stagnation recovery [Po2/Pinf]")
        plt.show()
        


if __name__ == "__main__":
    """This is where the code actually runs"""

    expander = normal_oblique("normal_oblique.json") # this line initializes the nozzle class 
    # expander.plotter()
    expander.run_normal_oblique()

    

