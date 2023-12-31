import json 
import numpy as np 
import math
import matplotlib.pyplot as plt 

class nozzle_solver:
    """This class contains important nozzle attributes that can assist in nozzle design"""
    def __init__(self, json_file):
        self.json_file = json_file 


    # First, here are some functions that allow us to take the information 
    # from the text file so that we can plot it. 
    def load_json(self):
        """This function loads the json file that references the text file"""
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.nozzle_points = input_vals['nozzle_points'] # in cm
            self.gamma = input_vals["gamma"] 
            self.P_e = input_vals["P_e"]*10**3 # in Pa
            self.P_o = input_vals["P_o"]*10**6 # in Pa
            self.T_o = input_vals["T_o"] # in K
            self.m_weight = input_vals["m_weight"] # in g/mol
            self.end_epsilon = input_vals["end_epsilon"]


    def pull_nozzle_points(self):
        """This function parses out the text file so it can be used in this code."""
        nozzle_points = self.nozzle_points
        with open(nozzle_points, 'r') as text_handle:
            points = [list(map(float, line.strip().split())) for line in text_handle]
        self.nozzle_profile = np.array(points)


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
    
    
    def calc_epsilon(self, m_j, expansion_ratio):
        """This function retrieves the value that can trigger the stopping condition in newton solver"""
        m_j_sq = m_j**2
        # first, second, third, and exponent are just portions of the equation for del_F/del_M. 
        # Splitting it up this way allows for a code that's easier to write and read.
        first = 1/m_j
        second = 2/(self.gamma + 1)
        third = 1 + ((self.gamma - 1)/2)*m_j_sq
        exponent = (self.gamma + 1)/(2*(self.gamma-1))
        epsilon = abs(((first * (second * third)**exponent)) - expansion_ratio)/expansion_ratio
        return epsilon


    def newton_solver(self, start_mach, expansion_ratio):
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
        return m_j
    

    # Now, here are some functions that help to find R_g, Thrust, isp, etc...
    def calc_R_g(self):
        """This function calculates the gas constant"""
        R_u = 8314.4612 #j/kg-kmol-k 
        self.R_g = R_u/self.m_weight # because g/mol = kg/kmol 


    def calc_exit_and_throat_areas(self):
        """This function calculates the area of the throat and the area of the exit"""
        radius_throat = self.nozzle_profile[3][1] # this is the throat location
        print("radius_throat", radius_throat)
        self.a_star = (np.pi * radius_throat**2)*(1/(100**2))
        radius_exit = self.nozzle_profile[len(self.nozzle_profile) - 1, 1] 
        print("exit_radius", radius_exit)
        self.a_e = (np.pi * radius_exit**2)*(1/(100**2))
        self.a_exit_over_a_star = self.a_e/self.a_star
        print("exit area over throat area:", self.a_exit_over_a_star)
        

    def calc_is_choked(self):
        """This function checks if the throat is choked using mass flow rates"""
        first = 2/(self.gamma - 1)
        second = (self.P_o/self.P_e)
        exponent = (self.gamma - 1)/self.gamma
        self.M_e = np.sqrt(first * ((second**exponent) - 1))
        self.T_e = self.T_o/(1 + (((self.gamma-1)/(2))*self.M_e**2))
        self.V_e = self.M_e*np.sqrt(self.gamma * self.R_g * self.T_e)
        self.rho_e = (self.P_e)/(self.R_g*self.T_e) 
        self.m_dot_actual = self.rho_e*self.a_e*self.V_e
        # self.choked = True
        # print("actual", self.m_dot_actual)
        self.m_dot_required = (self.P_o/np.sqrt(self.T_o))*self.a_star*np.sqrt((self.gamma/self.R_g)*(2/(self.gamma + 1))**((self.gamma + 1)/(self.gamma - 1)))
        # print("required", self.m_dot_required)
        if self.m_dot_actual < self.m_dot_required: 
            self.choked = False
        elif self.m_dot_actual >= self.m_dot_required: 
            self.choked = True 


    def calc_thrust(self):
        """This function calculates the thrust of the nozzle"""
        self.thrust = self.m_dot_actual*self.V_e + self.P_e*self.a_e


    def calc_ISP(self):
        """This function calculates the ISP of the nozzle based on the thrust"""
        self.gravity = 9.81 #m/s^2
        self.ISP = (1/self.gravity)*(self.thrust/self.m_dot_actual)


    # This part does the calculations across the length of the nozzle
    def calc_area_ratio(self):
        """This function calculates the ratio of the desired area over the throat area"""
        area_vec = []
        for i in range(len(self.nozzle_profile)):
            radius = self.nozzle_profile[i, 1]
            # print("radius", radius)
            a_e = (np.pi*(radius**2))*(1/(100**2))
            ratio = a_e/self.a_star
            # print("ratio:", ratio)
            area_vec.append(ratio)
        self.areas = np.array(area_vec)


    def calc_machs(self):
        """This function calculates the mach number as the gas moves down the throat"""
        mach_val = []
        for i in range(len(self.areas)):
            if i < 3:
                self.start_mach = 0.5
                mach_number = self.newton_solver(self.start_mach, self.areas[i])
            elif i >= 3:
                self.start_mach = 3.0
                mach_number = self.newton_solver(self.start_mach, self.areas[i])
            mach_val.append(mach_number)
        self.machs = np.array(mach_val)
        
    
    def calc_temps(self):
        """This function plots the temperature profile over the nozzle"""
        temp_val = []
        for i in range(len(self.machs)):
            temp = self.T_o/(1 + ((self.gamma - 1)/self.gamma) * self.machs[i]**2)
            temp_val.append(temp)
        self.temps = np.array(temp_val)

    
    def calc_pressures(self):
        """This function plots the pressure profile over the nozzle in kPa"""
        pressure_val = []
        for i in range(len(self.machs)):
            pressure = self.P_o*10**-6/((1 + ((self.gamma - 1)/self.gamma) * self.machs[i]**2)**((self.gamma)/(self.gamma - 1)))
            pressure_val.append(pressure)
        self.pressures = np.array(pressure_val)


    # Here are the plots 
    def plot_nozzle(self):
        plt.figure("Nozzle")
        plt.plot(self.nozzle_profile[:,0], self.nozzle_profile[:,1], color = "red")
        plt.plot(self.nozzle_profile[:,0], self.nozzle_profile[:,2], color = "red")
    

    def plot_mach_numbers(self):
        plt.figure("Mach Numbers")
        plt.plot(self.nozzle_profile[:,0], self.machs[:], color = "blue")


    def plot_temps(self):
        plt.figure("Temperatures [k]")
        plt.plot(self.nozzle_profile[:,0], self.temps[:], color = "green")


    def plot_pressures(self):
        plt.figure("Pressures [MPa]")
        plt.plot(self.nozzle_profile[:,0], self.pressures[:], color = "green")


    def run_solver(self):
        """When called, this function runs all the needed functions to complete HW 4, problem 2"""
        self.load_json()
        self.pull_nozzle_points()
        self.calc_R_g()
        self.calc_exit_and_throat_areas()
        self.calc_is_choked()
        self.calc_area_ratio()
        self.calc_machs()
        self.calc_temps()
        self.calc_pressures()
        self.calc_thrust()
        self.calc_ISP()
        print("T_exit [k]", self.T_e)
        print("T_o [k]",self.T_o)
        print("m_dot [kg/s]", self.m_dot_actual)
        print("R_g[J/kg-k]:", self.R_g)
        print("Exit_Mach:", self.M_e)
        print("Exit_vel[m/s]:", self.V_e)
        print("thrust [N]", self.thrust)
        print("ISP[s]:", self.ISP)
        
        nozzle.plot_nozzle()
        nozzle.plot_mach_numbers()
        nozzle.plot_temps()
        nozzle.plot_pressures()


if __name__ == "__main__":
    """This is where the code actually runs"""
    nozzle = nozzle_solver("nozzle_info.json") # this line initializes the nozzle class 
    run = nozzle.run_solver()
    plt.show()

