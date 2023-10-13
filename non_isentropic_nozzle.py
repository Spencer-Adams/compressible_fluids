import json 
import numpy as np 
import math
import matplotlib.pyplot as plt 

class isentropic_nozzle_solver:
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
            self.one_ATM = input_vals["one_ATM"] 
            self.P_o_1 = input_vals["P_o_1"]*self.one_ATM # in Pa
            # self.P_o_2 = input_vals["P_o_2"]*self.one_ATM # in Pa
            self.P_o_2 = self.P_o_1 # in Pa
            # self.P_o_1 = self.P_o_2 # in Pa
            self.T_o_1 = input_vals["T_o_1"] # in K
            self.m_weight = input_vals["m_weight"] # in g/mol=kg/kmol
            self.R_g = 8314.41/self.m_weight 
            self.A_t = input_vals["A_t"] # in m^2
            self.A_e = input_vals["A_e"] # in m^2
            self.end_epsilon = input_vals["end_epsilon"]


    def calc_mass_flow_rate(self): # pt a
        """This function calculates the non-isentropic flow rate based on the throat area"""
        first = ((self.A_t * self.P_o_1)/np.sqrt(self.T_o_1))
        second = self.gamma/self.R_g
        third = 2/(self.gamma + 1)
        exponent = (self.gamma + 1)/(self.gamma - 1)
        self.m_dot_isentropic = first*np.sqrt(second*((third)**exponent))


    # Now, here are the functions associated with the newton solver (HW 4, problem 1)
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
    

    def calc_expansion_ratio(self):
        """This function defines the stagnation pressure ratio at the exit over the throat"""
        self.expansion_ratio = self.A_e/self.A_t


    def calc_epsilon(self, m_j): #### could be a non class-specific function. 
        """This function retrieves the value that can trigger the stopping condition in newton solver"""
        epsilon = abs(self.calc_f_m_j(m_j))/self.expansion_ratio
        return epsilon
    

    def isentropic_newton_solver(self, start_mach):
        """This function numerically solves for the mach number given the"""
        m_j = start_mach
        epsilon = 10 #arbitrary placeholder
        while epsilon >= self.end_epsilon:
            # print("mach_number: ",m_j)
            f_m_j = self.calc_f_m_j(m_j)
            del_f_del_M = self.calc_del_f_del_M(m_j)
            M_j_plus_one = m_j - (f_m_j/del_f_del_M)
            epsilon = self.calc_epsilon(m_j)
            m_j = M_j_plus_one
        self.m_e = m_j
        return m_j
    

    def calc_exit_pressure(self):  #### adjust to test P_o_2 as well 
        """This function plots the pressure profile over the nozzle in kPa"""
        self.P_e = self.P_o_2/((1 + ((self.gamma - 1)/2) * self.m_e**2)**((self.gamma)/(self.gamma - 1)))


    def calc_exit_temperature(self): 
        """This function plots the temperature profile over the nozzle"""
        self.T_o_2 = self.T_o_1
        self.T_e = self.T_o_2/(1 + ((self.gamma - 1)/2) * self.m_e**2)


    def calc_rho_e(self):
        """This function calculates the flow density at the nozzle exit"""
        self.rho_e = (self.P_e)/(self.R_g*self.T_e) 


    def calc_V_e(self):
        """This function calculates the exit velocity based on exit mach number"""
        self.V_e = self.m_e*np.sqrt(self.gamma*self.R_g*self.T_e)


    def calc_m_dot_e(self): #### maybe delete?
        """This function calculates the mass flow at the nozzle exit"""
        self.m_dot_e = self.rho_e*self.A_e*self.V_e


    def calc_m_dot_e_required(self):
        """This function calculates the required mass flow for choked flow at the exit"""
        self.m_dot_required = (self.P_o_2/np.sqrt(self.T_o_2))*self.A_t*np.sqrt((self.gamma/self.R_g)*(2/(self.gamma + 1))**((self.gamma + 1)/(self.gamma - 1)))


    def is_choked(self):
        """this function compares the m_dot_required to m_dot_actual"""
        if self.m_dot_e < self.m_dot_required: 
            self.choked = False
        elif self.m_dot_e >= self.m_dot_required: 
            self.choked = True 


    def calc_thrust_momentum(self):
        """This function calculates the thrust do only to momentum (same for vacuum and sea level)"""
        self.thrust_momentum = self.m_dot_e*self.V_e


    def calc_total_thrust(self):
        """This function adds the pressure parts of the thrust equation to the momentum part"""
        self.total_thrust = self.thrust_momentum + (self.P_e*self.A_e - self.one_ATM*self.A_e)
    

    def calc_ISP(self):
        """This function calculates the ISP of the nozzle"""
        self.ISP = (1/9.807)*(self.total_thrust/self.m_dot_e)


    def calc_vacuum_thrust_momentum(self):
        """This function calculates the thrust do only to momentum (same for vacuum and sea level)"""
        self.vacuum_thrust_momentum = self.m_dot_e*self.V_e


    def calc_vacuum_total_thrust(self):
        """This function adds the pressure parts of the vacuum thrust equation to the momentum part"""
        self.vacuum_total_thrust = self.vacuum_thrust_momentum + self.P_e*self.A_e


    def calc_vacuum_ISP(self):
        """This function calculates the vacuum ISP of the nozzle"""
        self.vacuum_ISP = (1/9.807)*(self.vacuum_total_thrust/self.m_dot_e)


    def run_isentropic_solver(self):
        """This function runs all of the necessary functions to run the isentropic solver"""
        print("Isentropic case:")
        self.load_project_one_json()
        self.calc_mass_flow_rate()
        print("mass flow:", self.m_dot_isentropic)
        self.calc_expansion_ratio()
        m_e = self.isentropic_newton_solver(3)
        print("m_e: ", m_e)
        self.calc_exit_pressure()
        print("P_e [Pa]: ", self.P_e)
        self.calc_exit_temperature()
        print("T_e [k]: ", self.T_e)
        self.calc_rho_e()
        print("rho_e [kg/m^3]:", self.rho_e)
        self.calc_V_e()
        print("V_e [m/s]:", self.V_e)
        self.calc_m_dot_e()
        print("m_dot_e [kg/s]:", self.m_dot_e)
        self.calc_m_dot_e_required()
        print("m_dot_e_required [kg/s]", self.m_dot_required)
        self.is_choked()
        print("Is it choked?", self.choked)
        self.calc_thrust_momentum()
        print("thrust (momentum) [N]:", self.thrust_momentum)
        self.calc_total_thrust()
        print("total thrust: [N]", self.total_thrust)
        self.calc_ISP()
        print("ISP: [s]", self.ISP)
        self.calc_vacuum_thrust_momentum()
        print("vacuum thrust (momentum) [N]:", self.vacuum_thrust_momentum)
        self.calc_vacuum_total_thrust()
        print("vacuum total thrust: [N]", self.vacuum_total_thrust)
        self.calc_vacuum_ISP()
        print("vacuum ISP: [s]", self.vacuum_ISP, "\n")


# Here's the non_isentropic stuff we need.
class non_isentropic_nozzle_solver:
    """This class contains important nozzle attributes that can assist in nozzle design"""
    def __init__(self, json_file):
        self.json_file = json_file 


    def load_project_one_json(self):
        """This function loads the json file that references the text file"""
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.gamma = input_vals["gamma"]
            self.one_ATM = input_vals["one_ATM"] 
            self.P_o_1 = input_vals["P_o_1"]*self.one_ATM # in Pa
            self.P_o_2 = input_vals["P_o_2"]*self.one_ATM # in Pa
            self.T_o_1 = input_vals["T_o_1"] # in K
            self.m_weight = input_vals["m_weight"] # in g/mol=kg/kmol
            self.R_g = 8314.41/self.m_weight 
            self.A_t = input_vals["A_t"] # in m^2
            self.A_e = input_vals["A_e"] # in m^2
            self.end_epsilon = input_vals["end_epsilon"]


    def calc_mass_flow_rate(self): # pt a
        """This function calculates the non-isentropic flow rate based on the throat area"""
        first = ((self.A_t * self.P_o_1)/np.sqrt(self.T_o_1))
        second = self.gamma/self.R_g
        third = 2/(self.gamma + 1)
        exponent = (self.gamma + 1)/(self.gamma - 1)
        self.m_dot_non_isentropic = first*np.sqrt(second*((third)**exponent))


    # Now, here are the functions associated with the non_isentropic newton solver 
    def calc_g_m_j(self, m_j):
        """This function returns a value for g(M_j)"""
        m_j_sq = m_j**2
        # first, second, third, and exponent are just portions of the equation for g(M_j). 
        # Splitting it up this way allows for a code that's easier to write and read.
        exponent_one = 1 / (self.gamma - 1)
        exponent_two = self.gamma/(self.gamma - 1)
        first_num = 2
        first_denom = (self.gamma + 1)*((self.gamma * m_j_sq)-((self.gamma-1)/2))**exponent_one 
        first = first_num/first_denom 
        second_num = (((self.gamma + 1)/2)*m_j)**2
        second_denom = (1 + ((self.gamma - 1)/2)*m_j_sq)
        second = second_num/second_denom 
        G_M_j = (first * ((second)**exponent_two)) - (self.P_o_2/self.P_o_1)
        return G_M_j
    

    def calc_del_g_del_M(self, m_j):
        """This function returns a value for the derivative of F with respect to M, del_F/del_M"""
        m_j_sq = m_j**2
        # first, second, third, and exponent are just portions of the equation for del_G/del_M. 
        # Splitting it up this way allows for a code that's easier to write and read.
        big_denom = (self.gamma + 1)*m_j*(2 + m_j_sq*(self.gamma - 1))*(1 + self.gamma*(2*m_j_sq - 1))
        exponent_one = 3 - ((2*self.gamma)/(self.gamma - 1))
        exponent_two = self.gamma/(self.gamma - 1)
        exponent_three = ((-1)/(self.gamma - 1))
        first = (2**(exponent_one))*self.gamma*((m_j_sq-1)**2)
        second = ((((self.gamma + 1)*m_j)**2)/(1 + ((self.gamma-1)/2)*m_j_sq))**exponent_two
        third = (0.5 + self.gamma*(m_j_sq - 0.5))**exponent_three
        del_G_DEL_M = -(first*second*third)/big_denom
        return del_G_DEL_M
    

    def calc_stag_pressure_ratio(self):
        """This function defines the stagnation pressure ratio at the exit over the throat"""
        self.stag_pressure_ratio = self.P_o_2/self.P_o_1


    def calc_M_1_epsilon(self, m_j): #### could be a non class-specific function. 
        """This function retrieves the value that can trigger the stopping condition in newton solver"""
        epsilon = abs(self.calc_g_m_j(m_j))/self.stag_pressure_ratio
        return epsilon
    

    def M_1_newton_solver(self, start_mach): #### this solves for M1
        """This function numerically solves for the mach number given the pressure ratio"""
        m_j = start_mach
        epsilon = 10 #arbitrary placeholder
        while epsilon >= self.end_epsilon:
            # print("mach_number: ",m_j)
            g_m_j = self.calc_g_m_j(m_j)
            del_g_del_M = self.calc_del_g_del_M(m_j)
            M_j_plus_one = m_j - (g_m_j/del_g_del_M)
            epsilon = self.calc_M_1_epsilon(m_j)
            m_j = M_j_plus_one
        self.M_1 = m_j
        return m_j
    

    def calc_A_ratio_at_shock(self, m_1):
        """This calculates the Area ratio at the shock"""
        m_1_sq = m_1**2
        first = 1/m_1
        second = 2/(self.gamma + 1)
        third = 1 + ((self.gamma-1)/2)*m_1_sq
        exponent = (self.gamma+1)/(2*(self.gamma - 1))
        self.a_shock_over_a_throat = first*(second*third)**exponent


    def calc_A_2(self):
        """This takes advantage of the fact that the pressure times the charact area on either side of the shock is equal"""
        self.A_two = (self.P_o_1*self.A_t)/self.P_o_2


    def calc_A_e_over_A_2(self):
        """This calculates A_e over A_1 to run the next (post-shock) newton solver"""
        self.A_e_over_A2 = self.A_e/self.A_two


        # Now, here are the functions associated with the newton solver after the shock 
    def calc_f_m_j(self, m_j):
        """This function returns a value for F(M_j)"""
        m_j_sq = m_j**2
        # first, second, third, and exponent are just portions of the equation for f(M_j). 
        # Splitting it up this way allows for a code that's easier to write and read.
        first = 1/(m_j)
        second = 2/(self.gamma + 1)
        third = 1+(((self.gamma - 1)/2)*m_j_sq)
        exponent = (self.gamma + 1) / (2*(self.gamma - 1))
        F_M_j = (first)*(second * third)**(exponent) - self.A_e_over_A2
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
    
    
    def calc_M_2_epsilon(self, m_j): #### could be a non class-specific function. 
        """This function retrieves the value that can trigger the stopping condition in newton solver"""
        epsilon = abs(self.calc_f_m_j(m_j))/self.A_e_over_A2
        return epsilon


    def M_2_newton_solver(self, start_mach):
        """This function numerically solves for the mach number given the"""
        m_j = start_mach
        epsilon = 10 #arbitrary placeholder
        while epsilon >= self.end_epsilon:
            f_m_j = self.calc_f_m_j(m_j)
            del_f_del_M = self.calc_del_f_del_M(m_j)
            M_j_plus_one = m_j - (f_m_j/del_f_del_M)
            epsilon = self.calc_M_2_epsilon(m_j)
            m_j = M_j_plus_one
        self.m_e = m_j
        return m_j
    

    def M_2(self):
        """This calculates M_2"""
        numerator = (1 + ((self.gamma - 1)/2)*self.M_1**2)
        denominator = (self.gamma*self.M_1**2 - (self.gamma - 1)/2)
        self.m_ee = np.sqrt(numerator/denominator)
        
    
    def calc_exit_pressure(self):  #### adjust to test P_o_2 as well 
        """This function plots the pressure profile over the nozzle in kPa"""
        self.P_e = self.P_o_2/((1 + ((self.gamma - 1)/2) * self.m_e**2)**((self.gamma)/(self.gamma - 1)))


    def calc_exit_temperature(self): 
        """This function plots the temperature profile over the nozzle"""
        self.T_o_2 = self.T_o_1
        self.T_e = self.T_o_2/(1 + ((self.gamma - 1)/2) * self.m_e**2)


    def calc_rho_e(self):
        """This function calculates the flow density at the nozzle exit"""
        self.rho_e = (self.P_e)/(self.R_g*self.T_e) 


    def calc_V_e(self):
        """This function calculates the exit velocity based on exit mach number"""
        self.V_e = self.m_e*np.sqrt(self.gamma*self.R_g*self.T_e)


    def calc_m_dot_e(self):
        """This function calculates the mass flow at the nozzle exit"""
        self.m_dot_e = self.rho_e*self.A_e*self.V_e


    def calc_thrust_momentum(self):
        """This function calculates the thrust do only to momentum (same for vacuum and sea level)"""
        self.thrust_momentum = self.m_dot_e*self.V_e


    def calc_total_thrust(self):
        """This function adds the pressure parts of the thrust equation to the momentum part"""
        self.total_thrust = self.thrust_momentum + (self.P_e*self.A_e - self.one_ATM*self.A_e)
    

    def calc_ISP(self):
        """This function calculates the ISP of the nozzle"""
        self.ISP = (1/9.807)*(self.total_thrust/self.m_dot_e)


    def calc_vacuum_thrust_momentum(self):
        """This function calculates the thrust do only to momentum (same for vacuum and sea level)"""
        self.vacuum_thrust_momentum = self.m_dot_e*self.V_e


    def calc_vacuum_total_thrust(self):
        """This function adds the pressure parts of the vacuum thrust equation to the momentum part"""
        self.vacuum_total_thrust = self.vacuum_thrust_momentum + self.P_e*self.A_e


    def calc_vacuum_ISP(self):
        """This function calculates the vacuum ISP of the nozzle"""
        self.vacuum_ISP = (1/9.807)*(self.vacuum_total_thrust/self.m_dot_e)


    def run_non_isentropic_solver(self):
        """This function runs all of the necessary functions to run the isentropic solver"""
        print("\n non_isentropic case:")
        self.load_project_one_json()
        self.calc_mass_flow_rate()
        print("mass flow:", self.m_dot_non_isentropic)
        self.calc_stag_pressure_ratio()
        m_1 = self.M_1_newton_solver(4) # 4 is an initial guess for M_1
        print("m_1: ", m_1)
        self.calc_A_ratio_at_shock(m_1)
        print("A_1/A*: ", self.a_shock_over_a_throat)
        self.calc_A_2()
        self.calc_A_e_over_A_2()
        m_e = self.M_2_newton_solver(0.1)
        print("m_e: ", m_e)
        self.M_2()
        print("other_m_e:", self.m_ee)
        self.calc_exit_pressure()
        print("P_e [Pa]: ", self.P_e)
        self.calc_exit_temperature()
        print("T_e [k]: ", self.T_e)
        self.calc_rho_e()
        print("rho_e [kg/m^3]:", self.rho_e)
        self.calc_V_e()
        print("V_e [m/s]:", self.V_e)
        self.calc_m_dot_e()
        print("m_dot_e [kg/s]:", self.m_dot_e)
        self.calc_thrust_momentum()
        print("thrust (momentum) [N]:", self.thrust_momentum)
        self.calc_total_thrust()
        print("total thrust: [N]", self.total_thrust)
        self.calc_ISP()
        print("ISP: [s]", self.ISP)
        self.calc_vacuum_thrust_momentum()
        print("vacuum thrust (momentum) [N]:", self.vacuum_thrust_momentum)
        self.calc_vacuum_total_thrust()
        print("vacuum total thrust: [N]", self.vacuum_total_thrust)
        self.calc_vacuum_ISP()
        print("vacuum ISP: [s]", self.vacuum_ISP, "\n")


if __name__ == "__main__":
    """This is where the code actually runs"""

    nozzle = non_isentropic_nozzle_solver("project_one_nozzle_info.json") # this line initializes the nozzle class 
    run = nozzle.run_non_isentropic_solver()


    nozzle_isentropic = isentropic_nozzle_solver("project_one_nozzle_info.json")
    run_isentropic = nozzle_isentropic.run_isentropic_solver()
    

