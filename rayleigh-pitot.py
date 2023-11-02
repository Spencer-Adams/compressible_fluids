import json 
import numpy as np 

class pitot_rayleigh:
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
            self.gamma = input_vals["gamma"] 
            self.P_inf = input_vals["P_inf"]
            self.P_o_2 = input_vals["P_o_2"]
            self.end_epsilon = input_vals["end_epsilon"]


    # Now, here are the functions associated with the newton solver (HW 6, problem 3, pt 1)
    def calc_f_m_j(self, m_j, pressure_ratio):
        """This function returns a value for F(M_j)"""
        m_j_sq = m_j**2
        # Splitting it up this way allows for a code that's easier to write and read.
        numerator = (((self.gamma+1)/2)*m_j_sq)**(self.gamma/(self.gamma-1))
        denominator = (((2*self.gamma)/(self.gamma + 1))*m_j_sq-((self.gamma-1)/(self.gamma+1)))**(1/(self.gamma-1))
        F_M_j = (numerator/denominator) - pressure_ratio
        return F_M_j


    def calc_del_f_del_M(self, m_j):
        """This function returns a value for the derivative of F with respect to M, del_F/del_M"""
        m_j_sq = m_j**2
        # Splitting it up this way allows for a code that's easier to write and read.
        numerator = self.gamma*m_j*(2*m_j_sq-1)*(m_j_sq*((self.gamma+1)/2))**(1/(self.gamma-1))
        denominator = (((2*self.gamma)/(self.gamma+1))*m_j_sq-((self.gamma-1)/(self.gamma+1)))**(self.gamma/(self.gamma-1))
        del_F_DEL_M = (numerator/denominator)
        return del_F_DEL_M
    

    def calc_is_critical(self, P_o_2, P_inf):
        """This function takes in pressures and then determines if the flow is critical or not"""
        is_critical = False
        ratio = (P_o_2/P_inf)
        gammas = ((self.gamma+1)/(2))**(self.gamma/(self.gamma-1))
        if ratio > gammas:
            is_critical = True 
        else:
            is_critical = False 
        return is_critical
    
    
    def calc_epsilon(self, m_j, pressure_ratio):
        """This function retrieves the value that can trigger the stopping condition in newton solver"""
        # Splitting it up this way allows for a code that's easier to write and read.
        f_m_j = self.calc_f_m_j(m_j, pressure_ratio)
        epsilon = abs(f_m_j)/pressure_ratio
        return epsilon


    def newton_solver(self, m_j, pressure_ratio):
        """This function numerically solves for the mach number given the"""
        epsilon = 10 #arbitrary placeholder            
        while epsilon >= self.end_epsilon:
            # print("mach_number: ",m_j)
            f_m_j = self.calc_f_m_j(m_j, pressure_ratio)
            del_f_del_M = self.calc_del_f_del_M(m_j)
            M_j_plus_one = m_j - (f_m_j/del_f_del_M)
            epsilon = self.calc_epsilon(m_j, pressure_ratio)
            m_j = M_j_plus_one
        return m_j
    

    def calc_machs(self):
        """This function calculates the mach number at the different pressure values"""
        mach_val = []
        for i in range(len(self.P_inf)):
            pressure_ratio = self.P_o_2[i]/self.P_inf[i]
            is_critical = self.calc_is_critical(self.P_o_2[i], self.P_inf[i])
            if is_critical:
                self.start_mach = 3.0
                mach_number = self.newton_solver(self.start_mach, pressure_ratio)
            else:
                self.start_mach = 0.1
                mach_number = self.newton_solver(self.start_mach, pressure_ratio)
            mach_val.append(mach_number)
        machs = np.array(mach_val)
        self.machs = np.real(machs)
        
    
    def run_solver(self):
        """When called, this function runs all the needed functions to complete HW 4, problem 2"""
        self.calc_machs()
        print("mach numbers:", self.machs)
        

if __name__ == "__main__":
    """This is where the code actually runs"""
    pitot = pitot_rayleigh("rayleigh-pitot.json") # this line initializes the nozzle class 
    run = pitot.run_solver()

