import json 
import numpy as np 
import sympy as sp
import matplotlib.pyplot as plt

class pitot_rayleigh: # only does supersonic relationship. 
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
            self.half_wedge = np.radians(input_vals["half_wedge[deg]"])
            self.wedge = 2*self.half_wedge
            self.alpha = np.radians(input_vals["alpha"])
            self.M1_start_guess = input_vals["M1_start_guess"]
            self.M_crit_start_guess = input_vals["M_crit_start_guess"]
            self.Xfoil_cpo = input_vals["Xfoil_cpo"]
            self.unswept_thickness = input_vals["unswept_thickness"]
            self.sweep = np.radians(input_vals["sweep[deg]"])
            self.alpha_2 = np.radians(input_vals["alpha_2[deg]"])
            self.alphas = np.array(input_vals["alphas"])
            self.CPos = np.array(input_vals["CPos"])


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
            print("true")
        else:
            is_critical = False
            print("false") 
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
        

    def calc_v_M4(self, M): # this is an angle value
        """This function calculates the angle as a function of the entrance mach"""
        self.v_M4 = ((self.gamma+1)/(self.gamma-1))**0.5*np.arctan((((self.gamma-1)/(self.gamma+1))*(M**2-1))**0.5)-np.arctan((M**2-1)**0.5)

    
    def calc_v_M2(self, M_2guess): # this is an angle value
        """This function calculates the angle as a function of the mach after the turn"""
        v_M2 = ((self.gamma+1)/(self.gamma-1))**0.5*np.arctan((((self.gamma-1)/(self.gamma+1))*(M_2guess**2-1))**0.5)-np.arctan((M_2guess**2-1)**0.5) 
        return v_M2


    def calc_del_v_M2(self, M_2guess):
        """This function is the derivative of the Prandtl-meyer equation"""
        del_v_M2 = (1/M_2guess)*(((M_2guess**2-1)**0.5)/(1+((self.gamma-1)/2)*(M_2guess**2)))
        return del_v_M2
    
    
    def calc_funcs_in_newton(self, m_j):
        """this helps the newton solver be more general"""
        del_function = self.calc_del_v_M2(m_j)
        function = (2*self.half_wedge + self.calc_v_M2(m_j) - self.v_M4)/del_function
        
        return function
        
           
    def backward_newton_solver(self, m_j):
        """This function numerically solves for the mach number given the"""
        epsilon = 10 #arbitrary placeholder
        while epsilon >= self.end_epsilon:
            delta = self.calc_funcs_in_newton(m_j) #((2*self.half_wedge + v_m_j - self.v_M4) /del_v_del_Mj) 
            m_j = m_j - delta
            # print("m_j",m_j)
            epsilon = abs(delta)
            # print("epsilon", epsilon)
            self.M_2 = m_j
        return m_j
    

        ### now, calculate the beta angle to do the shock part of the question (section 2 to 3)
    def calc_a_b_c(self, M_1):
        """the a part of the cubic function"""
        self.a = (1 + ((self.gamma-1)/2)*M_1**2)*np.tan(self.half_wedge + self.alpha) #### theta only works again because it's the same for the second station
        self.b = (M_1**2-1)
        self.c = (1 + ((self.gamma+1)/2)*M_1**2)*np.tan(self.half_wedge + self.alpha) #### theta only works again because it's the same for the second station

    
    ### calculate beta explicitly
    def calc_lambda(self):
        """the lambda coefficient used for solving beta"""
        self.lamb = (self.b**2- 3*self.a*self.c)**0.5

    
    def calc_chi(self, M_1):
        """the chi coefficient used for solving beta"""
        self.chi = ((self.b**3 - 9*self.a*(self.a+((self.gamma+1)/4)*M_1**4*np.tan(self.half_wedge + self.alpha)))/(self.lamb**3))
    

    def calc_beta(self):
        """beta value used for calculating M3"""
        delta = 1
        numerator = self.b + 2*self.lamb*np.cos((4*np.pi*delta + np.arccos(self.chi))/(3))
        denominator = 3*self.a
        self.beta = np.arctan(numerator/denominator)

    
    def calc_Mn2(self, M_1):
        """This function is used to calculate properties as they cross normal to the oblique shock"""
        numerator = 1 + ((self.gamma-1)/2)*(M_1*np.sin(self.beta))**2
        denominator = self.gamma*(M_1*np.sin(self.beta))**2 - ((self.gamma-1)/2)
        self.Mn2 = np.sqrt((numerator)/(denominator))
        
    
    def calc_M_2_guess(self):
        """this function calculates M_3 using Mn3"""
        self.M_2_guess = (self.Mn2)/(np.sin(self.beta-(self.half_wedge + self.alpha)))

    
    def calc_M_1(self, M_1):
        """This function uses explicit solution for M to find M_1 using a,b,c using lambda"""
        #### M_1 is a starting guess
        epsilon = 10 # starting value
        while epsilon >= self.end_epsilon:
        # print("HERE")
            self.calc_a_b_c(M_1)
            self.calc_lambda()
            self.calc_chi(M_1)
            self.calc_beta()
            self.calc_Mn2(M_1)
            self.calc_M_2_guess()
            M_1 = M_1 - 0.0001
            epsilon = abs(((self.M_2_guess-self.M_2)))
        return M_1
    

    def calc_Cp_cr(self, m_inf):
        """Helps find critical mach"""
        first = 2/(self.gamma*m_inf**2)
        second_num = 1 + ((self.gamma-1)/2)*m_inf**2
        second_denom = (self.gamma+1)/2
        exponent = self.gamma/(self.gamma-1)
        second_total = (second_num/second_denom)**exponent - 1
        Cp_cr = first*second_total
        return Cp_cr
    

    def calc_Cp_not_crit(self, m_inf):
        """Helps find critical mach"""
        Cp = self.Xfoil_cpo/(np.sqrt(1-m_inf**2))
        return Cp
    

    def get_prandtl_gluaert_derivative(self):
        M = sp.symbols("M")
        gamma = sp.symbols("gamma")
        CPo = sp.symbols("CPo")
        expression = (2/(gamma*M**2))*(((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1) - CPo/((1-M**2)**0.5)
        diff = sp.diff(expression, M)
        return diff
    

    def newton_M_from_CP_crit_and_prandtl(self, M, gamma, CPo):

        epsilon = 10 #arbitrary placeholder
        while epsilon >= 0.0000001:
            func = (2/(gamma*M**2))*(((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1) - CPo/((1-M**2)**0.5) #### might be wrong. 
            diff = -1.0*CPo*M/(1 - M**2)**1.5 + 4*((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1))*(gamma/2 - 1/2)/(M*(gamma - 1)*(M**2*(gamma/2 - 1/2) + 1)) - 4*(((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1)) - 1)/(M**3*gamma)
            delta = func/diff
            M = M - delta
            epsilon = abs((delta))
        return M


    def calc_M_crit_swept(self,M_unswept):
        
        M_crit_swept = M_unswept/(np.sqrt(1 - (np.sin(self.sweep))**2*(np.cos(self.alpha_2))**2))
        return M_crit_swept
    

    def calc_swept_fineness_ratio(self):

        swept_thickness = self.unswept_thickness*np.cos(self.sweep)
        return swept_thickness
    

    def get_Karmen_tsien_derivative(self):
        
        M = sp.symbols("M")
        gamma = sp.symbols("gamma")
        CPo = sp.symbols("CPo")
        expression = (2/(gamma*M**2))*(((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1) - CPo/((1-M**2)**0.5 + (M**2/(1+(1-M**2)**0.5))*(CPo/2))
        diff = sp.diff(expression, M)
        return diff
    

    def get_Laitone_derivative(self):
        
        M = sp.symbols("M")
        gamma = sp.symbols("gamma")
        CPo = sp.symbols("CPo")
        expression = (2/(gamma*M**2))*(((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1) - CPo/((1-M**2)**0.5 + ((M**2*(1+((gamma-1)/2)*M**2))/(1+(1-M**2)**0.5))*(CPo/2))
        diff = sp.diff(expression, M)
        return diff
    

    def newton_M_from_CP_crit_and_karmen(self, M, gamma, CPo):

        epsilon = 10 #arbitrary placeholder
        while epsilon >= 0.0000001:
            karmen_func = (2/(gamma*M**2))*(((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1) - CPo/((1-M**2)**0.5 + (M**2/(1+(1-M**2)**0.5))*(CPo/2)) #### might be wrong. 
            karmen_diff = -CPo*(-0.5*CPo*M**3/((1 - M**2)**0.5*((1 - M**2)**0.5 + 1)**2) - CPo*M/((1 - M**2)**0.5 + 1) + 1.0*M/(1 - M**2)**0.5)/(CPo*M**2/(2*((1 - M**2)**0.5 + 1)) + (1 - M**2)**0.5)**2 + 4*((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1))*(gamma/2 - 1/2)/(M*(gamma - 1)*(M**2*(gamma/2 - 1/2) + 1)) - 4*(((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1)) - 1)/(M**3*gamma)
            delta = karmen_func/karmen_diff
            M = M - delta
            epsilon = abs((delta))
        return M
    

    def newton_M_from_CP_crit_and_laitone(self, M, gamma, CPo):

        epsilon = 10 #arbitrary placeholder
        while epsilon >= 0.0000001:
            laitone_func = (2/(gamma*M**2))*(((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1) - CPo/((1-M**2)**0.5 + ((M**2*(1+((gamma-1)/2)*M**2))/(1+(1-M**2)**0.5))*(CPo/2)) #### might be wrong. 
            laitone_diff = -CPo*(-0.5*CPo*M**3*(M**2*(gamma/2 - 1/2) + 1)/((1 - M**2)**0.5*((1 - M**2)**0.5 + 1)**2) - CPo*M**3*(gamma/2 - 1/2)/((1 - M**2)**0.5 + 1) - CPo*M*(M**2*(gamma/2 - 1/2) + 1)/((1 - M**2)**0.5 + 1) + 1.0*M/(1 - M**2)**0.5)/(CPo*M**2*(M**2*(gamma/2 - 1/2) + 1)/(2*((1 - M**2)**0.5 + 1)) + (1 - M**2)**0.5)**2 + 4*((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1))*(gamma/2 - 1/2)/(M*(gamma - 1)*(M**2*(gamma/2 - 1/2) + 1)) - 4*(((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1)) - 1)/(M**3*gamma)
            delta = laitone_func/laitone_diff
            M = M - delta
            epsilon = abs((delta))
        return M


    def plot_M_at_diff_angles(self, M_glauert, M_karmen, M_laitone):
        
        glauert_graph = np.zeros(len(self.CPos))
        karmen_graph = np.zeros(len(self.CPos))
        laitone_graph = np.zeros(len(self.CPos))

        for i in range(len(self.CPos)):
            glauert_graph[i] = self.newton_M_from_CP_crit_and_prandtl(M_glauert,self.gamma, self.CPos[i])

        for j in range(len(self.CPos)):
            karmen_graph[j] = self.newton_M_from_CP_crit_and_karmen(M_karmen, self.gamma, self.CPos[j])

        for k in range(len(self.CPos)):
            laitone_graph[k] = self.newton_M_from_CP_crit_and_laitone(M_laitone, self.gamma, self.CPos[k])
        
        print("glauert_correction critical Machs:\n", glauert_graph)
        print("karmen_correction critical Machs:\n", karmen_graph)
        print("laitone_correction critical Machs:\n", laitone_graph)

        plt.plot(self.alphas, glauert_graph, label = "Prandtl-Glauert Correction")
        plt.plot(self.alphas, karmen_graph, label = "Karman_Tsien Correction")
        plt.plot(self.alphas, laitone_graph, label = "Laitone Correction")
        plt.xlabel("alpha[deg]")
        plt.ylabel("Mcrit")
        plt.legend()
        plt.show()


    def run_solver(self):
        """When called, this function runs all the needed functions to complete HW 8"""
        self.calc_machs()
        for i in range(len(self.machs)):
            print("\n")
            print(" P_o_2:", self.P_o_2[i])
            print(" P_inf:", self.P_inf[i])
            print(" M_4:", self.machs[i])

        self.calc_v_M4(self.machs[0])
        M_2 = self.backward_newton_solver(self.machs[0])
        print(" M_2: ", M_2, "\n")
        M_1 = self.calc_M_1(self.M1_start_guess)
        print(" M1:", M_1, "\n")
        prandtl_glauert_derivative = self.get_prandtl_gluaert_derivative()
        # print("prandtl_gluaert_derivative", prandtl_glauert_derivative, "\n")
        M_crit = self.newton_M_from_CP_crit_and_prandtl(self.M_crit_start_guess, self.gamma, self.Xfoil_cpo)
        print("M_crit:", M_crit)
        M_crit_unswept = self.calc_M_crit_swept(M_crit)
        print("Swept M_crit", M_crit_unswept)
        swept_fineness_ratio = self.calc_swept_fineness_ratio()
        print("swept fineness ratio:", swept_fineness_ratio)
        print("unswept fineness ratio:", self.unswept_thickness)
        karmen_tsien_derivative = self.get_Karmen_tsien_derivative()
        # print("karmen_tsien_derivative:", karmen_tsien_derivative)
        laitone_derivative = self.get_Laitone_derivative()
        # print("laitone_derivative:", laitone_derivative)
        # Mcrit_prandtl = self.newton_M_from_CP_crit_and_prandtl_2(self.M_crit_start_guess, self.gamma, self.Xfoil_cpo)
        self.plot_M_at_diff_angles(0.5,0.5,0.5) #0.5 is a starting guess
        

if __name__ == "__main__":
    """This is where the code actually runs"""
    pitot = pitot_rayleigh("hw8.json") # this line initializes the nozzle class 
    run = pitot.run_solver()

