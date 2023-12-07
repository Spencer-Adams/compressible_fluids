import json 
import numpy as np 
import matplotlib.pyplot as plt 


class Double_wedge:
    """This class contains important nozzle attributes that can assist in double wedge analysis"""
    def __init__(self, json_file):
        self.json_file = json_file
        self.load_json() 

    # First, here are some functions that allow us to take the information 
    # from the text file so that we can plot it. 
    def load_json(self):
        """This function loads the json file that references the text file"""
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.standard_atm = self.pull_nodes(input_vals["Standard_Atm"])
            self.Chord = input_vals["chord[m]"] 
            self.half_wedge = np.deg2rad(input_vals["half_wedge"]) # value in input should be given in degrees 
            self.end_epsilon = input_vals["end_epsilon"]
            

    def pull_nodes(self, standard_atm):
        """This function grabs the text files inside the json and turns them into arrays"""
        with open(standard_atm, 'r') as text_handle:
            # Skip the first two lines
            lines = text_handle.readlines()[2:]
            # Process the remaining lines
            points = [list(map(float, line.strip().split())) for line in lines]
        standard_atm = np.array(points)
        return standard_atm


    # HERE IS SOME EXPANSION FAN STUFF!!!
    def exp_fan_calc_mu_1(self, M1):
        """This function calculates the forward mach angle of an expansion fan 
        
        params
        --------
        M1: mach number before the expansion fan"""
        mu_1 = np.arcsin(1/M1)
        return mu_1


    def exp_fan_calc_v_M1(self, M1, gamma): # this is an angle value
        """This function calculates the angle as a function of the entrance mach"""
        v_M1 = ((gamma+1)/(gamma-1))**0.5*np.arctan((((gamma-1)/(gamma+1))*(M1**2-1))**0.5)-np.arctan((M1**2-1)**0.5)
        return v_M1
    

    def exp_fan_calc_v_M2(self, M_2guess, gamma): # this is an angle value
        """This function calculates the angle as a function of the mach after the turn"""
        v_M2 = ((gamma+1)/(gamma-1))**0.5*np.arctan((((gamma-1)/(gamma+1))*(M_2guess**2-1))**0.5)-np.arctan((M_2guess**2-1)**0.5) 
        return v_M2


    def exp_fan_calc_del_v_M2(self, M_2guess, gamma):
        """This function is the derivative of the Prandtl-meyer equation"""
        del_v_M2 = (1/M_2guess)*(((M_2guess**2-1)**0.5)/(1+((gamma-1)/2)*(M_2guess**2)))
        return del_v_M2


    def exp_fan_M_2_newton_solver(self, m_j, M1, theta, gamma):
        """This function numerically solves for the mach number given the mach number before the expansion fan
        
        params
        --------
        m_j: starting guess mach_2
        M1: Mach number going into entrance"""
        v_M1 = self.exp_fan_calc_v_M1(M1, gamma)
        epsilon = 10 #arbitrary placeholder
        while epsilon >= self.end_epsilon:
            v_m_j = self.exp_fan_calc_v_M2(m_j, gamma)
            del_v_del_Mj = self.exp_fan_calc_del_v_M2(m_j, gamma)
            delta = ((theta + v_M1 - v_m_j) /del_v_del_Mj) 
            m_j = m_j + delta
            epsilon = abs(delta)
        self.M_2 = m_j
        return m_j
    
    
    def exp_fan_calc_mu_next(self, Mnext):
        """mu angle after the expansion fan 
        
        params
        --------
        m_j: starting guess mach_2
        Mnext: Mach number after expansion fan"""
        mu_next = np.arcsin(1/Mnext)
        return mu_next

    
    def exp_fan_calc_exit_mach_wave_angle(self, theta_in, mu_next):
        """This function calculates the exit mach wave angle
        
        params
        --------
        theta_in: wedge angle going into the expansion fan
        mu_next mu angle after expansion fan"""
        exit_mach_wave_angle = mu_next - theta_in
        return exit_mach_wave_angle


    def calc_P_o(self, P, M, gamma):
        """Solves the stagnation pressure given the static pressure and mach number at a given point"""
        P_o = P*((1 + ((gamma-1)/2)*M**2))**(gamma/(gamma-1))
        return P_o


    def calc_T_o(self, T,M, gamma):
        """Solves the stagnation temperature given the static temperature and mach number at a given point"""
        T_o = T*((1 + ((gamma-1)/2)*M**2))
        return T_o

    
    def exp_fan_calc_P_stat_next(self, P_o_next, M_next, gamma):
        """Calculates the static pressure based on the stagnation pressure and mach number after the expansion fan
        
        params
        --------
        P_o_next: stagnation pressure after expansion fan
        M_next: Mach number after expansion fan"""
        P_stat_next = P_o_next/(((1 + ((gamma-1)/2)*M_next**2))**(gamma/(gamma-1)))
        return P_stat_next
    

    def exp_fan_calc_T_stat_next(self, T_o_next, M_next, gamma):
        """Calculates the static pressure based on the stagnation pressure and mach number after the expansion fan
        
        params
        --------
        T_o_next: stagnation temperature after expansion fan
        M_next: Mach number after expansion fan"""
        T_stat_next = T_o_next/((1 + ((gamma-1)/2)*M_next**2))
        return T_stat_next
    

    def exp_fan_find_vals_after_fan(self, Mbefore, theta_in, Pstat_before, Tstat_before, gamma):
        """This function provides values for pressure, mach number, and temperature after flow passes through an expansion fan
        
        params
        --------
        Mbefore: Mach number before the expansion fan
        theta_in: wedge angle going into the expansion fan
        Pstat_before: static pressure before the expansion fan
        Tstat_before: static temperature before the expansion fan
        """
        mu_before = self.exp_fan_calc_mu_1(Mbefore)
        Mafter = self.exp_fan_M_2_newton_solver(Mbefore, Mbefore, theta_in, gamma)
        mu_after = self.exp_fan_calc_mu_next(Mafter)
        wave_angle_after = self.exp_fan_calc_exit_mach_wave_angle(theta_in, mu_after)
        Po_before = self.calc_P_o(Pstat_before, Mbefore, gamma)
        To_before = self.calc_T_o(Tstat_before, Mbefore, gamma)
        To_after = To_before # expansion fans are isentropic
        Po_after = Po_before # expansion fans are isentropic
        Tstat_after = self.exp_fan_calc_P_stat_next(To_after, Mafter, gamma)
        Pstat_after = self.exp_fan_calc_P_stat_next(Po_after, Mafter, gamma)
        value = np.array([Mafter,  Pstat_after, Tstat_after])
        return value
    

    ### now, calculate the beta angle to do the shock part of the question (section 2 to 3)
    def oblique_calc_a_b_c(self, Mbefore, theta, gamma):
        """the a part of the cubic function
        
        params
        --------
        Mbefore: mach before the oblique shockwave
        theta: total wedge ange"""
        a = (1 + ((gamma-1)/2)*Mbefore**2)*np.tan(theta) #### theta only works again because it's the same for the second station
        b = (Mbefore**2-1)
        c = (1 + ((gamma+1)/2)*Mbefore**2)*np.tan(theta) #### theta only works again because it's the same for the second station
        abc = np.array([a, b, c])
        return abc
    

    ### calculate beta explicitly
    def oblique_calc_lambda(self, Mbefore, theta, gamma):
        """the lambda coefficient used for solving beta
        
        params
        --------
        Mbefore: mach before the oblique shockwave
        theta: total wedge ange"""

        abc = self.oblique_calc_a_b_c(Mbefore, theta, gamma)
        lamb = (abc[1]**2- 3*abc[0]*abc[2])**0.5
        return lamb

    
    def oblique_calc_chi(self, Mbefore, theta, lamb, gamma):
        """the chi coefficient used for solving beta
        
        params
        --------
        Mbefore: mach before the oblique shockwave
        theta: total wedge ange
        lamb: lambda coefficient"""
        abc = self.oblique_calc_a_b_c(Mbefore, theta, gamma)
        chi = (abc[1]**3 - 9*abc[0]*(abc[0]+((gamma+1)/4)*Mbefore**4*np.tan(theta)))/(lamb**3)
        return chi
    

    def oblique_calc_beta(self, Mbefore, theta, lamb, chi, gamma):
        """beta used for calculating M after oblique shock
        
        params
        --------
        Mbefore: mach before the oblique shockwave
        theta: total wedge ange
        chi: chi coefficient
        lamb: lambda coefficient"""

        if theta == 0.0:
            beta = np.arcsin(1/Mbefore)
        else:
            abc = self.oblique_calc_a_b_c(Mbefore, theta, gamma)
            delta = 1
            numerator = abc[1] + 2*lamb*np.cos((4*np.pi*delta + np.arccos(chi))/(3))
            denominator = 3*abc[0]
            beta = np.arctan(numerator/denominator)
        return beta


    ## now, calculate Mn2 
    def oblique_calc_Mn_before(self, Mbefore, beta):
        """This function is used to calculate properties as they cross normal to the oblique shock
        
        params
        --------
        Mbefore: mach before the oblique shockwave
        beta: shock angle"""
        Mnbefore = Mbefore*np.sin(beta)
        return Mnbefore

    
    def oblique_calc_Mn_after(self, Mbefore, beta, gamma):
        """This function is used to calculate properties as they cross normal to the oblique shock
        
        params
        --------
        Mbefore: mach before the oblique shockwave
        beta: shock angle"""
        numerator = 1 + ((gamma-1)/2)*(Mbefore*np.sin(beta))**2
        denominator = gamma*(Mbefore*np.sin(beta))**2 - ((gamma-1)/2)
        Mnafter = np.sqrt((numerator)/(denominator))
        return Mnafter
        
    
    def oblique_calc_M_after(self, Mnafter, beta, theta):
        """this function calculates Mafter using Mnafter
        
        params
        --------
        Mnafter: normal Mach number after the oblique shock
        beta: shock angle"""
        Mafter = (Mnafter)/(np.sin(beta-theta))
        return Mafter

    
    def oblique_calc_Pstat_after(self, Pstat_before, Mnbefore, gamma):
        """this function calculates P_stat_after using Mnbefore
        
        params
        --------
        Pstat_before: static pressure before oblique shock
        Mnbefore: mach number normal to shock before shock"""
        Pstat_after = Pstat_before*(1+(2*gamma/(gamma+1))*(Mnbefore**2-1))
        return Pstat_after


    def oblique_calc_Tstat_after(self, Tstat_before, Mnbefore, gamma):
        """this function calculates P_stat_after using Mnbefore
        
        params
        --------
        Tstat_before: static temp before oblique shock
        Mnbefore: mach number normal to shock before shock"""
        Tstat_after = Tstat_before*(1+(2*gamma/(gamma+1))*(Mnbefore**2-1))*((2+(gamma-1)*Mnbefore**2)/((gamma+1)*Mnbefore**2))
        return Tstat_after
    

    def oblique_calc_Po_after(self, Pstat_after, Mafter, gamma):
        """this function calculates P_stag_after using M after and Pstat_after
        
        params
        --------
        Pstat_after: static pressure after oblique shock
        Mafter: mach number after shock"""
        Po_after = Pstat_after*((1 + ((gamma-1)/2)*Mafter**2))**(gamma/(gamma-1))


    def oblique_calc_To_after(self, Tstat_after, Mafter, gamma):
        """this function calculates T_stag_after using Mafter and Tstat_after
        
        params
        --------
        Pstat_after: static pressure after oblique shock
        Mafter: mach number after shock"""
        To_after = Tstat_after*((1 + ((gamma-1)/2)*Mafter**2))
        return To_after
    

    def oblique_find_vals_after_oblique(self, Mbefore, theta, Pstat_before, Tstat_before, gamma):

        lamb = self.oblique_calc_lambda(Mbefore, theta, gamma)
        chi = self.oblique_calc_chi(Mbefore, theta, lamb, gamma)
        beta = self.oblique_calc_beta(Mbefore, theta, lamb, chi, gamma)
        Mn_before = self.oblique_calc_Mn_before(Mbefore, beta)
        Mn_after = self.oblique_calc_Mn_after(Mbefore, beta, gamma)
        Mafter = self.oblique_calc_M_after(Mn_after, beta, theta)
        Pstat_after = self.oblique_calc_Pstat_after(Pstat_before, Mn_before, gamma)
        Tstat_after = self.oblique_calc_Tstat_after(Tstat_before, Mn_before, gamma)
        Po_after = self.oblique_calc_Po_after(Pstat_after,Mafter, gamma)
        To_after = self.oblique_calc_To_after(Tstat_after, Mafter, gamma)
        value = np.array([Mafter, Pstat_after, Tstat_after])
        return value
    

    def calc_CD_inviscid(self, half_wedge, alpha, M_inf, P1, P2, P4, P3, P5, gamma):
        
        coeff = (1/(2*np.cos(gamma)))*(1/((gamma/2)*M_inf**2))
        first = ((P2/P1)-(P5/P1))*np.sin(half_wedge+alpha)
        second = ((P3/P1)-(P4/P1))*np.sin(half_wedge-alpha)
        CD = coeff*(first + second)
        return CD
    

    def calc_CL_inviscid(self, half_wedge, alpha, M_inf, P1, P2, P4, P3, P5, gamma):

        coeff = (1/(2*np.cos(gamma)))*(1/((gamma/2)*M_inf**2))
        first = ((P2/P1)-(P5/P1))*np.cos(half_wedge+alpha)
        second = ((P4/P1)-(P3/P1))*np.cos(half_wedge-alpha)
        CL = coeff*(first + second)
        return CL
    

    def calc_reynolds(self, rho, v, L, mu, alpha, half_wedge):

        Rey = rho*v*L/mu
        Re = (Rey*np.cos(alpha))/np.cos(half_wedge)
        return Re


    def run_double_wedge_solver(self, alpha, Minf, half_wedge, Pstat_before, Tstat_before, gamma, type, Re):
        """This function runs all of the necessary functions to run the isentropic solver"""
        if (half_wedge) - (alpha) < 0:
                theta_real = (half_wedge) - (alpha)
                vals3 = self.exp_fan_find_vals_after_fan(Minf, -theta_real, Pstat_before, Tstat_before, gamma)
            
        elif (half_wedge) - (alpha) >= 0:
            theta_real = (half_wedge) - (alpha)
            vals3 = self.oblique_find_vals_after_oblique(Minf, theta_real, Pstat_before, Tstat_before, gamma)
            
        if (half_wedge) + (alpha) < 0:
            theta_real = (half_wedge) + (alpha)
            vals2 = self.exp_fan_find_vals_after_fan(Minf, -theta_real, Pstat_before, Tstat_before, gamma)
            
        elif (half_wedge) + (alpha) >= 0:
            theta_real = (half_wedge) + (alpha)
            vals2 = self.oblique_find_vals_after_oblique(Minf, theta_real, Pstat_before, Tstat_before, gamma)

        second_theta_real = (2*half_wedge)
        vals5 = self.exp_fan_find_vals_after_fan(vals3[0], second_theta_real, vals3[1], vals3[2], gamma)
        vals4 = self.exp_fan_find_vals_after_fan(vals2[0], second_theta_real, vals2[1], vals2[2], gamma)

        pressures_1_2_4_3_5 = np.array([Pstat_before, vals2[1], vals4[1], vals3[1], vals5[1]])


        inviscid_CD_wing = self.calc_CD_inviscid(half_wedge, alpha, Minf, pressures_1_2_4_3_5[0], pressures_1_2_4_3_5[1], pressures_1_2_4_3_5[2], pressures_1_2_4_3_5[3], pressures_1_2_4_3_5[4], gamma)
            # print("\ninviscid_CD:\n", inviscid_CD_wing)
        inviscid_CL_wing = self.calc_CL_inviscid(half_wedge, alpha, Minf, pressures_1_2_4_3_5[0], pressures_1_2_4_3_5[1], pressures_1_2_4_3_5[2], pressures_1_2_4_3_5[3], pressures_1_2_4_3_5[4], gamma)
        

        if type == "inviscid":
            L_over_D = inviscid_CL_wing/inviscid_CD_wing

        elif type == "turbulent":
            # Cdfric = 7/(225*Re**(1/7))
            Cdfric = 14/(225*np.cos(self.half_wedge)*(Re**(1/7)))
            Tavg = Tstat_before*(1+(2/9)*(((gamma - 1)/2)*Minf**2))
            correction = 1/((Tstat_before/Tavg)**(5/2)*((Tavg+120)/(Tstat_before+120)))**(1/7)
            L_over_D = inviscid_CL_wing/(inviscid_CD_wing + Cdfric*correction)

        elif type == "transitional":
            # Cdfric = 7/(225*Re**(1/7))
            Cdfric = 14/(225*np.cos(self.half_wedge)*(Re**(1/7)))
            Tavg = Tstat_before*(1+(2/9)*(((gamma - 1)/2)*Minf**2))
            correction = 1/((Tstat_before/Tavg)**(5/2)*((Tavg+120)/(Tstat_before+120)))**(1/7)
            L_over_D = inviscid_CL_wing/(inviscid_CD_wing + Cdfric*correction)

        elif type == "laminar":
            Cdfric = 1.328/np.sqrt(Re)
            Tavg = Tstat_before*(1+(1-(8/15))*((gamma-1)/2)*Minf**2)
            correction = 1/np.sqrt((Tstat_before/Tavg)**(5/2)/((Tstat_before + 120)/(Tavg + 120)))
            L_over_D = inviscid_CL_wing/(inviscid_CD_wing + Cdfric*correction)


        return L_over_D
    

if __name__ == "__main__":
    """This is where the code actually runs"""
    # alpha = np.radians(1.0)
    wedge = Double_wedge("Final_project.json") # this line initializes the nozzle class
    
    ## Problem 2, part a
    Machs_part_2_a = np.array([1.5,2,4,8]) 
    alphas = np.radians(np.linspace(0,15,45))
    L_over_D_alpha = []
    alpha = []
    gamma_prob_2_pt_a = ((wedge.standard_atm[0][5]**2)*(wedge.standard_atm[0][3]))/(wedge.standard_atm[0][1])
    # print(gamma_prob_2_pt_a)
    for i in range(len(Machs_part_2_a)):
        L_over_D_alpha = []
        alpha = []
        for j in range(len(alphas)):
            vinf = Machs_part_2_a[i]*wedge.standard_atm[0][5]
            Reynolds = wedge.calc_reynolds(wedge.standard_atm[0][3],vinf, wedge.Chord, wedge.standard_atm[0][4]*10**-3, alphas[j], wedge.half_wedge)
            L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[0][1], wedge.standard_atm[0][2], gamma_prob_2_pt_a, "inviscid", Reynolds)) 
            alpha.append(np.degrees(alphas[j]))
        plt.plot(alpha, L_over_D_alpha, label = Machs_part_2_a[i])
    L_over_D_final = np.array(L_over_D_alpha)
    alpha_final = np.array(alpha)

    # wedge.plotter(np.degrees(alpha_final), L_over_D_final)
    plt.xlabel("alpha[deg]")
    plt.ylabel("L/D")
    plt.title("P2 pt.a) Inviscid L/D w.r.t alpha")
    plt.legend()
    plt.figure()


    ## Problem 2, part b
    altitudes = [10, 20, 30, 40, 50]
    Machs_part_2_a = np.linspace(1,10,20)
    alphas = np.radians(np.linspace(0,9,30))
    L_over_D_alpha = []
    L_over_D_max = []
    alpha = []
    gamma_prob_2_pt_a = ((wedge.standard_atm[0][5]**2)*(wedge.standard_atm[0][3]))/(wedge.standard_atm[0][1])
    # print(gamma_prob_2_pt_a)
    for k in range(len(altitudes)):
        index = np.where(wedge.standard_atm[:,0] == altitudes[k])[0][0]
        gamma_prob_2_pt_a = ((wedge.standard_atm[index][5]**2)*(wedge.standard_atm[index][3]))/(wedge.standard_atm[index][1])
        L_over_D_max = []
        for i in range(len(Machs_part_2_a)):
            L_over_D_alpha = []
            alpha = []
            for j in range(len(alphas)):
                vinf = Machs_part_2_a[i]*wedge.standard_atm[index][5]
                Reynolds = wedge.calc_reynolds(wedge.standard_atm[index][3],vinf, wedge.Chord, wedge.standard_atm[index][4]*10**-3, alphas[j], wedge.half_wedge)
                L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a, "inviscid", Reynolds)) 
                alpha.append(np.degrees(alphas[j]))
            L_over_D_max.append(max(L_over_D_alpha))
        # print(wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a)
        plt.plot(Machs_part_2_a, L_over_D_max, label = str(altitudes[k]) + ' km')
    plt.xlabel("machs")
    plt.ylabel("max L/D")
    plt.ylim(0, 15)
    # plt.gca().set_aspect('equal')
    plt.title("P2 pt.b) Inviscid L/D Max w.r.t mach number")
    plt.legend()
    plt.figure()



    ## Problem 3, part a
    ### Here is the turbulent part 
    Machs_part_2_a = np.array([2,5,10]) 
    alphas = np.radians(np.linspace(0,15,45))
    L_over_D_alpha = []
    alpha = []
    gamma_prob_2_pt_a = ((wedge.standard_atm[0][5]**2)*(wedge.standard_atm[0][3]))/(wedge.standard_atm[0][1])
    # print(gamma_prob_2_pt_a)
    for i in range(len(Machs_part_2_a)):
        L_over_D_alpha = []
        alpha = []
        for j in range(len(alphas)):
            vinf = Machs_part_2_a[i]*wedge.standard_atm[0][5]
            Reynolds = wedge.calc_reynolds(wedge.standard_atm[0][3],vinf, wedge.Chord, wedge.standard_atm[0][4]*10**-3,alphas[j], wedge.half_wedge)
            L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[0][1], wedge.standard_atm[0][2], gamma_prob_2_pt_a, "turbulent", Reynolds)) 
            alpha.append(np.degrees(alphas[j]))
        plt.plot(alpha, L_over_D_alpha, label = "Viscous " + str(Machs_part_2_a[i]))
    L_over_D_final = np.array(L_over_D_alpha)
    alpha_final = np.array(alpha)


    # Here is the inviscid part for comparison
    Machs_part_2_a = np.array([2,5,10]) 
    alphas = np.radians(np.linspace(0,15,45))
    L_over_D_alpha = []
    alpha = []
    gamma_prob_2_pt_a = ((wedge.standard_atm[0][5]**2)*(wedge.standard_atm[0][3]))/(wedge.standard_atm[0][1])
    # print(gamma_prob_2_pt_a)
    for i in range(len(Machs_part_2_a)):
        L_over_D_alpha = []
        alpha = []
        for j in range(len(alphas)):
            vinf = Machs_part_2_a[i]*wedge.standard_atm[0][5]
            Reynolds = wedge.calc_reynolds(wedge.standard_atm[0][3],vinf, wedge.Chord, wedge.standard_atm[0][4]*10**-3, alphas[j], wedge.half_wedge)
            L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[0][1], wedge.standard_atm[0][2], gamma_prob_2_pt_a, "inviscid", Reynolds)) 
            alpha.append(np.degrees(alphas[j]))
        plt.plot(alpha, L_over_D_alpha, label = "Inviscid " + str(Machs_part_2_a[i]))
    L_over_D_final = np.array(L_over_D_alpha)
    alpha_final = np.array(alpha)

    # wedge.plotter(np.degrees(alpha_final), L_over_D_final)
    plt.xlabel("alpha[deg]")
    plt.ylabel("L/D")
    plt.title("P3 pt.a) Inviscid Compared to Viscous L/D w.r.t alpha")
    plt.legend()
    plt.figure()


    ## Problem 3, part b
    altitudes = [10, 20, 30, 40, 50]
    Machs_part_2_a = np.linspace(1,10,20)
    alphas = np.radians(np.linspace(0,9,30))
    L_over_D_alpha = []
    L_over_D_max = []
    alpha = []
    gamma_prob_2_pt_a = ((wedge.standard_atm[0][5]**2)*(wedge.standard_atm[0][3]))/(wedge.standard_atm[0][1])
    # print(gamma_prob_2_pt_a)
    for k in range(len(altitudes)):
        index = np.where(wedge.standard_atm[:,0] == altitudes[k])[0][0]
        gamma_prob_2_pt_a = ((wedge.standard_atm[index][5]**2)*(wedge.standard_atm[index][3]))/(wedge.standard_atm[index][1])
        L_over_D_max = []
        for i in range(len(Machs_part_2_a)):
            L_over_D_alpha = []
            alpha = []
            for j in range(len(alphas)):
                vinf = Machs_part_2_a[i]*wedge.standard_atm[index][5]
                Reynolds = wedge.calc_reynolds(wedge.standard_atm[index][3],vinf, wedge.Chord, wedge.standard_atm[index][4]*10**-3, alphas[j], wedge.half_wedge)
                L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a, "turbulent", Reynolds)) 
                alpha.append(np.degrees(alphas[j]))
            L_over_D_max.append(max(L_over_D_alpha))
        # print(wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a)
        plt.plot(Machs_part_2_a, L_over_D_max, label = str(altitudes[k]) + ' km')
    plt.xlabel("machs")
    plt.ylabel("max L/D")
    plt.ylim(0, 15)
    # plt.gca().set_aspect('equal')
    plt.title("P3 pt.b) Viscous L/D Max w.r.t mach number")
    plt.legend()
    plt.figure()


    ## Problem 4, part a
    altitudes = [10, 20, 30, 40, 50]
    Machs_part_2_a = np.linspace(1.25,10,20)
    alphas = np.radians(np.linspace(0,9,30))
    L_over_D_alpha = []
    L_over_D_max = []
    alpha = []
    gamma_prob_2_pt_a = ((wedge.standard_atm[0][5]**2)*(wedge.standard_atm[0][3]))/(wedge.standard_atm[0][1])
    # print(gamma_prob_2_pt_a)
    for k in range(len(altitudes)):
        index = np.where(wedge.standard_atm[:,0] == altitudes[k])[0][0]
        gamma_prob_2_pt_a = ((wedge.standard_atm[index][5]**2)*(wedge.standard_atm[index][3]))/(wedge.standard_atm[index][1])
        L_over_D_max = []
        for i in range(len(Machs_part_2_a)):
            L_over_D_alpha = []
            alpha = []
            for j in range(len(alphas)):
                vinf = Machs_part_2_a[i]*wedge.standard_atm[index][5]
                Reynolds = wedge.calc_reynolds(wedge.standard_atm[index][3],vinf, wedge.Chord, wedge.standard_atm[index][4]*10**-3, alphas[j], wedge.half_wedge)
                if Reynolds <= 500000:
                    L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a, "laminar", Reynolds))
                else:
                    L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a, "turbulent", Reynolds))
                alpha.append(np.degrees(alphas[j]))
            L_over_D_max.append(max(L_over_D_alpha))
        # print(wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a)
        plt.plot(Machs_part_2_a, L_over_D_max, label = str(altitudes[k]) + ' km')
    plt.xlabel("machs")
    plt.ylabel("max L/D")
    plt.ylim(0, 15)
    plt.title("P4 pt.a) L/D Max accounting for flow regime w.r.t mach number")
    plt.legend()
    plt.figure()


    ## Problem 4, part b
    altitudes = [10, 20, 30, 40, 50]
    Machs_part_2_a = np.linspace(1,10,20)
    alphas = np.radians(np.linspace(0,9,30))
    L_over_D_alpha = []
    L_over_D_max = []
    Reynolds_list = np.zeros(len(Machs_part_2_a))
    Reynolds_alpha_list = []
    alpha = []
    gamma_prob_2_pt_a = ((wedge.standard_atm[0][5]**2)*(wedge.standard_atm[0][3]))/(wedge.standard_atm[0][1])
    # print(gamma_prob_2_pt_a)
    for k in range(len(altitudes)):
        index = np.where(wedge.standard_atm[:,0] == altitudes[k])[0][0]
        gamma_prob_2_pt_a = ((wedge.standard_atm[index][5]**2)*(wedge.standard_atm[index][3]))/(wedge.standard_atm[index][1])
        L_over_D_max = []
        for i in range(len(Machs_part_2_a)):
            L_over_D_alpha = []
            Reynolds_alpha_list = []
            alpha = []
            # Reynolds_list = []
            for j in range(len(alphas)):
                vinf = Machs_part_2_a[i]*wedge.standard_atm[index][5]
                Reynolds = wedge.calc_reynolds(wedge.standard_atm[index][3],vinf, wedge.Chord, wedge.standard_atm[index][4]*10**-3, alphas[j], wedge.half_wedge)
                Reynolds_alpha_list.append(Reynolds)
                if Reynolds <= 500000:
                    L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a, "laminar", Reynolds))
                else:
                    L_over_D_alpha.append(wedge.run_double_wedge_solver(alphas[j], Machs_part_2_a[i], wedge.half_wedge, wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a, "turbulent", Reynolds)) 
                alpha.append(np.degrees(alphas[j]))
            L_over_D_max.append(max(L_over_D_alpha))
            index_max = np.argmax(L_over_D_alpha)
            Reynolds_list[i] = Reynolds_alpha_list[index_max]

        # print(wedge.standard_atm[index][1], wedge.standard_atm[index][2], gamma_prob_2_pt_a)
        plt.plot(Reynolds_list, L_over_D_max, label = str(altitudes[k]) + ' km')
    plt.xlabel("reynolds Number")
    plt.xscale('log')
    plt.ylabel("Max L/D")
    plt.ylim(0, 15)
    # plt.gca().set_aspect('equal')
    plt.title("P4 pt.b) L/D Max w.r.t mach number")
    plt.legend()

    plt.show()


    

