"""
 *
 * @authors  
 * Source code originated by J. Berry & R. Ritter(jasminab@usc.edu, ritterr@usc.edu)
 * University of Southern California, BME 504 
 *
 * Spindle Model incorporating the functionality of the muscle's spindle. 
 * 
"""


import math
import numpy as np

class SpindleModule():
    
    #CLASS VARIABLES
    gamma_dynamic = gamma_static = L = L_dot = L_ddot = 0
    
    
    #CONSTANT VARIABLES
    dt = 1
    K_SR = 10.4649 #sensory region spring constant [FU/L0]
    K_PR = .15 #polar region spring constant [FU/L0]
    M = 0.0002 # intrafusal fiber mass
    
    X = 0.7 #percentage of the secondary afferent on sensory region
    LN_SR = 0.0423 #Sensory region threshold length (L0)
    LN_PR = 0.89 #Polar region threshold length (L0)
    R = 0.46 #Fascicle length below which force production is zero (L0)
    
    L0_SR = 0.04 #sensory region rest length (L0)
    L0_PR = 0.76 #polar region rest length (L0)
    L_secondary = 0.04 #secondary afferent rest length (L0)
    
    f_dynamic = 0
    f_static = 0
    T_ddot_bag1 = 0
    T_dot_bag1 = 0
    T_bag1 = 0
    T_ddot_bag2 = 0
    T_dot_bag2 = 0
    T_bag2 = 0
    T_ddot_chain = 0
    T_dot_chain = 0
    T_chain = 0

    S = 0.156 #amount of partial occlusion that occurs in primary afferent
    
    def __init__(self):
        print("In the SPINDLE: ")
    
    def chain_func(self, L, L_dot, L_ddot):
        
        freq_chain = 90 #fusimotor frequency to activation
        
        beta0 = 0.0822 #Passive damping coefficient [FU/(L0/s)]
        beta2_chain = -0.069 #Coef. of damping due to static fusimotor input [FU/(L0/s)]
        Gamma2_chain = 0.0954 #Coef. of force generation due to stat. fusimotor input [FU]
        
        if L_dot >= 0:
            C = 1; # polar region (L)engthening | C is constant describing the experimentally observed asymmetric effect of velocity on force production during lengthening and shortening
        else:
            C = 0.4200 #polar region (S)hortening
            
        G = 10000 # (7250) term relating the sensory region's stretch to afferent firing 
        
        #converts fusimotor frequency(gamma_static) to an equivalent activation level (f_static)
        f_static_chain = pow(self.gamma_static,self.p) /pow(self.gamma_static^self.p) + pow(freq_chain,self.p)
    
        beta = beta0 + beta2_chain * f_static_chain
        Gamma = Gamma2_chain * self.f_static
        
        self.T_ddot_chain = (self.K_SR/self.M) * [C * beta * np.sign(L_dot-self.T_dot_chain/self.K_SR) * abs(pow(L_dot-self.T_dot_chain/self.K_SR,self.a)) * (L-self.L0_SR - self.T_chain/self.K_SR -self.R) + self.K_PR *(L - self.L0_SR - self.T_chain/self.K_SR - self.L0_PR) + self.M * L_ddot + Gamma - self.T_chain]
        self.T_dot_chain = self.T_ddot_chain*self.dt/1000 + self.T_dot_chain
        self.T_chain = self.T_dot_chain * self.dt/1000 + self.T_chain
        
        AP_primary_chain =  G * [self.T_dot_chain/self.K_SR - (self.LN_SR-self.L0_SR)]# Afferent Potential 
        AP_secondary_chain = G * (self.X * self.L_secondary/self.L0_SR * [self.T_chain/self.K_SR - (self.LN_SR - self.L0_SR)] + (1-self.X) * self.L_secondary/self.L0_PR * (L - self.T_chain/self.K_SR - self.L0_SR - self.LN_PR) )
        
    def bagONE_func(self, L, L_dot, L_ddot):
        
        tau_bag1 = 0.149
        freq_bag1 = 60
        
        beta0 = 0.0605
        beta1 = 0.2592
        Gamma1 = 0.0289
        
        G = 20000
        
        if L_dot >= 0:
            C = 1
        else:
            C = 0.42
        
        df_dynamic = (pow(self.gamma_dynamic, self.p)/(pow(self.gamma_dynamic, self.p)+pow(freq_bag1, self.p))-self.f_dynamic)/tau_bag1;
        self.f_dynamic = self.dt/1000*df_dynamic + self.f_dynamic;
        
        beta = beta0 + beta1 * self.f_dynamic
        Gamma = Gamma1 * self.f_dynamic
        
        self.T_ddot_bag1 = self.K_SR/self.M * (C * beta * np.sign(L_dot-self.T_dot_bag1/self.K_SR)*abs(pow(L_dot-self.T_dot_bag1/self.K_SR, self.a))*(L-self.L0_SR-self.T_bag1/self.K_SR-self.R)+self.K_PR*(L-self.L0_SR-self.T_bag1/self.K_SR-self.L0_PR)+self.M*L_ddot+Gamma-self.T_bag1);
        self.T_dot_bag1 = self.T_ddot_bag1*self.dt/1000 + self.T_dot_bag1
        self.T_bag1 = self.T_dot_bag1*self.dt/1000 + self.T_bag1
        
        AP_bag1 = G * (self.T_bag1/self.K_SR-(self.LN_SR-self.L0_SR))

        
        
    def bagTWO_func(self, L, L_dot, L_ddot):
        
        tau_bag2 = 0.205
        freq_bag2 = 60
        
        beta0 = 0.0822
        beta2 = -0.046
        Gamma2 = 0.0636


        G = 10000 # (7250)
        
        if L_dot >= 0:
            C = 1 
        else:
            C = 0.42
        
        df_static = (pow(self.gamma_static, self.p)/(pow(self.gamma_static, self.p)+pow(self.freq_bag2, self.p))-self.f_static)/tau_bag2
        self.f_static = self.dt/1000*df_static + self.f_static;
        
        beta = beta0 + beta2 * self.f_static
        Gamma = Gamma2 * self.f_static;
        
        T_ddot_bag2 = self.K_SR/self.M * (C * beta * np.sign(L_dot-self.T_dot_bag2/self.K_SR)*abs(pow(L_dot-self.T_dot_bag2/self.K_SR,self.a))*(L-self.L0_SR-self.T_bag2/self.K_SR-self.R)+self.K_PR*(L-self.L0_SR-self.T_bag2/self.K_SR-self.L0_PR)+self.M*L_ddot+Gamma-self.T_bag2)
        T_dot_bag2 = T_ddot_bag2*self.dt/1000 + self.T_dot_bag2;
        T_bag2 = T_dot_bag2*self.dt/1000 + self.T_bag2;
        
        AP_primary_bag2 = G*(T_bag2/self.K_SR-(self.LN_SR-self.L0_SR));
        AP_secondary_bag2 = G*(self.X*self.L_secondary/self.L0_SR*(T_bag2/self.K_SR-(self.LN_SR-self.L0_SR))+(1-self.X)*self.L_secondary/self.L0_PR*(L-T_bag2/self.K_SR-self.L0_SR-self.LN_PR));

        
# Afferent firing model
S = 0.156; #represents the amount of partial occlusion(blockage or closing of blood vessel) that occurs in primary afferent firing also originated from direct experimental measurements (Fallon)


#AP_bag1 = bag1_model(L,L_dot,L_ddot);
#[AP_primary_bag2,AP_secondary_bag2] = bag2_model(L,L_dot,L_ddot);
#[AP_primary_chain,AP_secondary_chain] = chain_model(L,L_dot,L_ddot);
"""
if AP_bag1 < 0
    AP_bag1 = 0;
end

if AP_primary_bag2 < 0
    AP_primary_bag2 = 0;
end

if AP_primary_chain < 0
    AP_primary_chain = 0;
end


if AP_secondary_bag2 < 0
    AP_secondary_bag2 = 0;
end

if AP_secondary_chain < 0
    AP_secondary_chain = 0;
end


if AP_bag1 > (AP_primary_bag2+AP_primary_chain)
    Larger = AP_bag1;
    Smaller = AP_primary_bag2+AP_primary_chain;
else
    Larger = AP_primary_bag2+AP_primary_chain;
    Smaller = AP_bag1;
end
Output_primary = Larger + S * Smaller;
Output_secondary = AP_secondary_bag2 + AP_secondary_chain;




Primary = Output_primary;
Secondary = Output_secondary;
   
PrimaryAfferent = SecondaryAfferent = 0 #values of the function    
"""