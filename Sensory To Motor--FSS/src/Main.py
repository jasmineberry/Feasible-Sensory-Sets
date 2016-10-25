"""
 *
 * @authors  
 * Source code originated by J. Berry & R. Ritter(jasminab@usc.edu, ritterr@usc.edu)
 * University of Southern California, BME 504 
 *
 * Starting Program.......
 * Demonstrating sensory changes as a varying of muscle fiber length changes, fiber velocities, and force
 * for 2 link, 2-DOF, (6 muscle) tendon-driven muscle arm model in various arm posture tasks: 
 * 1) circular hand motion 
 * 2) 
 * 
"""

import numpy as np 
from scipy.io import loadmat
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import SpindleModule
 
class Main:  
    
    l1 = .35 # length of first arm, meters, 13.7795 in 
    l2 = .27 # length of second arm, meters, 10.6299 in 
    q1range = 130 * (np.pi/180) # limited range of movement, in radians, 130 deg range
    q2range = 150 * (np.pi/180) # 150 deg

    def __init__(self):
        
        print ("Starting simulation...")
        #self.createCircTrajectory()
        self.allPossPoints()
        #self.superimposedCircPoints()
        x, y, q1, q2 = self.geometricLoc(Main.l1, Main.l2, Main.q1range, Main.q2range)
        
        """
        file = open("valuesDUMP.txt", "w")
        for i in range(0, len(x)):
            print (str(x[i])+"," ) #, y[i])
            file.write(str(x[i])+",\n" )


        file.close()
        """
        
        print ("Sim complete!")

    def createCircTrajectory(self): #creating simple circle path; modify circle parameters

        ChangeInTime = 0.001
        
        r = 11 # circle radius
        x = r * np.cos(np.pi/6 + (np.pi/2) * np.arange(0, 2*np.pi, ChangeInTime))
        y = r * np.sin(np.pi/6 + (np.pi/2) * np.arange(0, 2*np.pi, ChangeInTime)) + 13
        
        # create the figure
        plt.figure()
        plt.title("Circular Path Trajectory")
        plt.plot(x, y)
        plt.xlabel('x position(m)')
        plt.ylabel('y position (m)')
        plt.grid()
    

    def allPossPoints(self): #defines all possible effector endpoint for a specific range. 
        
        theta1 = np.arange(0, 130 * (np.pi/180), 0.0174533) # all possible theta1 values
        theta2 = np.arange(0, 150 * (np.pi/180), 0.0174533) # all possible theta2 values
        
        [T1, T2] = np.meshgrid(theta1, theta2) # generate a grid of theta1 and theta2 values
        
        #Using Forward Kinematics Formula we find all possible values of x coordinate positions of the limb for theta1 and theta 2
        x = Main.l1 * np.cos(T1) + Main.l2 * np.cos(T1 + T2) # compute x coordinates
        y = Main.l1 * np.sin(T1) + Main.l2 * np.sin(T1 + T2) # compute y coordinates
        
        plt.figure(2)
        plt.plot(x, y, 'r.');
        #axis equal;
        plt.xlabel('x position (m)')
        plt.ylabel('y position (m)')
        plt.title('All possible x-y coordinates for theta1 & theta2 limb degree ranges') 
        plt.grid()
    
    def superimposedCircPoints(self): #creates graph that superimposes a graph of all possible ranges of theta 1&2 on the circle path trajectory
        
        #Circle Path
        ChangeInTime = 0.001
        
        r = .2794 # circle radius in meters, 11 in
        xCirc = r * np.cos(np.pi/6 + (np.pi/2) * np.arange(0, 2*np.pi, ChangeInTime))
        yCirc = r * np.sin(np.pi/6 + (np.pi/2) * np.arange(0, 2*np.pi, ChangeInTime)) + .3302 #.3302m = 13 inches
        
        #Possible endpoint ranges
        theta1 = np.arange(0, 130 * (np.pi/180), 0.1) # all possible theta1 values
        theta2 = np.arange(0, 150 * (np.pi/180), 0.1) # all possible theta2 values

        [T1, T2] = np.meshgrid(theta1, theta2) # generate a grid of theta1 and theta2 values

        x1 = Main.l1 * np.cos(T1)
        y1 = Main.l1 * np.sin(T1)
        #print(np.nd(theta1))
        x2 = []
        y2 = []
        
        
        for i in theta1:
            for j in theta2:

                tempX = Main.l1 * np.cos(i) + Main.l2 * np.cos(i + j)
                tempY = Main.l1 * np.sin(i) + Main.l2 * np.sin(i + j)
                if tempX <= max(xCirc) and (tempX) >= min(xCirc) and (tempY) <= max(yCirc) and (tempY) >= min(yCirc):
                    x2.append(tempX)
                    y2.append(tempY)

        
        #x2 = Main.l1 * np.cos(T1) + Main.l2 * np.cos(T1 + T2)
        #y2 = Main.l1 * np.sin(T1) + Main.l2 * np.sin(T1 + T2)

        #plot        
        plt.figure(3)
        plt.plot(xCirc, yCirc, 'b')
        plt.plot(x1, y1, 'g.') #For link/joint #1
        plt.plot(x2, y2, 'r.') #For link/joint #2

        plt.xlabel('x position (m)')
        plt.ylabel('y position (m)')
        plt.title('All possible x-y coordinates for theta1 & theta2 limb degree ranges') 
        plt.grid()
    
    
    def my_range(self, start, end, step):
        while start <= end:
            yield start
            start += step    
    
    def geometricLoc(self, link1, link2, qRange1, qRange2): #finds all locations of endpoints based on link lengths and join angles; acc. to FNpg.28
        
        x = y = q1 = q2 = [] 
        count = 0
        
        fileX = open("valuesDUMPX.txt", "w")
        fileY = open("valuesDUMPY.txt", "w")
        for i in self.my_range (0, qRange1, 0.0174533): #iterates by one degree
            for j in self.my_range(0,qRange2, 0.0174533):
                tempX = link2 * np.cos(i + j) + link1 * np.cos(i)
                tempY = link2 * np.sin(i + j) + link1 * np.sin(i)
                x.append(tempX)
                y.append(tempY)
                q1.append(i)
                q2.append(j)
                """
                #DATA DUMP of coordinates for plotly graph
                #if (i <= (43 * (np.pi/180)) and j <= (50*(np.pi/180)) ):
                #if (i >= (0 * (np.pi/180)) and i <= (43 *(np.pi/180)) and j >= (51*(np.pi/180))  and j <= (100*(np.pi/180)) ):
                #if (i >= (0 * (np.pi/180)) and i <= (43 *(np.pi/180)) and j >= (101*(np.pi/180))  and j <= (150*(np.pi/180)) ):
                #if (i >= (44 * (np.pi/180)) and i <= (87 *(np.pi/180)) and j >= (0*(np.pi/180))  and j <= (50*(np.pi/180)) ):
                #if (i >= (44 * (np.pi/180)) and i <= (87 *(np.pi/180)) and j >= (51*(np.pi/180))  and j <= (100*(np.pi/180)) ):
                #if (i >= (44 * (np.pi/180)) and i <= (87 *(np.pi/180)) and j >= (101*(np.pi/180))  and j <= (150*(np.pi/180)) ):
                #if (i >= (88 * (np.pi/180)) and i <= (130 *(np.pi/180)) and j >= (0*(np.pi/180))  and j <= (50*(np.pi/180)) ):
                #if (i >= (88 * (np.pi/180)) and i <= (130 *(np.pi/180)) and j >= (51*(np.pi/180))  and j <= (100*(np.pi/180)) ):
                #if (i >= (88 * (np.pi/180)) and i <= (130 *(np.pi/180)) and j >= (101*(np.pi/180))  and j <= (150*(np.pi/180)) ):
                    fileX.write(str(tempX)+", " )
                    fileY.write(str(tempY)+", " )
                    count+=1
                """
                
        fileX.close()  
        fileY.close() 
        return x, y, q1, q2
    


# Moment Arms, muscle-number = [ m1 ,  m2 ,  m3 ,  m4 ,  m5 ,  m6]
# m1 = Deltoid Anterior
# m2 = Deltoid Posterior
# m3 = Bicep brachii rotates forearm and flexes elbow
# m4 = Tricep brachii rotates forearm and flexes elbow 
# m5 = brachialis 
# m6 = anconeus
# q1 = shoulder movement flexion (-) extension (+)
# q2 = elbow flexion and extension 

rSM = [1.9, -0.8, 3.6, 2.1, 0,       0];
rEM = [0,    0,    3.6, 2.1, 1.8, -1.2];



Main()
plt.show()  


