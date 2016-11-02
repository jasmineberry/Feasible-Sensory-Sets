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
import math
import csv
from mpl_toolkits.mplot3d import Axes3D

import SpindleModule
from bokeh.core.properties import AngleSpec
 
class Main:  
    
    l1 = .35 # length of first arm, meters, 13.7795 in 
    l2 = .27 # length of second arm, meters, 10.6299 in 
    q1range = 130 * (np.pi/180) # limited range of movement, in radians, 130 deg range
    q2range = 150 * (np.pi/180) # 150 deg

    def __init__(self):
        
        print ("Starting simulation...")
        
        """Creating Trajectories"""
        AnglesCirc1, AnglesCirc2 = self.createCircTrajectory()
        AnglesLine1a, AnglesLine1b, AnglesLine2a, AnglesLine2b, AnglesLine3a, AnglesLine3b, AnglesLine4a, AnglesLine4b, AnglesLine5a, AnglesLine5b  = self.createLineTrajectory()
        AnglesBF1, AnglesBF2 = self.createBackForthTrajectory()
        #self.allPossPoints()
        #self.superimposedCircPoints()
        x, y = self.geometricLoc(Main.l1, Main.l2, Main.q1range, Main.q2range)
        self.plotConfigSpaceStatic(Main.q1range, Main.q2range)
        
        
        """Establishing Muscle Lengths and Moment Arms"""
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
        
        R = np.array([rSM, rEM])
        
        # Optimal length (tendon) for base ReferenceAngle at 0degs for each joint, in meters
        OL1 = 0.098
        OL2 = 0.137
        OL3 = 0.124
        OL4 = 0.121
        OL5 = 0.086
        OL6 = 0.027
        OptimalLengths = np.matrix([OL1,OL2,OL3,OL4,OL5,OL6])
        
        ReferenceAngleJ1 = ReferenceAngleJ2 = 0 #shoulder, elbow in degrees
        
        
        #Circle Trajectory 
        AnglesCirc1 = np.array(AnglesCirc1)
        AnglesCirc2 = np.array(AnglesCirc2)
        AngChange_Circ1 = AnglesCirc1 - ReferenceAngleJ1 
        AngChange_Circ2 = AnglesCirc2 - ReferenceAngleJ2

        TotalAngleChange_Circ = np.matrix([AngChange_Circ1, AngChange_Circ2])

        #Line Trajectory
        AnglesLine1a = np.array(AnglesLine1a)
        AnglesLine1b = np.array(AnglesLine1b) 
        AnglesLine2a = np.array(AnglesLine2a)
        AnglesLine2b = np.array(AnglesLine2b)
        AnglesLine3a = np.array(AnglesLine3a)
        AnglesLine3b = np.array(AnglesLine3b)
        AnglesLine4a = np.array(AnglesLine4a)
        AnglesLine4b = np.array(AnglesLine4b)
        AnglesLine5a = np.array(AnglesLine5a)
        AnglesLine5b = np.array(AnglesLine5b)
        AngChange_Line1a = AnglesLine1a - ReferenceAngleJ1 
        AngChange_Line1b = AnglesLine1b - ReferenceAngleJ2
        AngChange_Line2a = AnglesLine2a - ReferenceAngleJ1 
        AngChange_Line2b = AnglesLine2b - ReferenceAngleJ2
        AngChange_Line3a = AnglesLine3a - ReferenceAngleJ1 
        AngChange_Line3b = AnglesLine3b - ReferenceAngleJ2
        AngChange_Line4a = AnglesLine4a - ReferenceAngleJ1 
        AngChange_Line4b = AnglesLine4b - ReferenceAngleJ2
        AngChange_Line5a = AnglesLine5a - ReferenceAngleJ1 
        AngChange_Line5b = AnglesLine5b - ReferenceAngleJ2
        
        
        TotalAngleChange_Line1 = np.matrix([AngChange_Line1a, AngChange_Line1b])
        TotalAngleChange_Line2 = np.matrix([AngChange_Line2a, AngChange_Line2b])
        TotalAngleChange_Line3 = np.matrix([AngChange_Line3a, AngChange_Line3b])
        TotalAngleChange_Line4 = np.matrix([AngChange_Line4a, AngChange_Line4b])
        TotalAngleChange_Line5 = np.matrix([AngChange_Line5a, AngChange_Line5b])

        #Back and Forth Trajectory 
        AnglesBF1 = np.array(AnglesBF1)
        AnglesBF2 = np.array(AnglesBF2)
        AngChange_BF1 = AnglesBF1 - ReferenceAngleJ1 
        AngChange_BF2 = AnglesBF2 - ReferenceAngleJ2

        TotalAngleChange_BF = np.matrix([AngChange_BF1, AngChange_BF2])

        
        # Calculate change in tendon length (ChangeInExcursion)
        ChangeInExcursion_Circ = -R.T * TotalAngleChange_Circ #also known as relative muscle length
        
        ChangeInExcursion_Line1 = -R.T * TotalAngleChange_Line1
        ChangeInExcursion_Line2 = -R.T * TotalAngleChange_Line2
        ChangeInExcursion_Line3 = -R.T * TotalAngleChange_Line3
        ChangeInExcursion_Line4 = -R.T * TotalAngleChange_Line4
        ChangeInExcursion_Line5 = -R.T * TotalAngleChange_Line5 
 
        ChangeInExcursion_BF = -R.T * TotalAngleChange_BF
 
        # Calculate muscle length for the different tasks
        MuscleLengths_Circ = ChangeInExcursion_Circ + OptimalLengths.T
        
        MuscleLengths_Line1 = ChangeInExcursion_Line1 + OptimalLengths.T
        MuscleLengths_Line2 = ChangeInExcursion_Line2 + OptimalLengths.T
        MuscleLengths_Line3 = ChangeInExcursion_Line3 + OptimalLengths.T
        MuscleLengths_Line4 = ChangeInExcursion_Line4 + OptimalLengths.T
        MuscleLengths_Line5 = ChangeInExcursion_Line5 + OptimalLengths.T

        MuscleLengths_BF = ChangeInExcursion_BF + OptimalLengths.T
        
        
        #Writing Data to CSV

        self.createCSVFile('data/CircleTraj.csv', MuscleLengths_Circ)
        
        self.createCSVFile('data/Line1Traj.csv', MuscleLengths_Line1)
        self.createCSVFile('data/Line2Traj.csv', MuscleLengths_Line2)
        self.createCSVFile('data/Line3Traj.csv', MuscleLengths_Line3)
        self.createCSVFile('data/Line4Traj.csv', MuscleLengths_Line4)
        self.createCSVFile('data/Line5Traj.csv', MuscleLengths_Line5)
        
        self.createCSVFile('data/BFTraj.csv', MuscleLengths_BF)

        print ("Sim complete!")

    def createCSVFile(self, filename, matrixData):
        
        csv.register_dialect(
            'mydialect',
            delimiter = ',',
            quotechar = '"',
            doublequote = True,
            skipinitialspace = True,
            lineterminator = '\n',
            quoting = csv.QUOTE_MINIMAL)
            
        arrayofdata=[[1,2,4,5,'something','spam',2.334], [3,1,6,3,'anything','spam',0]]
        
        #print())
        #postureNum = list(np.array(range(1,len(matrixData.T)+1)).T)
        
        postureNum = (range(1,len(matrixData.T)+1))
        #matrixData = list(np.array(matrixData.T))
        
        matrixData = np.concatenate((np.array(postureNum,ndmin=2), \
                                    np.array(matrixData,ndmin=2)), axis = 0)
        
        matrixData = list(matrixData.T)
        #matrixData = matrixData.T
        with open(filename, 'w') as mycsvfile:
            thedatawriter = csv.writer(mycsvfile, dialect='mydialect')
            titleLabels = ['Posture Number', 'Deltoid Anterior', 'Deltoid Posterior', 'Bicep', 'Tricep', 'Brachialis', 'Anconeus']
            
            thedatawriter.writerow(titleLabels)
            for row in matrixData:
                thedatawriter.writerow(row)
                
                
    def createCircTrajectory(self): #creating simple circle path; modify circle parameters; plots cart. space, config space, 

        ChangeInStep = 0.001
        
        r = .15 # circle radius in meters
        x = r * np.cos(np.pi/6 + (np.pi/2) * np.arange(0, 2*np.pi, ChangeInStep))
        y = r * np.sin(np.pi/6 + (np.pi/2) * np.arange(0, 2*np.pi, ChangeInStep)) +.4
        
        # create the figures
        plt.figure()
        plt.title("Circular Path Trajectory-Cartesian Space")
        plt.plot(x, y)
        plt.xlim(-0.3, 0.3)
        plt.ylim(0, .6)
        plt.xlabel('x position (m)')
        plt.ylabel('y position (m)')
        plt.grid()
        
        q1angles, q2angles = self.inverseKin(x, y)
        q1anglesP = [a * (180/np.pi) for a in q1angles] #converting from radians to degrees
        q2anglesP = [b * (180/np.pi) for b in q2angles]
        
        plt.figure()
        plt.title("Circular Path Trajectory-Configuration Space")
        plt.plot(q1anglesP, q2anglesP, 'g.') 
        plt.xlabel('q1 in degrees')
        plt.ylabel('q2 in degrees')
        plt.grid()
        
        return (q1angles, q2angles)   

    def createLineTrajectory(self):
        
        ChangeInStep = 0.004 #0.001 creates 400 postures, 0.004 creates 100 postures
        
        x = np.arange(-.2, .2, ChangeInStep)
        y1 = .5 * x +.45
        y2 = .4 * x +.45
        y3 = .3 * x +.45
        y4 = .2 * x +.45
        y5 = .1 * x +.45
        
        # create the figures
        plt.figure()
        plt.title("Line Path Trajectory-Cartesian Space")
        plt.plot(x, y1)
        plt.plot(x, y2)
        plt.plot(x, y3)
        plt.plot(x, y4)
        plt.plot(x, y5)
        plt.xlim(-0.4, 0.4)
        plt.ylim(0, 0.6)
        plt.xlabel('x position(m)')
        plt.ylabel('y position (m)')
        plt.grid()
        
        
        q1anglesLine1, q2anglesLine1 = self.inverseKin(x, y1)
        q1anglesLine2, q2anglesLine2 = self.inverseKin(x, y2)
        q1anglesLine3, q2anglesLine3 = self.inverseKin(x, y3)
        q1anglesLine4, q2anglesLine4 = self.inverseKin(x, y4)
        q1anglesLine5, q2anglesLine5 = self.inverseKin(x, y5)

        q1anglesLineP1 = [a * (180/np.pi) for a in q1anglesLine1] #converting from radians to degrees
        q2anglesLineP1 = [b * (180/np.pi) for b in q2anglesLine1]
        q1anglesLineP2 = [a * (180/np.pi) for a in q1anglesLine2] #converting from radians to degrees
        q2anglesLineP2 = [b * (180/np.pi) for b in q2anglesLine2]
        q1anglesLineP3 = [a * (180/np.pi) for a in q1anglesLine3] #converting from radians to degrees
        q2anglesLineP3 = [b * (180/np.pi) for b in q2anglesLine3]
        q1anglesLineP4 = [a * (180/np.pi) for a in q1anglesLine4] #converting from radians to degrees
        q2anglesLineP4 = [b * (180/np.pi) for b in q2anglesLine4]
        q1anglesLineP5 = [a * (180/np.pi) for a in q1anglesLine5] #converting from radians to degrees
        q2anglesLineP5 = [b * (180/np.pi) for b in q2anglesLine5]
        
        plt.figure()
        plt.title("Line Path Trajectory-Configuration Space")
        plt.plot(q1anglesLineP1, q2anglesLineP1)
        plt.plot(q1anglesLineP2, q2anglesLineP2)
        plt.plot(q1anglesLineP3, q2anglesLineP3)
        plt.plot(q1anglesLineP4, q2anglesLineP4)
        plt.plot(q1anglesLineP5, q2anglesLineP5)

        plt.xlabel('q1 in degrees')
        plt.ylabel('q2 in degrees')
        plt.grid()
        
        print(q1anglesLine1, q2anglesLine1)
        return q1anglesLine1, q2anglesLine1, q1anglesLine2, q2anglesLine2, q1anglesLine3, q2anglesLine3, q1anglesLine4, q2anglesLine4, q1anglesLine5, q2anglesLine5 
    
    
    def createBackForthTrajectory(self): #sinusoidal trajectory
        
        ChangeInStep = 0.001
        freq = 6 * np.pi
        amplitude = .05
        yShift = .35
        x = np.arange(-.5, .5, ChangeInStep)
        y = np.sin(freq * np.arange(-.5, .5, ChangeInStep)) * amplitude + yShift

        # create the figures
        plt.figure()
        plt.title("Back and Forth Path Trajectory-Cartesian Space")
        plt.plot(x, y)
        plt.xlim(-0.6, 0.6)
        plt.ylim(0, .6)
        plt.xlabel('x position(m)')
        plt.ylabel('y position (m)')
        plt.grid()
        
        plt.figure()
        plt.title("Back and Forth Path Trajectory-Configuration Space")
        
        q1angles, q2angles = self.inverseKin(x, y)
        
        q1anglesP = [a * (180/np.pi) for a in q1angles] #converting from radians to degrees
        q2anglesP = [b * (180/np.pi) for b in q2angles]
        
        plt.plot(q1anglesP, q2anglesP, 'g')
        plt.xlabel('q1 in degrees')
        plt.ylabel('q2 in degrees')
        plt.grid()
        
        return q1angles, q2angles
        
    
    def inverseKin(self, x, y): #solving inverse kinematics for x,y cartesian location
        q1angles = []
        q2angles = []
    
        for i in range(len(x)):
            c = (math.pow(x[i],2) + math.pow(y[i],2) - math.pow(Main.l1,2) - math.pow(Main.l2,2)) / (2 * Main.l1 * Main.l2)
            s = np.sqrt(1 - math.pow(c,2))
            q1angles.append(np.arcsin((y[i] * (Main.l1 + Main.l2 * c) - x[i] * Main.l2 * s) / (math.pow(x[i],2) + math.pow(y[i],2)))) 
            q2angles.append(np.arccos((math.pow(x[i],2) + math.pow(y[i],2) - math.pow(Main.l1,2) - math.pow(Main.l2,2)) / (2 * Main.l1 * Main.l2))) 
        
        return (q1angles, q2angles)
            

    


    def allPossPoints(self): #defines all possible effector endpoint for a specific range. 
        
        theta1 = np.arange(0, 130 * (np.pi/180), 0.0174533) # all possible theta1 values
        theta2 = np.arange(0, 150 * (np.pi/180), 0.0174533) # all possible theta2 values
        
        [T1, T2] = np.meshgrid(theta1, theta2) # generate a grid of theta1 and theta2 values
        
        #Using Forward Kinematics Formula we find all possible values of x coordinate positions of the limb for theta1 and theta 2
        x = Main.l1 * np.cos(T1) + Main.l2 * np.cos(T1 + T2) # compute x coordinates
        y = Main.l1 * np.sin(T1) + Main.l2 * np.sin(T1 + T2) # compute y coordinates
        
        plt.figure()
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
        plt.figure()
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
        
        file = open("valuesDUMP.txt", "w")
        fileX = open("valuesDUMPX.txt", "w")
        fileY = open("valuesDUMPY.txt", "w")
        for i in self.my_range (0, qRange1, 0.0174533): #iterates by one degree
            for j in self.my_range(0,qRange2, 0.0174533):
                tempX = link2 * np.cos(i + j) + link1 * np.cos(i)
                tempY = link2 * np.sin(i + j) + link1 * np.sin(i)
                x.append(tempX)
                y.append(tempY)                
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
        file.close()
        fileX.close()  
        fileY.close() 
        return x, y
    
    def plotConfigSpaceStatic(self, qRange1, qRange2): #plots q1 angles versus q2 angles for a given task. 
        
        plotq1 = []
        plotq2 = []    
        plt.figure()
        for i in self.my_range (0, qRange1, 0.0174533): #iterates by one degree
            for j in self.my_range(0,qRange2, 0.0174533):
                plotq1.append(i)
                plotq2.append(j)
        
        plotq1 = [a * (180/np.pi) for a in plotq1]
        plotq2 = [b * (180/np.pi) for b in plotq2]
        plt.plot(plotq1, plotq2, 'g.') 

        plt.xlabel('q1 in degrees')
        plt.ylabel('q2 in degrees')
        plt.title('Static Case: Degree Angles in Configuration Space. ') 
        plt.grid()

Main()
plt.show()  


