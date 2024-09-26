import numpy as np
import matplotlib.pyplot as plt
import argparse
import celluloid
def parse():
    """Argument Parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument("initThetaOne",help="Initial Theta Pendulum 1",type=float)
    parser.add_argument("initThetaTwo",help="Initial Theta Pendulum 2",type=float)

    parser.add_argument("initThetaDot_One",help="Initial Theta Dot Pendulum 1",type=float)
    parser.add_argument("initThetaDot_Two",help="Initial Theta Dot Pendulum 2",type=float)

    parser.add_argument("Mass1",help="Mass of pendulum 1",type=float)
    parser.add_argument("Mass2",help="Mass of pendulum 2",type=float)

    parser.add_argument("g",help="acceleration due to gravity",type=float)
    parser.add_argument("L1",help="Length of Pendulum 1",type=float)
    parser.add_argument("L2",help="Length of Pendulum 2",type=float)
    
    return parser.parse_args()


def getThetaDoubleDotPendulum1(theta1 :type=float,theta2:type=float,theta_dot_1: type=float,theta_dot_2: type=float ,mass1: type=float,mass2: type=float,g: type=float, L_1: type=float, L_2: type=float):
    """Returns the value of the second time derivative of theta for the first pendulum"""
    return ((-g*(2*(mass1 + mass2))*np.sin(theta1)) - (mass2*g*np.sin(theta1-(2*theta2))) -(2*np.sin(theta1-theta2)*mass2*(((theta_dot_2**2)*L_2)+((theta_dot_1**2)*L_1*np.cos(theta1-theta2))))) / (L_1*((2*mass1)+mass2-(mass2*np.cos(2*(theta1-theta2)))))

def getThetaDoubleDotPendulum2(theta1,theta2,theta_dot_1,theta_dot_2,mass1,mass2,g,L_1,L_2):
    """Returns the value of the second time derivative of theta for the second pendulum"""
    return (2*np.sin(theta1-theta2)*(((theta_dot_1**2)*L_1*(mass1+mass2))+(g*(mass1+mass2)*np.cos(theta1))+(theta_dot_2**2)*L_2*mass2*np.cos(theta1-theta2))) / (L_2*(2*mass1)+mass2-(mass2*np.cos((2*(theta1-theta2)))))


def calculation(args,time_step,final_time):
    """Gets the value of theta at specified intervals in a specified time duration"""
    theta_one = args.initThetaOne
    theta_two = args.initThetaTwo

    theta_dot_one = args.initThetaDot_One
    theta_dot_two = args.initThetaDot_Two

    pendulum1_mass = args.Mass1
    pendulum2_mass = args.Mass2

    g = args.g
    L_1 = args.L1
    L_2 = args.L2


    theta1 = np.zeros(len(np.arange(0,final_time,time_step)))
    theta2 = np.zeros(len(np.arange(0,final_time,time_step)))

    thetadot_one = np.zeros(len(np.arange(0,final_time,time_step)))
    thetadot_two = np.zeros(len(np.arange(0,final_time,time_step)))

    for i in range(len(np.arange(0,final_time,time_step))):
        #Assigning values 
        theta1[i] = theta_one
        theta2[i] = theta_two
        thetadot_one[i] = theta_dot_one
        thetadot_two[i] = theta_dot_two

        theta_double_dot_1 = getThetaDoubleDotPendulum1(
            theta_one,theta_two,theta_dot_one, theta_dot_two,pendulum1_mass, pendulum2_mass, g, L_1, L_2 
        )
        #Gets the value of the second time derivative of theta for the first pendulum


        theta_double_dot_2 = getThetaDoubleDotPendulum2(
            theta_one,theta_two,theta_dot_one, theta_dot_two,pendulum1_mass, pendulum2_mass, g, L_1, L_2
        )
         #Gets the value of the second time derivative of theta for the second pendulum

        theta_one += (theta_dot_one * time_step)       
        theta_two += (theta_dot_two * time_step)
        theta_dot_one += (theta_double_dot_1 * time_step)
        theta_dot_two += (theta_double_dot_2 * time_step)
        #The change in each value should be equal to their rates of change multiplied by the time step
        
    return np.vstack((theta1,theta2)),np.vstack((thetadot_one,thetadot_two)),np.arange(0,final_time,time_step)

def animate_double_pendulum(t, theta,args,trail_size, fps_mul=1.0):
    """Animate the solution found by THETA()."""


    tet1 = theta[0] 
    tet2 = theta[1]

    dt = t[0]

    el = args.L1
    el2 = args.L2       #Lengths of the pendulums

    #x,y coordinates for the pendulum bobs
    x1 = el*np.sin(tet1)      
    y1 = -el*np.cos(tet1)

    x2 = x1 + el2*np.sin(tet2)
    y2 = y1 - el2*np.cos(tet2)

    fig = plt.figure(figsize=(8,6))
    plt.xlim((-2.1*el, 2.1*el))
    plt.axis('equal')

    cam = celluloid.Camera(fig)

    trail = np.ones((trail_size,2)) #Makes an array to hold the values for the trail
    trail *= np.NaN #Makes it so that values which weren't added aren't plotted

    
    for k in range(t.size):
        plt.plot([0,x1[k]],[0,y1[k]],'b-',linewidth=4)
        plt.plot([x1[k],x2[k]],[y1[k],y2[k]],'r-',linewidth=4)


        trail = np.roll(trail,(1,1))    #Rolls the values in the array so that the oldest value is replaced
        trail[0,0] = x2[k]
        trail[0,1] = y2[k]

        plt.plot(trail[:,0],trail[:,1],color='grey',zorder=-999999) #Plots the trail

        cam.snap()
    anim = cam.animate(interval=1,repeat=False)

    anim.save("Double Pendulum.gif") #Saves the animation as a gif

    plt.show()
    return anim


def main():
    args = parse()
    pendTheta,thetaDot,time = calculation(args,0.02,15)
    animate_double_pendulum(time,pendTheta,args,25)

if __name__ == "__main__":
    main()