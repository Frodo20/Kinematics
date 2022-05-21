import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation  
import imageio
from PIL import Image, ImageSequence


def Animate(q,save_plot=1):
        
        #字体定义
        font = {'family': 'Times New Roman',
        'color':  'Black',
        'weight': '900',
        'size': 18,
        }
        
        #q = self.Y[:,0:9]
        num = q.shape[0]
        #print(num)
        
        theta = np.linspace(0,2*np.pi,num)
        
         # create a figure, axis and plot element
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
        dpi_range = np.linspace(-1,1,100)
        
        
        #杆1
        xdata_1 = [0, 0]
        ydata_1 = [0, 0.25]
        line1, = ax.plot(xdata_1, ydata_1, lw=2) #线宽lw=2
        
        #杆2
        xdata_2 = [0,0.6]
        ydata_2 = [0.25,0]
        line2, = ax.plot(xdata_2,ydata_2,lw=2)
        
        #约束铰
        circle_line, = ax.plot([],[],lw=2)
        
        #滑移铰
        line3, = ax.plot([],[],lw=2)
        line4, = ax.plot([],[],lw=2)
        
        
        # initialization function
        def init():
            # creating an empty plot/frame
            ax.set_ylim(-0.4,0.6)
            ax.set_xlim(-0.1,0.9)
            del xdata_1[:]
            del ydata_1[:]
            del xdata_2[:]
            del ydata_2[:]
            ''' del xdata_3[:]
            del ydata_3[:]
            del xdata_4[:]
            del ydata_4[:] '''
            line1.set_data(xdata_1,ydata_1)
            line2.set_data(xdata_2,ydata_2)
            ''' line3.set_data(xdata_3,ydata_3)
            line4.set_data(xdata_4,ydata_4) '''

            return line1, line2
        
        
        ax.grid() #网格线
        plt.axis('equal') #横纵坐标等值
        
        
        def run(i):
            
            #计算转换矩阵
            A1_i = np.array([[np.cos(q[i,2]),-np.sin(q[i,2])],[np.sin(q[i,2]),np.cos(q[i,2])]])
            A2_i = np.array([[np.cos(q[i,5]),-np.sin(q[i,5])],[np.sin(q[i,5]),np.cos(q[i,5])]])
            S1 = np.array([[0.25],[0]])
            S2 = np.array([[0.65],[0]])
            
            #根据i更新1杆的位置
            r1 = q[i,0:2]
            r2 = q[i,0:2]+(np.dot(A1_i,S1)).T
            #print(r1)
            #print(r2)
            xdata_1 = [r1[0],r2[0,0]]
            #print (xdata_1)
            ydata_1 = [r1[0],r2[0,1]]
            #print(ydata_1)
            line1.set_data(xdata_1,ydata_1)
            
            
            #根据i更新2杆的位置
            r3 = q[i,3:5]+(np.dot(A2_i,S2)).T
            xdata_2 = [r2[0,0],r3[0,0]]
            #print(xdata_2)
            ydata_2 = [r2[0,1],r3[0,1]]
            #print(ydata_2)
            line2.set_data(xdata_2,ydata_2)
            
            #铰链
            circle_x = [0.01*np.cos(theta[j])+r2[0,0] for j in range(len(theta))]
            circle_y = [0.01*np.sin(theta[j])+r2[0,1] for j in range(len(theta))]
            circle_line.set_data(circle_x,circle_y)
            
            #滑移铰(上半部分)
            xdata_3 = [0.6-np.sin(q[i,8])*0.03-np.cos(q[i,8])*0.04, 0.6-np.sin(q[i,8])*0.03+np.cos(q[i,8])*0.04]
            ydata_3 = [-np.sin(q[i,8])*0.04+np.cos(q[i,8])*0.03, np.sin(q[i,8])*0.04+np.cos(q[i,8])*0.03]
            line3.set_data(xdata_3, ydata_3)
            
            #滑移跤(下半部分)
            xdata_4 = [0.6-np.sin(q[i,8])*0.03-np.cos(q[i,8])*0.04, 0.6-np.sin(q[i,8])*0.03+np.cos(q[i,8])*0.04]
            ydata_4 = [-np.sin(q[i,8])*0.04-np.cos(q[i,8])*0.03,np.sin(q[i,8])*0.04-np.cos(q[i,8])*0.03]
            line4.set_data(xdata_4, ydata_4)
            return line1, line2, circle_line, line3, line4, 
        
        #(0,0)处的铰
        ax.plot(dpi_range/100,np.sqrt(1-dpi_range**2)/100,'b',lw=2)
        ax.plot(dpi_range/100,-np.sqrt(1-dpi_range**2)/100,'b',lw=2)
        #(0.6,0)处的铰链
        ax.plot(dpi_range/100+0.6,np.sqrt(1-dpi_range**2)/100,'b',lw=2)
        ax.plot(dpi_range/100+0.6,-np.sqrt(1-dpi_range**2)/100,'b',lw=2)
        
        #绘制地平线
        line5, =ax.plot([-0.05,0.75],[0,0],  'r--', lw=1)
        plt.title('Animation',fontdict=font)
        
        #call the animator
        ani1 = animation.FuncAnimation(fig, run, frames=num-1, blit=False, interval=50, repeat=True, init_func=init)
            
        plt.show()

        if save_plot==1:

            ani2 = animation.FuncAnimation(fig, run, frames=num-1, blit=False, interval=50, repeat=False, init_func=init)

            ani2.save("../roll.gif",writer='imagemagick', fps=100)