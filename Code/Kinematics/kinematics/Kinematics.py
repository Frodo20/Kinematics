import numpy as np
import matplotlib.pyplot as plt 
from kinematics.Constraints import *
from kinematics.Drives import *

class Kinematics:
    def __init__(self,paras):
        self.T = 0
        self.Y = 0
        self.all_constraints = []
        for para in paras:
            '''
            para:[type, j, S_j, v_j, i, S_i, v_i, c, c_bar, C]
            '''
            type = para[0]
            
            #读取表格参数,创建约束矩阵phi
            if type =='ax':
                ax1 = ax(i=para[4],S_i=para[5],c=para[7])
                self.all_constraints.append(ax1)
            elif type == 'ay':
                ay1 = ay(i=para[4],S_i=para[5],c=para[7])
                self.all_constraints.append(ay1)
            elif type == 'aphi':
                aphi1 = aphi(i=para[4],c=para[7])
                self.all_constraints.append(aphi1)
            elif type =='ad':
                ad1 = ad(i=para[4],S_i=para[5],c=para[7],C=para[9])
                self.all_constraints.append(ad1)
            elif type == 'rx':
                rx1 = rx(i=para[4],S_i=para[5],j=para[1],S_j=para[2],c=para[7])
                self.all_constraints.append(rx1)
            elif type =='ry':
                ry1 = ry(i=para[4],S_i=para[5],j=para[1],S_j=para[2],c=para[7])
                self.all_constraints.append(ry1)
            elif type == 'rphi':
                rphi1 = rphi(i=para[4],j=para[1],c=para[7])
                self.all_constraints.append(rphi1)
            elif type=='rd':
                rd1 = rd(i=para[4],S_i=para[5],j=para[1],S_j=para[2],c=para[7])
                self.all_constraints.append(rd1)
            elif type=='r':
                r1 = r(i=para[4],S_i=para[5],j=para[1],S_j=para[2])
                self.all_constraints.append(r1)
            elif type=='t':
                t1 = t(i=para[4],S_i=para[5],v_i=para[6],j=para[1],S_j=para[2],v_j=para[3])
                self.all_constraints.append(t1)
            elif type=='rt':
                rt1 = t(i=para[4],S_i=para[5],v_i=para[6],j=para[1],S_j=para[2],c=para[7])
                self.all_constraints.append(rt1)

            elif type =='axd':
                axd1 = ax(i=para[4],S_i=para[5],c_bar=para[8])
                self.all_constraints.append(axd1)
            elif type == 'ayd':
                ayd1 = ay(i=para[4],S_i=para[5],c_bar=para[8])
                self.all_constraints.append(ayd1)
            elif type =='aphid':
                aphid1 = aphid(i=para[4],c_bar=para[8])
                self.all_constraints.append(aphid1)
            elif type =='add':
                add1 = ad(i=para[4],S_i=para[5],c=para[7],C=para[9],c_bar=para[8])
                self.all_constraints.append(add1)
            elif type == 'rphid':
                rphid1 = rphi(i=para[4],j=para[1],c_bar=para[8])
                self.all_constraints.append(rphid1)
            elif type=='rdd':
                rdd1 = rd(i=para[4],S_i=para[5],j=para[1],S_j=para[2],c_bar=para[8])
                self.all_constraints.append(rdd1)
            elif type=='tdd': #i、j互换
                tdd1 = tdd(i=para[1],S_i=para[2],v_i=para[3],j=para[4],S_j=para[5])
                self.all_constraints.append(tdd1)
            elif type=='tdd8': #i、j互换
                tdd1 = tdd8(i=para[1],S_i=para[2],v_i=para[3],j=para[4],S_j=para[5])
                self.all_constraints.append(tdd1)
            elif type=='tdd12': #i、j互换
                tdd1 = tdd12(i=para[1],S_i=para[2],v_i=para[3],j=para[4],S_j=para[5])
                self.all_constraints.append(tdd1)
    
    #构建约束矩阵            
    def constraints(self,q,t):
        i = 1
        for constraint in self.all_constraints:
            if i == 1:
                phi = constraint.constraints(q,t)
            else:
                phi = np.vstack((phi,constraint.constraints(q,t)))
            #print(i)
            #print(phi)
            i+=1
        #print (i)       
        return phi
    
    #构建雅克比矩阵
    def Jacobi(self,q,t,l):
        i =1
        for constraint in self.all_constraints:
            if i == 1:
                J = constraint.Jacobi(q,t,l)
            else:
                J = np.vstack((J,constraint.Jacobi(q,t,l)))
            #print(i)
            #print(constraint.Jacobi(q,t,l))
            #print(J)
            i+=1            
        return J

    #构建速度方程右项            
    def v(self,q,t):
        i=1
        for constraint in self.all_constraints:
            if i==1:
                v = constraint.v(q, t)
            else:
                v = np.vstack((v,constraint.v(q, t)))
            #print(i)
            #print(constraint.v(q, t))
            i+=1
        #print(v)
        return v

    #构建加速度右项        
    def gamma(self,q,dq,t):
        i=1
        for constraint in self.all_constraints:
            #print(constraint.gamma(q, dq, t))
            if i==1:
                gamma = constraint.gamma(q,dq,t)
                #print (gamma)
            else:
                gamma = np.vstack((gamma,constraint.gamma(q,dq,t)))
            #print(i)
            #print(gamma)
            i+=1
                
        return gamma
    
    #运动学求解
    def solve(self,q0=np.array([[0]]),t0=0,te=2.5,num=100,eps1=1e-10,eps2=1e-10,iter_max=15):
        #设定初值
        q = q0
        l = len(q)
        #运动学输出结果
        Y = np.zeros((num+1,3*l)) 

        
        dt = (te-t0)/num
        T = np.arange(t0,te+dt,dt)
        
        for i in range(num+1):
            t = t0+i*dt
            phi = self.constraints(q, t)
            #print(phi)   
            delta1 = np.linalg.norm(phi)
            iter_num = 0
            
            while delta1 > eps1:
                Jacobi = self.Jacobi(q, t, l)
                #print(Jacobi)
                if abs(np.linalg.det(Jacobi))<eps2:
                    print('Improper initial value,1')
                    
                
                dq = -np.dot(np.linalg.inv(Jacobi),phi)
                q = q+dq
                phi = self.constraints(q, t)
                delta1 = np.linalg.norm(phi)
                
                iter_num +=1
                if iter_num > iter_max:
                    print('Improper initial value,2')
                    return 
                
            Jacobi = self.Jacobi(q,t,l)
            v =self.v(q, t)
            if np.abs(np.linalg.det(Jacobi))<eps2:
                print('Sigular configuration')
                return
            
            dq = np.dot(np.linalg.inv(Jacobi),v)
            gamma = self.gamma(q, dq, t)
            ddq = np.dot(np.linalg.inv(Jacobi),gamma)
            Y[i,...] = np.c_[q.T,dq.T,ddq.T] #Y=[q,dq,ddq]
            
            q = q+dq*dt+ddq*dt**2/2
        
        self.Y=Y
        self.T=T
        #print(Y[100,:])
        return T,Y
    
    #绘制运动图像
    def plot(self,plot_x=1,plot_y=0,plot_phi=1,plot_dx=1,plot_dy=0,plot_dphi=1,plot_ddx=1,plot_ddy=0,plot_ddphi=1):
        
        font = {'family': 'Times New Roman',
        'color':  'darkred',
        'weight': '900',
        'size': 13,
        }
        
        #设置图的样式
        plt.style.use('seaborn-paper')
        #设置画布大小
        plt.figure(figsize=(12, 8),dpi=120)
        #plot q2
        
        if plot_x==1:
            #plt.subplot(1,2,1)
            plt.subplot(2,3,1)
            plt.plot(self.T, self.Y[...,3],linewidth=1.75)
            
            #plt.xlim((0,2.5))
            #plt.ylim((-0.05,0.25))
            
            #plt.xticks(np.arange(0,2.5,0.5))
            #plt.yticks(np.arange(-0.05,0.25,0.05))
            
            plt.grid()
            plt.title("$x_2-t$",fontdict=font)
            plt.xlabel("$t(s)$",fontsize=8,loc='right')
            plt.ylabel("$x_2(m)$",fontsize=8)
            #plt.show()
            
        if plot_y==1:
            plt.subplot(3,1,2)
            plt.plot(self.T, self.Y[...,4])
            plt.grid()
            plt.title("$Y_c-t$")
            plt.xlabel("$t(s)$")
            plt.ylabel("$Y_c(m)$")
            #plt.show()
            
        if plot_phi==1:
            #plt.subplot(1,2,2)
            plt.subplot(2,3,4)
            plt.plot(self.T, self.Y[...,5],linewidth=1.75)
            plt.grid()
            plt.title("$\phi_2-t$",fontdict=font)
            plt.xlabel("$t(s)$",fontsize=8,loc='right')
            plt.ylabel("$\phi_2(rad)$",fontsize=8)
            
            #plt.subplots_adjust(wspace=1,hspace=1)
            #plt.show()
        
        #plot dq2
         
        if plot_dx==1:
            #plt.subplot(1,2,1)
            plt.subplot(2,3,2)
            plt.plot(self.T, self.Y[...,12],linewidth=1.75)
            plt.grid()
            plt.title("$u_2-t$",fontdict=font)
            plt.xlabel("$t(s)$",fontsize=8,loc='right')
            plt.ylabel("$u_2(m/s)$",fontsize=8)
            #plt.show()
            
        if plot_dy==1:
            plt.subplot(3,1,2)
            plt.plot(self.T, self.Y[...,13])
            plt.grid()
            plt.title("$dY_c-t$")
            plt.xlabel("$t(s)$")
            plt.ylabel("$dY_c(m)$")
            #plt.show()
            
        if plot_dphi==1:
            #plt.subplot(1,2,2)
            plt.subplot(2,3,5)
            plt.plot(self.T, self.Y[...,14],linewidth=1.75)
            plt.grid()
            plt.title("$w_2-t$",fontdict=font)
            plt.xlabel("$t(s)$",fontsize=8,loc='right')
            plt.ylabel("$w_2(rad/s)$",fontsize=8)
            
            #plt.subplots_adjust(wspace=1,hspace=1)
            #plt.show()
            
        #plot ddq2
        
        if plot_ddx==1:
            #plt.subplot(1,2,1)
            plt.subplot(2,3,3)
            plt.plot(self.T, self.Y[...,21],linewidth=1.75)
            plt.ylim(-0.02,0)
            plt.grid()
            plt.title("$a_{2x}-t$",fontdict=font)
            plt.xlabel("$t(s)$",fontsize=8,loc='right')
            plt.ylabel("$a_{2x}(m/s^2)$",fontsize=8)
            #plt.show()
            
        if plot_ddy==1:
            plt.subplot(3,1,2)
            plt.plot(self.T, self.Y[...,22])
            plt.grid()
            plt.title("$ddY_c-t$")
            plt.xlabel("$t(s)$")
            plt.ylabel("$ddY_c(m)$")
            #plt.show()
            
        if plot_ddphi==1:
            #plt.subplot(1,2,2)
            plt.subplot(2,3,6)
            plt.plot(self.T, self.Y[...,23],linewidth=1.75)
            plt.grid()
            plt.title("$\\alpha_2-t$",fontdict=font)
            plt.xlabel("$t(s)$",fontsize=8,loc='right')
            plt.ylabel("$\\alpha_2(rad/s^2)$",fontsize=8)
            
        
        
        plt.subplots_adjust(wspace=0.4,hspace=0.5)
        plt.show()
    
    ''' def Animate(self):
        
        font = {'family': 'Times New Roman',
        'color':  'Black',
        'weight': '900',
        'size': 18,
        }
        
        q = self.Y[:,0:9]
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
            line1.set_data(xdata_1,ydata_1)
            line2.set_data(xdata_2,ydata_2)
            return line1, line2
        
        
        ax.grid()
        plt.axis('equal')
        
        
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
            
            #滑移铰
            xdata_3 = [0.6-np.sin(q[i,8])*0.03-np.cos(q[i,8])*0.05, 0.6-np.sin(q[i,8])*0.03+np.cos(q[i,8])*0.05]
            ydata_3 = [-np.sin(q[i,8])*0.05+np.cos(q[i,8])*0.03, np.sin(q[i,8])*0.05+np.cos(q[i,8])*0.03]
            line3.set_data(xdata_3, ydata_3)
            
            xdata_4 = [0.6-np.sin(q[i,8])*0.03-np.cos(q[i,8])*0.05, 0.6-np.sin(q[i,8])*0.03+np.cos(q[i,8])*0.05]
            ydata_4 = [-np.sin(q[i,8])*0.05-np.cos(q[i,8])*0.03,np.sin(q[i,8])*0.05-np.cos(q[i,8])*0.03]
            line4.set_data(xdata_4, ydata_4)
            return line1, line2, circle_line, line3, line4, 
        
        #(0,0)处的铰
        ax.plot(dpi_range/100,np.sqrt(1-dpi_range**2)/100,'b',lw=2)
        ax.plot(dpi_range/100,-np.sqrt(1-dpi_range**2)/100,'b',lw=2)
        #(0.6,0)处的铰链
        ax.plot(dpi_range/100+0.6,np.sqrt(1-dpi_range**2)/100,'b',lw=2)
        ax.plot(dpi_range/100+0.6,-np.sqrt(1-dpi_range**2)/100,'b',lw=2)
        
        
        line5, =ax.plot([-0.05,0.75],[0,0],  'r--', lw=1)
        plt.title('Animation',fontdict=font)
        
        #call the animator
        ani = animation.FuncAnimation(fig, run, frames=num-1, blit=False, interval=50, repeat=True, init_func=init)
            
        plt.show() '''
            
        
        
    
        
        
            
        
            
