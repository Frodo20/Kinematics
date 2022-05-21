
#-------------约束库-------------#
from math import gamma
import numpy as np

R=np.array([[0,-1],[1,0]])#旋转矩阵

#Constraints between a body and ground(Absolute constraints)

#绝对X约束
class ax:
    def __init__(self,i=0,S_i=np.zeros((2,1)),c=0):  #构造函数
        self.i=i
        self.S_i=S_i
        self.c=c
    
    def constraints(self,q,t):  #构造约束
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        x = q[3*self.i-3,0]
        phi_i = q[3*self.i+2-3,0]
        #print(self.c)
        phi = x+self.S_i[0,0]*np.cos(phi_i)-self.S_i[1,0]*np.sin(phi_i)-self.c
        return np.array([[phi]])
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[i*3+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def Jacobi(self,q,t,l): #计算雅可比矩阵 l=3*n,即为自由坐标数
        J = np.zeros((1,3))
        B = np.dot(R,self.A(q,self.i))
        J[0,0] = 1
        J[0,1] = 0
        J[0,2] = np.dot(np.dot(np.array([1,0]),B),self.S_i)
        if(l>3):
            if (self.i==1):
                J = np.hstack((J,np.zeros((1,l-3)))) 
            else:
                J = np.hstack((np.zeros((1,3*self.i-3)),J))
                J = np.hstack((J,np.zeros((1,l-3*self.i))))
        return J
    
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def gamma(self,q,dq,t): #计算加速度右项
        w = dq[3*self.i-3+2,0]
        gamma = np.dot(np.dot(np.array([1,0]),self.A(q,self.i)),self.S_i)*w**2
        return np.array([gamma])
        
#绝对Y约束
class ay:
    def __init__(self,i=0,S_i=np.zeros((2,1)),c=0):  #构造函数
        self.i=i
        self.S_i=S_i
        self.c=c
    
    def constraints(self,q,t):  #构造约束
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        y = q[3*self.i+1-3,0]
        phi_i = q[3*self.i-3+2,0]
        phi = y+self.S_i[0,0]*np.sin(phi_i)+self.S_i[1,0]*np.cos(phi_i)-self.c
        return np.array([[phi]])
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def Jacobi(self,q,t,l): #计算雅可比矩阵
        J = np.zeros((1,3))
        B = np.dot(R,self.A(q,self.i))
        J[0,0] = 0
        J[0,1] = 1
        J[0,2] = np.dot(np.dot(np.array([0,1]),B),self.S_i)
        if(l>3):
            if (self.i==1):
                J = np.hstack((J,np.zeros((1,l-3)))) 
            else:
                J = np.hstack((np.zeros((1,3*self.i-3)),J))
                J = np.hstack((J,np.zeros((1,l-3*self.i))))
        return J
    
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def gamma(self,q,dq,t): #计算加速度右项
        w = dq[3*self.i-3+2,0]
        gamma = np.dot(np.dot(np.array([0,1]),self.A(q,self.i)),self.S_i)*w**2
        return np.array([gamma])

#绝对角度约束
class aphi: 
    def __init__(self,i=0,c=0): #构造函数
        self.i=i
        self.c=c
        pass
    
    def constraints(self,q,t): #构造约束矩阵
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        phi_i = q[3*self.i+2-3,0]
        phi = phi_i - self.c
        return np.array([[phi]])
    
    def Jacobi(self,q,t,l): #构造雅可比矩阵
        J = np.array([[0,0,1]])
        if(l>3):
            if (self.i==1):
                J = np.hstack((J,np.zeros((1,l-3)))) 
            else:
                J = np.hstack((np.zeros((1,3*self.i-3)),J))
                J = np.hstack((J,np.zeros((1,l-3*self.i))))
        return J
    
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def gamma(self,q,dq,t): #计算加速度右项
        return np.array([[0]])

#绝对距离约束
class ad: 
    def __init__(self, i=0, S_i=np.zeros((2,1)),C=np.zeros((2,1)),c=0):
        self.i=i
        self.S_i=S_i
        self.C=C
        self.c=c
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def h(self,q,i): #计算h矩阵
        r = q[3*i-3:3*i+2-3]
        h = r + np.dot(self.A(q,self.i),self.S_i)-self.C
        return h
        
    def constraints(self,q,t):  #构造约束
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        phi = np.dot((self.h(q,self.i)).T,self.h(q,self.i))-self.c**2
        
        return np.array([[phi]])
    
    def Jacobi(self,q,t,l): #计算雅可比矩阵
        J = np.zeros((1,3))
        B = np.dot(R,self.A(q,self.i))
        J[0,0] = 1
        J[0,1] = 1
        J[0,2] = np.dot(B,self.S_i)
        J = 2*(self.h(q,self.i)).T*J
        if(l>3):
            if (self.i==1):
                J = np.hstack((J,np.zeros((1,l-3)))) 
            else:
                J = np.hstack((np.zeros((1,3*self.i-3)),J))
                J = np.hstack((J,np.zeros((1,l-3*self.i))))
        return J
    
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def h1(self,q,dq,i):
        w = dq[3*i-3+2,0]
        r1 = dq[3*i-3:3*i+2-3]
        B = np.dot(R,self.A(q,self.i))
        h1 = r1 + w*np.dot(B,self.S_i)
        return h1
    
    def gamma(self,q,dq,t): #计算加速度右项
        w = dq[3*self.i-3+2,0]
        gamma = 2*w*np.dot((self.h(q,self.i)).T,np.dot(self.A(q,self.i),self.S_i))-2*np.dot((self.h1(q, dq,self.i)).T,self.h1(q, dq,self.i))
        return np.array([gamma])

#Constraints between two bodies(Relative constraints)

#相对X约束
class rx:
    def __init__(self,i=0,S_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1)),c=0):  #构造函数
        self.i=i
        self.S_i=S_i
        self.j = j
        self.S_j = S_j
        self.c=c
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def constraints(self,q,t):  #构造约束
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        xi = q[self.i*3-3,0]
        xj = q[self.j*3-3,0]
        phi_i = q[self.i*3+2-3,0]
        phi_j = q[self.j*3+2-3,0]
        phi = xj+np.dot(np.array([1,0]),np.dot(self.A(q,self.j),self.S_j))-xi-np.dot(np.array([1,0]),np.dot(self.A(q,self.i),self.S_i))-self.c
        return np.array([[phi]])
    
    
    def Jacobi(self,q,t,l): #计算雅可比矩阵
        J = np.zeros((1,l))
        J1 = np.zeros((1,6))
        Bi = np.dot(R,self.A(q,self.i))
        Bj = np.dot(R,self.A(q,self.j))
        J1[0,0] = -1
        J1[0,1] = 0
        J1[0,2] = -np.dot(np.dot(np.array([1,0]),Bi),self.S_i)
        J1[0,3] = 1
        J1[0,4] = 0
        J1[0,5] = np.dot(np.dot(np.array([1,0]),Bj),self.S_j)
        
        J[0,3*self.i-3]=J1[0,0]
        J[0,3*self.i-3+1]=J1[0,1]
        J[0,3*self.i-3+2]=J1[0,2]
        J[0,3*self.j-3]=J1[0,3]
        J[0,3*self.j-3+1]=J1[0,4]
        J[0,3*self.j-3+2]=J1[0,5]
        return J
    
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def gamma(self,q,dq,t): #计算加速度右项
        wi = dq[3*self.i+2-3,0]
        wj = dq[3*self.j+2-3,0]
        gamma = np.dot(np.array([1,0]),np.dot(self.A(q,self.j),self.S_j)*wj**2-np.dot(self.A(q,self.i),self.S_i)*wi**2)
        return np.array([gamma])
    
#相对Y约束
class ry:
    def __init__(self,i=0,S_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1)),c=0):  #构造函数
        self.i=i
        self.S_i=S_i
        self.j = j
        self.S_j = S_j
        self.c=c
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def constraints(self,q,t):  #构造约束
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        yi = q[3*self.i+1-3,0]
        yj = q[3*self.j+1-3,0]
        phi_i = q[3*self.i+2-3,0]
        phi_j = q[3*self.j+2-3,0]
        phi = yj+np.dot(np.array([0,1]),np.dot(self.A(q,self.j),self.S_j))-yi-np.dot(np.array([0,1]),np.dot(self.A(q,self.i),self.S_i))-self.c
        return np.array([[phi]])
    
    def Jacobi(self,q,t,l): #计算雅可比矩阵
        J = np.zeros((1,6))
        Bi = np.dot(R,self.A(q,self.i))
        Bj = np.dot(R,self.A(q,self.j))
        J[0,0] = 0
        J[0,1] = -1
        J[0,2] = -np.dot(np.dot(np.array([0,1]),Bi),self.S_i)
        J[0,3] = 0
        J[0,4] = 1
        J[0,5] = np.dot(np.dot(np.array([0,1]),Bj),self.S_j)
        
        J1 = np.zeros((1,l))
        J1[0,3*self.i-3]=J[0,0]
        J1[0,3*self.i-3+1]=J[0,1]
        J1[0,3*self.i-3+2]=J[0,2]
        J1[0,3*self.j-3]=J[0,3]
        J1[0,3*self.j-3+1]=J[0,4]
        J1[0,3*self.j-3+2]=J[0,5]
        return J1
    
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def gamma(self,q,dq,t): #计算加速度右项
        wi = dq[3*self.i+2-3,0]
        wj = dq[3*self.j+2-3,0]
        gamma = np.dot(np.array([0,1]),np.dot(self.A(q,self.j),self.S_j)*wj**2-np.dot(self.A(q,self.i),self.S_i)*wi**2)
        return np.array([gamma])
    
#相对角度约束
class rphi: 
    def __init__(self,i=0,j=0,c=0): #构造函数
        self.i=i
        self.j=j
        self.c=c
    
    def constraints(self,q,t): #构造约束矩阵
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        phi_i = q[3*self.i+2-3,0]
        phi_j = q[3*self.j+2-3,0]
        phi = phi_j - phi_i - self.c
        return np.array([[phi]])
    
    def Jacobi(self,q,t,l): #构造雅可比矩阵
        J = np.zeros((1,l))
        J[0,3*self.i-3+2]=-1
        J[0,3*self.j-3+2]=1
        return J
    
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def gamma(self,q,dq,t): #计算加速度右项
        return np.array([[0]])

#相对距离约束    
class rd:
    def __init__(self,i=0,S_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1)),c=0):
        self.i = i
        self.S_i = S_i
        self.j = j
        self.S_j = S_j
        self.c = c

    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i

    def h(self,q): #计算h矩阵
        rj = q[3*self.j-3:3*self.i+2-3]
        ri = q[3*self.i-3:3*self.i+2-3]
        h = rj + np.dot(self.A(q,self.j),self.S_j)-ri-np.dot(self.A(q,self.i),self.S_i)
        return h

    def constraints(self,q,t):#构造约束
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        phi = np.dot((self.h(q)).T,self.h(q))-self.c**2

        return np.array([[phi]])

    def Jacobi(self,q,t,l): #计算雅可比矩阵
        J = np.zeros((1,6))
        Bi = np.dot(R,self.A(q,self.i))
        Bj = np.dot(R,self.A(q,self.j))
        J[0,0:2] = -2*(self.h(q)).T
        J[0,2] = -2*np.dot((self.h(q)).T,np.dot(Bi,self.S_i))
        J[0,3:5] = 2*self.h(q)
        J[0,5] = 2*np.dot((self.h(q)).T,np.dot(Bj,self.S_j))

        J1=np.zeros((1,l))
        J1[0,3*self.i-3:3*self.i-3+2] = J[0,0:2]
        J1[0,3*self.i-3+2] = J[0,2]
        J1[0,3*self.j-3:3*self.j-3+2] = J[0,3:5]
        J1[0,3*self.j-3+2] = J[0,5]

        return J1
    
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def h1(self,q,dq):
        wi = dq[3*self.i-3+2,0]
        wj = dq[3*self.j-3+2,0]
        dqi = dq[3*self.i-3:3*self.i-3+2,0]
        dqj = dq[3*self.j-3:3*self.j-3+2,0]
        Bi = np.dot(R,self.A(q,self.i))
        Bj = np.dot(R,self.A(q,self.j))
        ri1 = dqi+np.dot(np.dot(Bi,self.S_i),wi)
        rj1 = dqj+np.dot(np.dot(Bj,self.S_j),wj)
        return rj1-ri1
    
    def gamma(self,q,dq,t): #计算加速度右项
        wi = dq[3*self.i-3+2,0]
        wj = dq[3*self.j-3+2,0]
        
        gamma = -2*np.dot((self.h1(q,dq)).T,self.h1(q,dq))+2*np.dot((self.h1(q,dq).T,np.dot(self.A(q,self.j),self.S_j)*wj**2-np.dot(self.A(q,self.i),self.S_i)*wi**2))
        return np.array([gamma])


#铰约束
class r:

    def __init__(self,i=0,S_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1))):
        self.i=i
        self.S_i = S_i
        self.j = j
        self.S_j = S_j
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def constraints(self,q,t):
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        ri = q[3*self.i-3:3*self.i-3+2]
        rj = q[3*self.j-3:3*self.j-3+2]
        #xi=q[3*self.i-3,0]
        #xj=q[3*self.j-3,0]
        #yi=q[3*self.i+1-3,0]
        #yj=q[3*self.j+1-3,0]
        #print(xj)
        #print(xi)
        #print(self.S_j)
        #print(self.S_i)
        #print(np.cos(np.pi/2))
        #print(np.dot(np.array([1,0]),np.dot(self.A(q,self.j),self.S_j)))
        #print(np.dot(np.array([1,0]),np.dot(self.A(q,self.i),self.S_i)))
        #print(self.A(q,self.i))
        #phi1 = xj+np.dot(np.array([1,0]),np.dot(self.A(q,self.j),self.S_j))-xi-np.dot(np.array([1,0]),np.dot(self.A(q,self.i),self.S_i))
        #print(phi1)
        #phi2 = yj+np.dot(np.array([0,1]),np.dot(self.A(q,self.j),self.S_j))-yi-np.dot(np.array([0,1]),np.dot(self.A(q,self.i),self.S_i))
        #print(phi2)
        phi = rj+np.dot(self.A(q,self.j),self.S_j)-ri-np.dot(self.A(q,self.i),self.S_i)
        #phi = np.vstack((phi1,phi2))
        #print(phi)
        return phi
    
    def Jacobi(self,q,t,l):
        J = np.zeros((2,6))
        Bi = np.dot(R,self.A(q,self.i))
        Bj = np.dot(R,self.A(q,self.j))
        J[0,0] = -1
        J[0,1] = 0
        J[0,2] = -np.dot(np.dot(np.array([1,0]),Bi),self.S_i)
        J[0,3] = 1
        J[0,4] = 0
        J[0,5] = np.dot(np.dot(np.array([1,0]),Bj),self.S_j)
        J[1,0] = 0
        J[1,1] = -1
        J[1,2] = -np.dot(np.dot(np.array([0,1]),Bi),self.S_i)
        J[1,3] = 0
        J[1,4] = 1
        J[1,5] = np.dot(np.dot(np.array([0,1]),Bj),self.S_j)
        
        J1 = np.zeros((2,l))
        J1[0,3*self.i-3]=J[0,0]
        J1[0,3*self.i-3+1]=J[0,1]
        J1[0,3*self.i-3+2]=J[0,2]
        J1[0,3*self.j-3]=J[0,3]
        J1[0,3*self.j-3+1]=J[0,4]
        J1[0,3*self.j-3+2]=J[0,5]
        
        J1[1,3*self.i-3]=J[1,0]
        J1[1,3*self.i-3+1]=J[1,1]
        J1[1,3*self.i-3+2]=J[1,2]
        J1[1,3*self.j-3]=J[1,3]
        J1[1,3*self.j-3+1]=J[1,4]
        J1[1,3*self.j-3+2]=J[1,5]
        return J1
    
    def v(self,q,t): #计算速度右项
        return np.array([[0],[0]])
    
    def gamma(self,q,dq,t):
        wi = dq[3*self.i+2-3,0]
        wj = dq[3*self.j+2-3,0]
        #print(wi)
        #print(wj)
        #print(np.dot(self.A(q,self.j),self.S_j))
        #print(np.dot(self.A(q,self.j),self.S_j)*wj**2-np.dot(self.A(q,self.i),self.S_i)*wi**2)
        #gamma1 = np.dot(np.array([1,0]),np.dot(self.A(q,self.j),self.S_j)*wj**2-np.dot(self.A(q,self.i),self.S_i)*wi**2)
        #print(gamma1)
        #gamma2 = np.dot(np.array([0,1]),np.dot(self.A(q,self.j),self.S_j)*wj**2-np.dot(self.A(q,self.i),self.S_i)*wi**2)
        #gamma = np.array([gamma1,gamma2])
        gamma = np.dot(self.A(q,self.j),self.S_j)*wj**2-np.dot(self.A(q,self.i),self.S_i)*wi**2
        #print(gamma)
        return gamma
    
class t: #滑移铰
    def __init__(self,i=0,S_i=np.zeros((2,1)),v_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1)),v_j=np.zeros((2,1))):
        self.i = i
        self.S_i = S_i
        self.v_i = v_i
        self.j = j
        self.S_j = S_j
        self.v_j = v_j
        
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def Aij(self,q):
        return np.dot((self.A(q,self.i)).T,self.A(q,self.j))
    
    def h(self,q):
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        #print(ri)
        #print(rj)
        h = rj+np.dot(self.A(q,self.j),self.S_j)-ri-np.dot(self.A(q,self.i),self.S_i)
        return h
        
    def constraints(self,q,t):
        Bi = np.dot(R,self.A(q,self.i))
        Bij = np.dot(R,self.Aij(q))
        phi1 = np.dot((np.dot(Bi,self.v_i)).T,self.h(q))
        phi2 = -np.dot((self.v_i).T,np.dot(Bij,self.v_j))
        phi = np.vstack((phi1,phi2))
        #phi = np.concatenate(phi1, phi2),axis=0)
        return phi
    
    def Jacobi(self,q,t,l):
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        
        J = np.zeros((2,6))
        #print(np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j)))
        #print((self.v_i).T)
        #print((self.A(q, self.i)).T)
        #print(rj-ri)
        J[0,0:2] = -np.dot((self.v_i).T,Bi.T)
        J[0,2] = -np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),rj-ri)-np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j))
        J[0,3:5] = np.dot((self.v_i).T,Bi.T)
        #print(self.i)
        #print(self.S_j)
        #print(self.Aij(q))
        J[0,5] = np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j))
        J[1,0] = 0
        J[1,1] = 0
        J[1,2] = -np.dot((self.v_i).T,np.dot(self.Aij(q),self.v_j))
        J[1,3] = 0
        J[1,4] = 0
        J[1,5] = np.dot((self.v_i).T,np.dot(self.Aij(q),self.v_j))
        
        J1 = np.zeros((2,l))
        J1[0,3*self.i-3:3*self.i+2-3]=J[0,0:2]
        J1[0,3*self.i-3+2]=J[0,2]
        J1[0,3*self.j-3:3*self.j+2-3]=J[0,3:5]
        J1[0,3*self.j+2-3]=J[0,5]
        
        J1[1,3*self.i-3]=J[1,0]
        J1[1,3*self.i-3+1]=J[1,1]
        J1[1,3*self.i-3+2]=J[1,2]
        J1[1,3*self.j-3]=J[1,3]
        J1[1,3*self.j-3+1]=J[1,4]
        J1[1,3*self.j-3+2]=J[1,5]
        return J1
              
    def v(self,q,t): #计算速度右项
        return np.array([[0],[0]])
    
    def gamma(self,q,dq,t):
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        dri = dq[3*self.i-3:3*self.i+2-3]
        drj = dq[3*self.j-3:3*self.j+2-3]
        wi = dq[3*self.i+2-3,0]
        #wj = dq[3*self.j+2-3,0]
        
        #print(np.dot(np.dot((self.v_i).T,Bi.T),rj-ri)*wi**2)
        gamma1 = np.dot(np.dot((self.v_i).T,Bi.T),rj-ri)*wi**2+2*np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),drj-dri)*wi
        #print(gamma1[0,0])
        gamma = np.array([[gamma1[0,0]],[0]])
        #print(gamma)
        return gamma    
    
#Composite joints

#Revolute-revolute joints is same as the rd, no need to code again

#Revolute-translational composite joint
class rt:
    def __init__(self,i=0,S_i=np.zeros((2,1)),v_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1)),c=0):
        self.i = i
        self.S_i = S_i
        self.v_i = v_i
        self.j = j
        self.S_j = S_j
        self.c = c
        
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def vi(self):
        v1 = self.v_i[0,0]
        v2 = self.v_i[1,0]
        vi1 = np.sqrt(v1**2+v2**2)
        return vi1
    
    def Aij(self,q):
        return np.dot((self.A(q,self.i)).T,self.A(q,self.j))
    
    def h(self,q):
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        #print(ri)
        #print(rj)
        h = rj+np.dot(self.A(q,self.j),self.S_j)-ri-np.dot(self.A(q,self.i),self.S_i)
        return h
        
    def constraints(self,q,t):
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        Bi = np.dot(R,self.A(q,self.i))
        Bij = np.dot(R,self.Aij(q))

        phi = (np.dot(np.dot((self.v_i).T,Bi.T),rj-ri)-np.dot((self.v_i).T,np.dot(Bij,self.S_j)-np.dot((self.v_i).T),np.dot(R.T,self.S_i)))/self.vi()-self.c
        return np.array([[phi]])
    
    def Jacobi(self,q,t,l):
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        
        J = np.zeros((1,6))
        #print(np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j)))
        #print((self.v_i).T)
        #print((self.A(q, self.i)).T)
        #print(rj-ri)
        J[0,0:2] = -np.dot((self.v_i).T,Bi.T)/self.vi()
        J[0,2] = (-np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),rj-ri)-np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j)))/self.vi()
        #print(self.i)
        #print(self.S_j)
        #print(self.Aij(q))
        J[0,3:5] = np.dot((self.v_i).T,Bi.T)/self.vi()
        J[0,5] = (np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j)))/self.vi()
        
        J1 = np.zeros((1,l))
        J1[0,3*self.i-3:3*self.i+2-3]=J[0,0:2]
        J1[0,3*self.i-3+2]=J[0,2]
        J1[0,3*self.j-3:3*self.j+2-3]=J[0,3:5]
        J1[0,3*self.j+2-3]=J[0,5]
        
        return J1
              
    def v(self,q,t): #计算速度右项
        return np.array([[0]])
    
    def gamma(self,q,dq,t):
        Bi = np.dot(R,self.A(q,self.i))
        Bij = np.dot(R,self.Aij(q))
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        dri = dq[3*self.i-3:3*self.i+2-3]
        drj = dq[3*self.j-3:3*self.j+2-3]
        wi = dq[3*self.i+2-3,0]
        wj = dq[3*self.j+2-3,0]
        
        gamma1 = (2*np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),(drj-dri)*wi)+np.dot(np.dot((self.v_i).T,Bi.T),(rj-ri)*wi**2)-np.dot((self.v_i).T,np.dot(Bij,self.S_j))*(wj-wi)**2)/self.vi()
        #print(gamma1)
        return np.array([gamma1])