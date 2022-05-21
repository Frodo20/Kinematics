import numpy as np
''' import sympy as sym

t = sym.symbols('t') #定义时间变量t为符号变量 '''

R=np.array([[0,-1],[1,0]])#旋转矩阵

#所有驱动的表达式均为0.65-0.1*t

#Absolute Drive
#绝对x驱动
class axd:
    def __init__(self,i=0,S_i=np.zeros((2,1)),c_bar=0):  #构造函数
        self.i=i
        self.S_i=S_i
        self.c=c_bar
    
    def constraints(self,q,t):  #构造约束
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        x = q[3*self.i-3,0]
        phi_i = q[2,0]
        #print(self.c)
        phi = x+self.S_i[0,0]*np.cos(phi_i)-self.S_i[1,0]*np.sin(phi_i)-(0.65-0.1*t) #c_bar=0.65-0.1*t
        return np.array([[phi]])
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[i*3+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def Jacobi(self,q,t,l): #计算雅可比矩阵
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
        return np.array([[-0.1]]) #\dot{c(t)}
    
    def gamma(self,q,dq,t): #计算加速度右项
        w = dq[3*self.i-3+2,0]
        gamma = np.dot(np.dot(np.array([1,0]),self.A(q,self.i)),self.S_i)*w**2
        return np.array([gamma]) 

#绝对y驱动   
class ayd:
    def __init__(self,i=0,S_i=np.zeros((2,1)),c_bar=0):  #构造函数
        self.i=i
        self.S_i=S_i
        self.c_bar=c_bar
    
    def constraints(self,q,t):  #构造约束
        '''
        :param q:[[x],[y],[phi_i]]
        :return: phi
        '''
        y = q[3*self.i+1-3,0]
        phi_i = q[2,0]
        phi = y+self.S_i[0,0]*np.sin(phi_i)+self.S_i[1,0]*np.cos(phi_i)-(0.65-0.1*t) #c_bar=0.65-0.1*t
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
        return np.array([[-0.1]])  #\dot{c(t)}
    
    def gamma(self,q,dq,t): #计算加速度右项
        w = dq[3*self.i+2-3,0]
        gamma = np.dot(np.dot(np.array([0,1]),self.A(q,self.i)),self.S_i)*w**2
        return np.array([gamma])

#绝对角度驱动
class aphid: 
    def __init__(self,i=0,c_bar=0):
        self.i=i

    def constraints(self,q,t):
        '''
        :param q: [[x],[y],[theta]]
        :return: phi
        '''
        theta = q[3*self.i+2-3, 0]
        phi=theta-0.65+0.1*t
        return np.array([[phi]])

    def Jacobi(self,q,t,l):
        J = np.array([[0, 0, 1]])
        if(l>3):
            if (self.i==1):
                J = np.hstack((J,np.zeros((1,l-3)))) 
            else:
                J = np.hstack((np.zeros((1,3*self.i-3)),J))
                J = np.hstack((J,np.zeros((1,l-3*self.i))))
        return J

    def v(self, q,t):
        return np.array([[-0.1]])

    def gamma(self, q, dq,t):
        return np.array([[0]])
    
#绝对距离驱动
class add: 
    def __init__(self, i=0, S_i=np.zeros((2,1)),C=np.zeros((2,1)),c_bar=0):
        self.i=i
        self.S_i=S_i
        self.C=C
        self.c_bar=c_bar
    
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
        return np.array([[2*(0.65-0.1*t)*(-0.1)]])
    
    def h1(self,q,dq,i):
        w = dq[3*i-3+2,0]
        r1 = dq[3*i-3:3*i+2-3]
        B = np.dot(R,self.A(q,self.i))
        h1 = r1 + w*np.dot(B,self.S_i)
        return h1
    
    def gamma(self,q,dq,t): #计算加速度右项
        w = dq[3*self.i+2-3,0]
        gamma = 2*w*np.dot((self.h(q,self.i)).T,np.dot(self.A(q,self.i),self.S_i))-2*np.dot((self.h1(q, dq,self.i)).T,self.h1(q, dq,self.i))+2*(0.1)**2
        return np.array([gamma])

#Relative Dirve

#相对角驱动
class rphid: 
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
        return np.array([[-0.1]])
    
    def gamma(self,q,dq,t): #计算加速度右项
        return np.array([[0]])

#相对旋转驱动
#rdd（与相对角驱动约束相同）

#相对距离约束
class rdd:
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
        return np.array([[2*(0.65-0.1*t)*(-0.1)]])
    
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
        
        gamma = -2*np.dot((self.h1(q,dq)).T,self.h1(q,dq))+2*np.dot((self.h1(q,dq).T,np.dot(self.A(q,self.j),self.S_j)*wj**2-np.dot(self.A(q,self.i),self.S_i)*wi**2))+2*(-0.1)**2
        return np.array([gamma])


#平行距离驱动
class tdd8: 
    def __init__(self,i=0,S_i=np.zeros((2,1)),v_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1)),c_bar=0):
        self.i=i
        self.S_i=S_i
        self.v_i=v_i
        self.j=j
        self.S_j=S_j
        self.c_bar=c_bar
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def Aij(self,q):
        #print(np.dot((self.A(q,self.i)).T,self.A(q,self.j)))
        return np.dot((self.A(q,self.i)).T,self.A(q,self.j))
    
    def vi(self,q):
        vi = np.dot(self.A(q,self.i),self.v_i)
        v1 = vi[0,0]
        v2 = vi[1,0]
        vi1 = np.sqrt(v1**2+v2**2)
        return vi1
    
    def constraints(self,q,t):
        '''
        :param q: [[x],[y],[theta]]
        :return: phi
        '''
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        #print(self.v_i)
        #print(self.Aij(q))
        #print(self.S_j)
        #A=np.dot((self.v_i).T,self.Aij(q))
        #print(np.dot(A,self.S_j))
        
        phi = np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),rj-ri)+np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j))-np.dot((self.v_i).T,self.S_i)-self.vi(q)*(0.65-0.08*t) #c_bar=0.65-0.1*t
        #phi = np.dot(np.dot((self.v_i).T,self.Aij(q)),self.S_j)
        #print(self.vi(q))
        print(phi)
        return np.array(phi)
    
    def Jacobi(self,q,t,l):
        J = np.zeros((1,6))
        
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[0:2]
        rj = q[0:2]
        
        J[0,0:2]=-np.dot((self.v_i).T,(self.A(q,self.i)).T)
        J[0,2] = np.dot(np.dot((self.v_i).T,Bi.T),rj-ri)
        J[0,3:5]=np.dot((self.v_i).T,(self.A(q,self.i)).T)
        J[0,5]=0
        
        J1 = np.zeros((1,l))
        J1[0,3*self.i-3:3*self.i+2-3]=J[0,0:2]
        J1[0,3*self.i-3+2]=J[0,2]
        J1[0,3*self.j-3:3*self.j+2-3]=J[0,3:5]
        J1[0,3*self.j+2-3]=J[0,5]
              
        return J1
        
    def v(self,q,t):
        v=-0.08*self.vi(q)
        #print(v)
        return np.array([[v]])
    
    def gamma(self,q,dq,t):
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[3*self.i-3:3*self.i-3+2]
        rj = q[3*self.j-3:3*self.j-3+2]
        dri = dq[3*self.i-3:3*self.i-3+2]
        drj = dq[3*self.j-3:3*self.j-3+2]
        wi = dq[3*self.i-3+2,0]
        
        gamma = np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),rj-ri)*wi**2-2*np.dot(np.dot((self.v_i).T,Bi.T),drj-dri)*wi   
        #print(gamma)
        return gamma

class tdd: 
    def __init__(self,i=0,S_i=np.zeros((2,1)),v_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1)),c_bar=0):
        self.i=i
        self.S_i=S_i
        self.v_i=v_i
        self.j=j
        self.S_j=S_j
        self.c_bar=c_bar
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def Aij(self,q):
        #print(np.dot((self.A(q,self.i)).T,self.A(q,self.j)))
        return np.dot((self.A(q,self.i)).T,self.A(q,self.j))
    
    def vi(self,q):
        vi = np.dot(self.A(q,self.i),self.v_i)
        v1 = vi[0,0]
        v2 = vi[1,0]
        vi1 = np.sqrt(v1**2+v2**2)
        return vi1
    
    def constraints(self,q,t):
        '''
        :param q: [[x],[y],[theta]]
        :return: phi
        '''
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        #print(self.v_i)
        #print(self.Aij(q))
        #print(self.S_j)
        #A=np.dot((self.v_i).T,self.Aij(q))
        #print(np.dot(A,self.S_j))
        
        phi = np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),rj-ri)+np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j))-np.dot((self.v_i).T,self.S_i)-self.vi(q)*(0.65-0.1*t) #c_bar=0.65-0.1*t
        #phi = np.dot(np.dot((self.v_i).T,self.Aij(q)),self.S_j)
        #print(self.vi(q))
        print(phi)
        return np.array(phi)
    
    def Jacobi(self,q,t,l):
        J = np.zeros((1,6))
        
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[0:2]
        rj = q[0:2]
        
        J[0,0:2]=-np.dot((self.v_i).T,(self.A(q,self.i)).T)
        J[0,2] = np.dot(np.dot((self.v_i).T,Bi.T),rj-ri)
        J[0,3:5]=np.dot((self.v_i).T,(self.A(q,self.i)).T)
        J[0,5]=0
        
        J1 = np.zeros((1,l))
        J1[0,3*self.i-3:3*self.i+2-3]=J[0,0:2]
        J1[0,3*self.i-3+2]=J[0,2]
        J1[0,3*self.j-3:3*self.j+2-3]=J[0,3:5]
        J1[0,3*self.j+2-3]=J[0,5]
              
        return J1
        
    def v(self,q,t):
        v=-0.1*self.vi(q)
        #print(v)
        return np.array([[v]])
    
    def gamma(self,q,dq,t):
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[3*self.i-3:3*self.i-3+2]
        rj = q[3*self.j-3:3*self.j-3+2]
        dri = dq[3*self.i-3:3*self.i-3+2]
        drj = dq[3*self.j-3:3*self.j-3+2]
        wi = dq[3*self.i-3+2,0]
        
        gamma = np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),rj-ri)*wi**2-2*np.dot(np.dot((self.v_i).T,Bi.T),drj-dri)*wi   
        #print(gamma)
        return gamma

class tdd12: 
    def __init__(self,i=0,S_i=np.zeros((2,1)),v_i=np.zeros((2,1)),j=0,S_j=np.zeros((2,1)),c_bar=0):
        self.i=i
        self.S_i=S_i
        self.v_i=v_i
        self.j=j
        self.S_j=S_j
        self.c_bar=c_bar
    
    def A(self,q,i): #计算坐标转换矩阵
        phi_i = q[3*i+2-3,0]
        A_i = np.array([[np.cos(phi_i),-np.sin(phi_i)],[np.sin(phi_i),np.cos(phi_i)]])
        return A_i
    
    def Aij(self,q):
        #print(np.dot((self.A(q,self.i)).T,self.A(q,self.j)))
        return np.dot((self.A(q,self.i)).T,self.A(q,self.j))
    
    def vi(self,q):
        vi = np.dot(self.A(q,self.i),self.v_i)
        v1 = vi[0,0]
        v2 = vi[1,0]
        vi1 = np.sqrt(v1**2+v2**2)
        return vi1
    
    def constraints(self,q,t):
        '''
        :param q: [[x],[y],[theta]]
        :return: phi
        '''
        ri = q[3*self.i-3:3*self.i+2-3]
        rj = q[3*self.j-3:3*self.j+2-3]
        #print(self.v_i)
        #print(self.Aij(q))
        #print(self.S_j)
        #A=np.dot((self.v_i).T,self.Aij(q))
        #print(np.dot(A,self.S_j))
        
        phi = np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),rj-ri)+np.dot((self.v_i).T,np.dot(self.Aij(q),self.S_j))-np.dot((self.v_i).T,self.S_i)-self.vi(q)*(0.65-0.12*t) #c_bar=0.65-0.1*t
        #phi = np.dot(np.dot((self.v_i).T,self.Aij(q)),self.S_j)
        #print(self.vi(q))
        print(phi)
        return np.array(phi)
    
    def Jacobi(self,q,t,l):
        J = np.zeros((1,6))
        
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[0:2]
        rj = q[0:2]
        
        J[0,0:2]=-np.dot((self.v_i).T,(self.A(q,self.i)).T)
        J[0,2] = np.dot(np.dot((self.v_i).T,Bi.T),rj-ri)
        J[0,3:5]=np.dot((self.v_i).T,(self.A(q,self.i)).T)
        J[0,5]=0
        
        J1 = np.zeros((1,l))
        J1[0,3*self.i-3:3*self.i+2-3]=J[0,0:2]
        J1[0,3*self.i-3+2]=J[0,2]
        J1[0,3*self.j-3:3*self.j+2-3]=J[0,3:5]
        J1[0,3*self.j+2-3]=J[0,5]
              
        return J1
        
    def v(self,q,t):
        v=-0.12*self.vi(q)
        #print(v)
        return np.array([[v]])
    
    def gamma(self,q,dq,t):
        Bi = np.dot(R,self.A(q,self.i))
        ri = q[3*self.i-3:3*self.i-3+2]
        rj = q[3*self.j-3:3*self.j-3+2]
        dri = dq[3*self.i-3:3*self.i-3+2]
        drj = dq[3*self.j-3:3*self.j-3+2]
        wi = dq[3*self.i-3+2,0]
        
        gamma = np.dot(np.dot((self.v_i).T,(self.A(q,self.i)).T),rj-ri)*wi**2-2*np.dot(np.dot((self.v_i).T,Bi.T),drj-dri)*wi   
        #print(gamma)
        return gamma