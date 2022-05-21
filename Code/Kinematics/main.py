import numpy as np
import matplotlib.pyplot as plt
from kinematics.Kinematics import Kinematics
from kinematics.Animate import Animate
''' import sympy as sym '''

#t = sym.symbols('t') #定义时间变量t为符号变量

if __name__ == "__main__":

    #表格参数输入
    '''
    para:[type, j, S_j, v_j, i, S_i, v_i, c, c_bar, C]
    '''
    
    ''' paras8=[\
    ['ax', 0, np.array([[0],[0]]), np.zeros((2,1)), 1, np.zeros((2,1)), np.zeros((2,1)), 0, 0], \
    ['ay', 0, np.array([[0], [0]]), np.zeros((2, 1)), 1, np.zeros((2, 1)), np.zeros((2, 1)), 0, 0], \
    ['r', 2, np.zeros((2,1)), np.zeros((2, 1)), 1, np.array([[0.25],[0]]), np.zeros((2, 1)), 0, 0], \
    ['t',2,np.zeros((2,1)),np.array([[1],[0]]),3,np.zeros((2,1)),np.array([[1],[0]]), 0, 0], \
    ['ax', 0, np.array([[0],[0]]), np.zeros((2,1)), 3, np.zeros((2,1)), np.zeros((2,1)), 0.6, 0], \
    ['ay', 0, np.array([[0],[0]]), np.zeros((2,1)), 3, np.zeros((2,1)), np.zeros((2,1)), 0, 0], \
    ['tdd8',2,np.zeros((2,1)),np.array([[1],[0]]),3,np.zeros((2,1)),np.zeros((2,1)), 0, 0] \
    ] '''

    paras10=[\
    ['ax', 0, np.array([[0],[0]]), np.zeros((2,1)), 1, np.zeros((2,1)), np.zeros((2,1)), 0, 0], \
    ['ay', 0, np.array([[0], [0]]), np.zeros((2, 1)), 1, np.zeros((2, 1)), np.zeros((2, 1)), 0, 0], \
    ['r', 2, np.zeros((2,1)), np.zeros((2, 1)), 1, np.array([[0.25],[0]]), np.zeros((2, 1)), 0, 0], \
    ['t',2,np.zeros((2,1)),np.array([[1],[0]]),3,np.zeros((2,1)),np.array([[1],[0]]), 0, 0], \
    ['ax', 0, np.array([[0],[0]]), np.zeros((2,1)), 3, np.zeros((2,1)), np.zeros((2,1)), 0.6, 0], \
    ['ay', 0, np.array([[0],[0]]), np.zeros((2,1)), 3, np.zeros((2,1)), np.zeros((2,1)), 0, 0], \
    ['tdd',2,np.zeros((2,1)),np.array([[1],[0]]),3,np.zeros((2,1)),np.zeros((2,1)), 0, 0] \
    ]

    ''' paras12=[\
    ['ax', 0, np.array([[0],[0]]), np.zeros((2,1)), 1, np.zeros((2,1)), np.zeros((2,1)), 0, 0], \
    ['ay', 0, np.array([[0], [0]]), np.zeros((2, 1)), 1, np.zeros((2, 1)), np.zeros((2, 1)), 0, 0], \
    ['r', 2, np.zeros((2,1)), np.zeros((2, 1)), 1, np.array([[0.25],[0]]), np.zeros((2, 1)), 0, 0], \
    ['t',2,np.zeros((2,1)),np.array([[1],[0]]),3,np.zeros((2,1)),np.array([[1],[0]]), 0, 0], \
    ['ax', 0, np.array([[0],[0]]), np.zeros((2,1)), 3, np.zeros((2,1)), np.zeros((2,1)), 0.6, 0], \
    ['ay', 0, np.array([[0],[0]]), np.zeros((2,1)), 3, np.zeros((2,1)), np.zeros((2,1)), 0, 0], \
    ['tdd12',2,np.zeros((2,1)),np.array([[1],[0]]),3,np.zeros((2,1)),np.zeros((2,1)), 0, 0] \
    ] '''

    ''' K8 = Kinematics(paras8) '''
    K10 = Kinematics(paras10)
    ''' K12 = Kinematics(paras12) '''
    
    '''
    q: [x1 y1 phi1 x2 y2 phi2 x3 y3 phi3]
    '''
    
    q0 = np.array([[0],[0],[np.pi/2],[0],[0.25],[2*np.pi-np.arctan(5/12)],[0.6],[0],[2*np.pi-np.arctan(5/12)]])
    
    
    ''' T8, Y8 = K8.solve(q0) '''
    T10, Y10 = K10.solve(q0)
    ''' T12, Y12 = K12.solve(q0) '''
    
    #绘制物体2的运动特性图
    ''' K8.plot() '''
    K10.plot()
    ''' K12.plot() '''
    #K.plotv()
    #K.plota()
    
    
    ''' #设置图的样式
    plt.style.use('seaborn-paper')
    #设置画布大小
    plt.figure(figsize=(7, 5),dpi=120)
    plt.plot(T8, Y8[...,21],linewidth=1.75,label='$d_1(t)$')
    plt.plot(T10, Y10[...,21],linewidth=1.75,label='$d_2(t)$')
    plt.plot(T12, Y12[...,21],linewidth=1.75,label='$d_3(t)$')
    #plt.ylim(-1,5)
    plt.xlabel("$t(s)$",fontsize=10,loc='right')
    plt.ylabel("$a_{2x}(m/s^2)$",fontsize=10)
    plt.grid()
    plt.legend()
    plt.show() '''
    
    q = Y10[:,0:9]

    #参数1表示要保存动图
    Animate(q,1)