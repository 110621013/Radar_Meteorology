# angle=1.5, observational range 25km, C band radar with 5 cm wavelength, PRF=1000
#a 繪製一個圖來顯示徑向風 V r （在垂直坐標上）相對於方位角（水平坐標）坐標的分佈。
#b 展開徑向風數據，繪製它，然後使用 VAD 方法計算 u ， v ，水平收斂/發散，剪切變形和拉伸變形。 （假設垂直速度 w 可以忽略）。
#c 使用原始折疊徑向風數據和梯度 VAD 技術計算 u、v，並將結果與 (b) 的結果進行比較

import numpy as np
import matplotlib.pyplot as plt
import math
pi = math.pi
ts = 0.001
lamb = 5*10**-2
max_vr = lamb/(4*ts)
alpha = 2*pi*1.5/360
r = 25000
w0 = 0

# read fort.18 data and return vr nparray(0~359)
def read_fort():
    with open('fort.18') as f: #, encoding='utf-8'
        lines = f.readlines()
    vr = []
    for line in lines:
        vr.append(float(line[-11:-1]))
    return np.array(vr)

def a_plot_vr(vr):
    plt.plot(vr)
    plt.title('VR-angle')
    plt.xlabel('angel')
    plt.ylabel('VR(m/s)')
    plt.savefig('hw6-a')
    #plt.show()
    plt.close()

class fourier():
    def __init__(self, data):
        self.data = data
        self.N = len(data)
    def cal_a0(self):
        total = 0.0
        for i in range(self.N):
            total += self.data[i]
        return total/self.N
    def cal_am(self, m):
        assert m >= 1
        assert m <= self.N/2 - 1
        total = 0.0
        for i in range(self.N):
            total += self.data[i]*np.cos(2*np.pi*m*i/self.N)
        return total*2/self.N
    def cal_bm(self, m):
        assert m >= 1
        assert m <= self.N/2 - 1
        total = 0.0
        for i in range(self.N):
            total += self.data[i]*np.sin(2*np.pi*m*i/self.N)
        return total*2/self.N

    def five_nums(self):
        return self.cal_a0(), self.cal_am(1), self.cal_bm(1), self.cal_am(2), self.cal_bm(2)

def b_unfold_plot_VAD_cal_var(vr): # assume vr[0] no folded, and fold once
    print('--b gogo--')
    vr = np.array(vr, copy=True)
    # unfold
    for i in range(1, len(vr)):
        if vr[i]-vr[i-1] > max_vr:
            vr[i] -= 2*max_vr
        elif vr[i]-vr[i-1] < -5.0:
            vr[i] += 2*max_vr
    
    # plot
    plt.plot(vr)
    plt.title('VR-angle unfolded')
    plt.xlabel('angel')
    plt.ylabel('unfolded VR(m/s)')
    plt.savefig('hw6-b')
    #plt.show()
    plt.close()

    # cal fourier
    a0, a1, b1, a2, b2 = fourier(vr).five_nums()
    print('fourier a0, a1, b1, a2, b2:', a0, a1, b1, a2, b2)
    # cal parameter
    u0 = b1/math.cos(alpha)
    v0 = a1/math.cos(alpha)
    div = 2*(a0-w0*math.sin(alpha))/(r*math.cos(alpha)**2)
    stretching_deformation = -2*a2/(r*math.cos(alpha)**2)
    shear_deformation = 2*b2/(r*math.cos(alpha)**2)
    print('u0, v0, div, stretching_deformation, shear_deformation:', u0, v0, div, stretching_deformation, shear_deformation)

def c_gradient_VAD_cal_uv(vr):
    vr = np.array(vr, copy=True)
    print('--c gogo--')
    # no unfold dvr/dbeta
    print('---> no unfold dvr/dbeta')
    dvr = [vr[i+1]-vr[i] for i in range(len(vr)-1)]
    plt.plot(dvr)
    plt.title('no unfold dvr/dbeta')
    plt.savefig('hw6-c-1')
    plt.close()
    # cal fourier
    _, a1, b1, _, _ = fourier(dvr).five_nums()
    print('fourier a1, b1:', a1, b1)
    u0 = a1/math.cos(alpha) * (len(dvr)/2/pi)
    v0 = -b1/math.cos(alpha) * (len(dvr)/2/pi)
    print('u0, v0:', u0, v0)
    
    # unfolded dvr/dbeta
    print('---> unfolded dvr/dbeta')
    for i in range(len(dvr)):
        if dvr[i] > max_vr:
            dvr[i] -= 2*max_vr
        elif dvr[i] < -max_vr:
            dvr[i] += 2*max_vr
    plt.plot(dvr)
    plt.title('unfolded dvr/dbeta')
    plt.savefig('hw6-c-2')
    plt.close()
    # cal fourier
    _, a1, b1, _, _ = fourier(dvr).five_nums()
    print('fourier a1, b1:', a1, b1)
    u0 = a1/math.cos(alpha) * ((len(dvr))/2/pi)
    v0 = -b1/math.cos(alpha) * (len(dvr)/2/pi)
    print('u0, v0:', u0, v0)

def RMSE_beta_sigma(way1, way2, beta_nums):
    if way1=='1' and way2=='1':
        return beta_nums

    total = 0.0
    if way1=='1':
        if way2=='cos':
            for beta in range(beta_nums):
                total += math.cos(beta*pi/180)
        if way2=='sin':
            for beta in range(beta_nums):
                total += math.sin(beta*pi/180)
        if way2=='cos2':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)
        if way2=='sin2':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)
    if way1=='cos':
        if way2=='1':
            for beta in range(beta_nums):
                total += math.cos(beta*pi/180)
        if way2=='cos':
            for beta in range(beta_nums):
                total += math.cos(beta*pi/180)**2
        if way2=='sin':
            for beta in range(beta_nums):
                total += math.sin(beta*pi/180)*math.cos(beta*pi/180)
        if way2=='cos2':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)*math.cos(beta*pi/180)
        if way2=='sin2':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)*math.cos(beta*pi/180)
    if way1=='sin':
        if way2=='1':
            for beta in range(beta_nums):
                total += math.sin(beta*pi/180)
        if way2=='cos':
            for beta in range(beta_nums):
                total += math.sin(beta*pi/180)*math.cos(beta*pi/180)
        if way2=='sin':
            for beta in range(beta_nums):
                total += math.sin(beta*pi/180)**2
        if way2=='cos2':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)*math.sin(beta*pi/180)
        if way2=='sin2':
            for beta in range(beta_nums):
                total += math.sin(2*beta*pi/180)*math.sin(beta*pi/180)
    if way1=='cos2':
        if way2=='1':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)
        if way2=='cos':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)*math.cos(beta*pi/180)
        if way2=='sin':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)*math.sin(beta*pi/180)
        if way2=='cos2':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)**2
        if way2=='sin2':
            for beta in range(beta_nums):
                total += math.cos(2*beta*pi/180)*math.sin(2*beta*pi/180)
    if way1=='sin2':
        if way2=='1':
            for beta in range(beta_nums):
                total += math.sin(2*beta*pi/180)
        if way2=='cos':
            for beta in range(beta_nums):
                total += math.sin(2*beta*pi/180)*math.cos(beta*pi/180)
        if way2=='sin':
            for beta in range(beta_nums):
                total += math.sin(2*beta*pi/180)*math.sin(beta*pi/180)
        if way2=='cos2':
            for beta in range(beta_nums):
                total += math.sin(2*beta*pi/180)*math.cos(2*beta*pi/180)
        if way2=='sin2':
            for beta in range(beta_nums):
                total += math.sin(2*beta*pi/180)**2
    return total

def RMSE_beta_vr_sigma(way, vr, beta_nums):
    total = 0.0
    if way=='1':
        for beta in range(beta_nums):
            total += vr[beta]
    if way=='cos':
        for beta in range(beta_nums):
            total += vr[beta]*math.cos(beta*pi/180)
    if way=='sin':
        for beta in range(beta_nums):
            total += vr[beta]*math.sin(beta*pi/180)
    if way=='cos2':
        for beta in range(beta_nums):
            total += vr[beta]*math.cos(2*beta*pi/180)
    if way=='sin2':
        for beta in range(beta_nums):
            total += vr[beta]*math.sin(2*beta*pi/180)
    return total

def new_b_RMSE_VAD(vr): 
    print('--new b gogo--')
    vr = np.array(vr, copy=True)
    # unfold
    for i in range(1, len(vr)):
        if vr[i]-vr[i-1] > max_vr:
            vr[i] -= 2*max_vr
        elif vr[i]-vr[i-1] < -5.0:
            vr[i] += 2*max_vr

    # 建立係數矩陣
    five_nums = np.zeros((5, 5), dtype=float)
    ans = np.zeros((5), dtype=float)
    way = ['1', 'cos', 'sin', 'cos2', 'sin2']
    for i in range(5):
        for j in range(5):
            five_nums[i, j] = RMSE_beta_sigma(way[i], way[j], 360)
        ans[i] = RMSE_beta_vr_sigma(way[i], vr, 360)
    #print('init', five_nums, ans)

    for i in range(4):
        for j in range(i+1, 5):
            scale = five_nums[j,i]/five_nums[i,i]
            five_nums[j] -= five_nums[i]*scale
            ans[j] -= ans[i]*scale
    #print('gause', five_nums, ans)

    for i in range(4, -1, -1):
        ans[i] /= five_nums[i,i]
        for j in range(0, i):
            ans[j] -= five_nums[j,i]*ans[i]
    #print('ok',five_nums, ans)

    print('RMSE a0, a1, b1, a2, b2:', ans)
    # cal parameter
    u0 = ans[2]/math.cos(alpha)
    v0 = ans[1]/math.cos(alpha)
    div = 2*(ans[0]-w0*math.sin(alpha))/(r*math.cos(alpha)**2)
    stretching_deformation = -2*ans[3]/(r*math.cos(alpha)**2)
    shear_deformation = 2*ans[4]/(r*math.cos(alpha)**2)
    print('u0, v0, div, stretching_deformation, shear_deformation:', u0, v0, div, stretching_deformation, shear_deformation)

def new_c_RMSE_gradient_VAD(vr):
    vr = np.array(vr, copy=True)
    print('--new c gogo--')
    # no unfold dvr/dbeta
    print('---> no unfold dvr/dbeta')
    dvr = np.array([vr[i+1]-vr[i] for i in range(len(vr)-1)])
    plt.plot(dvr)
    plt.title('newc folded dvr/dbeta')
    plt.savefig('hw6-newc-1')
    plt.close()
    
    # 建立係數矩陣
    five_nums = np.zeros((5, 5), dtype=float)
    ans = np.zeros((5), dtype=float)
    way = ['1', 'cos', 'sin', 'cos2', 'sin2']
    for i in range(5):
        for j in range(5):
            five_nums[i, j] = RMSE_beta_sigma(way[i], way[j], 359)
        ans[i] = RMSE_beta_vr_sigma(way[i], dvr, 359)
    #print('init', five_nums, ans)

    for i in range(4):
        for j in range(i+1, 5):
            scale = five_nums[j,i]/five_nums[i,i]
            five_nums[j] -= five_nums[i]*scale
            ans[j] -= ans[i]*scale
    #print('gause', five_nums, ans)

    for i in range(4, -1, -1):
        ans[i] /= five_nums[i,i]
        for j in range(0, i):
            ans[j] -= five_nums[j,i]*ans[i]
    #print('ok',five_nums, ans)

    print('RMSE a0, a1, b1, a2, b2:', ans)
    u0 = ans[1]/math.cos(alpha) * (dvr.shape[0]/2/pi)
    v0 = -ans[2]/math.cos(alpha) * (dvr.shape[0]/2/pi)
    shear_deformation = ans[3]/math.cos(alpha)**2/r  * (dvr.shape[0]/2/pi)
    stretching_deformation = -ans[4]/math.cos(alpha)**2/r  * (dvr.shape[0]/2/pi)
    print('u0, v0, shear_deformation, stretching_deformation:', u0, v0, shear_deformation, stretching_deformation)
    

    # unfolded dvr/dbeta
    print('---> unfolded dvr/dbeta')
    for i in range(dvr.shape[0]):
        if dvr[i] > max_vr:
            dvr[i] -= 2*max_vr
        elif dvr[i] < -max_vr:
            dvr[i] += 2*max_vr
    plt.plot(dvr)
    plt.title('newc unfolded dvr/dbeta')
    plt.savefig('hw6-newc-2')
    plt.close()

    # 建立係數矩陣
    five_nums = np.zeros((5, 5), dtype=float)
    ans = np.zeros((5), dtype=float)
    way = ['1', 'cos', 'sin', 'cos2', 'sin2']
    for i in range(5):
        for j in range(5):
            five_nums[i, j] = RMSE_beta_sigma(way[i], way[j], 359)
        ans[i] = RMSE_beta_vr_sigma(way[i], dvr, 359)
    #print('init', five_nums, ans)

    for i in range(4):
        for j in range(i+1, 5):
            scale = five_nums[j,i]/five_nums[i,i]
            five_nums[j] -= five_nums[i]*scale
            ans[j] -= ans[i]*scale
    #print('gause', five_nums, ans)

    for i in range(4, -1, -1):
        ans[i] /= five_nums[i,i]
        for j in range(0, i):
            ans[j] -= five_nums[j,i]*ans[i]
    #print('ok',five_nums, ans)

    print('RMSE a0, a1, b1, a2, b2:', ans)
    u0 = ans[1]/math.cos(alpha) * (dvr.shape[0]/2/pi)
    v0 = -ans[2]/math.cos(alpha) * (dvr.shape[0]/2/pi)
    shear_deformation = ans[3]/math.cos(alpha)**2/r  * (dvr.shape[0]/2/pi)
    stretching_deformation = -ans[4]/math.cos(alpha)**2/r  * (dvr.shape[0]/2/pi)
    print('u0, v0, shear_deformation, stretching_deformation:', u0, v0, shear_deformation, stretching_deformation)

if __name__ == '__main__':
    vr = read_fort()
    #a_plot_vr(vr)
    #b_unfold_plot_VAD_cal_var(vr)
    #c_gradient_VAD_cal_uv(vr)
    new_b_RMSE_VAD(vr)
    new_c_RMSE_gradient_VAD(vr)