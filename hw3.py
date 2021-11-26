import numpy as np
import matplotlib.pyplot as plt

n = [i for i in range(16)]
t = [100.00, 129.02, 127.14, 102.45, 89.27, 101.75, 117.30, 110.59, 89.41, 82.70, 98.25, 110.73, 97.55, 72.86, 70.98, 82.88]
N = 16

def q1_plot_t_n():
    plt.plot(n, t, label='T')
    plt.legend()
    plt.title('T-n plot')
    plt.xlabel('n')
    plt.ylabel('T')
    plt.show()
    plt.close()

def cal_a0():
    total = 0.0
    for i in range(N):
        total += t[i]
    return total/N
def cal_am(m):
    assert m >= 1
    assert m <= N/2 - 1
    total = 0.0
    for i in range(N):
        total += t[i]*np.cos(2*np.pi*m*i/N)
    return total*2/N
def cal_bm(m):
    assert m >= 1
    assert m <= N/2 - 1
    total = 0.0
    for i in range(N):
        total += t[i]*np.sin(2*np.pi*m*i/N)
    return total*2/N
def cal_ahalfN():
    total = 0.0
    for i in range(N):
        total += t[i]*np.cos(np.pi*i)
    return total/N
def cal_sm2(m):
    assert m >= 1
    assert m <= N/2 - 1
    am = cal_am(m)
    bm = cal_bm(m)
    return (am**2 + bm**2)/2
def cal_shalfN2():
    return (cal_ahalfN())**2
def cal_s2():
    total = 0.0
    for i in range(1, int(N/2)):
        total += cal_sm2(i)
    total += cal_shalfN2()
    return total

def q2_find_variance_contribution():
    sm2 = []
    for i in range(1, int(N/2)):
        sm2.append(cal_sm2(i))
    sm2.append(cal_shalfN2())
    print(sm2)

    s2 = cal_s2()
    contri = sm2[:]
    for i in range(len(sm2)):
        contri[i] /= s2
    print(contri)

def q3_plot_t_n_top5():
    # from q2, top 5:n=[3,1,2,4,5]
    a0 = cal_a0()
    a8 = cal_ahalfN()
    for num in [3, 1, 2, 4, 5]:
        am = cal_am(num)
        bm = cal_bm(num)
        new_t = []
        for i in range(N):
            new_t.append(a0+am*np.cos(2*np.pi*num*i/N)+bm*np.sin(2*np.pi*num*i/N)+a8/np.cos(np.pi*i))
        plt.plot(n, new_t, label='n={}'.format(str(num)))
    plt.legend()
    plt.title('new T-n plot')
    plt.xlabel('n')
    plt.ylabel('new T')
    plt.show()
    plt.close()

def q4_plot_t_n_top5_gradually():
    a0 = cal_a0()
    a8 = cal_ahalfN()
    superimpose_t = [a0+a8/np.cos(np.pi*i) for i in range(N)]
    label_s = ''
    for num in [3, 1, 2, 4, 5]:
        am = cal_am(num)
        bm = cal_bm(num)
        new_t = []
        # n eave create
        for i in range(N):
            new_t.append(am*np.cos(2*np.pi*num*i/N)+bm*np.sin(2*np.pi*num*i/N))
        
        # superimpose
        for i in range(N):
            superimpose_t[i] += new_t[i]
        
        label_s += '+'+str(num)
        plt.plot(n, superimpose_t, label='n={}'.format(label_s))
        plt.plot(n, t, label='original T')
        plt.legend()
        plt.title('superimpose T-n plot')
        plt.xlabel('n')
        plt.ylabel('superimpose T')
        plt.show()
        plt.close()

def q5_plot_t_n_all():
    a0 = cal_a0()
    a8 = cal_ahalfN()
    superimpose_t = [a0+a8/np.cos(np.pi*i) for i in range(N)]
    label_s = ''
    for num in range(1, 8):
        am = cal_am(num)
        bm = cal_bm(num)
        new_t = []
        # n eave create
        for i in range(N):
            new_t.append(am*np.cos(2*np.pi*num*i/N)+bm*np.sin(2*np.pi*num*i/N))
        
        # superimpose
        for i in range(N):
            superimpose_t[i] += new_t[i]
    delta_t = np.array(superimpose_t)-np.array(t)

    plt.plot(n, superimpose_t, label='n=all(1~8)')
    plt.plot(n, t, label='original T')
    plt.legend()
    plt.title('all-superimpose T-n plot')
    plt.xlabel('n')
    plt.ylabel('all-superimpose T')
    plt.show()
    plt.close()

    plt.plot(n, delta_t, label='delta t')
    plt.legend()
    plt.title('delta t-n plot')
    plt.xlabel('n')
    plt.ylabel('delta T')
    plt.show()
    plt.close()


if __name__ == '__main__':
    #q1_plot_t_n()
    #q2_find_variance_contribution()
    #q3_plot_t_n_top5()
    q4_plot_t_n_top5_gradually()
    #q5_plot_t_n_all()