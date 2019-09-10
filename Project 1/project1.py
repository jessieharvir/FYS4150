import numpy as np 
import matplotlib.pyplot as plt

def b(n):
    x = np.genfromtxt('prob1_b_n{}.txt'.format(n), usecols=0)
    u = np.genfromtxt('prob1_b_n{}.txt'.format(n), usecols=1)
    v = np.genfromtxt('prob1_b_n{}.txt'.format(n), usecols=2)

    #plot of analytical
    plt.plot(x, u, label='Analytical solution')
    #plot of numerical
    plt.plot(x, v, label = 'Numerical solution')
    plt.title('Comparison of analytical and numerical solution of Poisson Equation, with n = {}'.format(len(x)))
    plt.xlabel('x-values')
    plt.legend(loc='best')
    plt.ylabel('Values of solutions')
    return x, u, v


def c(n):
    x = np.genfromtxt('prob1_c_n{}.txt'.format(n), usecols=0)
    u = np.genfromtxt('prob1_c_n{}.txt'.format(n), usecols=1)
    v = np.genfromtxt('prob1_c_n{}.txt'.format(n), usecols=2)

    #plot of analytical
    plt.plot(x, u, label='Analytical solution')
    #plot of numerical
    plt.plot(x, v, label = 'Numerical solution')
    plt.title('Comparison of analytical and numerical solution of Poisson Equation, with n = {}'.format(len(x)))
    plt.xlabel('x-values')
    plt.legend(loc='best')
    plt.ylabel('Values of solutions v(x)')
    return x, u, v

n_values = [10, 100, 1000, 10000, 100000, 1000000, 10000000]
#x, u, v = c(n_values[2])

def relative_error(u, v):
    return (np.log10(np.abs((v-u)/u)))

def relative_error_many():
    error_max = []
    h_lst = []
    for i in range(len(n_values)):
        error = np.genfromtxt('prob1_c_n{}.txt'.format(n_values[i]), usecols=3)
        error_max.append(np.max(error))
        h = 1/(n_values[i] + 1)
        h_lst.append(h)
    error_max_array = np.asarray(error_max)
    h_array = np.asarray(np.log10(h_lst))
    plt.xlabel('log10(h)')
    plt.ylabel('max(log10 error)')
    plt.title('Plot of maximal error and step size')
    #plt.yscale('log')
    #plt.xscale('log')
    plt.plot(h_array, error_max_array, 'ro')

def plot_mange():
    for i in range(len(n_values)):
        x = np.genfromtxt('prob1_c_n{}.txt'.format(n_values[i]), usecols=0)
        u = np.genfromtxt('prob1_c_n{}.txt'.format(n_values[i]), usecols=1)
        v = np.genfromtxt('prob1_c_n{}.txt'.format(n_values[i]), usecols=2)
        plt.plot(x, v, label = 'Numerical solution, n = {}'.format(n_values[i]))
    
    plt.plot(x, u, label='Analytical solution')
    plt.title('Comparison of analytical and numerical solution of Poisson Equation')
    plt.xlabel('x-values')
    plt.legend(loc='best')
    plt.ylabel('Values of solutions v(x)')


#plot_mange()
relative_error_many()
plt.show()
