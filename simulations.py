import control
import matplotlib.pyplot as plt
from sympy import N, I, Symbol, solve, re, expand, sqrt, sin, exp
import sympy as sym

###Transfer Function Calculation

m = 0.462
g = 9.81
d = 0.42
delta = 0.65
R = 2200
L0 = 125 * 10 ** -3
L1 = 24.1 * 10 ** -3
alpha = 1.2
c = 6.811
k = 1885
b = 10.4
phi = 41
x1 = 0.4


def remove_imaginary():
    if str(N(G)).__contains__('*I'):
        tf = str(N(G)).replace('*I', '')
    else:
        tf = str(N(G))
    return tf


G, s, A, B, C, D, E, F = sym.symbols('G s A B C D E F')
G = (F * C) / expand(((s ** 2 - A - B * s) * (s - E)))
print(f'G = {G}')

x3 = sqrt((k * x1 - k * d - m * g * sin(phi)) * (delta - x1) ** 2 / c)
A = 5 * (2 * c * x3 ** 2 / (delta - x1) ** 3 - k) / (7 * m)
B = -5 * b / (7 * m)
C = 10 * c * x3 / (7 * m * (delta - x1) ** 2)
E = -R / (L0 + L1 * exp(-alpha * (delta - x1)))
F = 1 / (L0 + L1 * exp(-alpha * (delta - x1)))
G = C*F/(A*E - A*s + B*E*s - B*s**2 - E*s**2 + s**3)
tf = remove_imaginary()
print(f'Transfer Function: G = {tf}')

###Sumulations###

s = control.TransferFunction.s
G = 1374.11696533797/(s**3 + 15416.4482468152*s**2 + 250996.733474466*s + 51926213.4568898)

t_imp, x_imp = control.impulse_response(G)
plt.grid()
plt.plot(t_imp, x_imp)
plt.title('Impulse Response Plot')
plt.xlabel("Time (s)")
plt.ylabel("Distance (m)")
plt.show()

t_stp, x_stp = control.step_response(G)
plt.grid()
plt.plot(t_stp, x_stp)
plt.title('Step Response Plot')
plt.xlabel("Time (s)")
plt.ylabel("Distance (m)")
plt.show()

Kp = 1200
Kd = 20
Ki = 0

controller = Kp + Ki/s + Kd*s
Gfin = control.feedback(G, controller)
t_imp, g_imp = control.impulse_response(Gfin)
t_imp1, g_imp1 = control.step_response(Gfin)

plt.plot(t_imp, g_imp*1000)
plt.plot(t_imp1, g_imp1*1000)
plt.grid()
plt.title('PID controller with both step and impulse')
plt.legend(['Impulse Respnse', 'Step Response'])
plt.xlabel("Time (s)")
plt.ylabel("Distance (mm)")
plt.show()

