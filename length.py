import pandas
from math import sin, pi, cos
from matplotlib.pyplot import plot, show

degrees = 15
alpha = degrees * pi / 180

R = 0.0032
L = [0.0064]
for i in range(38):
    L.append(0.0005*(i+13))

Lb = 0.9
m = 0.0010355
sigma = 0.21

M = []
for i in range(len(L)):
    M.append(sigma * (2*pi*R*(L[i] - 2*R[i]) + 4*pi*R**2))

meu_1 = 0.03
meu_2 = 0.3


H = Lb * sin(alpha)
B = H * sin(1/2 * pi - alpha) / sin(alpha)

I_ball = 2/5 * m * R ** 2
g = 9.8

dt = 0.000001   # 1 micro sec (time step size)

"""
A_cyl = 2 * pi * R * (L - 2 * R)
A_ss = 2 * pi * R ** 2
sigma = M / (A_cyl + 2 * A_ss)
I_cyl = sigma * A_cyl * (1/4 * R ** 2 + 1/3 * (L - 2 * R)**2)
I_lss = sigma * A_ss * (1/12 * R ** 2 + (1/2 * R + (L - 2 * R))**2)
I_rss = sigma * A_ss * (1/3 * R ** 2)
I_cap = I_cyl + I_lss + I_rss
I_sys = I_ball + I_cap
d = (M * L/2 + m * R)/(M + m)
"""

A_cyl = []
for i in range(len(L)):
    A_cyl.append(2 * pi * R * (L[i] - 2 * R))
A_ss = 2 * pi * R ** 2
I_cyl = []
I_lss = []
for i in range(len(L)):
    I_cyl = sigma * A_cyl[i] * (1 / 4 * R ** 2 + 1 / 3 * (L[i] - 2 * R)**2)
    I_lss = sigma * A_ss * (1 / 12 * R ** 2 + (1/2 * R + (L[i] - 2 * R)) ** 2)
I_rss = sigma * A_ss * (1/3 * R ** 2)
I_cap = I_cyl + I_lss + I_rss
I_sys = I_ball + I_cap
d = []
for i in range(len(L)):
    d.append((M[i] * L[i] / 2 + m * R)/(M[i] + m))

time = [0] * len(L)

theta1 = [0] * len(L)
theta1_dot = [0] * len(L)
ref_pos = [0] * len(L)
ball_pos = [0] * len(L)
theta2_dot = [0] * len(L)
theta2 = [0] * len(L)
E1prime = [0] * len(L)
E2prime = [0] * len(L)
meu_1alpha = [0] * len(L)
meu_2alpha = [0] * len(L)

for i in range(len(L)):
    # print("i")
    while True:
        # print("two")
        meu_1alpha[i] = meu_1 / sin(alpha)
        meu_2alpha[i] = meu_2 / sin(alpha)
        if (ball_pos[i] - ref_pos[i]) >= L[i] - 2 * R:
            E1prime[i] = 1 / 2 * (I_ball + m * R ** 2) * theta1_dot[i] ** 2 + m * g * (
                    H - R * sin(alpha) - R * theta1[i] * sin(alpha)) + M[i] * g * (
                              H - L[i] / 2 * sin(alpha)) - meu_1alpha[i] * cos(alpha) * m * g * R * theta1[i]
            H = H - (L[i] - 2 * R) * sin(alpha)
            theta2_dot_square = (E1prime[i] - (m + M[i]) * g * (H - R * sin(alpha)) - M[i] * g * (
                    1 / 2 * (L[i] - 2 * R) * sin(alpha))) / (1 / 2 * (I_sys + (m + M[i]) * R ** 2))
            theta2_dot[i] = theta2_dot_square ** 1 / 2
            while 0 <= theta2[i] < pi:
                dtheta2 = theta2[i]
                theta2[i] += theta2_dot[i] * dt
                if theta2[i] > pi:
                    theta2[i] = pi
                dtheta2 = theta2[i] - dtheta2
                theta2_dot[i] += (((m + M[i]) * g * R * sin(alpha) - M[i] * g * (
                            1 / 2 * (L[i] - 2 * R) * cos(alpha + theta2[i])) - meu_2alpha[i] * cos(alpha) * (
                                               m + M[i]) * g * R) / (I_sys + (m + M[i]) * R ** 2)) * dt
                time[i] += dt
                ball_pos[i] += R * dtheta2
                if ball_pos[i] >= Lb:
                    print(time[i])
                    print(L[i])
                    print(degrees)
                    print(meu_1)
                    print(meu_2)
                    break
            if ball_pos[i] >= Lb:
                print(time[i])
                print(L[i])
                print(degrees)
                print(meu_1)
                print(meu_2)
                break
            theta1[i] = 0
            H = H - R * pi * sin(alpha)
            """
            E2prime[i] = 1 / 2 * (I_sys + (m + M) * R ** 2) * theta2_dot[i] ** 2 + (m + M) * g * (
                        H[i] - R * sin(alpha[i]) - (m + M) * g * R * pi * sin(alpha[i]) + M * g * 1 / 2 * (
                            L - 2 * R) * (sin(alpha[i] + pi))) - meu_2 * (m + M) * R * pi
            theta1_dot_square = (E2prime[i] - m * g * (H[i] - R * sin(alpha[i])) - M * g * (
                    H[i] - L / 2 * sin(alpha[i]))) / (1 / 2 * (I_ball + m * R ** 2))
            theta1_dot[i] = theta1_dot_square ** 1 / 2
            or
            theta1_dot[i] = theta2_dot[i]
            or
            theta1_dot[i] = 0
            ?   ~ all above mentioned assumptions are wrong
            """
            ref_pos[i] = ball_pos[i]

        if ball_pos[i] >= Lb:
            print(time[i])
            print(L[i])
            print(degrees)
            print(meu_1)
            print(meu_2)
            break
        # print("one")
        dtheta1 = theta1[i]
        theta1[i] += theta1_dot[i] * dt
        dtheta1 = theta1[i] - dtheta1
        # ref_theta1_dot = theta1_dot[i]
        theta1_dot[i] += (m * g * R * (sin(alpha) - meu_1alpha[i] * cos(alpha))) / (I_ball + m * R ** 2) * dt
        ball_pos[i] += R * dtheta1
        time[i] += dt

print("degrees = ", degrees)
print("L = ", L)
for i in range(len(L)):
    time[i] = round(time[i], 4)
print("time = ", time)
plot(L, time)
show()

df = pandas.DataFrame(time)
df.to_excel('length15degreef.xlsx')
