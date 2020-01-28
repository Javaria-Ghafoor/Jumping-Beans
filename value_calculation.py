import pandas
from math import sin, pi, cos
from matplotlib.pyplot import plot, show

alpha = [0.1960, 0.2218, 0.2474, 0.2733, 0.2994, 0.3256, 0.3522, 0.3791, 0.4061, 0.4334, 0.4611, 0.4892, 0.5178, 0.5468]

L = 0.0126       # length of capsule
R = 0.0024     # radius of ball bearing
Lb = 0.9         # length of board
m = 0.0004541    # mass of ball
M = 0.0000624    # mass of capsule
meu_1 = 0.032   # put meu_1 and meu_2 = 0 for no friction / ideal case solution
meu_2 = 0.34

H = []
for i in range(len(alpha)):
    H.append(Lb * sin(alpha[i]))

B = []
for i in range(len(alpha)):
    B.append(H[i] * sin(1/2 * pi - alpha[i]) / sin(alpha[i]))

I_ball = 2/5 * m * R ** 2
g = 9.8

dt = 0.000001   # 1 micro sec (time step size)

A_cyl = 2 * pi * R * (L - 2 * R)
A_ss = 2 * pi * R ** 2
sigma = M / (A_cyl + 2 * A_ss)
I_cyl = sigma * A_cyl * (1/4 * R ** 2 + 1/3 * (L - 2 * R)**2)
I_lss = sigma * A_ss * (1/12 * R ** 2 + (1/2 * R + (L - 2 * R))**2)
I_rss = sigma * A_ss * (1/3 * R ** 2)
I_cap = I_cyl + I_lss + I_rss
I_sys = I_ball + I_cap
d = (M * L/2 + m * R)/(M + m)

time = [0] * len(alpha)

theta1 = [0] * len(alpha)
theta1_dot = [0] * len(alpha)
ref_pos = [0] * len(alpha)
ball_pos = [0] * len(alpha)
theta2_dot = [0] * len(alpha)
theta2 = [0] * len(alpha)
E1prime = [0] * len(alpha)
E2prime = [0] * len(alpha)
meu_1alpha = [0] * len(alpha)
meu_2alpha = [0] * len(alpha)
degrees = [0] * len(alpha)

# comment line 23, 24, 26, 27 and uncomment line 29, 74, 79  to use the set of radian angles in Mehreen's data

for i in range(len(alpha)):
    degrees[i] = alpha[i] * 180/pi
    # print("i")
    while True:
        # print("two")
        meu_1alpha[i] = meu_1 / sin(alpha[i])
        meu_2alpha[i] = meu_2 / sin(alpha[i])
        if (ball_pos[i] - ref_pos[i]) >= L - 2 * R:
            E1prime[i] = 1 / 2 * (I_ball + m * R ** 2) * theta1_dot[i] ** 2 + m * g * (
                    H[i] - R * sin(alpha[i]) - R * theta1[i] * sin(alpha[i])) + M * g * (
                              H[i] - L / 2 * sin(alpha[i])) - meu_1alpha[i] * cos(alpha[i]) * m * g * R * theta1[i]
            H[i] = H[i] - (L - 2 * R) * sin(alpha[i])
            theta2_dot_square = (E1prime[i] - (m + M) * g * (H[i] - R * sin(alpha[i])) - M * g * (
                    1 / 2 * (L - 2 * R) * sin(alpha[i]))) / (1/2 * (I_sys + (m + M) * R ** 2))
            theta2_dot[i] = theta2_dot_square ** 1 / 2
            while 0 <= theta2[i] < pi:
                dtheta2 = theta2[i]
                theta2[i] += theta2_dot[i] * dt
                if theta2[i] > pi:
                    theta2[i] = pi
                dtheta2 = theta2[i] - dtheta2
                theta2_dot[i] += (((m + M) * g * R * sin(alpha[i]) - M * g * (
                            1 / 2 * (L - 2 * R) * cos(alpha[i] + theta2[i])) - meu_2alpha[i] * cos(alpha[i]) * (
                                               m + M) * g * R) / (I_sys + (m + M) * R ** 2)) * dt
                time[i] += dt
                ball_pos[i] += R * dtheta2
                if ball_pos[i] >= Lb:
                    print(time[i])
                    print(alpha[i])
                    print(degrees[i])
                    print(meu_1)
                    print(meu_2)
                    break
            if ball_pos[i] >= Lb:
                print(time[i])
                print(alpha[i])
                print(degrees[i])
                print(meu_1)
                print(meu_2)
                break
            theta1[i] = 0
            theta2[i] = 0
            theta1_dot[i] = 0
            H[i] = H[i] - R * pi * sin(alpha[i])
            """
            if (alpha[i] < arctan(meu_2)):
                theta1_dot[i] = .9/(cos(alpha[i])**2) * theta1_dot[i]
            E2prime[i] = 1 / 2 * (I_sys + (m + M) * R ** 2) * theta2_dot[i] ** 2 + (m + M) * g * (
                        H[i] - R * sin(alpha[i]) - (m + M) * g * R * pi * sin(alpha[i]) + M * g * 1 / 2 * (
                            L - 2 * R) * (sin(alpha[i] + pi))) - meu_2 * (m + M) * R * pi
            theta1_dot_square = (E1prime[i] - E2prime[i] - m * g * (H[i] - R * sin(alph
            theta1_dot[i] = theta1_dot_square ** 1 / 2a[i])) - M * g * (
                    H[i] - L / 2 * sin(alpha[i]))) / (1 / 2 * (I_ball + m * R ** 2))
            theta1_dot[i] = theta2_dot[i]
            or
            theta1_dot[i] = 0
            ?   ~ all above mentioned assumptions are wrong
            """
            ref_pos[i] = ball_pos[i]

        if ball_pos[i] >= Lb:
            print(time[i])
            print(alpha[i])
            print(degrees[i])
            print(meu_1)
            print(meu_2)
            break
        # print("one")
        dtheta1 = theta1[i]
        theta1[i] += theta1_dot[i] * dt
        dtheta1 = theta1[i] - dtheta1
        # ref_theta1_dot = theta1_dot[i]
        theta1_dot[i] += (m * g * R * (sin(alpha[i]) - meu_1alpha[i] * cos(alpha[i]))) / (I_ball + m * R ** 2) * dt
        ball_pos[i] += R * dtheta1
        time[i] += dt

print("degrees = ", degrees)
print("alpha = ", alpha)
print("time = ", time)
for i in range(len(alpha)):
    time[i] = round(time[i], 4)
plot(degrees, time)
show()

df = pandas.DataFrame(time)
df.to_excel('theta1dot0.xlsx')
