import pandas
from math import sin, pi, cos
from matplotlib.pyplot import plot, show


def getMinMax(arr):
    n = len(arr)

    # mx index
    j = 0

    # mn index
    k = 0

    # If array has even number of elements then
    # initialize the first two elements as minimum
    # and maximum
    if (n % 2 == 0):
        mx = max(arr[0], arr[1])
        mn = min(arr[0], arr[1])

        # set the starting index for loop
        i = 2

    # If array has odd number of elements then
    # initialize the first element as minimum
    # and maximum
    else:
        mx = mn = arr[0]

        # set the starting index for loop
        i = 1

    # In the while loop, pick elements in pair and
    # compare the pair with max and min so far
    while (i < n - 1):
        if arr[i] < arr[i + 1]:
            mx = max(mx, arr[i + 1])
            if mx == arr[i+1]:
                j = i+1
            mn = min(mn, arr[i])
            if mn == arr[i]:
                k = i
        else:
            mx = max(mx, arr[i])
            if mx == arr[i]:
                j = i
            mn = min(mn, arr[i + 1])
            if mn == arr[i+1]:
                k = i+1

            # Increment the index by 2 as two
        # elements are processed in loop
        i += 2
    maximum = [j, mx]
    minimum = [k, mn]
    return maximum, minimum


def getMin(arr):
    n = len(arr)
    # mn index
    k = 0
    # If array has even number of elements then
    # initialize the first two elements as minimum
    if (n % 2 == 0):
        mn = min(arr[0][1], arr[1][1])

        # set the starting index for loop
        i = 2

    # If array has odd number of elements then
    # initialize the first element as minimum
    else:
        mn = arr[0][1]

        # set the starting index for loop
        i = 1
    # In the while loop, pick elements in pair and
    # compare the pair with max and min so far
    while (i < n - 1):
        if arr[i][1] < arr[i + 1][1]:
            mn = min(mn, arr[i][1])
            if mn == arr[i][1]:
                k = i
        else:
            mn = min(mn, arr[i + 1][1])
            if mn == arr[i+1][1]:
                k = i+1

            # Increment the index by 2 as two
        # elements are processed in loop
        i += 2
    minimum = [k, mn]
    return minimum


def getMax(arr):
    n = len(arr)
    # mx index
    j = 0
    # If array has even number of elements then
    # initialize the first two elements as maximum
    if (n % 2 == 0):
        mx = max(arr[0][1], arr[1][1])

        # set the starting index for loop
        i = 2

    # If array has odd number of elements then
    # initialize the first element as maximum
    else:
        mx = arr[0][1]

        # set the starting index for loop
        i = 1

    # In the while loop, pick elements in pair and
    # compare the pair with max and min so far
    while (i < n - 1):
        if arr[i][1] < arr[i + 1][1]:
            mx = max(mx, arr[i + 1][1])
            if mx == arr[i + 1][1]:
                j = i + 1
        else:
            mx = max(mx, arr[i][1])
            if mx == arr[i][1]:
                j = i

            # Increment the index by 2 as two
        # elements are processed in loop
        i += 2
    maximum = [j, mx]
    return maximum




degrees = 30
alpha = degrees * pi / 180

Lb = 0.9
rho = 7725    # density of ball
sigma = 0.21  # density of capsule
meu_1 = 0.03
meu_2 = 0.3
meu_1alpha = meu_1 / sin(alpha)
meu_2alpha = meu_2 / sin(alpha)
H = Lb * sin(alpha)
B = H * sin(1/2 * pi - alpha) / sin(alpha)
g = 9.8
dt = 0.000001   # 1 micro sec (time step size)

R = []
n = 0
while True:
    R.append(0.03*(n+1))
    if 0.03*(n+1) >= 0.03:
        break
    n += 1

L = []           # j, i
m = []           # j
M = []           # j, i
I_ball = []      # j
A_cyl = []       # j, i
A_ss = []        # j
I_cyl = []       # j, i
I_lss = []       # j, i
I_rss = []       # j
I_cap = []       # j, i
I_sys = []       # j, i
d = []           # j, i
time = []        # j, i
theta1 = []      # j, i
theta1_dot = []  # j, i
ref_pos = []     # j, i
ball_pos = []    # j, i
theta2_dot = []  # j, i
theta2 = []      # j, i
E1prime = []     # j, i
E2prime = []     # j, i

max_bin = []
min_bin = []

for j in range(len(R)):
    L_R = []
    k = 0
    while True:
        L_R.append(0.0005*(k+3*R[j]/0.0005))
        if 0.0005*(k+3*R[j]/0.0005) >= 0.12:
            break
        k += 1
    L.append(L_R)
    m.append(rho * (4/3*pi*R[j]**3))
    M_R = []
    for i in range(len(L[j])):
        M_R.append(sigma * (2*pi*R[j]*(L[j][i] - 2*R[j]) + 4*pi*R[j]**2))
    M.append(M_R)
    I_ball.append(2/5 * m[j] * R[j] ** 2)
    A_cyl_R = []
    for i in range(len(L[j])):
        A_cyl_R.append(2 * pi * R[j] * (L[j][i] - 2 * R[j]))
    A_cyl.append(A_cyl_R)
    A_ss.append(2 * pi * R[j] ** 2)
    I_cyl_R = []
    I_lss_R = []
    for i in range(len(L[j])):
        I_cyl_R.append(sigma * A_cyl[j][i] * (1 / 4 * R[j] ** 2 + 1 / 3 * (L[j][i] - 2 * R[j]) ** 2))
        I_lss_R.append(sigma * A_ss[j] * (1 / 12 * R[j] ** 2 + (1 / 2 * R[j] + (L[j][i] - 2 * R[j])) ** 2))
    I_cyl.append(I_cyl_R)
    I_lss.append(I_lss_R)
    I_rss.append(sigma * A_ss[j] * (1 / 3 * R[j] ** 2))
    I_cap_R = []
    I_sys_R = []
    for i in range(len(L[j])):
        I_cap_R.append(I_cyl[j][i] + I_lss[j][i] + I_rss[j])
    I_cap.append(I_cap_R)
    for i in range(len(L[j])):
        I_sys_R.append(I_ball[j] + I_cap[j][i])
    I_sys.append(I_sys_R)
    d_R = []
    for i in range(len(L[j])):
        d_R.append((M[j][i] * L[j][i] / 2 + m[j] * R[j]) / (M[j][i] + m[j]))
    d.append(d_R)

    time.append([0] * len(L[j]))

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

    theta1.append([0] * len(L[j]))
    theta1_dot.append([0] * len(L[j]))
    ref_pos.append([0] * len(L[j]))
    ball_pos.append([0] * len(L[j]))
    theta2_dot.append([0] * len(L[j]))
    theta2.append([0] * len(L[j]))
    E1prime.append([0] * len(L[j]))
    E2prime.append([0] * len(L[j]))

    for i in range(len(L[j])):
        # print("i")
        while True:
            # print("two")
            if (ball_pos[j][i] - ref_pos[j][i]) >= L[j][i] - 2 * R[j]:
                E1prime[j][i] = 1 / 2 * (I_ball[j] + m[j] * R[j] ** 2) * theta1_dot[j][i] ** 2 + m[j] * g * (
                            H - R[j] * sin(alpha) - R[j] * theta1[j][i] * sin(alpha)) + M[j][i] * g * (
                                            H - L[j][i] / 2 * sin(alpha)) - meu_1alpha * cos(alpha) * m[j] * g * R[j] * \
                                theta1[j][i]
                H = H - (L[j][i] - 2 * R[j]) * sin(alpha)
                theta2_dot_square = (E1prime[j][i] - (m[j] + M[j][i]) * g * (H - R[j] * sin(alpha)) - M[j][i] * g * (
                    1 / 2 * (L[j][i] - 2 * R[j]) * sin(alpha))) / (1 / 2 * (I_sys[j][i] + (m[j] + M[j][i]) * R[j] ** 2))
                theta2_dot[j][i] = theta2_dot_square ** 1 / 2
                while 0 <= theta2[j][i] < pi:
                    dtheta2 = theta2[j][i]
                    theta2[j][i] += theta2_dot[j][i] * dt
                    if theta2[j][i] > pi:
                        theta2[j][i] = pi
                    dtheta2 = theta2[j][i] - dtheta2
                    theta2_dot[j][i] += (((m[j] + M[j][i]) * g * R[j] * sin(alpha) - M[j][i] * g * (
                            1 / 2 * (L[j][i] - 2 * R[j]) * cos(alpha + theta2[j][i])) - meu_2alpha * cos(alpha) * (
                                    m[j] + M[j][i]) * g * R[j]) / (I_sys[j][i] + (m[j] + M[j][i]) * R[j] ** 2)) * dt
                    time[j][i] += dt
                    ball_pos[j][i] += R[j] * dtheta2
                    if ball_pos[j][i] >= Lb:
                        print(time[j][i])
                        print(L[j][i])
                        print(degrees)
                        print(meu_1)
                        print(meu_2)
                        break
                if ball_pos[j][i] >= Lb:
                    print(time[j][i])
                    print(L[j][i])
                    print(degrees)
                    print(meu_1)
                    print(meu_2)
                    break
                theta1[j][i] = 0
                H = H - R[j] * pi * sin(alpha)
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
                ref_pos[j][i] = ball_pos[j][i]

            if ball_pos[j][i] >= Lb:
                print(time[j][i])
                print(L[j][i])
                print(degrees)
                print(meu_1)
                print(meu_2)
                break
            # print("one")
            dtheta1 = theta1[j][i]
            theta1[j][i] += theta1_dot[j][i] * dt
            dtheta1 = theta1[j][i] - dtheta1
            # ref_theta1_dot = theta1_dot[i]
            theta1_dot[j][i] += (m[j] * g * R[j] * (sin(alpha) - meu_1alpha * cos(alpha))) / (
                        I_ball[j] + m[j] * R[j] ** 2) * dt
            ball_pos[j][i] += R[j] * dtheta1
            time[j][i] += dt
    print("degrees = ", degrees)
    print("L = ", L[j])
    # for i in range(len(L[j])):
    #     time[j][i] = round(time[j][i], 4)
    print("time = ", time[j])
    plot(L[j], time[j])
    show()

    mx, mn = getMinMax(time[j])
    mx[0] = L[j][int(mx[0])]
    mn[0] = L[j][int(mn[0])]
    print("Local Minimum [Length(m), Time(s)] for radius ", R[j], "m is ", mn)
    print("Local Maximum [Length(m), Time(s)] for radius ", R[j], "m is ", mx)
    max_bin.append(mx)
    min_bin.append(mn)

mx = getMax(max_bin)
mx[0] = max_bin[int(mx[0])][0]
mn = getMin(min_bin)
mn[0] = min_bin[int(mn[0])][0]
print("Minimum [Length(m), Time(s)] is ", mn)
print("Maximum [Length(m), Time(s)] is ", mx)

# df = pandas.DataFrame(time)
# df.to_excel('optimizedseries5.xlsx')
