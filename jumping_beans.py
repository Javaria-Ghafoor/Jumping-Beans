"""
This is a physical simulation. To run the program, ensure that you have VPython 7 installed and configured on your IDE
"""

from math import sin, pi, cos
from vpython import triangle, vertex, vec, quad, compound, sphere, cylinder, color, canvas, mag, diff_angle, rate

scene = canvas(width=1200, height=600, center=vec(0, 0, 0), background=color.black)


def make_ramp(alpha, H):
    parts = []
    B = H * sin(0.5 * pi - alpha) / sin(alpha)
    front = triangle(
        v0=vertex(pos=vec(-0.5 * B, 0, 0.25 * H)),
        v1=vertex(pos=vec(-0.5 * B, H, 0.25 * H)),
        v2=vertex(pos=vec(0.5 * B, 0, 0.25 * H))
    )
    back = triangle(
        v0=vertex(pos=vec(-0.5 * B, 0, - 0.25 * H)),
        v1=vertex(pos=vec(-0.5 * B, H, - 0.25 * H)),
        v2=vertex(pos=vec(0.5 * B, 0, - 0.25 * H))
    )
    top = quad(
        v0=front.v1,
        v1=back.v1,
        v2=back.v2,
        v3=front.v2,
    )
    bottom = quad(
        v0=front.v0,
        v1=back.v0,
        v2=back.v2,
        v3=front.v2
    )
    side = quad(
        v0=front.v0,
        v1=back.v0,
        v2=back.v1,
        v3=front.v1
    )
    parts.append(front)
    parts.append(back)
    parts.append(top)
    parts.append(bottom)
    parts.append(side)
    obj = compound(parts)
    obj.height = H
    obj.base = B
    obj.angle = alpha
    # obj.color = color.cyan
    # obj.opacity = 0.1

    return obj


def make_ball(R):
    obj = sphere(radius=R, color=color.red)
    return obj


def make_capsule(R, L):
    parts = []
    sph1 = sphere(pos=vec(R, 0, 0), radius=R, color=color.cyan, opacity=0.3)
    cyl = cylinder(pos=vec(R, 0, 0), axis=vec(L - 2 * R, 0, 0), radius=R, color=color.cyan, opacity=0.3)
    sph2 = sphere(pos=sph1.pos+cyl.axis, radius=R, color=color.cyan, opacity=0.3)
    parts.append(sph1)
    parts.append(cyl)
    parts.append(sph2)
    obj = compound(parts)
    return obj

"""
variations of alpha in experiment: (in radians)
0.1960
0.2218
0.2474
0.2733
0.2994
0.3256
0.3522
0.3791
0.4061
0.4334
0.4611
0.4892
0.5178
0.5468
"""

alpha = 0.1960
L = 0.0235
R = 0.003175   # radius of ball and capsule assumed the same throughout
Lb = 0.9
m = 0.0010355
M = 0.0001266
meu_1 = 0  # simulation constructed for meu = constant case
meu_2 = 0

# make system:

H = Lb * sin(alpha)
B = H * sin(1/2 * pi - alpha) / sin(alpha)
ramp = make_ramp(alpha, H)
capsule = make_capsule(R, L)
# print(capsule.pos)
capsule.pos = vec(-1/2 * B + 1/2 * (L - 2 * R) + R, H + R, 0)
P = vec(capsule.pos.x, capsule.pos.y, capsule.pos.z)
origin = vec(capsule.pos.x - (1/2 * (L - 2 * R) + R), capsule.pos.y, capsule.pos.z)
capsule.rotate(angle=alpha, axis=vec(0, 0, -1), origin=origin)  # (angle, position, centre)
ball = make_ball(R)
ball.pos = vec(capsule.pos.x - 1/2 * (L - 2 * R), capsule.pos.y, capsule.pos.z)
ball.rotate(angle=alpha, axis=vec(0, 0, -1), origin=capsule.pos)

I_ball = 2/5 * m * R ** 2
g = 9.8

# A = ball.pos-capsule.pos
# B = capsule.pos-ball.pos
# print(mag(A))
# print(mag(B))
# print(0.5 * (L-2*R))
# print(norm(A))
# print(norm(-A))
# print(norm(B))

A = ball.pos-capsule.pos
dt = 0.000001   # 1 micro sec (time step size)
theta1 = 0
theta1_dot = 0
ref_pos = vec(ball.pos.x, ball.pos.y, ball.pos.z)

A_cyl = 2 * pi * R * (L - 2 * R)
A_ss = 2 * pi * R ** 2
sigma = M / (A_cyl + 2 * A_ss)
I_cyl = sigma * A_cyl * (1/4 * R ** 2 + 1/3 * (L - 2 * R)**2)
I_lss = sigma * A_ss * (1/12 * R ** 2 + (1/2 * R + (L - 2 * R))**2)
I_rss = sigma * A_ss * (1/3 * R ** 2)
I_cap = I_cyl + I_lss + I_rss
I_sys = I_ball + I_cap
d = (M * L/2 + m * R)/(M + m)

theta2 = 0
theta2_dot = 0

time = 0
# print(1//dt)

while True:
    # rate(50)
    if mag(ball.pos - ref_pos) >= L - 2 * R:
        # the following line added to keep simulation outlook smooth
        ball.pos = vec(ref_pos.x+(L-2*R)*cos(alpha), ref_pos.y-(L-2*R)*sin(alpha), ref_pos.z)
        # print('W')
        # updating via energy transfer on impact:
        E1prime = 1 / 2 * (I_ball + m * R ** 2) * theta1_dot ** 2 + m * g * (
                H - R * sin(alpha) - R * theta1 * sin(alpha)) + M * g * (
                             H - L / 2 * sin(alpha)) - meu_1 * cos(alpha) * m * g * R * theta1
        H = H - (L - 2 * R) * sin(alpha)
        theta2_dot_square = (E1prime - (m + M) * g * (H - R * sin(alpha)) - M * g * (
                1 / 2 * (L - 2 * R) * sin(alpha))) / (1 / 2 * (I_sys + (m + M) * R ** 2))
        theta2_dot = theta2_dot_square ** 1 / 2
        cap_ref = vec(capsule.pos.x, capsule.pos.y, capsule.pos.z)
        while 0 <= theta2 < pi:
            # rate(50)
            dtheta2 = theta2
            theta2 += theta2_dot * dt
            if theta2 > pi:
                theta2 = pi
            dtheta2 = theta2 - dtheta2
            theta2_dot += (((m + M) * g * R * sin(alpha) - M * g * (
                    1 / 2 * (L - 2 * R) * cos(alpha + theta2)) - meu_2 * cos(alpha) * (
                                       m + M) * g * R) / (I_sys + (m + M) * R ** 2)) * dt
            time += dt
            capsule.pos += vec(R * dtheta2 * cos(alpha), - R * dtheta2 * sin(alpha), 0)
            origin = vec(capsule.pos.x + (0.5 * L - R) * cos(alpha + theta2),
                         capsule.pos.y - (0.5 * L - R) * sin(alpha + theta2),
                         capsule.pos.z)
            # white_mark = sphere(pos=origin, radius=0.5, color=color.white)
            capsule.rotate(angle=dtheta2, axis=vec(0, 0, -1), origin=origin)
            ball.pos += vec(R * dtheta2 * cos(alpha), - R * dtheta2 * sin(alpha), 0)
            if ball.pos.y < 0:
                break
        if ball.pos.y < 0:
            break
        check_vec = capsule.pos - vec(0.5 * B, 0, 0)
        if diff_angle(check_vec, - vec(0.5 * B, 0, 0)) != alpha:
            print("Simulation and upgradation of positions and numerical values with dt = 1 micro sec ",
                  "resulted in an angle of inclination of the simulated capsule to be: ", diff_angle(check_vec, - vec(0.5 * B, 0, 0)))
            print("The true angle of inclination is: ", alpha)
            print("Adjusting the position of capsule for the simulation to visually seem smooth...")
            """
            This is python's numerical computation fault :(. Further refer to:
            https://www-uxsup.csx.cam.ac.uk/courses/moved.NumericalPython/paper_1.pdf
            to see how python creates this gap for numerical calculation errors.
            """
            # the following lines added to keep simulation outlook smooth
            capsule.pos = vec(cap_ref.x + R * pi * cos(alpha), cap_ref.y - R * pi * sin(alpha), cap_ref.z)
            origin = vec(capsule.pos.x + (0.5 * L - R) * cos(alpha),
                         capsule.pos.y - (0.5 * L - R) * sin(alpha),
                         capsule.pos.z)
            capsule.rotate(angle=pi, axis=vec(0, 0, -1), origin=origin)
        # print('J')

        ref_pos = vec(ball.pos.x, ball.pos.y, ball.pos.z)
        H = H - R * pi * sin(alpha)
        theta1 = 0
        theta2 = 0
    if ball.pos.y < 0:
        break
    dtheta1 = theta1
    theta1 += theta1_dot * dt
    dtheta1 = theta1 - dtheta1
    theta1_dot += (m * g * R * (sin(alpha) - meu_1 * cos(alpha))) / (I_ball + m * R ** 2) * dt
    ball.pos += vec(R * dtheta1 * cos(alpha), - R * dtheta1 * sin(alpha), 0)
    time += dt
    # print('O')

print("time taken at angle of inclination ", alpha, " is ", time, "s")
