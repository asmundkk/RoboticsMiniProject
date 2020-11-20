from numpy import pi, array, arctan2, sqrt, sin, cos, cross, dot
from numpy.linalg import norm

class link():

    def __init__(self, a, alpha, d):
        self.a = a                              # smallest distance z axis, also called r
        self.alpha = alpha                      # angle between z axis
        self.d = d                              # "vertical" distance between rotational axis
        self.th = float                         # joint angle

class vec():

    def __init__(self, v):
        self.x = v[0]
        self.y = v[1]
        self.z = v[2]
        self.v = array([self.x, self.y, self.z]).T            # column vec by default

        self.print()

    def print(self):
        print("vec: \n", self.v)

def normalize_theta(theta_list):
    """
    sets the joint values fo the robot in the range [-pi, pi]

    theta_list: list of joint angles
    """

    for i in range(len(theta_list)):
        theta_list[i] = 2 * pi - (abs(theta_list[i]) % 2 * pi)

        if theta_list[i] > pi and theta_list[i] < 2 * pi:
            theta_list[i] = 2 * pi - theta_list[i]
    return theta_list

def analytical_inverse_kinematics(T, S = 1, E = 1, W = 1):
    """
    calculates the inverse kinematics for UR 5 robot, given the T0e transformation

    T: 4x4 transformation matrix from home to end effector configuration
    S: shoulder up or down  1 | -1
    E: elbow up or down     1 | -1
    W: wrist up or down     1 | -1

    theta_list: list of all the joint angles
    """

    ne = T[0:3, 0]
    se = T[0:3, 1]
    ae = T[0:3, 2]
    pe = T[0:3, 3]

    print(ne)
    print(se)
    print(ae)
    print(pe)

    ne = vec(ne)
    se = vec(se)
    ae = vec(ae)
    pe = vec(pe)

    l1 = link( 0,       pi/2,   0.089)
    l2 = link(-0.425,   0,      0)
    l3 = link(-0.392,   0,      0)
    l4 = link( 0,       pi/2,   0.109)
    l5 = link( 0,      -pi/2,   0.095)
    l6 = link( 0,       0,      0.082)

    print("p5")
    p5 = vec(pe.v - l6.d * ae.v)
    ph1 = arctan2(S * p5.y, S * p5.x)
    de1 = arctan2(l4.d, sqrt(p5.x**2 + p5.y**2 - l4.d**2))
    l1.th = ph1 + S * de1

    print("ph1", ph1)
    print("de1", de1)
    print("th1", l1.th)

    x1 = vec([cos(l1.th), sin(l1.th), 0])
    y1 = vec([0, 0, 1])
    z1 = vec([sin(l1.th), -cos(l1.th), 0])

    cross4 = cross(z1.v, ae.v)
    z4 = vec(W * S *(cross4 / (norm(cross4))))
    x4 = vec(cross(z1.v, z4.v))

    y4 = z1
    l5.th = arctan2(dot(-ae.v, x4.v), dot(ae.v, y4.v))
    l6.th = arctan2(dot(-z4.v, ne.v), dot(-z4.v, se.v))

    p3 = vec(p5.v - l5.d * z4.v - l4.d * z1.v)
    ph = S * sqrt(p3.x**2 + p3.y**2)
    pv = p3.z -l1.d

    #  ci and si means cos(thi) and sin(thi)
    c3 = (ph**2 + pv**2 - l2.a**2 - l3.a**2) / (2 * l2.a * l3.a)
    s3 = E * sqrt(1 - c3**2)
    l3.th = arctan2(s3, c3)

    c2 = (ph * (l2.a + l3.a * c3) + pv * l3.a * s3) / (ph**2 + pv**2)
    s2 = (pv * (l2.a + l3.a * c3) - ph * l3.a * s3) / (ph**2 + pv**2)
    l2.th = arctan2(s2, c2)

    th234 = arctan2(dot(z4.v, x1.v), dot(z4.v, y1.v))
    l4.th = th234 - (l2.th + l3.th)

    theta_list = normalize_theta(array([l1.th, l2.th, l3.th, l4.th, l5.th, l6.th]))
    return theta_list

if __name__ == '__main__':

    T = array([[0, 0, 1, 0.4671],
               [1, 0, 0, 0.1099],
               [0, 1, 0, 0.4195],
               [0, 0, 0, 1]])

    theta_list = analytical_inverse_kinematics(T)

    print("\n", theta_list)
