import numpy as np
import matplotlib.pyplot as plt


class FourBarLinkage:
    '''
    Four bar linkage simulation. 1st joint (p1) located at (0,0). 4th joint at distance d along x-axis
    '''
    def __init__(self, d12=1, d23=1, d34=1, d14=1):
        self.d12 = d12  # distance between p1 and p2
        self.d23 = d23
        self.d34 = d34
        self.d14 = d14

        self.p1 = np.array([[0], [0]])
        self.p4 = np.array([[self.d14], [0]])
        self.p2, self.p3 = None, None
        self.set_position(np.pi/2)  # set p2 and p3


    def set_position(self, theta):
        # From theta (orientation of p1->p2) compute p2 and p3
        self.p2 = self.p1 + np.array(self.d12*[[np.cos(theta)], [np.sin(theta)]])
        self.p3 = self.bilateration_method(self.p2, self.p4, self.d23, self.d34)

    def bilateration_method(self, p1, p2, d2, d3):
        #  Compute p3 from known p1 and p2 and distances between all 3 points
        d1 = np.linalg.norm(p2-p1)
        s1, s2, s3 = d1 ** 2, d2 ** 2, d3 ** 2
        sign = 1  # 1 = p3 to the left, -1 = to the right
        area = (sign / 4) * np.sqrt((s1 + s2 + s3) ** 2 - 2 * (s1**2+s2**2+s3**2))
        z1 = s1 + s2 - s3
        z2 = 4 * area
        Z = (1/(2*s1))*np.array([[z1, -z2], [z2, z1]])  # bilateration matrix
        v12 = p2 - p1  # vector from p1 to p2
        p3 = p1 + np.matmul(Z, v12)

        return p3

    def plot(self):
        plt.cla()
        plt.axis([-1, 2, -1.5, 1.5])
        plt.plot([self.p1[0], self.p2[0]], [self.p1[1], self.p2[1]])
        plt.plot([self.p2[0], self.p3[0]], [self.p2[1], self.p3[1]])
        plt.plot([self.p3[0], self.p4[0]], [self.p3[1], self.p4[1]])
        plt.draw()

    def run_trajectory(self):
        plt.axis([-1, 2, -1.5, 1.5])
        plt.ion()
        plt.show()

        for theta in np.arange(0, np.pi, np.pi/100):
            self.set_position(theta)
            self.plot()
            plt.pause(0.001)


def main():
    print('Simulate four bar linkage')
    fbl = FourBarLinkage()
    fbl.set_position(np.pi/4)
    print('p1 at: ', fbl.p1[0], fbl.p1[1])
    print('p2 at: ', fbl.p2[0], fbl.p2[1])
    print('p3 at: ', fbl.p3[0], fbl.p3[1])
    print('p4 at: ', fbl.p4[0], fbl.p4[1])
    fbl.run_trajectory()



if __name__ == '__main__':
    main()
