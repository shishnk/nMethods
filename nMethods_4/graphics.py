import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def firstTest():  # две окружности, которые не пересекаются
    angle = np.linspace(0, 2 * np.pi, 100)

    radius = 2

    x1 = radius * np.cos(angle) + 2
    y1 = radius * np.sin(angle)

    x2 = radius * np.cos(angle) - 3
    y2 = radius * np.sin(angle)

    figure, axes = plt.subplots(1)

    axes.plot(x1, y1, label='F1', linewidth=2, color='g')
    axes.plot(x2, y2, label='F2', linewidth=2, color='y')

    axes.legend()

    mesh1 = np.arange(-5, 4.1, 0.1)
    mesh2 = np.arange(-3, 3.1, 0.1)
    meshX, meshY = np.meshgrid(mesh1, mesh2)

    z = (meshX-2)**2 + (meshY**2) - 4 + (meshX+3)**2 + meshY**2 - 4

    contourf_ = axes.contourf(meshX, meshY, z, levels=12)

    return angle, x1, y1, x2, y2, axes, figure, mesh1, mesh2, meshX, meshY, z, contourf_


def secondTest():  # две окружности, которые пересекаются в одной точке
    angle = np.linspace(0, 2 * np.pi, 100)

    radius = 2

    x1 = radius * np.cos(angle) - 3
    y1 = radius * np.sin(angle)

    x2 = radius * np.cos(angle) + 1
    y2 = radius * np.sin(angle)

    figure, axes = plt.subplots(1)

    mesh1 = np.arange(-5, 3.1, 0.1)
    mesh2 = np.arange(-3, 3.1, 0.1)
    meshX, meshY = np.meshgrid(mesh1, mesh2)

    z = (meshX+3)**2 + (meshY**2) - 4 + (meshX-1)**2 + meshY**2 - 4

    contourf_ = axes.contourf(meshX, meshY, z, levels=12)

    return angle, x1, y1, x2, y2, axes, figure, mesh1, mesh2, meshX, meshY, z, contourf_


def thirdTest():  # две окружности, которые пересекаются в двух точках
    angle = np.linspace(0, 2 * np.pi, 100)

    radius = 2

    x1 = radius * np.cos(angle)
    y1 = radius * np.sin(angle) - 2

    x2 = radius * np.cos(angle)
    y2 = radius * np.sin(angle) + 1

    figure, axes = plt.subplots(1)

    mesh1 = np.arange(-3, 3.1, 0.1)
    mesh2 = np.arange(-4, 3.1, 0.1)
    meshX, meshY = np.meshgrid(mesh1, mesh2)

    z = (meshX)**2 + (meshY + 2)**2 - 4 + (meshX)**2 + (meshY - 1)**2 - 4

    contourf_ = axes.contourf(meshX, meshY, z, levels=12)

    return angle, x1, y1, x2, y2, axes, figure, meshX, meshY, z, contourf_


def main():
    x = []
    y = []

    with open("coords.txt") as file:
        for line in file:
            xC, yC = line.split()
            x.append((float)(xC))
            y.append((float)(yC))

    # angle, x1, y1, x2, y2, axes, figure, mesh1, mesh2, meshX, meshY, z, contourf_ = firstTest()

    # angle, x1, y1, x2, y2, axes, figure, mesh1, mesh2, meshX, meshY, z, contourf_ = secondTest()

    angle, x1, y1, x2, y2, axes, figure, meshX, meshY, z, contourf_ = thirdTest()

    # fourthTest()

    # fifthTest()

    # sixthTest()

    axes.plot(x1, y1, label='F1', linewidth=2, color='g')
    axes.plot(x2, y2, label='F2', linewidth=2, color='y')

    axes.legend()

    plt.xlim(-7, 6)
    plt.ylim(-5, 4)

    axes.xaxis.set_major_locator(ticker.MultipleLocator(1))
    axes.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    axes.yaxis.set_major_locator(ticker.MultipleLocator(1))
    axes.yaxis.set_minor_locator(ticker.MultipleLocator(1))

    plt.xlabel("x")
    plt.ylabel("y")

    plt.plot(x, y, '-o', markersize=6, color='c')
    plt.plot(x[-1], y[-1], 'o', markersize=9, color='r')

    axes.contourf(meshX, meshY, z, levels=12)

    figure.colorbar(contourf_, shrink=0.93)

    plt.grid()
    axes.set_aspect(1)
    plt.show()

if __name__ == "__main__":
    main()
