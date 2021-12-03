import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# две окружности, которые не пересекаются
def firstTest(forLevels):
    angle = np.linspace(0, 2 * np.pi, 100)

    radius = 2

    x1 = radius * np.cos(angle) + 2
    y1 = radius * np.sin(angle)

    x2 = radius * np.cos(angle) - 3
    y2 = radius * np.sin(angle)

    figure, axes = plt.subplots(1)

    mesh1 = np.arange(-5, 4.1, 0.1)
    mesh2 = np.arange(-3, 3.1, 0.1)
    meshX, meshY = np.meshgrid(mesh1, mesh2)

    z = (meshX-2)**2 + (meshY**2) - 4 + (meshX+3)**2 + meshY**2 - 4

    contourf_ = axes.contourf(meshX, meshY, z, levels=forLevels)

    return angle, x1, y1, x2, y2, axes, figure, mesh1, mesh2, meshX, meshY, z, contourf_

# две окружности, которые пересекаются в одной точке
def secondTest(forLevels):
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

    contourf_ = axes.contourf(meshX, meshY, z, levels=forLevels)

    return angle, x1, y1, x2, y2, axes, figure, mesh1, mesh2, meshX, meshY, z, contourf_

# две окружности, которые пересекаются в двух точках
def thirdTest(forLevels):
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

    z = ((meshX)**2 + (meshY + 2)**2 - 4)**2 + ((meshX)**2 + (meshY - 1)**2 - 4)**2

    contourf_ = axes.contourf(meshX, meshY, z, levels=forLevels)

    return angle, x1, y1, x2, y2, axes, figure, meshX, meshY, z, contourf_

# две окружности, которые пересекаются в двух точках + прямая
def fourthTest(forLevels):
    angle = np.linspace(0, 2 * np.pi, 100)

    radius = 2

    x1 = radius * np.cos(angle)
    y1 = radius * np.sin(angle) - 2

    x2 = radius * np.cos(angle)
    y2 = radius * np.sin(angle) + 1

    lineX = np.linspace(-2.9, 1.76, 10)
    lineY = 500.0 / 441 * lineX + 1

    figure, axes = plt.subplots(1)

    mesh1 = np.arange(-3, 3.1, 0.1)
    mesh2 = np.arange(-4, 3.1, 0.1)
    meshX, meshY = np.meshgrid(mesh1, mesh2)

    z = ((meshX)**2 + (meshY + 2)**2 - 4)**2 + ((meshX)**2 + (meshY - 1)**2 - 4)**2 + (500.0/441 * meshX + 1 - meshY)**2

    contourf_ = axes.contourf(meshX, meshY, z, levels=forLevels)

    return angle, x1, y1, x2, y2, lineX, lineY, axes, figure, meshX, meshY, z, contourf_

# три попарно пересекающиеся
def fifthTest(forLevels):
    lineX1 = np.linspace(-2.9, 1.9, 10)
    lineY1 = lineX1 + 1

    lineX2 = np.linspace(-2.9, 2.9, 10)
    lineY2 = 1.0 / 10 * lineX2

    lineX3 = np.linspace(-0.9, 2.9, 10)
    lineY3 = -lineX3 + 2

    figure, axes = plt.subplots(1)

    mesh1 = np.arange(-3, 3.1, 0.1)
    mesh2 = np.arange(-4, 3.1, 0.1)
    meshX, meshY = np.meshgrid(mesh1, mesh2)

    z = (meshX + 1 - meshY)**2 + (1.0/10 * meshX - meshY)**2 + (-meshX + 2 - meshY)**2

    contourf_ = axes.contourf(meshX, meshY, z, levels=forLevels)

    return lineX1, lineY1, lineX2, lineY2, lineX3, lineY3, axes, figure, meshX, meshY, z, contourf_

# прямая, которая пересекает синусоиду
def sixthTest(forLevels):
    argX = np.linspace(-10, 10, 1000)
    func = 2 + 4 * np.sin(2 * argX + 1)

    lineX = np.linspace(-2, 4, 10)
    lineY = lineX

    figure, axes = plt.subplots(1)

    mesh1 = np.arange(-7, 7, 0.1)
    mesh2 = np.arange(-2.1, 5.3, 0.1)
    meshX, meshY = np.meshgrid(mesh1, mesh2)

    z = (2 + 4 * np.sin(2 * meshX + 1) - meshY)**2 + (meshX - meshY)**2

    contourf_ = axes.contourf(meshX, meshY, z, levels=forLevels)

    return argX, func, lineX, lineY, figure, axes, meshX, meshY, z, contourf_

def main():
    x = []
    y = []

    with open("coords.txt") as file:
        for line in file:
            xC, yC = line.split()
            x.append(float(xC))
            y.append(float(yC))

    forLevels = np.linspace(0, 10, 20)

    angle, x1, y1, x2, y2, axes, figure, mesh1, mesh2, meshX, meshY, z, contourf_ = firstTest(forLevels)
    #angle, x1, y1, x2, y2, axes, figure, mesh1, mesh2, meshX, meshY, z, contourf_ = secondTest(forLevels)
    #angle, x1, y1, x2, y2, axes, figure, meshX, meshY, z, contourf_ = thirdTest(forLevels)
    #angle, x1, y1, x2, y2, lineX, lineY, axes, figure, meshX, meshY, z, contourf_ = fourthTest(forLevels)
    #lineX1, lineY1, lineX2, lineY2, lineX3, lineY3, axes, figure, meshX, meshY, z, contourf_ = fifthTest(forLevels)
    #argX, func, lineX, lineY, figure, axes, meshX, meshY, z, contourf_ = sixthTest(forLevels)

    # для 1-3 теста
    axes.plot(x1, y1, label='F1', linewidth=2, color='g')
    axes.plot(x2, y2, label='F2', linewidth=2, color='y')
    # axes.plot(lineX, lineY, linewidth=2, color='m')  # для 4-го теста

    # для 5-го теста
    # axes.plot(lineX1, lineY1, label='F1', linewidth=2, color='g')
    # axes.plot(lineX2, lineY2, label='F2', linewidth=2, color='y')
    # axes.plot(lineX3, lineY3, label='F2', linewidth=2, color='m')

    # для 6-го теста
    #axes.plot(argX, func, label='F1', linewidth=2, color='g')
    #axes.plot(lineX, lineY, label='F2', linewidth=4, color='m')

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

    axes.contourf(meshX, meshY, z, levels=forLevels)

    figure.colorbar(contourf_, shrink=1)

    plt.grid()
    axes.set_aspect(1)
    plt.show()

if __name__ == "__main__":
    main()