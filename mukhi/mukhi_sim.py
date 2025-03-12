import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.path import Path

# Параметры симуляции
n1 = 5            # число частиц первой группы
n2 = 10            # число частиц второй группы
n = n1 + n2        # общее число частиц
R = 1              # "радиус" банки (используется для генерации многоугольника)
r = 0.05           # радиус частицы (для расчёта столкновений)
V = 8              # масштаб начальной скорости
dt = 0.005          # временной шаг
k = 1              # коэффициент отражения перпендикулярной составляющей (стенка)
k_pr = 0.98        # коэффициент для параллельной составляющей (стенка)
k_m = 1            # коэффициент обмена для перпендикулярных компонент (частицы)
k_m_pr = 1         # коэффициент обмена для параллельных компонент (частицы)

# Количество сторон многоугольника (n-угольник)
poly_n = 4  # можно менять для получения разного числа сторон

# Генерация вершин многоугольника
poly_angles = np.linspace(0, 2*np.pi, poly_n, endpoint=False)
poly_x = R * np.cos(poly_angles)
poly_y = R * np.sin(poly_angles)
poly_verts = np.column_stack((poly_x, poly_y))
# Создаём объект для проверки принадлежности точек многоугольнику
polygon_path = Path(poly_verts)

# Инициализация начальных положений частиц (используем полярные координаты)
Rs = R * np.random.rand(n)
Phis1 = np.pi * np.random.rand(n1)
Phis2 = np.pi * (1 + np.random.rand(n2))
Phis = np.concatenate([Phis1, Phis2])
X = Rs * np.cos(Phis)
Y = Rs * np.sin(Phis)

# Инициализация скоростей
VX = V * (2 * np.random.rand(n) - 1)
VY = V * (2 * np.random.rand(n) - 1)

# Настройка визуализации
fig = plt.figure(figsize=[15, 7])
ax = fig.add_subplot(1, 1, 1)
ax.axis('equal')

# Отрисовка частиц двух групп
Mukhi1 = ax.plot(X[0:n1], Y[0:n1], 'bo')[0]
Mukhi2 = ax.plot(X[n1:n], Y[n1:n], 'ro')[0]

# Отрисовка траекторий для первых трёх частиц
Trace1X = np.array([X[0]])
Trace1Y = np.array([Y[0]])
Trace1 = ax.plot(Trace1X, Trace1Y, '--')[0]
Trace2X = np.array([X[1]])
Trace2Y = np.array([Y[1]])
Trace2 = ax.plot(Trace2X, Trace2Y, '--')[0]
Trace3X = np.array([X[2]])
Trace3Y = np.array([Y[2]])
Trace3 = ax.plot(Trace3X, Trace3Y, '--')[0]

# Отрисовка полигона (банки) с помощью его вершин
polygon_patch = plt.Polygon(poly_verts, closed=True, fill=False, edgecolor='k')
ax.add_patch(polygon_patch)
plt.draw()

def kadr(frame):
    global X, Y, VX, VY, dt, Trace1X, Trace1Y, Trace2X, Trace2Y, Trace3X, Trace3Y

    # Вычисляем предварительные новые координаты
    Xn = X + VX * dt
    Yn = Y + VY * dt

    # Проверяем столкновения с границей многоугольника для каждой частицы
    for i in range(n):
        if not polygon_path.contains_point([Xn[i], Yn[i]]):
            # Частица оказалась вне полигона, ищем ближайший отрезок границы
            closest_dist = np.inf
            closest_edge = None
            closest_point = None
            for j in range(poly_n):
                x1, y1 = poly_x[j], poly_y[j]
                x2, y2 = poly_x[(j + 1) % poly_n], poly_y[(j + 1) % poly_n]
                edge_vec = np.array([x2 - x1, y2 - y1])
                p_vec = np.array([Xn[i] - x1, Yn[i] - y1])
                edge_len_sq = np.dot(edge_vec, edge_vec)
                t_val = np.clip(np.dot(p_vec, edge_vec) / edge_len_sq, 0, 1)
                projection = np.array([x1, y1]) + t_val * edge_vec
                dist = np.linalg.norm(np.array([Xn[i], Yn[i]]) - projection)
                if dist < closest_dist:
                    closest_dist = dist
                    closest_edge = (np.array([x1, y1]), np.array([x2, y2]))
                    closest_point = projection
            # Вычисляем внешнюю нормаль для данного ребра.
            # Для выпуклого многоугольника нормаль можно получить через точку-середину ребра:
            midpoint = (closest_edge[0] + closest_edge[1]) / 2
            norm_mid = np.linalg.norm(midpoint)
            if norm_mid != 0:
                normal = midpoint / norm_mid
            else:
                normal = np.array([0, 0])
            # Разлагаем скорость на перпендикулярную и тангенциальную составляющие
            v = np.array([VX[i], VY[i]])
            v_perp = np.dot(v, normal) * normal
            v_tan = v - v_perp
            # Отражаем скорость: перпендикулярная инвертируется с коэффициентом k, тангенциальная умножается на k_pr
            v_new = -k * v_perp + k_pr * v_tan
            VX[i] = v_new[0]
            VY[i] = v_new[1]
            # Корректируем позицию, проецируя её на ближайшую точку на ребре
            Xn[i] = closest_point[0]
            Yn[i] = closest_point[1]

    # Столкновения между частицами (аналогично исходной реализации)
    for i in range(n - 1):
        for j in range(i + 1, n):
            if (Xn[i] - Xn[j])**2 + (Yn[i] - Yn[j])**2 < 4 * r**2:
                dx = Xn[i] - Xn[j]
                dy = Yn[i] - Yn[j]
                distance = np.sqrt(dx**2 + dy**2)
                if distance == 0:
                    continue
                # Нормаль по линии, соединяющей центры частиц
                n_coll = np.array([dx, dy]) / distance

                # Разложение скоростей для частицы i
                v_i = np.array([VX[i], VY[i]])
                v_perp_i = np.dot(v_i, n_coll) * n_coll
                v_tan_i = v_i - v_perp_i
                # Разложение скоростей для частицы j
                v_j = np.array([VX[j], VY[j]])
                v_perp_j = np.dot(v_j, n_coll) * n_coll
                v_tan_j = v_j - v_perp_j

                # Обмен перпендикулярными компонентами с коэффициентом k_m
                v_perp_i_new = v_perp_j * k_m
                v_perp_j_new = v_perp_i * k_m
                # Тангенциальные компоненты умножаются на k_m_pr
                v_tan_i_new = v_tan_i * k_m_pr
                v_tan_j_new = v_tan_j * k_m_pr

                v_i_new = v_perp_i_new + v_tan_i_new
                v_j_new = v_perp_j_new + v_tan_j_new

                VX[i], VY[i] = v_i_new[0], v_i_new[1]
                VX[j], VY[j] = v_j_new[0], v_j_new[1]
                print('Столкновение частиц', i, j)

    # Обновляем положения частиц
    X = Xn
    Y = Yn

    # Обновляем траектории для первых трёх частиц
    Trace1X = np.append(Trace1X, X[0])
    Trace1Y = np.append(Trace1Y, Y[0])
    Trace2X = np.append(Trace2X, X[1])
    Trace2Y = np.append(Trace2Y, Y[1])
    Trace3X = np.append(Trace3X, X[2])
    Trace3Y = np.append(Trace3Y, Y[2])

    # Обновляем графические данные
    Mukhi1.set_data(X[0:n1], Y[0:n1])
    Mukhi2.set_data(X[n1:n], Y[n1:n])
    Trace1.set_data(Trace1X, Trace1Y)
    Trace2.set_data(Trace2X, Trace2Y)
    Trace3.set_data(Trace3X, Trace3Y)

    return [Mukhi1, Mukhi2, Trace1, Trace2, Trace3]

kino = FuncAnimation(fig, kadr, interval=dt*1000, blit=True)
plt.show()
