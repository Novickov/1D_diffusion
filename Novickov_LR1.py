# Для рассчетов и построения графика необходимы дополнительные библиотеки, импортируем их:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Исходные данные (вариант 17):

''' Моделируемый процесс - диффузия
    Подложка - Si <100> p-типа
    Концентрация P: C_P = 6*10E+16 см^-3
    Режим загонки: T_zg = 1150 C; t_zg = 15 мин
    Режим разгонки: T_rzg = 950 C; t_rzg = 15 мин
    Глубина p-n: x = 0,45 мкм
'''

# Передаем начельные условия и физические константы переменным:

C_P = 6*10E+16
T_zg = 1150+273 # K
t_zg = 15*60 # c
T_rzg = 950+273 # K
t_rzg = 15*60 # c
x = float(input('Введите значение глубины расчета в мкм: '))
C0 = 1E+21 # Предел растворимости P в Si, см^-3
k = 8.617E-5 # Константа Больцмана, эВ
x_pn = 0.45

# Точность расчетов:

n = 100 # Число рассчитываемых точек
dx = x * 1E-4 / n # Шаг по глубине, см
dt = 1 # Шаг по времени, с 

# Приступаем к рассчетам:

def Eg(T):
    return 1.17 - 4.73E-4 * (T ** 2) / (T + 636)
def Nc(T):
    return 6.2E+15 * (T ** (3/2))
def Nv(T):
    return 3.5E+15 * (T ** (3/2))
def ni(T):
    return ((Nc(T) * Nv(T)) ** 0.5) * np.exp(-Eg(T) / (2 * k * T))
def D(C, T):
    return 3.85*np.exp(-3.66/(k*T))+4.44*(C/ni(T))*np.exp(-4/(k*T))+44.2*(C/ni(T))**2*np.exp(-4.37/(k*T))

# Создаем нулевые массивы длиной n

C_list = [0] * n
delta_list = [0] * n
lyambda_list = [0] * n
a_list = [0] * n
b_list = [0] * n
r_list = [0] * n
d_list = [0] * n
x_list = [0] * n
pn_list = [0] * n

# Процесс загонки:

# Задаем граничные условия

C_list[0] = C0
b_list[0] = 0
d_list[0] = 0
a_list[0] = 1
r_list[0] = C0
delta_list[0] = - d_list[0] / a_list[0]
lyambda_list[0] = r_list[0] / a_list[0]

C_list[n-1] = 0
b_list[n-1] = 0
d_list[n-1] = 0
a_list[n-1] = 1
r_list[n-1] = 0

# Заполняем массив значений x:

for i in range(n):
    x_list[i] = dx * i

# Запускаем основной алгоритм для загонки:

for j in range(t_zg):
    for i in range(1, n - 1): 
        b_list[i] = 1
        d_list[i] = 1
        a_list[i] = - (2 + (dx ** 2) / (D(C_list[i], T_zg) * dt))
        r_list[i] = - (dx ** 2) / (D(C_list[i], T_zg) * dt) * C_list[i]
    for i in range(1, n):
        delta_list[i] = - d_list[i] / (a_list[i] + b_list[i] * delta_list[i-1])
        lyambda_list[i] = (r_list[i] - b_list[i] * lyambda_list[i-1]) / (a_list[i] + b_list[i] * delta_list[i-1])
    C_list[n-1] = lyambda_list[n-1]
    for i in range(n-2, -1, -1):
        C_list[i] = C_list[i+1] * delta_list[i] + lyambda_list[i]

# Вычисляем индекс, по которому будем в дальнейшем определять значение концентрации (для задания 3):
i_pn = int(x_pn * n / x) 

# Вычисляем значение коэффициента диффузии для процесса загонки (для задания 3):
D_zgpn = D(C_list[i_pn],T_zg)

# Передаем функции для построения графиков массивы значений по оси x и y для загонки:

fig, axes = plt.subplots()
axes.plot(x_list, C_list)

# Процесс разгонки:

# Граничные условия:

b_list[0] = 0
d_list[0] = 1
a_list[0] = -1
r_list[0] = 0
delta_list[0] = - d_list[0] / a_list[0]
lyambda_list[0] = r_list[0] / a_list[0]

b_list[n-1] = 0
d_list[n-1] = 0
a_list[n-1] = 1
r_list[n-1] = 0

# Запускаем основной алгоритм для разгонки:

for j in range (0, t_rzg):
    for i in range(1, n - 1):
        b_list[i] = 1
        d_list[i] = 1
        a_list[i] = - (2 + (dx ** 2) / (D(C_list[i], T_rzg) * dt))
        r_list[i] = - dx ** 2 / (D(C_list[i], T_rzg) * dt) * C_list[i]
    for i in range (1, n):
        delta_list[i] = - d_list[i] / (a_list[i] + b_list[i] * delta_list[i-1])
        lyambda_list[i] = (r_list[i] - b_list[i] * lyambda_list[i-1]) / (a_list[i] + b_list[i] * delta_list[i-1])
    C_list[n-1] = lyambda_list[n-1]
    for i in range(n-2, -1, -1):
        C_list[i] = C_list[i+1] * delta_list[i] + lyambda_list[i]
    for i in range(0, n):
        pn_list[i] = np.absolute(C_P - C_list[i])

# Рассчитываем глубину залегания p-n-перехода:

min_pn = np.argmin(pn_list)
print('Глубина залегания p-n перехода =', round(x_list[min_pn] * 1E+4, 3), 'мкм')

# Определяем режимы диффузионных процессов 

# Режим разгонки:

t_rzgpn = (x_pn**2/(1E+8))/(4*D(C_list[i_pn],T_rzg)*np.log(C0/C_P))
print('Время разгонки =', t_rzgpn, ', с')

# Режим загонки:

t_zgpn = D(C_list[i_pn],T_rzg)/D_zgpn*(3.14*C0/(2*C_list[i_pn]))**2*t_rzgpn
print('Время загонки =', t_zgpn, ', с')


# Передаем функции для построения графиков массивы значений по оси x и y для разгонки:

axes.plot(x_list, C_list)

# Ограничиваем длину осей на графике:

axes.set_xlim(0, x * 1E-4)
axes.set_ylim(0, C0)

# Назначаем название для графика, подписываем оси:

plt.title('Распеделение концентрации по глубине залегания', pad = 20)
plt.xlabel('х, см')
plt.ylabel('С, см⁻³')

# Добавляем сетки на график:

axes.grid(which='major', color = '#666666')
axes.minorticks_on()
axes.grid(which='minor', color = 'gray', linestyle = ':')

# Выводим график на экран:

plt.show()