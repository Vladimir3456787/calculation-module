import copy
import math
import time

import numpy as np
from PyQt5.QtWidgets import QInputDialog, QMessageBox
from matplotlib import patches
from shapely import Point, LineString

from Data import Data
from Visualizer import DrawFigure
from MongoConnetion import MongoConnection, uri


class Calculate:
    def __init__(self):
        self.data = Data()
        self.triangle_coordinates = None
        self.info_layers = None
        self.main_window = None
        self.draw_graph = DrawFigure()
        self.mongo = MongoConnection(uri)

    def set_main_window(self, main_window):
        self.main_window = main_window

    def distance_between_points(self, x1, y1, x2, y2):

        return ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5

    def get_data(self):
        self.triangle_coordinates = self.data.get_triangle_coordinates()
        self.info_layers = self.data.get_area_info()

    # def calculate(self, draw, main_window=None):
    #     if self.main_window is None:
    #         self.set_main_window(main_window)
    #     areas = self.mongo.get_areas()
    #     while len(areas) == 0:
    #         time.sleep(1)
    #         areas = self.mongo.get_areas()
    #     areas_number = self.mongo.get_number_areas()
    #     for index, area in enumerate(areas_number, start=0):
    #         if area == self.main_window.combobox_area.currentText():
    #             self.data.change_triangle_coordinates(areas[index])
    #     self.get_data()
    #     if self.triangle_coordinates:
    #         cells = []
    #         for layer in self.info_layers:
    #             x0_layer = layer['x0']
    #             x1_layer = layer['x1']
    #             y0_layer = layer['y0']
    #             y1_layer = layer['y1']
    #             characteristics = layer['characteristics']
    #             m = 1000
    #             S = abs(abs(y1_layer) - abs(y0_layer)) * abs(abs(x1_layer) - abs(x0_layer))
    #             for index, coordinates in enumerate(self.triangle_coordinates):
    #                 if coordinates['x0'] >= x0_layer and coordinates['y0'] >= y0_layer and coordinates[
    #                     'x1'] <= x1_layer and coordinates['y1'] <= y1_layer and x0_layer <= coordinates[
    #                     'x2'] <= x1_layer and y0_layer <= \
    #                         coordinates['y2'] <= y1_layer:
    #                     a = ((coordinates['x1'] - coordinates['x0']) ** 2 + (
    #                                 coordinates['y1'] - coordinates['y0']) ** 2) ** 0.5
    #                     b = ((coordinates['x2'] - coordinates['x1']) ** 2 + (
    #                                 coordinates['y2'] - coordinates['y1']) ** 2) ** 0.5
    #                     c = ((coordinates['x2'] - coordinates['x0']) ** 2 + (
    #                                 coordinates['y2'] - coordinates['y0']) ** 2) ** 0.5
    #                     s = (a + b + c) / 2
    #                     # Площадь треугольника по формуле Герона
    #                     S_triangle = (s * (s - a) * (s - b) * (s - c)) ** 0.5
    #                     V_triangle = S_triangle
    #                     m = float(layer['density']) * V_triangle
    #                     p_total = m / V_triangle
    #                     air_density = 1.225  # Плотность воздуха при нормальных условиях в кг/м^3
    #                     air_v = V_triangle * (1 - float(layer['density']) / air_density)
    #                     pd = m / (V_triangle + air_v)
    #                     ps = m / V_triangle
    #                     n = 1 - pd / ps
    #                     E = n / m
    #                     fi = 0
    #                     c = 0
    #                     for ind, characteristic in enumerate(characteristics):
    #                         if E <= characteristic[0]:
    #                             fi = characteristic[1]
    #                             c = characteristic[2]
    #                     if fi == 0:
    #                         fi = characteristics[-1][1]
    #                         c = characteristics[-1][2]
    #                     lambd = pd * 9.81
    #                     h = 2 * S_triangle / a
    #                     sigma = lambd * h
    #                     sigma_rad = math.radians(sigma)
    #                     fi_rad = math.radians(fi)
    #                     u_x = sigma_rad * math.cos(fi_rad) * h / E  # смещение по x
    #                     u_y = sigma_rad * math.sin(fi_rad) * h / E  # смещение по y
    #                     cell_info = {
    #                         'name': layer['name'],
    #                         'color': layer['color'],
    #                         'x0': coordinates['x0'],
    #                         'x1': coordinates['x1'],
    #                         'x2': coordinates['x2'],
    #                         'y0': coordinates['y0'],
    #                         'y1': coordinates['y1'],
    #                         'y2': coordinates['y2'],
    #                         'density': layer['density'],
    #                         'middleX': coordinates['middleX'],
    #                         'middleY': coordinates['middleY'],
    #                         'number': coordinates['number'],
    #                         'napr': sigma,
    #                         'destruction_threshould': ((sigma_rad + c * math.degrees(math.atan(fi_rad))) * (
    #                                     math.tan(math.radians(45 + fi / 2)) ** 2) - c * math.degrees(math.atan(fi_rad)))
    #
    #                     }
    #                     cells.append(cell_info)
    #         # return cells
    #         if draw == 1:
    #             self.data.change_check_calculate(True)
    #             self.draw_graph.draw_subsidence(self.main_window, cells)
    #         else:
    #             calc = Calculate()
    #             calc.calcluate_stress(cells)
    def calculate_stress_cells(self, main_window=None):
        try:

            cells = []
            self.get_data()
            self.set_main_window(main_window)
            a = 30
            vertical_lines = self.data.get_vertical_lines()
            self.get_data()
            self.set_main_window(main_window)
            m = 0  # мощность угля
            m_all = []

            for index, layer in enumerate(self.data.get_info_layers(), start=0):
                if ('уголь') in layer['name'] and layer['number'] == 1:
                    m = layer['height']
                if not ('уголь') in layer['name'] and layer['number'] == 1:
                    m_all.append(layer['height'])
            h = sum(m_all)
            f = []
            qcop = 9.8

            b = 1  # решили, что берем 1 метр по ширине
            ysum = 0

            limitska = {
                (0, 20): 0.70,
                (20, 30): 0.6,
                (30, 40): 0.45,
                (40, 50): 0.25,
                (50, 60): 0.2,
                (60, 70): 0.15,
                (70, 1000): 0.15,
            }
            limitskq = {
                (0, 20): 0.55,
                (20, 30): 0.8,
                (30, 40): 0.95,
                (40, 50): 0.95,
                (50, 60): 0.95,
                (60, 70): 0.95,
                (70, 1000): 0.95,
            }
            kv = 1  # для одиночной выработки, когда есть только один пласт угля
            ks = 1  # для одиночной выработки, когда есть только один пласт угля

            # при расчетном сопротивление
            limitskL = {
                (0, 30): 1.8,
                (30, 60): 1.5,
                (60, 90): 1.2,
                (90, 120): 1,
                (120, 10000): 1,
            }

            collection_areas = self.data.get_areas_info()
            if len(collection_areas) == 0:
                print("No collection")
                return None
            areas_number = self.data.get_areas_number()
            areas = self.data.get_areas_info()
            for index, area in enumerate(areas_number, start=0):
                if area == self.main_window.combobox_area.currentText():
                    self.data.change_triangle_coordinates(areas[index])
                    areas = areas[index]
            power = 10 ** 7

            q2 = 40000  # тоже из бд, это пока пусть так будет, так как я еще не знаю как это связать с крепью
            collection_rocks = self.mongo.get_collection_density()
            H = None
            L1 = None
            nds_triangles = []
            p_float = None
            p_v = None
            characteristics = None
            for index, coordinates in enumerate(areas, start=0):
                if not 'уголь' in coordinates:
                    for index_name, rock in enumerate(collection_rocks, start=0):
                        if coordinates['name'] == rock['Породы']:
                            f = rock[
                                'f']  # крепости пород из мощностей,нужно сделать привязку к выбору точек (то есть в каком находится пласт, ту и крепость берем)
                            # расчет удельного веса пород (она должна зависит от нахождения точек, то есть плотность будет разная)
                            p_float = float(rock['p'])
                            p_v = float(rock['p_v'])
                            ysum = p_float * qcop  # расчет удельного веса пород
                            break
                    area_info = self.data.get_area_info()

                    if area_info:
                        for layer in area_info:
                            if layer['name'] == coordinates['name']:
                                characteristics = layer['characteristics']
                                break

                    for limit_range, index_lim in limitska.items():
                        lower_limit, upper_limit = limit_range
                        if lower_limit < a <= upper_limit:
                            ka = index_lim
                            break

                    for limit_range, index_lim in limitskq.items():
                        lower_limit, upper_limit = limit_range
                        if lower_limit < a <= upper_limit:
                            kq = index_lim
                            break

                    for v_line in vertical_lines:
                        if v_line['x0'] == coordinates['x0'] or v_line['x0'] == coordinates['x1'] or v_line['x0'] == \
                                coordinates['x2']:
                            if v_line['y0'] > v_line['y1']:
                                H = abs(abs(v_line['y0']) - abs(
                                    coordinates['y' + str(0)]))  # мощность всех слоев сверху
                                L1 = abs(coordinates['y' + str(0)]) + abs(
                                    v_line['y1'])  # мощность от точки до начала угля
                                H = round(H, 4)
                                L1 = round(L1, 4)
                            if v_line['y0'] < v_line['y1']:
                                H = abs(abs(v_line['y1']) - abs(
                                    coordinates['y' + str(0)]))  # мощность всех слоев сверху
                                L1 = abs(coordinates['y' + str(0)]) + abs(
                                    v_line['y0'])  # мощность от точки до начала угля
                                H = round(H, 4)
                                L1 = round(L1, 4)

                    Q_f = f * power  # значение Q для текущего индекса f
                    Ut = (63.12 * m + 76.92 * math.exp(
                        (4.2 * ysum * round(H, 4)) / Q_f) - 9.8 * b + 33.12 * math.log(
                        round(L1, 4)) - 102) / 1000

                    Rc = f * 10
                    for limit_range, index in limitskL.items():
                        lower_limit, upper_limit = limit_range
                        if lower_limit < Rc <= upper_limit:
                            kL = index
                            break

                    U = Ut * ka * kq * ks * kv * kL
                    coordinates['y' + str(0)] = coordinates['y' + str(0)] - U

                    ####Напряжение####
                    sigma = ysum * H

                    ####Давление####
                    if self.main_window is None:
                        self.set_main_window(main_window)
                    self.get_data()
                    areas_number = self.data.get_areas_number()
                    areas = self.data.get_areas_info()
                    areas = copy.deepcopy(areas)
                    for index, area in enumerate(areas_number, start=0):
                        if area == self.main_window.combobox_area.currentText():
                            # self.data.change_triangle_coordinates(areas[index])
                            areas = areas[index]
                            break

                    sigma_rad = math.radians(sigma)
                    n = 1 - p_float / p_v
                    if H != 0:
                        E = n / H
                    else:
                        E = n
                        # raise ValueError("H is zero, division by zero encountered.")
                    fi = 0
                    c = 0
                    for ind, characteristic in enumerate(characteristics):
                        if E <= characteristic[0]:
                            fi = characteristic[1]
                            c = characteristic[2]
                    if fi == 0:
                        fi = characteristics[-1][1]
                        c = characteristics[-1][2]
                    fi_rad = math.radians(fi)
                    sigma_dav = ((sigma_rad + c * math.degrees(math.atan(fi_rad))) * abs((
                                                                                                 math.tan(math.radians(
                                                                                                     45 + fi / 2)) ** 2) - c * math.degrees(
                        math.atan(fi_rad))))
                    cell_info = {
                        'name': layer['name'],
                        'color': layer['color'],
                        'x0': coordinates['x0'],
                        'x1': coordinates['x1'],
                        'x2': coordinates['x2'],
                        'y0': coordinates['y0'],
                        'y1': coordinates['y1'],
                        'y2': coordinates['y2'],
                        'density': layer['density'],
                        'middleX': coordinates['middleX'],
                        'middleY': coordinates['middleY'],
                        'number': coordinates['number'],
                        'napr': round(sigma_dav / 1000, 3),
                        'destruction_threshould': ((sigma_rad + c * math.degrees(math.atan(fi_rad))) * (
                                math.tan(math.radians(45 + fi / 2)) ** 2) - c * math.degrees(math.atan(fi_rad)))

                    }
                    # cells.append(cell_info)
                    self.data.add_rock_stress(cell_info)
            self.draw_graph.draw_gradient(main_window)


            # self.calcluate_stress(cells, main_window)
        except Exception as e:
            QMessageBox.critical(None, "Ошибка", str(e))

    def check_condition(self, elem, current_x0, current_y0, current_x1, current_y1, current_x2, current_y2):
        x0_in_region = elem['x0'] in [current_x1, current_x2] and ((elem['y0'] > current_y0 or elem['y0'] > current_y1 or elem['y0'] > current_y2) or (elem['y1'] > current_y0 or elem['y1'] > current_y1 or elem['y1'] > current_y2)or (elem['y2'] > current_y0 or elem['y2'] > current_y1 or elem['y2'] > current_y2))
        # x1_in_region = elem['x1'] in [current_x0, current_x2] and ((elem['y0'] > current_y0 or elem['y0'] > current_y1 or elem['y0'] > current_y2) or (elem['y1'] > current_y0 or elem['y1'] > current_y1 or elem['y1'] > current_y2)or (elem['y2'] > current_y0 or elem['y2'] > current_y1 or elem['y2'] > current_y2))
        # x2_in_region = elem['x2'] in [current_x0, current_x1] and ((elem['y0'] > current_y0 or elem['y0'] > current_y1 or elem['y0'] > current_y2) or (elem['y1'] > current_y0 or elem['y1'] > current_y1 or elem['y1'] > current_y2)or (elem['y2'] > current_y0 or elem['y2'] > current_y1 or elem['y2'] > current_y2))
        # x3_in_region = elem['x0'] in [current_x1, current_x2] and (
        #             (elem['y0'] > current_y0 or elem['y0'] > current_y1 or elem['y0'] > current_y2) or (
        #                 elem['y1'] > current_y0 or elem['y1'] > current_y1 or elem['y1'] > current_y2) or (
        #                         elem['y2'] > current_y0 or elem['y2'] > current_y1 or elem['y2'] > current_y2))
        # x4_in_region = elem['x1'] in [current_x0, current_x2] and (
        #             (elem['y0'] > current_y0 or elem['y0'] > current_y1 or elem['y0'] > current_y2) or (
        #                 elem['y1'] > current_y0 or elem['y1'] > current_y1 or elem['y1'] > current_y2) or (
        #                         elem['y2'] > current_y0 or elem['y2'] > current_y1 or elem['y2'] > current_y2))
        # x5_in_region = elem['x2'] in [current_x0, current_x1] and (
        #             (elem['y0'] > current_y0 or elem['y0'] > current_y1 or elem['y0'] > current_y2) or (
        #                 elem['y1'] > current_y0 or elem['y1'] > current_y1 or elem['y1'] > current_y2) or (
        #                         elem['y2'] > current_y0 or elem['y2'] > current_y1 or elem['y2'] > current_y2))

        return x0_in_region
    def calcluate_stress(self, cells, main_window=None):
        try:

            cells_sorted = sorted(cells, key=lambda x: x['y2'])
            for index, cell_info in enumerate(cells_sorted, start=0):

                current_x0 = cell_info['x0']
                current_y0 = cell_info['y0']
                current_x1 = cell_info['x1']
                current_y1 = cell_info['y1']
                current_x2 = cell_info['x2']
                current_y2 = cell_info['y2']
                total_napr = 0



                for i, elem in enumerate(cells_sorted, start=0):
                    is_hits = self.check_condition(elem, current_x0, current_y0, current_x1, current_y1, current_x2, current_y2)
                    if is_hits:
                        total_napr += elem['napr']

                cell_for_gradient = {
                    'name': cell_info['name'],
                    'color': cell_info['color'],
                    'x0': cell_info['x0'],
                    'x1': cell_info['x1'],
                    'x2': cell_info['x2'],
                    'y0': cell_info['y0'],
                    'y1': cell_info['y1'],
                    'y2': cell_info['y2'],
                    'density': cell_info['density'],
                    'middleX': cell_info['middleX'],
                    'middleY': cell_info['middleY'],
                    'number': cell_info['number'],
                    'napr': cell_info['napr'],
                    'destruction_threshould': cell_info['destruction_threshould'],
                    'total_napr': total_napr
                }
                self.data.add_rock_stress(cell_for_gradient)
            self.draw_graph.draw_gradient(main_window)
        except Exception as e:
            QMessageBox.critical(None, "Ошибка", str(e))

    def predict_coordinates(self, x, y):

        a = np.sum((x - np.mean(x)) * (y - np.mean(y))) / np.sum((x - np.mean(x)) ** 2)
        b = np.mean(y) - a * np.mean(x)

        # Прогнозируем значения y для каждой координаты x
        y_predicted = [a * x_coord + b for x_coord in x]

        # Создаем массив прогнозируемых координат (x, y_predicted)
        column_info = {"X": x, "Y": y_predicted}

        return column_info

    def calculate_analytycs(self, stress=False, number=1, q=0.0, a=0.0):
        try:
            L= 0
            collection_analytics = self.mongo.get_collection_analytics()
            # q = 0.75  # заменить на ввод из переменной q пользователя первичная или повторная
            area_info = self.data.get_area_info()
            for layer in area_info:
                if layer['number'] == number:
                    L = layer['x0']+layer['x1']
                    break
            m = 0  # Добавить переменную m1 из ввода угля
            # a = 0.3  # Добавить переменную a для аналитической формулы (ввод пользователя)
            m1 = []

            for index, layer in enumerate(self.data.get_info_layers(), start=0):
                if ('уголь') in layer['name'] and layer['number'] == number:
                    m = layer['height']
                else:
                    m1.append(layer['height'])

            if stress:
                f = [6, 8, 9, 10]
                L1 = 0
                nsk = []
                Nosdtk = []
                max_Nosdk = float('-inf')
                H = m + sum(m1)
            else:
                H = m + sum(m1)  # Сумма всех значений массива m1 (то есть мощности когда в массив будешь выводить)
                L1 = 0
                ns = []
                Nosdt = []
                max_Nosd = float('-inf')

            limits = {
                (0, 0.05): 0,
                (0.05, 0.1): 1,
                (0.1, 0.2): 2,
                (0.2, 0.3): 3,
                (0.3, 0.4): 4,
                (0.4, 0.5): 5,
                (0.5, 0.6): 6,
                (0.6, 0.7): 7,
                (0.7, 0.8): 8,
                (0.8, 0.9): 9,
                (0.9, 1): 10,
                (1, 1.1): 11,
                (1.1, 100): 12,
            }
            min_y = []
            for doc in collection_analytics.find({}):
                N1_N2_array = doc.get('N1_N2', [])
                S = doc.get('S1(z)', [])

                for i in range(int(L)):
                    L1 += 1
                    N1_N2 = L1 / H

                    if stress:
                        P = 1 / H * sum([m * f[index] for index, m in enumerate(m1)])

                    for limit_range, index in limits.items():
                        lower_limit, upper_limit = limit_range
                        if lower_limit < N1_N2 <= upper_limit:
                            if stress:
                                if index <= len(N1_N2_array) - 1:
                                    nk = 1908 * q * (m * math.cos(a)) ** (0.2) * (N1_N2_array[index] * N1_N2_array[index]) ** (
                                        0.3) * P ** (-4)
                                    nsk.append(nk)
                                    for j in range(1):
                                        Nosdk = nk * S[0]  # Changed S to S[j])
                                        max_Nosdk = max(max_Nosdk, Nosdk)
                                        print(f"оседание {j + 1}: {Nosdk}")
                                        Nosdtk.append(max_Nosdk)
                                        min_y.append(max_Nosdk)
                                    break
                            else:
                                if index <= len(N1_N2_array)-1:
                                    n = q * m * math.cos(a) * N1_N2_array[index] * N1_N2_array[index]
                                    ns.append(n)
                                    for j in range(1):
                                        Nosd = n * S[0]  # Changed S to S[j])
                                        max_Nosd = max(max_Nosd, Nosd)
                                        print(f"оседание {j + 1}: {Nosd}")
                                        Nosdt.append(max_Nosd)
                                        min_y.append(max_Nosd)
                                    break
            if min_y:
                max_subsidence = min(min_y)
            else:
                max_subsidence = 0
            if stress:
                return Nosdtk, max_subsidence
            else:
                return Nosdt, max_subsidence
        except Exception as e:
            QMessageBox.critical(None, "Ошибка", str(e))

    def subsidence_analytics_calculate(self, main_window, stress=False):
        try:
            q, ok_pressed = QInputDialog.getText(self.main_window, "Разработка", "Введите значение разработки\n(первичная "
                                                                                 "0.75/вторичная 0.85):")
            if ok_pressed:
                try:
                    q = float(q)
                except ValueError:
                    msg = QMessageBox()
                    msg.setWindowTitle("Внимание")
                    msg.setText("Возникли проблемы!")
                    msg.setIcon(QMessageBox.Warning)
                    msg.exec_()
                    return

            a, ok_pressed = QInputDialog.getText(self.main_window, "Угол наклона", "Введите угол наклона в градусах:")
            if ok_pressed and a.isdigit():
                a = float(a) / 100
            else:
                msg = QMessageBox()
                msg.setWindowTitle("Внимание")
                msg.setText("Значение не введено!")
                msg.setIcon(QMessageBox.Warning)
                msg.exec_()
                return
            columns = self.data.get_current_column_1()
            subsidences = []
            max_subsidence = 0
            if stress:
                for column in columns:
                    result = self.calculate_analytycs(True, column['Number'], q, a)
                    max_subsidence = result[1]
                    subsidences.append(result[0])
            else:
                for column in columns:
                    result =self.calculate_analytycs(False, column['Number'], q, a)
                    max_subsidence = result[1]
                    subsidences.append(result[0])

            point_subsidence_x = []
            point_subsidence_y = []
            max_subsidence = 0
            y_max_1 = 0
            y_max_2 = 0
            x = -1
            area_info = self.data.get_area_info()
            for column in columns:
                for layer in area_info:
                    if layer['number'] == column['Number'] and not ('уголь') in layer['name']:
                        y_max_1 = y_max_1 + layer['height']
                        y_max_2 = y_max_2 + layer['height']
                        break
            upper_lines = self.data.get_high_horizontal_lines()

            for index, column in enumerate(columns):
                if index > len(upper_lines)-1:
                    break
                x = upper_lines[index]['x0']
                for index_0, subsidence in enumerate(subsidences[index]):
                    x += 1
                    line_for_intersections = [(x, 0), (x,9000)]
                    line_coords = [(upper_lines[index]['x0'], upper_lines[index]['y0']), (upper_lines[index]['x1'], upper_lines[index]['y1'])]
                    line1 = LineString(line_for_intersections)
                    line2 = LineString(line_coords)

                    intersection = line1.intersection(line2)

                    # Проверяем, пересекается ли точка с линией
                    if max_subsidence < subsidence:
                        max_subsidence = subsidence
                    if intersection:
                        x_intersect, y_intersect = intersection.xy
                        y_max = y_intersect[0]
                        point_subsidence_x.append(x_intersect)
                        point_subsidence_y.append(y_max-max_subsidence)
                        continue


            self.draw_graph.draw_analytics_calculate(main_window, point_subsidence_x, point_subsidence_y, max_subsidence)
        except Exception as e:
            QMessageBox.critical(None, "Ошибка", str(e))

    def calculate_NDS(self, main_window, max_x=0.0):
        try:
            vertical_lines = self.data.get_vertical_lines()
            self.get_data()
            self.set_main_window(main_window)
            a = self.data.get_angle()
            m = 0  # мощность угля
            m_all = []

            for index, layer in enumerate(self.data.get_info_layers(), start=0):
                if ('уголь') in layer['name'] and layer['number'] == 1:
                    m = layer['height']
                if not ('уголь') in layer['name'] and layer['number'] == 1:
                    m_all.append(layer['height'])
            h = sum(m_all)
            f = []
            qcop = 9.8

            b = 1  # решили, что берем 1 метр по ширине
            ysum = 0

            limitska = {
                (0, 20): 0.70,
                (20, 30): 0.6,
                (30, 40): 0.45,
                (40, 50): 0.25,
                (50, 60): 0.2,
                (60, 70): 0.15,
                (70, 1000): 0.15,
            }
            limitskq = {
                (0, 20): 0.55,
                (20, 30): 0.8,
                (30, 40): 0.95,
                (40, 50): 0.95,
                (50, 60): 0.95,
                (60, 70): 0.95,
                (70, 1000): 0.95,
            }
            kv = 1  # для одиночной выработки, когда есть только один пласт угля
            ks = 1  # для одиночной выработки, когда есть только один пласт угля

            # при расчетном сопротивление
            limitskL = {
                (0, 30): 1.8,
                (30, 60): 1.5,
                (60, 90): 1.2,
                (90, 120): 1,
                (120, 10000): 1,
            }
            max_ut = 0
            array_ut = []
            collection_areas = self.data.get_areas_info()
            if len(collection_areas) == 0:
                print("No collection")
                return None
            areas_number = self.data.get_areas_number()
            areas = self.data.get_areas_info()
            areas = copy.deepcopy(areas)
            for index, area in enumerate(areas_number, start=0):
                if area == self.main_window.combobox_area.currentText():
                    # self.data.change_triangle_coordinates(areas[index])
                    areas = areas[index]
            power = 10 ** 7

            q2 = 40000  # тоже из бд, это пока пусть так будет, так как я еще не знаю как это связать с крепью
            collection_rocks = self.mongo.get_collection_density()
            H = None
            L1 = None
            nds_triangles = []
            for index, coordinates in enumerate(areas, start=0):
                if not 'уголь' in coordinates['name']:
                    if coordinates['x2'] <= max_x and coordinates['x1'] <= max_x and coordinates['x0'] <= max_x:
                        for index_name, rock in enumerate(collection_rocks, start=0):
                            if coordinates['name'] == rock['Породы']:
                                f = rock[
                                    'f']  # крепости пород из мощностей,нужно сделать привязку к выбору точек (то есть в каком находится пласт, ту и крепость берем)
                                # расчет удельного веса пород (она должна зависит от нахождения точек, то есть плотность будет разная)
                                p_float = float(rock['p'])
                                ysum = p_float * qcop  # расчет удельного веса пород
                                break

                        for limit_range, index in limitska.items():
                            lower_limit, upper_limit = limit_range
                            if lower_limit < a <= upper_limit:
                                ka = index
                                break

                        for limit_range, index in limitskq.items():
                            lower_limit, upper_limit = limit_range
                            if lower_limit < a <= upper_limit:
                                kq = index
                                break

                        for i in range(3):  # три раза потому что 3 Y
                            for v_line in vertical_lines:
                                if v_line['x0'] == coordinates['x0'] or v_line['x0'] == coordinates['x1'] or v_line['x0'] == \
                                        coordinates['x2']:
                                    if v_line['y0'] > v_line['y1']:
                                        H = abs(abs(v_line['y0']) - abs(
                                            coordinates['y' + str(i)]))  # мощность всех слоев сверху
                                        L1 = abs(coordinates['y' + str(i)]) + abs(
                                            v_line['y1'])  # мощность от точки до начала угля
                                        H = round(H, 4)
                                        L1 = round(L1, 4)
                                    if v_line['y0'] < v_line['y1']:
                                        H = abs(abs(v_line['y1']) - abs(
                                            coordinates['y' + str(i)]))  # мощность всех слоев сверху
                                        L1 = abs(coordinates['y' + str(i)]) + abs(
                                            v_line['y0'])  # мощность от точки до начала угля
                                        H = round(H, 4)
                                        L1 = round(L1, 4)

                            Q_f = f * power  # значение Q для текущего индекса f
                            Ut = (63.12 * m + 76.92 * math.exp(
                                (4.2 * ysum * round(H, 4)) / Q_f) - 9.8 * b + 33.12 * math.log(
                                round(L1, 4)) - 102) / 1000

                            Rc = f * 10
                            for limit_range, index in limitskL.items():
                                lower_limit, upper_limit = limit_range
                                if lower_limit < Rc <= upper_limit:
                                    kL = index
                                    break

                            U = Ut * ka * kq * ks * kv * kL
                            coordinates['y' + str(i)] = coordinates['y' + str(i)] - U
                            array_ut.append(Ut)
                            print(Ut)

                        nds_triangles.append(coordinates)
                    else:
                        nds_triangles.append(coordinates)
                        #####это на будующие с учетом крепей
                        # Ukrep = (57+36.43*m**3+5352/q2+1128*(ysum*H/Qsum)**2-4230/L1)/1000
                        # Ukreps.append(Ukrep)
                        # print(Ukreps)

            if array_ut:
                max_ut = max(array_ut)
            else:
                max_ut = 0
            return nds_triangles, max_ut
        except Exception as e:
            QMessageBox.critical(None, "Ошибка", str(e))

    def calculate_NDS_analytyc(self, main_window):

        a, ok_pressed = QInputDialog.getText(self.main_window, "Угол наклона", "Введите угол наклона в градусах:")
        if ok_pressed and a.isdigit():
            a = float(a)
        else:
            msg = QMessageBox()
            msg.setWindowTitle("Внимание")
            msg.setText("Значение не введено!")
            msg.setIcon(QMessageBox.Warning)
            msg.exec_()
            return
        vertical_lines = self.data.get_vertical_lines()
        self.get_data()
        self.set_main_window(main_window)
        m = 0  # мощность угля
        m_all = []

        for index, layer in enumerate(self.data.get_info_layers(), start=0):
            if ('уголь') in layer['name'] and layer['number'] == 1:
                m = layer['height']
            if not ('уголь') in layer['name'] and layer['number'] == 1:
                m_all.append(layer['height'])
        h = sum(m_all)
        f = []
        qcop = 9.8

        b = 1  # решили, что берем 1 метр по ширине
        ysum = 0

        limitska = {
            (0, 20): 0.70,
            (20, 30): 0.6,
            (30, 40): 0.45,
            (40, 50): 0.25,
            (50, 60): 0.2,
            (60, 70): 0.15,
            (70, 1000): 0.15,
        }
        limitskq = {
            (0, 20): 0.55,
            (20, 30): 0.8,
            (30, 40): 0.95,
            (40, 50): 0.95,
            (50, 60): 0.95,
            (60, 70): 0.95,
            (70, 1000): 0.95,
        }
        kv = 1  # для одиночной выработки, когда есть только один пласт угля
        ks = 1  # для одиночной выработки, когда есть только один пласт угля

        # при расчетном сопротивление
        limitskL = {
            (0, 30): 1.8,
            (30, 60): 1.5,
            (60, 90): 1.2,
            (90, 120): 1,
            (120, 10000): 1,
        }

        collection_areas = self.data.get_areas_info()
        if len(collection_areas) == 0:
            print("No collection")
            return None
        areas_number = self.data.get_areas_number()
        areas = self.data.get_areas_info()
        for index, area in enumerate(areas_number, start=0):
            if area == self.main_window.combobox_area.currentText():
                self.data.change_triangle_coordinates(areas[index])
                areas = areas[index]
        power = 10 ** 7

        q2 = 40000  # тоже из бд, это пока пусть так будет, так как я еще не знаю как это связать с крепью
        collection_rocks = self.mongo.get_collection_density()
        H = None
        L1 = None
        nds_triangles = []
        for index, coordinates in enumerate(areas, start=0):

            if not 'уголь' in coordinates['name']:

                for index_name, rock in enumerate(collection_rocks, start=0):
                    if coordinates['name'] == rock['Породы']:
                        f = rock[
                            'f']  # крепости пород из мощностей,нужно сделать привязку к выбору точек (то есть в каком находится пласт, ту и крепость берем)
                        # расчет удельного веса пород (она должна зависит от нахождения точек, то есть плотность будет разная)
                        p_float = float(rock['p'])
                        ysum = p_float * qcop  # расчет удельного веса пород
                        break

                for limit_range, index in limitska.items():
                    lower_limit, upper_limit = limit_range
                    if lower_limit < a <= upper_limit:
                        ka = index
                        break

                for limit_range, index in limitskq.items():
                    lower_limit, upper_limit = limit_range
                    if lower_limit < a <= upper_limit:
                        kq = index
                        break

                for i in range(3):  # три раза потому что 3 Y
                    for v_line in vertical_lines:
                        if v_line['x0'] == coordinates['x0'] or v_line['x0'] == coordinates['x1'] or v_line['x0'] == coordinates['x2']:
                            if v_line['y0'] > v_line['y1']:
                                H = abs(abs(v_line['y0']) - abs(coordinates['y' + str(i)])) # мощность всех слоев сверху
                                L1 = abs(coordinates['y' + str(i)]) + abs(v_line['y1'])  # мощность от точки до начала угля
                                H = round(H, 4)
                                L1 = round(L1, 4)
                            if v_line['y0'] < v_line['y1']:
                                H = abs(abs(v_line['y1']) - abs(coordinates['y' + str(i)])) # мощность всех слоев сверху
                                L1 = abs(coordinates['y' + str(i)]) + abs(v_line['y0'])  # мощность от точки до начала угля
                                H = round(H, 4)
                                L1 = round(L1, 4)


                    Q_f = f * power  # значение Q для текущего индекса f
                    Ut = (63.12 * m + 76.92 * math.exp((4.2 * ysum * round(H,4)) / Q_f) - 9.8 * b + 33.12 * math.log(
                        round(L1,4)) - 102) / 1000

                    Rc = f * 10
                    for limit_range, index in limitskL.items():
                        lower_limit, upper_limit = limit_range
                        if lower_limit < Rc <= upper_limit:
                            kL = index
                            break

                    U = Ut * ka * kq * ks * kv * kL
                    coordinates['y' + str(i)] = coordinates['y' + str(i)] - U

                    print(Ut)

                nds_triangles.append(coordinates)
                self.data.add_changed_triangles(coordinates)
                #####это на будующие с учетом крепей
                # Ukrep = (57+36.43*m**3+5352/q2+1128*(ysum*H/Qsum)**2-4230/L1)/1000
                # Ukreps.append(Ukrep)
                # print(Ukreps)

        return nds_triangles
