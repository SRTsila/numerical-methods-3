import math
from copy import deepcopy

matrix_const = [[0.10, 1.51, -0.2],
                [-0.1, -0.10, 1.00],
                [-0.50, -0.30, 0.20]]
values_const = [1.41, 1.60, -1.40]

matrix_for_3_and_4_tasks = [[-0.50, -0.30, 0.20, -1.4],
                            [0.10, 1.51, -0.2, 1.41],
                            [-0.1, -0.10, 1.00, 1.6]]


def gauss_method(matrix, value_col):
    for i in range(len(matrix)):
        a = matrix[i][i]
        for j in range(i, len(matrix)):
            matrix[i][j] /= a
        value_col[i] /= a
        for j in range(i + 1, len(matrix)):
            c = -matrix[j][i]
            for z in range(i, len(matrix)):
                matrix[j][z] += c * matrix[i][z]
            value_col[j] += c * value_col[i]

    for i in range(len(matrix) - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            c = -matrix[j][i]
            matrix[j][i] += c * matrix[i][i]
            value_col[j] += c * value_col[i]
    return {"x1": value_col[0], "x2": value_col[1], "x3": value_col[2]}


def gauss_with_main_element_method(matrix):
    order = {0: None, 1: None, 2: None}
    res = []
    max = 0
    for i in range(3):
        for j in range(3):
            if abs(matrix[i][j]) > max:
                max = abs(matrix[i][j])
                order[0] = i
                res.append(i)
    max = 0
    for i in range(3):
        if i == order[0]:
            continue
        for j in range(3):
            if abs(matrix[i][j]) > max:
                max = abs(matrix[i][j])
                order[1] = i
                res.append(i)
    res = [order[0], order[1]]
    for i in range(3):
        if i not in res:
            order[2] = i
            break
    new_matrix = [matrix[order[0]], matrix[order[1]], matrix[order[2]]]
    return gauss_method(new_matrix, deepcopy(values_const))


def jacobi_method(matrix):
    solution = tuple((1 for _ in range(len(matrix))))

    counter = 1
    accuracy = 0.5 * 10e-4
    while accuracy >= 0.5 * 10e-4:
        new_solution = []
        for i in range(len(matrix)):
            line = matrix[i]
            line_sum = line[len(line) - 1]
            for j in range(len(line) - 1):
                if i != j:
                    line_sum -= line[j] * solution[j]
            new_solution.append(line_sum / line[i])

        accuracy = get_distance(new_solution, solution)
        solution = new_solution
        counter += 1
    print('Число итераций: ', counter)
    return {"x1": new_solution[0], "x2": new_solution[1], "x3": new_solution[2]}


def gauss_seiedel_method(matrix):
    current_res = tuple((0 for _ in range(len(matrix))))
    accuracy = 0.5 * 10e-4
    counter = 0
    while accuracy >= 0.5 * 10e-4:
        new_solution = [None for _ in range(len(matrix))]
        for i in range(len(matrix)):
            line = matrix[i]

            element = line[-1]
            for j in range(len(line) - 1):
                if i != j:
                    element -= (line[j] * (current_res[j] if new_solution[j] is None else new_solution[j]))

            new_solution[i] = (element / line[i])
        accuracy = get_distance(new_solution, current_res)
        current_res = new_solution
        counter += 1
    print("Число итераций до достижения точности: ", counter)
    return {"x1": new_solution[0], "x2": new_solution[1], "x3": new_solution[2]}


def get_distance(x, y):
    return math.sqrt(sum([(x[i] - y[i]) ** 2 for i in range(len(x))]))


def main():
    matrix = deepcopy(matrix_const)
    res = deepcopy(values_const)
    accuracy_result = [3.0, 1.0, 2.0]

    gauss_method_res = gauss_method(matrix, res)
    matrix = deepcopy(matrix_const)
    res = deepcopy(values_const)
    gaus_method_with_main_element = gauss_with_main_element_method(matrix)
    print("Метод Гаусса:")
    print(gauss_method_res)
    print("Погрешность: ", get_distance(accuracy_result, list(gauss_method_res.values())))
    print()
    print("Метод Гаусса c выбором главного элемента:")
    print(gaus_method_with_main_element)
    print("Погрешность: ", get_distance(accuracy_result, list(gaus_method_with_main_element.values())))
    print()
    print("Метод Якоби")
    print(jacobi_method(deepcopy(matrix_for_3_and_4_tasks)))
    print()
    print("Метод Гаусса-Зейделя")
    print(gauss_seiedel_method(matrix_for_3_and_4_tasks))


if __name__ == "__main__":
    main()
