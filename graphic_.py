from tkinter import *
import re
from fractions import Fraction
from math import gcd
from functools import reduce

# цвета интерфейса (фон, текст, результат)
colors = ['olive', 'lightblue', 'yellow']

# создание окна
root = Tk()
root.title('Решение химических уравнений')
root.geometry('430x210')
root['bg'] = colors[0]


# -------- РАЗБОР ФОРМУЛЫ --------
def parse_formula(formula):
    formula = formula.strip()  # исходная строка формулы

    def parse_inner(s, pos):
        atoms_count = {}  # словарь: элемент -> количество атомов

        while pos < len(s):
            symbol = s[pos]  # текущий символ строки

            if symbol == '(':
                # разбираем выражение в скобках
                inside, pos = parse_inner(s, pos + 1)

                # читаем число после скобки
                num_str = ''
                while pos < len(s) and s[pos].isdigit():
                    num_str += s[pos]
                    pos += 1
                multiplier = int(num_str) if num_str else 1  # множитель

                # умножаем элементы внутри скобок
                for el in inside:
                    atoms_count[el] = atoms_count.get(el, 0) + inside[el] * multiplier

            elif symbol == ')':
                return atoms_count, pos + 1

            elif symbol.isupper():
                # читаем название элемента (например Fe, Ca)
                element = symbol
                pos += 1

                while pos < len(s) and s[pos].islower():
                    element += s[pos]
                    pos += 1

                # читаем число атомов
                num_str = ''
                while pos < len(s) and s[pos].isdigit():
                    num_str += s[pos]
                    pos += 1

                count = int(num_str) if num_str else 1  # количество атомов
                atoms_count[element] = atoms_count.get(element, 0) + count

            else:
                raise ValueError("Ошибка в формуле")

        return atoms_count, pos

    result, _ = parse_inner(formula, 0)
    return result


# -------- РАЗДЕЛЕНИЕ УРАВНЕНИЯ --------
def split_equation(eq):
    if '->' not in eq:
        raise ValueError("Нет стрелки ->")

    left_part, right_part = eq.split('->')

    # убираем возможные коэффициенты перед формулами
    left = [re.sub(r'^\d+', '', x).strip() for x in left_part.split('+')]
    right = [re.sub(r'^\d+', '', x).strip() for x in right_part.split('+')]

    return left, right


# -------- РЕШЕНИЕ СИСТЕМЫ --------
def solve_matrix(matrix):
    matrix = [row[:] for row in matrix]  # копия матрицы
    rows_count = len(matrix)             # количество строк (элементов)
    cols_count = len(matrix[0])          # количество веществ

    pivot_columns = []  # список ведущих столбцов
    current_row = 0     # текущая строка

    for col in range(cols_count):
        pivot = None

        # ищем строку с ненулевым элементом
        for r in range(current_row, rows_count):
            if matrix[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        # меняем строки местами
        matrix[current_row], matrix[pivot] = matrix[pivot], matrix[current_row]

        # делаем ведущий элемент равным 1
        pivot_value = matrix[current_row][col]
        matrix[current_row] = [x / pivot_value for x in matrix[current_row]]

        # обнуляем остальные строки
        for r in range(rows_count):
            if r != current_row and matrix[r][col] != 0:
                factor = matrix[r][col]
                matrix[r] = [
                    matrix[r][i] - factor * matrix[current_row][i]
                    for i in range(cols_count)
                ]

        pivot_columns.append(col)
        current_row += 1

    # ищем свободные переменные
    free_columns = [c for c in range(cols_count) if c not in pivot_columns]

    if not free_columns:
        return None

    solution = [Fraction(0)] * cols_count  # итоговый вектор коэффициентов
    solution[free_columns[0]] = Fraction(1)

    # выражаем остальные через свободную переменную
    for i, pc in enumerate(pivot_columns):
        solution[pc] = -matrix[i][free_columns[0]]

    return solution


# -------- ПЕРЕВОД В ЦЕЛЫЕ --------
def to_ints(fracs):
    common_denom = 1  # общий знаменатель

    for f in fracs:
        common_denom = common_denom * f.denominator // gcd(common_denom, f.denominator)

    nums = [int(f * common_denom) for f in fracs]

    g = reduce(gcd, [abs(x) for x in nums if x != 0])
    return [x // g for x in nums]


# -------- ОСНОВНАЯ ЛОГИКА --------
def balance(eq):
    left, right = split_equation(eq)
    all_parts = left + right

    atoms_list = []   # список словарей атомов для каждого вещества
    elements = []     # список всех элементов

    for formula in all_parts:
        parsed = parse_formula(formula)
        atoms_list.append(parsed)

        for el in parsed:
            if el not in elements:
                elements.append(el)

    matrix = []

    for el in elements:
        row = []

        for i, atoms in enumerate(atoms_list):
            value = Fraction(atoms.get(el, 0))  # количество атомов элемента

            if i >= len(left):
                value = -value  # продукты со знаком минус

            row.append(value)

        matrix.append(row)

    solution = solve_matrix(matrix)

    if solution is None:
        raise ValueError("Не получилось решить")

    coeffs = to_ints(solution)

    if any(c <= 0 for c in coeffs):
        raise ValueError("Некорректный результат")

    # сбор строки результата
    def format_part(c, f):
        return f if c == 1 else str(c) + f

    left_str = " + ".join(format_part(coeffs[i], left[i]) for i in range(len(left)))
    right_str = " + ".join(format_part(coeffs[i + len(left)], right[i]) for i in range(len(right)))

    return left_str + " -> " + right_str


# -------- КНОПКА --------
def solve():
    eq = equation_input.get()  # текст из поля ввода

    try:
        result = balance(eq)
        result_label.config(text=result)
    except Exception as e:
        result_label.config(text="Ошибка: " + str(e))


# -------- ИНТЕРФЕЙС --------
title = Label(root, text='Решение химических уравнений', bg=colors[0], fg=colors[1])
title.pack()

info = Label(root, text='Формат: A + B -> C + D', bg=colors[0], fg=colors[1])
info.place(x=5, y=25)

equation_input = Entry(root, width=50)  # поле ввода уравнения
equation_input.place(x=5, y=50)

btn = Button(root, text='Уравнять', command=solve)
btn.place(x=350, y=46)

result_label = Label(root, text='', bg=colors[0], fg=colors[2])  # вывод результата
result_label.place(x=5, y=70)

examples = Text(root, height=6, width=40)
examples.place(x=5, y=100)
examples.insert('1.0',
                'Примеры:\n'
                'N2 + H2 -> NH3\n'
                'Fe + O2 -> Fe2O3\n'
                'C3H8 + O2 -> CO2 + H2O')
examples.config(state="disabled")

root.mainloop()