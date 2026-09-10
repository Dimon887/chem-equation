# Copyright (c) 30.01.2026-03.07.2026 Бесшабашнов Дмитрий
from tkinter import *
import re
from fractions import Fraction

colors = ['olive', 'lightblue', 'yellow']

root = Tk()
root.title('Решение уравнений реакций')
root.geometry('430x210')
root['bg'] = colors[0]

# Алгоритм

# ─────────────────────────────────────────────
#  Разбор формулы
# ─────────────────────────────────────────────

def parse_formula(formula: str) -> dict[str, int]:
    """Рекурсивно разбирает химическую формулу, возвращает {элемент: количество}."""
    formula = formula.strip()

    def parse(s: str, pos: int) -> tuple[dict, int]:
        counts: dict[str, int] = {}

        while pos < len(s):
            c = s[pos]

            if c == '(':
                sub, pos = parse(s, pos + 1)
                # читаем множитель после ')'
                num_start = pos
                while pos < len(s) and s[pos].isdigit():
                    pos += 1
                mult = int(s[num_start:pos]) if pos > num_start else 1
                for el, cnt in sub.items():
                    counts[el] = counts.get(el, 0) + cnt * mult

            elif c == ')':
                return counts, pos + 1

            elif c.isupper():
                # читаем название элемента
                el = c
                pos += 1
                while pos < len(s) and s[pos].islower():
                    el += s[pos]
                    pos += 1
                # читаем коэффициент
                num_start = pos
                while pos < len(s) and s[pos].isdigit():
                    pos += 1
                num = int(s[num_start:pos]) if pos > num_start else 1
                counts[el] = counts.get(el, 0) + num

            else:
                #raise ValueError(f"Неожиданный символ '{c}' в формуле «{formula}»")
                p = ValueError(f"Неожиданный символ '{c}' в формуле «{formula}»")
                resultlabel.configure(text=p)

        return counts, pos

    result, _ = parse(formula, 0)
    return result


def split_equation(equation: str) -> tuple[list[str], list[str]]:
    """Разделяет уравнение на левую и правую части, возвращает списки формул."""
    arrow_match = re.search(r'->', equation)
    if not arrow_match:
        #raise ValueError("Уравнение должно содержать стрелку '->'")
        p = ValueError("Уравнение должно содержать стрелку '->'")
        resultlabel.configure(text=p)
    left_str, right_str = equation[:arrow_match.start()], equation[arrow_match.end():]
    left  = [f.strip() for f in left_str.split('+')]
    right = [f.strip() for f in right_str.split('+')]
    # убираем числовые коэффициенты, если пользователь их ввёл
    def strip_coeff(f: str) -> str:
        return re.sub(r'^\d+', '', f).strip()
    return [strip_coeff(f) for f in left], [strip_coeff(f) for f in right]


# ─────────────────────────────────────────────
#  Решение системы линейных уравнений
#  методом Гаусса над полем дробей
# ─────────────────────────────────────────────

def gauss_null_space(matrix: list[list[Fraction]]) -> list[Fraction]:
    """
    Находит ненулевое решение однородной системы Ax = 0
    методом Гаусса–Жордана над Q (дробями).
    Возвращает вектор x или None, если решения нет.
    """
    rows = [row[:] for row in matrix]   # копия
    m = len(rows)
    n = len(rows[0])

    pivot_cols: list[int] = []
    row_idx = 0

    for col in range(n):
        # ищем ведущий элемент
        pivot = None
        for r in range(row_idx, m):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue

        rows[row_idx], rows[pivot] = rows[pivot], rows[row_idx]
        # нормируем строку
        pv = rows[row_idx][col]
        rows[row_idx] = [x / pv for x in rows[row_idx]]
        # обнуляем остальные строки в этом столбце
        for r in range(m):
            if r != row_idx and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [rows[r][j] - factor * rows[row_idx][j] for j in range(n)]

        pivot_cols.append(col)
        row_idx += 1

    # свободные переменные
    free_cols = [c for c in range(n) if c not in pivot_cols]
    if not free_cols:
        return None  # только тривиальное решение

    # полагаем первую свободную переменную = 1
    free_col = free_cols[0]
    solution = [Fraction(0)] * n
    solution[free_col] = Fraction(1)

    for i, pc in enumerate(pivot_cols):
        if i < len(rows):
            solution[pc] = -rows[i][free_col]

    return solution


def lcm(a: int, b: int) -> int:
    from math import gcd
    return a * b // gcd(a, b)


def to_integers(fractions: list[Fraction]) -> list[int]:
    """Переводит список дробей в минимальные целые числа."""
    # общий знаменатель
    denom = 1
    for f in fractions:
        denom = lcm(denom, f.denominator)
    ints = [int(f * denom) for f in fractions]
    # делим на НОД
    from math import gcd
    from functools import reduce
    g = reduce(gcd, [abs(x) for x in ints if x != 0])
    return [x // g for x in ints]


# ─────────────────────────────────────────────
#  Основная логика
# ─────────────────────────────────────────────

def balance(equation: str) -> str:
    left_formulas, right_formulas = split_equation(equation)
    all_formulas = left_formulas + right_formulas
    n = len(all_formulas)   # количество коэффициентов

    # множество всех элементов
    elements: list[str] = []
    formula_atoms: list[dict[str, int]] = []
    for f in all_formulas:
        atoms = parse_formula(f)
        formula_atoms.append(atoms)
        for el in atoms:
            if el not in elements:
                elements.append(el)

    m = len(elements)  # количество уравнений (по числу элементов)

    # строим матрицу: левые вещества со знаком +, правые со знаком −
    left_count = len(left_formulas)
    matrix: list[list[Fraction]] = []
    for el in elements:
        row = []
        for i, atoms in enumerate(formula_atoms):
            coeff = Fraction(atoms.get(el, 0))
            if i >= left_count:
                coeff = -coeff   # правая часть → меняем знак
            row.append(coeff)
        matrix.append(row)

    solution = gauss_null_space(matrix)
    if solution is None:
        #raise ValueError("Не удалось найти решение — уравнение, возможно, не балансируется.")
        p = ValueError("Не удалось найти решение — уравнение, возможно, не балансируется.")
        resultlabel.configure(text=p)

    coeffs = to_integers(solution)

    if any(c <= 0 for c in coeffs):
        '''raise ValueError("Получены нулевые или отрицательные коэффициенты — "
                         "проверьте правильность уравнения.")'''
        p = ValueError("Получены нулевые или отрицательные коэффициенты — проверьте правильность уравнения.")
        resultlabel.configure(text=p)

    # собираем результат
    def fmt(coeff: int, formula: str) -> str:
        return formula if coeff == 1 else f"{coeff}{formula}"

    left_part  = " + ".join(fmt(coeffs[i], left_formulas[i])
                             for i in range(left_count))
    right_part = " + ".join(fmt(coeffs[left_count + i], right_formulas[i])
                             for i in range(len(right_formulas)))
    return f"{left_part} -> {right_part}"

def alg():
    eq = str(equation.get())
    res = balance(eq)
    resultlabel.configure(text=res)

# Interface

q = Label(root, text='Решение химических уравнений реакций', bg=colors[0], fg=colors[1])
q.pack()

q0 = Label(root, text='Формат: вещество + вещество -> продукт + продукт', bg=colors[0], fg=colors[1])
q0.place(x=5, y=25)

equation = Entry(root, width=50)
equation.place(x=5, y=50)

btn = Button(root, text='Уравнять', bg=colors[0], fg=colors[1], command=alg)
btn.place(x=350,y=46)

resultlabel = Label(root, text='', bg=colors[0], fg=colors[2])
resultlabel.place(x=5, y=70)

q1 = Text(root, height=6, width=40)
q1.place(x=5, y=100)
q1.insert('1.0', 'Примеры:\nN2 + H2 -> NH3\nFe + O2 -> Fe2O3\nC3H8 + O2 -> CO2 + H2O\nKMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2\nAl + H2SO4 -> Al2(SO4)3 + H2\nCa3(PO4)2 + H2SO4 -> CaSO4 + H3PO4')
q1.config(state="disabled")
root.mainloop()