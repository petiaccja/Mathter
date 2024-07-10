import sys
import copy
import pathlib
import math

def generate_accessors(elements, count):
    dimension = len(elements)
    total = math.prod([dimension] * count)

    accessors = []

    for code in range(total):
        name = ""
        indices = []
        my_code = copy.copy(code)
        for _ in range(count):
            index = my_code % dimension
            my_code //= dimension
            name = elements[index] + name
            indices = [index, *indices]
        accessors.append((name, indices))
    return accessors


def generate_accessors_for_dim(dimension):
    if dimension == 1:
        elements = ["x"]
    if dimension == 2:
        elements = ["x", "y"]
    if dimension == 3:
        elements = ["x", "y", "z"]
    if dimension == 4:
        elements = ["x", "y", "z", "w"]
    return [
        *generate_accessors(elements, 1),
        *generate_accessors(elements, 2),
        *generate_accessors(elements, 3),
        *generate_accessors(elements, 4),
    ]


template = """Swizzle<T, Dim, Packed, {indices}> {name};"""

dimension = int(sys.argv[1])
file = sys.argv[2]

accessors = generate_accessors_for_dim(dimension)
code = "\n".join([template.format(name=name, indices=", ".join([str(i) for i in indices])) for name, indices in accessors])

pathlib.Path(file).write_text(code, encoding="utf8")