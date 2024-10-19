import math
import random

if __name__ == "__main__":
    random.seed(0)

    XY = 10.0
    n_obstacles = 5000
    FIRST_MAX_N_AREAS = 1000
    SECOND_MAX_N_AREAS = 2000

    print(XY, n_obstacles, FIRST_MAX_N_AREAS, SECOND_MAX_N_AREAS)
    print(1.0, 3.0, 1.2, 7.0)

    for _ in range(n_obstacles // 2):
        x = random.uniform(0.0001, 9.0)
        y = -2.0 * math.cos((x - 5.0) / 0.4) - math.sqrt(x) / 5.0 + 3.0
        print(x, y)

    for _ in range(n_obstacles // 2):
        x = random.uniform(0.0001, 9.0)
        # y = 1.5 * math.cos((x - 2.5) / 0.6) + 7.0 + math.sqrt(x) / 3.0
        y = -2.0 * math.cos((x - 5.0) / 0.4) - math.sqrt(x) / 5.0 + 6.0
        print(x, y)
