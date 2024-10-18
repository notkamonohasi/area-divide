import math
import random

if __name__ == "__main__":
    random.seed(0)

    XY = 10.0
    n_obstacles = 5000
    N_MAX_AREAS = 1000

    print(XY, n_obstacles, N_MAX_AREAS)

    for _ in range(n_obstacles // 2):
        x = random.uniform(2.0, 8.0)
        y = -2.0 * math.cos((x - 5.0) / 0.4) + 3.0
        print(x, y)

    for _ in range(n_obstacles // 2):
        x = random.uniform(1.5, 6.5)
        y = 1.5 * math.cos((x - 2.5) / 0.6) + 7.0 + math.sqrt(x) / 3.0
        print(x, y)
