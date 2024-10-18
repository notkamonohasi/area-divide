from pathlib import Path

from Vis.const import ROOT_DIR
from Vis.vis import draw

if __name__ == "__main__":
    sim_dir = ROOT_DIR.parent.joinpath("Sim")
    in_path = sim_dir.joinpath("in.txt")
    out_path = sim_dir.joinpath("out.txt")

    obstacles: list[tuple[float, float]] = []
    coordinates: list[tuple[tuple[float, float], tuple[float, float]]] = []
    path: list[tuple[float, float]] = []

    with open(in_path, "r") as f:
        _, n, _ = map(float, f.readline().split())
        for _ in range(int(n)):
            y, x = map(float, f.readline().split())
            obstacles.append((y, x))

    with open(out_path, "r") as f:
        n = int(f.readline())
        for _ in range(n):
            y, x, dy, dx = map(float, f.readline().split())
            coordinates.append(((y, x), (dy, dx)))
        m = int(f.readline())
        for _ in range(m):
            y, x = map(float, f.readline().split())
            path.append((y, x))

    draw(obstacles, coordinates, path, Path("./hoge.png"))
