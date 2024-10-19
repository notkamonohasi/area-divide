from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt


def draw(
    obstacles: List[tuple[float, float]],
    coordinates: List[Tuple[Tuple[float, float], Tuple[float, float]]],
    path: Optional[List[tuple[float, float]]],
    save_path: Path,
) -> None:
    """
    左上と右下の座標のリストを受け取り、複数の正方形を描画する。

    Parameters:
        coordinates (List[Tuple[Tuple[float, float], Tuple[float, float]]]):
            各正方形の左上と右下の座標のリスト。
    """
    fig, ax = plt.subplots()

    ys: list[float] = []
    xs: list[float] = []
    for ob in obstacles:
        ys.append(ob[0])
        xs.append(ob[1])
    ax.scatter(ys, xs, s=0.1)

    for top_left, bottom_right in coordinates:
        # 幅と高さを計算
        width = bottom_right[0] - top_left[0]
        height = top_left[1] - bottom_right[1]

        # 正方形を描画
        square = patches.Rectangle(
            top_left, width, -height, linewidth=0.2, edgecolor="gray", facecolor="none"
        )
        ax.add_patch(square)

    if path is not None:
        path_y: list[float] = []
        path_x: list[float] = []
        for p in path:
            path_y.append(p[0])
            path_x.append(p[1])
        plt.plot(path_y, path_x, linewidth=3.0, color="red")
        # plt.scatter(path_y, path_x)

    # 軸の範囲を設定
    all_x = [
        x
        for top_left, bottom_right in coordinates
        for x in [top_left[0], bottom_right[0]]
    ]
    all_y = [
        y
        for top_left, bottom_right in coordinates
        for y in [top_left[1], bottom_right[1]]
    ]
    ax.set_xlim(min(all_x), max(all_x))
    ax.set_ylim(min(all_y), max(all_y))
    ax.set_aspect("equal")  # アスペクト比を1:1に固定

    plt.savefig(save_path)
    plt.clf()
