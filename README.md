# TrajectoryHomework
智能移动机器人作业2
## 运行方式

clone本仓库代码到ros工作空间中

编译

```bash
catkin_make
source devel/setup.sh
```

运行

```bash
roslaunch astar_path_planner astar_planner.launch
```

RVIZ左侧可切换path和trajectory，trajectory会等待5s后开始生成
