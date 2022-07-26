
Note about  problem-16-22106-pre.txt:
dataset    https://grail.cs.washington.edu/projects/bal/

Header: 16 22106 83718bal problem file loaded...
bal problem have 16 cameras and 22106 points.
Forming 83718 observations.



structure : problem-16-22106-pre.txt  ???

#CameraObs #IdMapPoint #2dNormalizedMapPointor3dto2dImageProjection?
example
12 0     -7.204300e+02 3.143400e+02    :  el map puntoO es observado por la posicion de camara 12, el punto tiene como coordenadas normalizadas   -7.204300e+02, 3.143400e+02  ??


0 to 15 pose cameras
0 to 22105 map points

Los 22016 map points fueron observados 83719 veces por todas las camaras 



WHAT IS THE FORMAT of problem-16-22106-pre.txt   ??
A partir de la linea 83720 hasta la 150181 que significa???  (unkonw values   150181-83719= 66462 ??)


Data Format
Each problem is provided as a bzip2 compressed text file in the following format.

TO CHECK

<num_cameras> <num_points> <num_observations>
<camera_index_1> <point_index_1> <x_1> <y_1>
...
<camera_index_num_observations> <point_index_num_observations> <x_num_observations> <y_num_observations>
<camera_1>
...
<camera_num_cameras>
<point_1>
...
<point_num_points>
