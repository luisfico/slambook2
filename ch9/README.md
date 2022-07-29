
Note about  problem-16-22106-pre.txt:
dataset    https://grail.cs.washington.edu/projects/bal/

Header: 16 22106 83718 bal problem file loaded...
bal problem have 16 cameras and 22106 points.
Forming 83718 observations.



structure : problem-16-22106-pre.txt  ???

#CameraObs #IdMapPoint #2dNormalizedMapPointor3dto2dImageProjection?
example
12 0     -7.204300e+02 3.143400e+02    :  el map puntoO es observado por la posicion de camara 12, el punto tiene como coordenadas normalizadas  o  es en pixels pero centradas -w/2??
 -7.204300e+02, 3.143400e+02  ??


0 to 15 pose cameras
0 to 22105 map points

Los 22016 map points fueron observados 83718 veces por todas las camaras 



WHAT IS THE FORMAT of problem-16-22106-pre.txt   ??
A partir de la linea 83720 hasta la 150181 que significa???  (unkonw values   150181-83719= 66462 ??)

A partir de esa linea se termino de observar los map points
<num_cameras> <num_points> <num_observations>
<camera_index_1> <point_index_1> <x_1> <y_1>
...

Luego viene





********************

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

A partir de la linea   150181  
se dan las posiciones de las 16 camaras
el pos inversa    camara esta bien pero las otras ko

camara1

r1  -1.6943983532198115e-02
r2  1.1171804676513932e-02
r3  2.4643508831711991e-03
tx  7.3030995682610689e-01
ty  -2.6490818471043420e-01
tz  -1.7127892627337182e+00
f   1.4300319432711681e+03
k1  -7.5572758535864072e-08
k2  3.2377569465570913e-14

osea hay 16 cameras con 9 params iniciales
termina hasta 83863          83863-83719= 144   = 16x9

que es el resto   la profundidad de los puntos map SI y estan en orden

pos 3d de los puntos iniciales
px
py
pz