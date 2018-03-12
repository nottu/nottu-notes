Compilar con
  gcc main.c vector.c matrix.c optimization.c matrix_factor.c -o ej1
  gcc ej2.c  vector.c matrix.c optimization.c matrix_factor.c -o ej2

Ejecutar con
 ej1 [rosenbrock | wood ] [StepBacktrack | StepInterpol] [filename] [maxiter] [tol_gradient] [tol_x] [tol_function] [initial_alpha]
 ej2 [trainY] [trainX] [testY] [testX] [num1] [num2]

Ejemplo
  ./ej1 rosenbrock StepBacktrack ros2.txt 1000000 1E-10 1E-50 1E-50 0.5
  ./ej2 trainY.csv trainX.csv testY.csv testX.csv 0 8