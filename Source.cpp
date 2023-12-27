#include <cstdio>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <iostream>
//#include <gtest/gtest.h>

//#pragma once;








//очистить выделенную память
void clear(double **arr, int n)
{
    for(int i = 0; i < n; i++)
        delete[] arr[i];
    delete[] arr;
}
//создать копию массива
double** clone(double **arr, int n)
{
    double **newArr = new double*[n];
    for(int row = 0; row < n; row++)
    {
        newArr[row] = new double[n];
        for(int col = 0; col < n; col++)
            newArr[row][col] = arr[row][col];
    }
    return newArr;
}
//print the matrix
void show(double **matrix, int n)
{
    for(int row = 0; row < n; row++){
            for(int col = 0; col < n; col++)
                printf("%lf\t", matrix[row][col]);
            printf("\n");
    }
    printf("\n");
}
//матричное умножение матриц
double** matrix_multi(double **A, double **B, int n)
{
	
    double **result = new double*[n];
    //заполнение нулями
    for(int row = 0; row < n; row++){
        result[row] = new double[n];
        for(int col = 0; col < n; col++){
            result[row][col] = 0;
        }
    }
	//std::cout<<"ABCWEWE"<<std::endl;
    for(int row = 0; row < n; row++){
        for(int col = 0; col < n; col++){
            for(int j = 0; j < n; j++){
					//std::cout<<"ABCWEWE"<<std::endl;

               result[row][col] +=  A[row][j] * B[j][col];
            }
        }
    }
    return result;
}
//умножение матрицы на число
void scalar_multi(double **m, int n, double a){
    for(int row = 0; row < n; row++)
        for(int col = 0; col < n; col++){
            m[row][col] *= a;
        }
}
//вычисление суммы двух квадратных матриц
void sum(double **A, double **B, int n)
{
    for(int row = 0; row < n; row++)
        for(int col = 0; col < n; col++)
            A[row][col] += B[row][col];
}
 
//вычисление определителя
double det(double **matrix, int n) //квадратная матрица размера n*n
{
    double **B = clone(matrix, n);
    //приведение матрицы к верхнетреугольному виду
    for(int step = 0; step < n - 1; step++)
        for(int row = step + 1; row < n; row++)
        {
            double coeff = -B[row][step] / B[step][step]; //метод Гаусса
            for(int col = step; col < n; col++)
                B[row][col] += B[step][col] * coeff;
        }
    //Рассчитать определитель как произведение элементов главной диагонали
    double Det = 1;
    for(int i = 0; i < n; i++)
        Det *= B[i][i];
    //Очистить память
    clear(B, n);
    return Det;
}

 /*TEST(TestInverse, Subtest)
 {
	ASSERT_TRUE(1==1); 
 }*/
 
 
 double** Inversion(double **arr, int n)
{
   double N1 = 0, Ninf = 0; //норма матрицы по столбцам и по строкам
    double **A0 = clone(arr, n);       //инициализация начального приближения
    for(int row = 0; row < n; row++){
        double colsum = 0, rowsum = 0;
        for(int col = 0; col < n; col++){
            rowsum += fabs(A0[row][col]);
            colsum += fabs(A0[col][row]);
        }
        N1 = std::max(colsum, N1);
        Ninf = std::max(rowsum, Ninf);
    }
    //транспонирование
    for(int row = 0; row < n - 1; row++){
        for(int col = row + 1; col < n; col++)
            std::swap(A0[col][row], A0[row][col]);
    }
    scalar_multi(A0, n, (1 / (N1 * Ninf))); //нормирование матрицы
    //инициализация удвоенной единичной матрицы нужного размера
    double **E2 = new double*[n];
    for(int row = 0; row < n; row++)
    {
        E2[row] = new double[n];
            for(int col = 0; col < n; col++){
                if(row == col)
                    E2[row][col] = 2;
                else
                    E2[row][col] = 0;
            }
    }
    double **inv = clone(A0, n); //A_{0}
    double EPS = 0.0000001;   //погрешность
    if(det(arr, n) != 0){ //если матрица не вырождена
        while(fabs(det(matrix_multi(arr, inv, n), n) - 1) >= EPS) //пока |det(A * A[k](^-1)) - 1| >= EPS
        {
            double **prev = clone(inv, n); //A[k-1]
            inv = matrix_multi(arr, prev, n);   //A.(A[k-1]^(-1))
            scalar_multi(inv, n, -1);         //-A.(A[k-1]^(-1))
            sum(inv, E2, n);                   //2E - A.(A[k-1]^(-1))
            inv = matrix_multi(prev, inv, n); //(A[k-1]^(-1)).(2E - A.(A[k-1]^(-1)))
            clear(prev, n);
        }
        
        //show(res, n);
        //вывод матрицы на экран
         
    

        
    }
   return inv;
  
    
    
}
 
 
/*TEST (matrixTests, inversion)
{
	
	int n = 3;
	
	double **testMatrix = new double*[n];
    for(int row = 0; row < n; row++)
    {
        testMatrix[row] = new double[n];

    }
	            
                testMatrix[0][0] = 2;
                testMatrix[0][1] = 5;
                testMatrix[0][2] = 7;
				testMatrix[1][0] = 6;
				testMatrix[1][1] = 3;
				testMatrix[1][2] = 4;
				testMatrix[2][0] = 5;
				testMatrix[2][1] = -2;
				testMatrix[2][2] = -3;
				
		
    double **invTest = Inversion(testMatrix, n);
    double **resTest = matrix_multi(testMatrix,invTest, n);

    ASSERT_LT(abs(resTest[0][0]-1),1e-6);
	ASSERT_LT(abs(resTest[0][1]),1e-6);
	
}
 
 
 TEST (matrixTests, inversion2)
{
	
	int n = 3;
	
	double **testMatrix = new double*[n];
    for(int row = 0; row < n; row++)
    {
        testMatrix[row] = new double[n];

    }
	            
               testMatrix[0][0] = 2;
                testMatrix[0][1] = 5;
                testMatrix[0][2] = 7;
				testMatrix[1][0] = 6;
				testMatrix[1][1] = 3;
				testMatrix[1][2] = 4;
				testMatrix[2][0] = 5;
				testMatrix[2][1] = -2;
				testMatrix[2][2] = -3;
				
		
    double **inversed = new double*[n];
    
    
    for(int row = 0; row < n; row++)
    {
        inversed[row] = new double[n];

    }
    
    inversed[0][0] = 1;
                inversed[0][1] = -1;
                inversed[0][2] = 1;
				inversed[1][0] = -38;
				inversed[1][1] = 41;
				inversed[1][2] = -34;
				inversed[2][0] = 27;
				inversed[2][1] = -29;
				inversed[2][2] = 24;
    
    double **invTest = new double*[n];
    
    
    for(int row = 0; row < n; row++)
    {
        invTest[row] = new double[n];

    }
    
    invTest = Inversion(testMatrix, n);


for(int i = 0; i < 3; i++)
{
	for(int j = 0; j < 3; j++)
	{
		ASSERT_LT(abs(inversed[i][j]-invTest[i][j]),1e-6);
		
	}
}
    
	

	
}
 TEST (matrixTests, inversion3)
{
	
	int n = 3;
	
	double **testMatrix = new double*[n];
    for(int row = 0; row < n; row++)
    {
        testMatrix[row] = new double[n];

    }
	            
                testMatrix[0][0] = 1;
                testMatrix[0][1] = 5;
                testMatrix[0][2] = 9;
				testMatrix[1][0] = 11;
				testMatrix[1][1] = 31;
				testMatrix[1][2] = 43;
				testMatrix[2][0] = 5;
				testMatrix[2][1] = -22;
				testMatrix[2][2] = 0;
				
		
    double **invTest = Inversion(testMatrix, n);
    double **resTest = matrix_multi(testMatrix,invTest, n);

    ASSERT_LT(abs(resTest[0][0]-1),1e-6);
	ASSERT_LT(abs(resTest[0][1]),1e-6);
	ASSERT_LT(abs(resTest[0][2]),1e-6);
	ASSERT_LT(abs(resTest[1][0]),1e-6);
	ASSERT_LT(abs(resTest[1][1]-1),1e-6);
	ASSERT_LT(abs(resTest[1][2]),1e-6);
	ASSERT_LT(abs(resTest[2][1]),1e-6);
	ASSERT_LT(abs(resTest[2][2]-1),1e-6);
	ASSERT_LT(abs(resTest[2][3]),1e-6);

	
}*/
int main()
{std::cout<<"MURad";
	srand(time(0));
    //Исходная матрица, динамический двухмерный массив
    int n=50;// scanf("%d", &n);
    double h = 1.0/n;
    
    double *x = new double(n+1);
    for(int i = 0; i < n; i++)
    {
        x[i] = i*h;
            
                
                //scanf("%lf", &A[row][col]);
    }
    
    double *A = new double(n*n);
    delete A;
    
    //double **B = Inversion(A,n);
        
       
    // show(B,  n);
      //clear(A,n);
     //   clear(B,n);

	//std::cout<<std::endl;

 //matrix_multi(A,A,n);
	//std::cout<<std::endl;
 
    /* Численное вычисление обратной матрицы по методу Ньютона-Шульца
        1. Записать начальное приближение [Pan, Schreiber]:
            1) Транспонировать данную матрицу
            2) Нормировать по столбцам и строкам
        2. Повторять процесс до достижения заданной точности.
    */
 
    
   

    std::cout<<std::endl;
    std::cout<<"TEST\n\n\n\n";
     //double **test ;
  // clone(inv, n);
         // double **test1 = clone(A, n);

   //  matrix_multi(inv,A,4);
   //::testing::InitGoogleTest();
   
   return 0;//RUN_ALL_TESTS();
   
    //return 0;
}
