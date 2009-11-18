/*
Lee Ourand
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

void load_matrix(int n, int m, double *matrix[], ifstream & input_file);
void print_matrix(int n, int m, double *matrix[]);
void swap_rows(int n, int m, double *matrix[], int row1, int row2);
void reduce_rows(int n, int m, double *matrix[], int column);
int find_pivot(int n, int m, double *matrix[], int column);
int gaussian_elimination(int n, int m, double *matrix[], double vector[]);
int lu_factorization(int n, double *matrix[]);
void lu_decomposition(int n, double *matrix[], double vector[]);
void load_vector(int n, double right_vector[], ifstream & input_file);
void print_vector(int n, double vector[]);
void forward_multiply(int n, double *matrix[], double vector[], int row);


int main()
{
	int type;
	int n; //square rows in matrix
	int m; //columns in matrix
	int e; // error code
	string filename;
	cout << "What type of matrix solver? (enter number)" << endl
	<< "1. Gaussian Elimination" << endl << "2. LU Decomposition" << endl
	<< "3. Jacobi" << endl << "4. Gauss-Seidel" << endl << endl;
	cin >> type;
	
	cout << endl;
	cout << "Data File Name: ";
	cin >> filename;
	
	ifstream input_file;
	input_file.open(filename.c_str() ,ios::in);
	
	input_file >> n;
	switch(type)
	{
		case 1:
			m = n * n +1;
			break;
		
		case 2 :
			m = n * n * 2;
			break;
	}
	
	double *matrix[m];
	double right_vector[n*n];
	switch(type)
	{
		case 1 : 
			load_vector(n*n, right_vector, input_file);
			load_matrix(n*n, m, matrix, input_file);
			print_matrix(n*n, m, matrix);
			e = gaussian_elimination(n*n, m, matrix, right_vector);
			if(!e)
			{
				cout << "Singular Matrix. Unsolvable" << endl;
			} else {
				print_matrix(n*n, m, matrix);
			}
			break;
		
		case 2 :
			load_vector(n*n, right_vector, input_file);
			load_matrix(n*n, m, matrix, input_file);
			print_matrix(n*n, m, matrix);
			e = lu_factorization(n*n, matrix);
			lu_decomposition(n*n, matrix, right_vector);
			print_matrix(n*n, m, matrix);					
			break;
			
		default :
			
			break;
	}		
	
	input_file.close();
	return 0;
}

/* Read in the matrix.
   Note: We're storing ROWS in the arrays. Each row (array) has a corresponding
   pointer in the "matrix" array above. This allows for easily swapping rows
   when performing row operations.
*/
void load_matrix(int n, int m, double *matrix[], ifstream & input_file)
{
	for(int i=0; i<n; i++)
	{
		matrix[i] = new double[m];
		// Load in the matrix
		for(int j=0; j<n; j++)
		{
			input_file >> matrix[i][j];
		}
	}
}

/* Print Matrix */
void print_matrix(int n, int m, double *matrix[])
{
	fstream output_file;
	output_file.open("results.csv", ios::out);
	for(int i=0; i<n; i++)
	{
		for(int j=0; j < m; j++)
		{
			output_file << matrix[i][j] << ", ";
		}
		output_file << endl;
	}
}

void load_vector(int n, double vector[], ifstream & input_file)
{
	for(int i=0; i<n; i++)
	{
		input_file >> vector[i];
	}	
}

void print_vector(int n, double vector[])
{
	for(int i=0; i<n; i++)
	{
		cout << vector[i] << endl;
	}
}

void forward_multiply(int n, double *matrix[], double vector[], int start_column)
{
	double x[n];
	//Set x to 0
	for(int i=0; i<n; i++)
	{
		x[i] = 0;
	}
	for(int i=0; i<n; i++)
	{
		for(int j=start_column; j<start_column + n; j++)
		{
			x[i] += matrix[i][j] * vector[j-start_column];
		}
		cout << endl;	
	}
	
	//swap x and vector
	for(int i=0; i<n; i++)
	{
		vector[i] = x[i];
	}
		
}

void swap_rows(int n, int m, double *matrix[], int row1, int row2)
{
	double vec[n];
	for(int p=0; p<m; p++)
	{
		vec[p] = matrix[row1][p];
		matrix[row1][p] = matrix[row2][p];
		matrix[row2][p] = vec[p];
	}
}

void reduce_rows(int n, int m, double *matrix[], int column)
{
	for(int i=column+1; i<n; i++)
	{
		double coeff = matrix[i][column]/matrix[column][column];
		for(int k=column; k<m; k++)
		{
			matrix[i][k] = matrix[i][k] - (coeff * matrix[column][k]);
			if(abs(matrix[i][k]) < 0.000001)
			{
				matrix[i][k] = 0;
			}
		}
	}
}

void back_solve(int n, int m, double *matrix[])
{
	double x[n];
	//initialize x to 0
	for(int i=0; i< n; i++)
	{
		x[i] = 0;
	}
	
	for(int j=n -1; j>=0; j--)
	{
		double sum = 0;
		for(int i=j+1; i<n; i++)
		{
			sum += matrix[j][i] * x[i];
		}
		x[j] = 1 / matrix[j][j] * (matrix[j][n] - sum);
		if(abs(x[j]) < 0.000001)
		{
			x[j] = 0;
		}
	}
	
	print_vector(n, x);	
}

int find_pivot(int n, int m, double *matrix[], int column)
{
	int position;
	double temp = 0.0;
	for(int i=column; i<n; i++)
	{
		if(abs(matrix[i][column]) > abs(temp))
		{
			temp = matrix[i][column];
			position = i;
		}
	}
	return position;
}

// BELOW ARE THE VARIOUS MATRIX SOLVER FUNCTIONS //

/* Gaussian Elimination */
int gaussian_elimination(int n, int m, double *matrix[], double vector[])
{
	int error = 1;
	int position = 0;
	// Add Right Hand Side As Augmented Matrix
	for(int i=0; i<n; i++)
	{
		matrix[i][m-1] = vector[i];
	}
	
	for(int j=0; j<n; j++)
	{
		position = find_pivot(n, n+1, matrix, j);
		// Check for 0 pivot (error)
		if(matrix[position][j] == 0.0)
		{
			error = 0;
		}
		swap_rows(n, n+1, matrix, position, j);
		reduce_rows(n, n+1, matrix, j);
	}
	back_solve(n, n+1, matrix);
	return error;
}

int lu_factorization(int n, double *matrix[])
{
	int m = n*2;
	int delta = 1;
	int position = 0;
	
	// Add Identity Matrix to Right Hand
	for(int i=0; i < n; i++)
	{
		for(int j=m-n; j< m; j++)
		{
			if(i==(j-n))
			{
				matrix[i][j] = 1;
			} else {
				matrix[i][j] = 0;	
			}
		}
	}
	
	// Factor A into LU
	for(int j=0; j<n;j++)
	{
		for(int k=j;k<n;k++)
		{
			double sum = 0;
			for(int i=0; i<j-1; i++)
			{
				sum += matrix[k][i] * matrix[i][j];
			}
			matrix[k][j] = matrix[k][j] - sum;
			if(matrix[k][j] < 0.0001)
			{
				matrix[k][j] = 0;
			}
		}
		
		position = find_pivot(n, n, matrix, j);
		if(matrix[position][j] == 0.0)
		{
			delta = 0;
		} else {
			if(position > j)
			{
				swap_rows(n, 2*n, matrix, position, j);
				delta = delta * -1;
			}
		}
		
		for(int k=j+1; k<n;k++)
		{
			double sum=0;
			for(int i=0; i<j-1;i++)
			{
				sum += matrix[j][i] * matrix[i][k];
			}
			matrix[j][k] = 1 / matrix[j][j] * (matrix[j][k] - sum);
			if(matrix[j][k] < 0.0001 )
			{
				matrix[j][k] = 0;
			}
		} 
		
		delta = matrix[j][j] * delta;		
		
	}
	return delta;
}

void lu_decomposition(int n, double *matrix[], double vector[])
{
	double y[n];
	double x[n];
	forward_multiply(n, matrix, vector, n-1);
	
	for(int i=0; i<n; i++)
		x[i] = 0;
	
	// Compute Y
	for(int k=0; k < n; k++)
	{
		double sum = 0;
		for(int i=0; i < k; i++)
		{
			sum += matrix[k][i] * y[i];
		}
		
		y[k] = 1/matrix[k][k] * (vector[k] - sum);
	}
	
	for(int k=n; k >=0; k--)
	{
		double sum = 0;
		for(int i=k+1; i<n; i++)
		{
			sum += matrix[k][i] * x[i];
		}
		x[k] = y[k] - sum;
	}
	
	print_vector(n, x);
}

