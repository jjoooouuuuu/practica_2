#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 512
float Mat[N][N];
float MatDD[N][N];
float V1[N];
float V2[N];
float V3[N];

void InitData(){
int i,j;
srand(334411);
for( i = 0; i < N; i++ )
 for( j = 0; j < N; j++ ){
 Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
 if ( (abs(i - j) <= 3) && (i != j))
 MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
 else if ( i == j )
 MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
 else MatDD[i][j] = 0.0;
 }
for( i = 0; i < N; i++ ){
 V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
}
}


void PrintVect( float vect[N], int from, int numel ){
	int i;
	for (i = from; i < from + numel; i++) {
		printf (" %f ",vect[i]);
	}
	printf("\n");
}


void PrintRow( float mat[N][N], int row, int from, int numel ){
	int i;
	for (i = from; i < from + numel; i++){
		printf (" %f ",mat[row][i]);
	}
	printf("\n");
}


void MultEscalar( float vect[N], float vectres[N], float alfa ){
	int i;
	for (i = 0; i < N; i++){
		vectres[i] = vect[i]*alfa;
	}
}


float Scalar( float vect1[N], float vect2[N] ){
	int i;
	float res = 0;
	for (i = 0; i < N; i++){
		res = res + (vect1[i]*vect2[i]);
	}
	return res;
}


float Magnitude( float vect[N] ){
	int i;
	float arrel = 0;
	float res = 0;
	for (i = 0; i < N; i++){
		arrel += vect[i]*vect[i];
	}
	res = sqrt(arrel);
	return res;
}


int Ortogonal( float vect1[N], float vect2[N] ){
	float res_scalar = Scalar(vect1,vect2);
	float res = 0;
	if (res_scalar == 0){
		res = 1;
	}
	else {
		res = 0;
	}
	return res;
}


void Projection( float vect1[N], float vect2[N], float vectres[N] ){
	int i;
	float res_scalar;
	float mag;
	float cof;
	mag = Magnitude(vect2);
	res_scalar = Scalar(vect1,vect2);
	cof = res_scalar/mag;
	for (i = 0; i < N; i++){
		vectres[i] = vect2[i] * cof;
	}
}


float Infininorm( float M[N][N] ){
	int i;
	int u;
	float sum = 0;
	float res = 0;
	for (i = 0; i < N; i++){
		for (u = 0; u < N; u++){
			sum += fabs(M[i][u]);
		}
		if (sum > res){
			res = sum;
		}
		sum = 0;
	}
	return res;
}


float Onenorm( float M[N][N] ){
	int i;
        int u;
        float sum = 0;
        float res = 0;
        for (u = 0; u < N; u++){
                for (i = 0; i < N; i++){
                        sum += fabs(M[i][u]);
		}
                if (sum > res){
                        res = sum;
                }
                sum = 0;
        }
        return res;
}


float NormFrobenius( float M[N][N] ){
	int i;
	int u;
	float forb = 0;
	float res;
	for (i = 0; i < N; i++){
		for (u = 0; u < N; u++){
			forb += M[i][u] * M[i][u];
		}
	}
	res = sqrt(forb);
	return res;
}


int DiagonalDom( float M[N][N] ){
	int i;
	int u;
	float sum = 0;
	int res;
	int dom = 0;
	for (i = 0; i < N; i++){
		for (u = 0; u < N; u++){
			sum += fabs(M[i][u]);
		}
		sum -= fabs(M[i][i]);
		if (M[i][i] > sum){
			dom = 1;
		}else{
			dom = 0;
			break;
		}
	}
	if (dom > 0){
		res = 1;
	}else{
		res = 0;
	}
	return res;
}


int main(){
	InitData();
	printf("Mostres vectors:");
	printf("\n");
	printf("Vector 1:");
	printf("\n");
	PrintVect(V1,0,10);
	PrintVect(V1,256,10);
	printf("Vector 2:");
	printf("\n");
	PrintVect(V2,0,10);
        PrintVect(V2,256,10);
	printf("Vector 3:");
	printf("\n");
	PrintVect(V3,0,10);
        PrintVect(V3,256,10);
	printf("Mostres matrius:");
	printf("\n");
	printf("Matriu Mat:");
	printf("\n");
	PrintRow(Mat,0,0,10);
	PrintRow(Mat,100,0,10);
	printf("Matriu MatDD:");
	printf("\n");
	PrintRow(MatDD,0,0,10);
	PrintRow(MatDD,100,95,10);
	printf("Multiplicació escalar de dos vectors:");
	float multscalar = Scalar(V2,V3);
	printf(" %f \n",multscalar);
	printf("Magnitud d'un vector:");
	float mag = Magnitude(V2);
	printf(" %f \n",mag);
	printf("Projecció:");
	float V4[N];
	Projection(V2,V3,V4);
	PrintVect(V4,0,10);
	float infi = Infininorm(Mat);
	printf ("Infininorma: %f \n",infi);
	float one = Onenorm(Mat);
	printf ("Norma-ú: %f \n",one);
	float frobenius = NormFrobenius(Mat);
	printf ("Norma de Frobenius: %f \n",frobenius);

}





