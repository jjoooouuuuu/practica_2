#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 512
float Mat[N][N];
float MatDD[N][N];
float V1[N];
float V2[N];
float V3[N];

// Funció d'inicialització
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

// 1. Funció per mostrar els elements d'un vector a partir d'una posició 'from'
void PrintVect( float vect[N], int from, int numel ){
	int i;
	for (i = from; i < from + numel; i++) {
		printf (" %f ",vect[i]);
	}
	printf("\n");
}

// 2. Funció per mostrar una fila d'una matriu a partir d'una posició 'from'
void PrintRow( float mat[N][N], int row, int from, int numel ){
	int i;
	for (i = from; i < from + numel; i++){
		printf (" %f ",mat[row][i]);
	}
	printf("\n");
}

// 3. Multiplicació d'un escalar per un vector
void MultEscalar( float vect[N], float vectres[N], float alfa ){
	int i;
	for (i = 0; i < N; i++){
		vectres[i] = vect[i]*alfa;
	}
}

// 4. Producte escalar entre dos vectors
float Scalar( float vect1[N], float vect2[N] ){
	int i;
	float res = 0;
	for (i = 0; i < N; i++){
		res = res + (vect1[i]*vect2[i]);
	}
	return res;
}

// 5. Càlcul de la magnitud d'un vector
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

// 6. Comprova si dos vectors són ortogonals
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

// 7. Càlcul de la projecció d'un vector 'u' sobre un vector 'v'
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

// 8. Càlcul de la infini-norma d'una matriu
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

// 9. Càlcul de la norma-ú d'una matriu
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

// 10. Càlcul de la norma de Frobenius d'una matriu
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

// 11. Comprova si una matriu és diagonal dominant
int DiagonalDom( float M[N][N] ){
	int i;
	int u;
	float diag;
	int res;
	int dom = 1;
	for (i = 0; i < N; i++){
		float sum = 0;
		for (u = 0; u < N; u++){
			sum += fabs(M[i][u]);
		}
		diag = fabs(M[i][i]);
		sum = sum - diag;
		if (diag < sum){
			dom = 0;
		}
	}
	if (dom == 1){
		res = 1;
	}else{
		res = 0;
	}
	return res;
}

// 12. Multiplicació d'una matriu per un vector
void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ){
	int i;
	int u;
	for (i = 0; i < N; i++){
		vectres[i] = 0;
		for (u = 0; u < N; u++){
			vectres[i] += M[i][u] * vect[u];
		}
	}
}

//Càlcul de la precissió del Mètode de Jacobí
//Funció càlcul de la norma del residu
float calcular_norma(float residu[N]) {
    float suma = 0.0;
    for (int i = 0; i < N; i++) {
        suma += residu[i] * residu[i];
    }
    return sqrt(suma);
}
//Funció càlcul del residu
void calcular_residu(float M[N][N], float vect[N], float x[N], float residu[N]) {
    for (int i = 0; i < N; i++) {
        residu[i] = vect[i];
        for (int j = 0; j < N; j++) {
            residu[i] -= M[i][j] * x[j];
        }
    }
}
// 13. Mètode de Jacobi per resoldre sistemes d'equacions
int Jacobi(float M[N][N], float vect[N], float vectres[N], unsigned iter) {
	int dom = DiagonalDom(M);
	if (dom == 0){ 
		printf("No es pot aplicar el mètode Jacobi amb aquesta matriu, ja que no és diagonal dominant\n");
	}else{
		float x_prev[N];
		int i;
		int j;
		int k;
			for (int i = 0; i < N; i++) vectres[i] = 0.0;
			for (unsigned k = 0; k < iter; k++) {
				for (int i = 0; i < N; i++) x_prev[i] = vectres[i];
				for (int i = 0; i < N; i++) {
					float sum = vect[i];
					for (int j = 0; j < N; j++) {
						if (j != i) sum -= M[i][j] * x_prev[j];
					}
					vectres[i] = sum / M[i][i];
				}
			}
		PrintVect(vectres,0,10);
	}
}


int main(){
	printf("PRÀCTICA 2: Operacions amb matrius i vectors en C\n");

	InitData();

	printf("\n");
	printf("A.- Mostres vectors:");
	printf("\n");

	printf("\n");
	printf("Vector 1:");
	printf("\n");
	PrintVect(V1,0,10);
	PrintVect(V1,256,10);

	printf("\n");
	printf("Vector 2:");
	printf("\n");
	PrintVect(V2,0,10);
        PrintVect(V2,256,10);

	printf("\n");
	printf("Vector 3:");
	printf("\n");
	PrintVect(V3,0,10);
        PrintVect(V3,256,10);

	printf("\n");
	printf("Mostres matrius:");
	printf("\n");

	printf("\n");
	printf("B.- Matriu Mat:");
	printf("\n");
	PrintRow(Mat,0,0,10);
	PrintRow(Mat,100,0,10);

	printf("\n");
	printf("C.- Matriu MatDD:");
	printf("\n");
	PrintRow(MatDD,0,0,10);
	PrintRow(MatDD,100,95,10);

	printf("\n");
	printf("D.- Propietats de les matrius (1 = SI, 0 = NO):\n");

	printf("\n");
	printf ("Mat:\n");
	float infi1 = Infininorm(Mat);
        printf ("Infininorma: %f \n",infi1);
        float one1 = Onenorm(Mat);
        printf ("Norma-ú: %f \n",one1);
        float frobenius1 = NormFrobenius(Mat);
        printf ("Norma de Frobenius: %f \n",frobenius1);
        int diagonal_dom1 = DiagonalDom(Mat);
        printf ("Diagonal dominant: %d \n",diagonal_dom1);

	printf("\n");
	printf("MatDD:\n");
	float infi2 = Infininorm(MatDD);
        printf ("Infininorma: %f \n",infi2);
        float one2 = Onenorm(MatDD);
        printf ("Norma-ú: %f \n",one2);
        float frobenius2 = NormFrobenius(MatDD);
        printf ("Norma de Frobenius: %f \n",frobenius2);
        int diagonal_dom2 = DiagonalDom(MatDD);
        printf ("Diagonal dominant: %d \n",diagonal_dom2);

	printf("\n");
	printf("E.- Producte escalar de dos vectors:\n");

	printf("\n");
	printf("V1 · V2: ");
	float multscalar1 = Scalar(V1,V2);
        printf(" %f \n",multscalar1);


	printf("V1 · V3: ");
	float multscalar2 = Scalar(V1,V3);
        printf(" %f \n",multscalar2);


	printf("V2 · V3: ");
	float multscalar3 = Scalar(V2,V3);
	printf(" %f \n",multscalar3);

	printf("\n");
	printf("F.- Magnitud d'un vector:\n");

	printf("\n");
	printf("V1:");
	float mag1 = Magnitude(V1);
	printf(" %f \n",mag1);

	printf("V2:");
        float mag2 = Magnitude(V2);
        printf(" %f \n",mag2);

	printf("V3:");
        float mag3 = Magnitude(V3);
        printf(" %f \n",mag3);

	printf("\n");
	printf("G.- Ortogonals (1 = SI, 0 = NO):\n");

	printf("\n");
	int ort1 = Ortogonal(V1,V2);
	printf("V1 i V2: %d\n",ort1);

	int ort2 = Ortogonal(V1,V3);
        printf("V1 i V3: %d\n",ort2);

	int ort3 = Ortogonal(V2,V3);
        printf("V2 i V3: %d\n",ort3);

	printf("\n");
	printf("H.- Multiplicació d'un vector per un escalar:\n");
	float V7[N];
	MultEscalar(V3,V7,2);
	PrintVect(V7,0,10);
	PrintVect(V7,256,10);

	printf("\n");
	printf("I.- Projecció:\n");

	printf("\n");
	printf("V2 sobre V3:\n");
	float V4[N];
	Projection(V2,V3,V4);
	PrintVect(V4,0,10);

	printf("\n");
	printf("V1 sobre V2:\n");
        float V8[N];
        Projection(V1,V2,V8);
        PrintVect(V8,0,10);

	printf("\n");
	printf("J.- Multiplicació d'una matriu amb un vector:\n");

	printf("\n");
	printf("Mat x V2:\n");
	float V5[N];
	Matriu_x_Vector(Mat,V2,V5);
	PrintVect(V5,0,10);

	printf("\n");
	printf("K.- Resolució sistemes d'equacions:\n");

	printf("\n");
	float V6[N];
	printf ("Resolució sistema: MatDD·X = V3 (1 iteració)\n");
	Jacobi(MatDD,V3,V6,1);

	printf("\n");
	float V9[N];
        printf ("Resolució sistema: MatDD·X = V3 (1000 iteracions)\n");
	Jacobi(MatDD,V3,V9,1000);

	printf("\n");
	float V10[N];
        printf ("Resolució sistema: Mat·X = V3\n");
	Jacobi(Mat,V3,V10,1);

	printf("\n");
	printf("Precissió del Mètode de Jacobi per la resolució del primer sistema:\n");
	float residu1[N];
	calcular_residu(MatDD,V3,V6,residu1);
	float norm_res1 = calcular_norma(residu1);
	printf(" %f \n",norm_res1);

	printf("\n");
        printf("Precissió del Mètode de Jacobi per la resolució del segon sistema:\n");
        float residu2[N];
        calcular_residu(MatDD,V3,V9,residu2);
        float norm_res2 = calcular_norma(residu2);
        printf(" %f \n",norm_res2);
}
