 #include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define K 100 // low rank approximation
#define ITER 25   // number of iterations basically how many tiem we are multiplying A^TA to the Q we are obtaining .

static float *repgm(const char *filename, int *column, int *row)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp) { fprintf(stderr,"ERROR : %s does not exist\n", filename); return NULL; }

    char fmt[3];
    if (fscanf(fp, "%2s", fmt) != 1 || fmt[0]!='P' || fmt[1]!='5') {
        printf("not valid PGM format\n");
        fclose(fp); return NULL;
    }

    int c;
    while ((c = fgetc(fp)) == '#') { 
        char b[256]; 
        fgets(b, sizeof(b), fp);
     }
    ungetc(c, fp);

    int w=0,h=0,maxv=0;
    if (fscanf(fp,"%d %d",&w,&h)!=2){fclose(fp);return NULL;}
    if (fscanf(fp,"%d",&maxv)!=1){fclose(fp);return NULL;}
    fgetc(fp); // eat newline

    int  size=w*h;
    unsigned char *buf = malloc(size);
    if(!buf){
        fclose(fp);
        return NULL;
    }
    fread(buf,1,size,fp);
    fclose(fp);

    float *img = malloc(size*sizeof(float));
    for(size_t i=0;i<size;i++){
         img[i] = buf[i]/255.0f;
    }
    free(buf);
    *column=w; *row=h;
    return img;
}
static int wrpgm(const char *filename, const float *img, int colm, int row)
{
    FILE *fp = fopen(filename,"wb");
    if(!fp){
        printf("Cannot write %s\n",filename);
        return 1;
    }
    fprintf(fp,"P5\n%d %d\n255\n",colm,row);
    int size=colm*row;
    for(int i=0;i<size;i++){
        float x=img[i];
        if(x<0){
            x=0;
        }
        if(x>1)
        {
            x=1;
        }
        unsigned char v=(unsigned char)(x*255.0f+0.5f);
        fwrite(&v,1,1,fp);
    }
    fclose(fp);
    return 0;
}
//  basically allocating memory 
double **a2d(int rows, int cols)
{
    double **mat = malloc(rows * sizeof(double*));
    for(int i=0;i<rows;i++)
        mat[i] = malloc(cols* sizeof(double));
    return mat;
}
// norm of columns as the function name sugested .
void normofcol(double **X, int N, int k)
{
    for(int j=0;j<k;j++){
        double norm=0.0;
        for(int i=0;i<N;i++) {
            norm+=X[i][j]*X[i][j];
        }
        norm=sqrt(norm);
        if(norm<1e-12) {
            norm=1e-12 ;
        }
        for(int i=0;i<N;i++){ 
        X[i][j]/=norm;
        }
    }
}
// evaluating A^TA
void ATA(double **A, double **C, int M, int N)
{
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++){
            double s=0.0;
            for(int k=0;k<M;k++)
                s=s+A[k][i]*A[k][j];
            C[i][j]=s;
        }
}
// evaluating A^TA*X
void ATAX(double **C, double **X, double **Y, int N, int k)
{
    for(int i=0;i<N;i++)
        for(int j=0;j<k;j++)
        {
            double sum=0.0;
            for(int t=0;t<N;t++){
                sum += C[i][t]*X[t][j];
                }
                    Y[i][j]=sum;
        }
}
// QR dioginalistion after every step of power iteration .
void QR(int m, int n, double **A, double **Q, double **R)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            R[i][j]=0.0;

    for(int j=0;j<n;j++)
    {
        for(int i=0;i<m;i++)
            Q[i][j]=A[i][j];

        for(int k=0;k<j;k++)
        {
            double dot=0.0;
            for(int i=0;i<m;i++)
                dot += Q[i][k]*A[i][j];
            R[k][j]=dot;
            for(int i=0;i<m;i++)
                Q[i][j]-=dot*Q[i][k];
        }

        double norm=0.0;
        for(int i=0;i<m;i++) norm+=Q[i][j]*Q[i][j];
        norm=sqrt(norm);
        // if norm is less then 1e-12 than in divide we lose numerical stability  .
        if(norm<1e-12){ 
            norm=1e-12;
        }
        R[j][j]=norm;
        for(int i=0;i<m;i++) Q[i][j]/=norm;
    }
}
// evaluating AV_i because we have to evaluate 
void AVi(double **A, double *v, double *out, int M, int N)
{
    for(int i=0;i<M;i++){
        double s=0.0;
        for(int j=0;j<N;j++) {
             s=s+A[i][j]*v[j];
        }
                out[i]=s;
    }
}

// Frobenius Norm
double fnorm(double **A, int M, int N) {
    double sum = 0.0;
    for(int i=0; i<M; i++)
        for(int j=0; j<N; j++)
            sum =sum+ A[i][j]*A[i][j];
    return sqrt(sum);
}

double foerror(double **A, float *A_k, int M, int N) {
    double sum = 0.0;
    for(int i=0; i<M; i++)
        for(int j=0; j<N; j++) {
            double diff = A[i][j] - (double)A_k[i*N + j];
            sum =sum+ diff*diff;
        }
    return sqrt(sum);
}
// free the storage 
void f2d(double **mat, int rows)
{
    for(int i=0;i<rows;i++)
        free(mat[i]);
    free(mat);
}

int main()
{
    const char *infile="image.pgm";
    const char *outfile="compressed120einstein.pgm";
    int w,h;
    float *img = repgm(infile,&w,&h);
    if(!img){printf("Failed to read input\n");return 1;}
    int M=h, N=w;

    printf("Image size: %dx%d\n",M,N);

    
    double **A = a2d(M,N);
    for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
            A[i][j] = (double)img[i*N+j];
    free(img);
    double **C = a2d(N,N);
    ATA(A,C,M,N);
    double **X = a2d(N,K);
    double **Y = a2d(N,K);
    double **Q = a2d(N,K);
    double **R = a2d(K,K);
    // making the random block matrix of (Nxk)
    srand(0);
    for(int i=0;i<N;i++)
        for(int j=0;j<K;j++)
            X[i][j]=(double)rand()/RAND_MAX;
    normofcol(X,N,K);
    // power iteration than compained by QR .
    for(int it=0;it<ITER;it++){
        ATAX(C,X,Y,N,K);
        QR(N,K,Y,Q,R);
        for(int i=0;i<N;i++)
            for(int j=0;j<K;j++){
                X[i][j]=Q[i][j];
            }
        printf("Iteration %d/%d done\n",it+1,ITER);
    }
    // evaluating sigma by ||AV_i||=sigma
    double sigma[K];
    double *tmp = malloc(M*sizeof(double));
    double *v = malloc(N*sizeof(double));
    for(int j=0;j<K;j++){
        for(int i=0;i<N;i++) {
            v[i]=X[i][j];
         }
          AVi(A,v,tmp,M,N);
        double s=0.0;
        for(int i=0;i<M;i++) s=s+tmp[i]*tmp[i];
        sigma[j]=sqrt(s);
    }
    printf("\nTop singular values (k=%d):\n",K);
    for(int j=0;j<K;j++)
        printf("sigma[%2d] = %.6f\n",j+1,sigma[j]);
       // evaluating  U_i
    double **U = a2d(M,K);
    for(int j=0;j<K;j++){
        for(int i=0;i<N;i++) v[i]=X[i][j];
        AVi(A,v,tmp,M,N);
        double s=sigma[j];
        for(int i=0;i<M;i++)
            U[i][j]=tmp[i]/s;
    }
    // A_k recontruction 
    float *out = calloc((size_t)M*N,sizeof(float));
    for(int r=0;r<K;r++){
        for(int i=0;i<M;i++){
            double u=U[i][r];
            for(int j=0;j<N;j++)
                out[i*N+j] =out[i*N+j]+ (float)(sigma[r]*u*X[j][r]);
        }
    }
    for(int i=0;i<M*N;i++){
        if(out[i]<0) out[i]=0;
        if(out[i]>1) out[i]=1;
    }
    wrpgm(outfile,out,N,M);
    printf("Wrote %s (k=%d)\n",outfile,K);

    // error Analysis
    double noA = fnorm(A, M, N);
    double noDiff = foerror(A, out, M, N);
    double erranal = noDiff / noA;

    printf("\nFrobenius Norm  A        : %.6f\n", noA);
    printf("Frobenius Norm  (A - A_k): %.6f\n", noDiff);
    printf(" Error analisis: %.6f\n", erranal);

    free(out);
    free(tmp);
    free(v);
    f2d(A,M);
    f2d(C,N);
    f2d(X,N);
    f2d(Y,N);
    f2d(Q,N);
    f2d(R,K);
    f2d(U,M);
    return 0;
}
