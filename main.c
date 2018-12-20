#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "eispack.h"

float ozdegerler[100];
float ozvektorler[100];
float dizi_schur[2];
int main()
{
    int kac_kere;
    printf("Carkin donme sayisini giriniz:");
    scanf("%d",&kac_kere);
    cark(kac_kere);
    system("PAUSE");
    return 0;
}
void kokbul(float a,float b,float c,float d){
    double f = (3 * c / a - b*b / (a*a)) / 3;
	double g = (2 * b*b*b / (a*a*a) - 9 * b*c / (a*a) + 27 * d / a) / 27;
	double h = g*g / 4 + f*f*f / 27;
	double i = sqrt(g*g / 4 - h);
	double j = pow(i, 1.0/3.0);
	double k = acos(-1 * (g /(2 * i)));
	double l = j*-1;
	double m = cos(k / 3);
	double n = sqrt(3) * sin(k / 3);
	double p = (b / (3 * a))*-1;
	double x1 = 2 * j*cos(k / 3) - (b / (3*a));
	double x2 = l*(m + n) + p;
	double x3 = l*(m - n) + p;
	ozdegerler[0]=x1;
    ozdegerler[1]=x2;
    ozdegerler[2]=x3;
}
void matris_eleman_gir(float *matris,int boyut)
{
    int counter=0,counter_2=0;
    for(counter=0; counter<boyut; counter++)
    {
        for(counter_2=0; counter_2<boyut; counter_2++)
        {
            printf("%d.Satir %d.Sutun Elemanini giriniz:",counter+1,counter_2+1);
            scanf("%f",&*(matris+counter*boyut+counter_2));
        }
    }
}
void schur (float *matris,int boyut){
 int i,j;
 double sinir=0;
 double ozdegerler_toplami=0;
 for(i=0;i<boyut;i++)
    {
        ozdegerler_toplami+=pow(ozdegerler[i],2);
        for(j=0;j<boyut;j++)
            {
                sinir+=pow(*(matris+i*boyut+j),2);
            }
    }
    dizi_schur[0]=ozdegerler_toplami;
    dizi_schur[1]=sinir;
}
void ozdegerbul(float *matris,int boyut)
{
    if(boyut==1)
    {
        ozdegerler[0]=*(matris+0*boyut+0);
    }
    if(boyut==2)
    {
        int i,j;
        float a,b,c,delta;
        float dizi1[2]= {1,1};
        float dizi2[2]= {1,1};
        float kok1,kok2;
        for(i=0; i<boyut; i++)
        {
            for(j=0; j<boyut; j++)
            {
                *(matris+i*boyut+j)*=-1;
            }
        }
        dizi1[1]=*(matris+0*boyut+0);
        dizi2[1]=*(matris+1*boyut+1);
        a=dizi1[0] * dizi2[0];
        b=(dizi1[1] * dizi2[0])+(dizi1[0] * dizi2[1]);
        c=dizi1[1] * dizi2[1];
        c -= *(matris+0*boyut+1)* *(matris+1*boyut+0);
        delta=(b*b)-(4*a*c);
        kok1=(-b-sqrt(delta))/(2*a);
        kok2=(-b+sqrt(delta))/(2*a);
        ozdegerler[0]=kok1;
        ozdegerler[1]=kok2;
        printf("k1=%.2f\nk2=%.2f\n",kok1,kok2);
    }
    if(boyut==3)
    {
        double a,b,c,d,e,f,k,l,m,katsayi,katsayi_2,katsayi_3,katsayi_4,kok_1,kok_2,kok_3;
        a=*(matris+0*boyut+0);
        b=*(matris+0*boyut+1);
        c=*(matris+0*boyut+2);
        d=*(matris+1*boyut+0);
        e=*(matris+1*boyut+1);
        f=*(matris+1*boyut+2);
        k=*(matris+2*boyut+0);
        l=*(matris+2*boyut+1);
        m=*(matris+2*boyut+2);
        katsayi=(double)-1;
        katsayi_2=(double)(a+m+e);
        katsayi_3=(double)(c*k)+(f*l)+(b*d)-(a*m)-(e*m)-(a*e);
        katsayi_4=(double)(d*l*c)+(f*b*k)+(a*e*m)-(c*e*k)-(f*a*l)-(m*b*d);
        kokbul(katsayi,katsayi_2,katsayi_3,katsayi_4);
    }
    /*if(boyut==4){

        double a,b,c,d,e,f,g,h,s,i,j,k,l,m,n,o,dort,dort_1,dort_2,dort_3,dort_4;

        a=*(matris+0*boyut+0);
        b=*(matris+0*boyut+1);
        c=*(matris+0*boyut+2);
        d=*(matris+0*boyut+3);
        e=*(matris+1*boyut+0);
        f=*(matris+1*boyut+1);
        g=*(matris+1*boyut+2);
        h=*(matris+1*boyut+3);
        s=*(matris+2*boyut+0);
        i=*(matris+2*boyut+1);
        j=*(matris+2*boyut+2);
        k=*(matris+2*boyut+3);
        l=*(matris+3*boyut+0);
        m=*(matris+3*boyut+1);
        n=*(matris+3*boyut+2);
        o=*(matris+3*boyut+3);

        dort=(double) 1;
        dort_1=(double) -1*(j+o+f+a);
        dort_2=(double) ((j+o+f+a)-(h*m)-(k*n)-(i*g)+(j*o)+(f*j)+(f*o)-(b*e)-(c*s)-(d*l));
        dort_3=(double) (h*a*m)+(k*a*n)+(i*g*a)-(j*o*a)-(f*j*a)-(f*o*a)+(h*j*m)+(k*n*f)+(i*g*o)-(f*j*o)-(i*n*h)-(m*g*k)+(g*s*b)+(h*l*b)-(j*e*b)-(o*e*b)+(o*s*c)+(f*s*c)-(e*i*c)-(l*k*c)+(f*l*d)+(l*j*d)-(s*n*d)-(e*m*d);
        dort_4=(double) (f*j*o*a)+(i*n*h*a)+(m*g*k*a)-(h*j*m*a)-(k*n*f*a)-(i*g*o*a)+(h*l*j*b)+(k*n*e*b)+(o*g*s*b)-(e*j*o*b)-(s*n*h*b)-(l*g*k*b)+(e*i*o*c)+(s*m*h*c)+(l*f*k*c)-(h*i*l*c)-(e*m*k*c)-(o*f*s*c)+(g*i*l*d)+(e*m*j*d)+(f*s*n*d)-(e*i*n*d)-(s*m*g*d)-(f*j*l*d);
        }*/
    /*void kokbul3(double *x1, double *x2, double *x3, double a, double b, double c, double d)
{
	double f = (3 * c / a - b*b / (a*a)) / 3;
	double g = (2 * b*b*b / (a*a*a) - 9 * b*c / (a*a) + 27 * d / a) / 27;
	double h = g*g / 4 + f*f*f / 27;
	double i = sqrt(g*g / 4 - h);
	double j = pow(i, 1.0 / 3.0);
	double k = acos(-1 * (g / (2 * i)));
	double l = j*-1;
	double m = cos(k / 3);
	double n = sqrt(3) * sin(k / 3);
	double p = (b / (3 * a))*-1;
	*x1 = 2 * j*cos(k / 3) - (b / (3 * a));
	*x2 = l*(m + n) + p;
	*x3 = l*(m - n) + p;
}

void kokbul4(double a, double b, double c, double d, double e)
{
	double f = c - (3 * b*b / 8);
	double g = d + (b*b*b / 8) - (b*c / 2);
	double h = e - (3 * b*b*b*b / 256) + (b*b*c / 16) - (b*d / 4);
	double a1 = 1;
	double a2 = (f / 2);
	double a3 = (f*f - 4 * h) / 16;
	double a4 = -1*g*g / 64;
	double y1,y2,y3;
	kokbul3(&y1,&y2,&y3,a1,a2,a3,a4);
	double km = -1 * pow(10, -7);
	double kmm = -1 * pow(10, -25);
	if (y1 > km && y1 < kmm){ y1 = 0; }
	if (y2 > km && y2 < kmm){ y2 = 0; }
	if (y3 > km && y3 < kmm){ y3 = 0; }
	double p, q;
	if (y1 == 0){
		p = sqrt(y2);
		q = sqrt(y3);
	}else if (y2 == 0){
		p = sqrt(y1);
		q = sqrt(y3);
	}else{
		p = sqrt(y1);
		q = sqrt(y2);
	}
	double r = -1 * g / (8 * p*q);
	double s = b / (4 * a);
	double x1 = p + q + r - s;
	double x2 = p - q - r - s;
	double x3 = -1 * p + q - r - s;
	double x4 = -1 * p - q + r - s;
	printf("x1 : %lf\nx2 : %lf\nx3 : %lf\nx4 : %lf\n", x1, x2, x3, x4);
}*/
    if (boyut>=4)
    {
    int matz = 1;
    int i;
    double* a = ( double * ) malloc (boyut * boyut * sizeof ( double ) );
    double* w = ( double * ) malloc (boyut * sizeof ( double ) );
    double* x = ( double * ) malloc (boyut * boyut * sizeof ( double ) );
    for(i = 0; i < boyut*boyut;i++)
    {
        double v;
        printf(" Matris[%d][%d] : ", i/boyut + 1, i%boyut + 1);
        scanf("%lf", &v);
        a[i] = v;
        *(matris+i)=v;
    }
    int ierr = rs ( boyut, a, w, matz, x );
    if ( ierr != 0 )
    {
        printf ( "\nHata :(\n" );
        return;
    }
    for ( i = 0; i < boyut; i++ )
    {
        ozdegerler[i]=w[i];
    }
    return 0;
    }
}
void ozvektorbul(float *matris,int boyut)
{
    int i;
    double* a = ( double * ) malloc (boyut * boyut * sizeof ( double ) );
    double* w = ( double * ) malloc (boyut * sizeof ( double ) );
    double* x = ( double * ) malloc (boyut * boyut * sizeof ( double ) );
    for(i = 0; i < boyut*boyut;i++)
    {
        double v;
        printf(" Matris[%d][%d] : ", i/boyut + 1, i%boyut + 1);
        scanf("%lf", &v);
        a[i] = v;
        *(matris+i)=v;
    }
	int matz = 1;
	int ierr = rs ( boyut, a, w, matz, x );
    for(i=0;i<boyut;i++)
    {
        ozdegerler[i]=w[i];
    }
    for ( i = 0; i < boyut; i++ )
    {
        ozvektorler[i]=x[i];
    }
}
int nilpotent(float *matris,int boyut)
{
    float yeni_matris[boyut][boyut];
    int i,j,k;
    float toplam=0;
    for(i=0; i<boyut; i++)
    {
        for(j=0; j<boyut; j++)
        {
            for(k=0; k<boyut; k++)
            {
                toplam += (*(matris+i*boyut+k)) * (*(matris+k*boyut+j));
            }
            yeni_matris[i][j]=toplam;
            toplam=0;
        }
    }
    int test=1;
    for(i=0; i<boyut; i++)
    {
        for(j=0; j<boyut; j++)
        {
            *(matris+i*boyut+j)=yeni_matris[i][j];
        }
    }
    for(i=0; i<boyut; i++)
    {
        for(j=0; j<boyut; j++)
        {
            if(yeni_matris[i][j]!=0)
            {
                test=0;
                break;
            }
        }
        if(test==0)
        {
            break;
        }
    }
    if(test)
        {
            printf("\nGirilen Matris Nilpotenttir.\n");
        }
        else
            printf("\nGirilen Matris Nilpotent Degildir.\n");
    return test;
}
void cark(int kac_kere)
{
    srand(time(NULL));
    int islem;
    int rastgele_sayi;
    int counter=0,counter_2=0;
    int boyut;
    for(counter=0; counter<kac_kere; counter++)
    {
        rastgele_sayi=rand()%241;
        printf("\n%d.sayi=%d\n",counter+1,rastgele_sayi);
        if(counter==kac_kere-1)
        islem=rastgele_sayi;
    }
    islem=islem%4;
    switch(islem)
    {
    case(0):
    {
        FILE *dosya;
        dosya=fopen("ozdeger.txt","a");
        printf("Ozdeger Bul Islemindesiniz.\nIslem yapilacak kare matrisin boyutunu giriniz:");
        scanf("%d",&boyut);
        float matris[boyut][boyut];
        if(boyut<=3){
        matris_eleman_gir(*matris,boyut);
        ozdegerbul(*matris,boyut);
        }
        if(boyut>=4){
            ozdegerbul(*matris,boyut);
        }
        fprintf(dosya,"Girilen Matris=\n");
        for(counter=0;counter<boyut;counter++)
        {
            fprintf(dosya,"|");
            for(counter_2=0;counter_2<boyut;counter_2++)
            {
                fprintf(dosya," %.1f ",matris[counter][counter_2]);
            }
            fprintf(dosya,"|\n");
        }
        fprintf(dosya,"Özdeðerler:");
        for(counter=0;counter<boyut;counter++)
        {
            fprintf(dosya,"%d.Özdeðer=%.3f   ",counter+1,ozdegerler[counter]);
        }
        fprintf(dosya,"\n");
        fclose(dosya);
        break;
    }
    case(1):
    {
        srand(time(NULL));
        FILE *dosya;
        dosya=fopen("nilpotent.txt","a");
        printf("Nilpotent Bul Islemindesiniz.\nIslem yapilacak kare matrisin boyutunu giriniz:");
        scanf("%d",&boyut);
        float matris[boyut][boyut];
        int i,j;
        for(i=0; i<boyut; i++)
        {
            for(j=0; j<boyut; j++)
            {
                matris[i][j]=(-15+rand()%30);
            }
        }
        fprintf(dosya,"Girilen Matris=\n");
        for(i=0;i<boyut;i++)
        {
            fprintf(dosya,"|");
            for(j=0;j<boyut;j++)
            {
                fprintf(dosya," %.2f ",matris[i][j]);
            }
            fprintf(dosya,"|\n");
        }
        int test=nilpotent(*matris,boyut);
        fprintf(dosya,"Sonuc Matrisi=\n");
        for(i=0;i<boyut;i++)
        {
            fprintf(dosya,"|");
            for(j=0;j<boyut;j++)
            {
                fprintf(dosya," %.2f ",matris[i][j]);
            }
            fprintf(dosya,"|\n");
        }
        if(test)
        {
            fprintf(dosya,"\nGirilen Matris Nilpotenttir.\n");
        }
        else
            fprintf(dosya,"\nGirilen Matris Nilpotent Degildir.\n");
        fclose(dosya);
        break;
    }
        case(2):
            {
            int i,j;
            FILE *dosya;
            dosya=fopen("schur.txt","a");
            printf("Schur Islemindesiniz.\nIslem yapilacak kare matrisin boyutunu giriniz:");
            scanf("%d",&boyut);
            float matris[boyut][boyut];
            if(boyut<=3)
            {
                matris_eleman_gir(*matris,boyut);
                ozdegerbul(*matris,boyut);
            }
            else
            {
                ozdegerbul(*matris,boyut);
            }
            fprintf(dosya,"Girilen Matris=\n");
            for(i=0;i<boyut;i++)
            {
                fprintf(dosya,"|");
                for(j=0;j<boyut;j++)
                {
                    fprintf(dosya," %.2f ",matris[i][j]);
                }
                fprintf(dosya,"|\n");
            }
            fprintf(dosya,"Özdeðerler:");
            for(i=0;i<boyut;i++)
                {
                    fprintf(dosya,"%d.Özdeðer=%.3f   ",i+1,ozdegerler[i]);
                }
            fprintf(dosya,"\n");
            schur(*matris,boyut);
            fprintf(dosya,"Sýnýr Deðer=%f\tSchur Deðeri=%f\n",dizi_schur[1],dizi_schur[0]);
            if(dizi_schur[0]<=dizi_schur[1])
            {
                fprintf(dosya,"%f<=%f ==>Schur Teoremi Saðlanýr.\n",dizi_schur[0],dizi_schur[1]);
            }
            else
            {
                 fprintf(dosya,"%f>=%f ==>Schur Teoremi Saðlanmaz.\n",dizi_schur[0],dizi_schur[1]);
            }
            fclose(dosya);
            break;
            }
            case(3):
            {
            int i,j;
            FILE *dosya;
            dosya=fopen("ozvektor.txt","a");
            printf("Ozvektor Bul Islemindesiniz\nIslem yapilacak kare matrisin boyutunu giriniz:");
            scanf("%d",&boyut);
            float matris[boyut][boyut];
            ozvektorbul(*matris,boyut);
            if(boyut<=3)
            {
                ozdegerbul(*matris,boyut);
                for(i=0;i<boyut*boyut;i++)
                {
                    fprintf(dosya,"\n%d.ozvektor:%lf\n",i%boyut+1,ozvektorler[i]);
                    if(i==boyut-1)
                        break;
                }
            }
            else
            {
                for(i=0;i<boyut*boyut;i++)
                {
                fprintf(dosya,"\n%d.ozvektor:%lf\n",i%boyut+1,ozvektorler[i]);
                if(i==boyut-1)
                    break;
                }
            }
            fprintf(dosya,"Girilen Matris=\n");
            for(i=0;i<boyut;i++)
            {
                fprintf(dosya,"|");
                for(j=0;j<boyut;j++)
                {
                    fprintf(dosya," %.2f ",matris[i][j]);
                }
            fprintf(dosya,"|\n");
            }
            fprintf(dosya,"Özdeðerler:");
            for(i=0;i<boyut;i++)
            {
                fprintf(dosya,"%d.Özdeðer=%.3f   ",i+1,ozdegerler[i]);
            }
            fclose(dosya);
            break;
            }
            }
    }

