#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <cstdlib>
using namespace std;
const float zero=1e-12; //stala przyblizenia zera

//Podstawowy////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool gauss_podstawowy(int n, float **macierz, float *wynik)
{
    float robocza,rob;
    int i,j,k;

//eliminacja wspolczynników
    for(i=0; i<n-1; i++)
    {
        for(j=i+1; j<n; j++)
        {
            if(fabs(macierz[i][i])<zero)//funkcja 'fabs' wartosc bezwgledna
                return false;

            robocza=((-macierz[j][i])/(macierz[i][i]));

            for(k=i+1; k<=n; k++)
            {
                macierz[j][k]=((macierz[j][k])+(robocza * macierz[i][k]));
            }
        }
    }

//obliczenie  niewiadomych i wyswietlenie
    for(i=n-1; i>=0; i--)
    {
        rob = macierz[i][n];
        for(j=n-1; j>=i+1; j--)
        {
            rob=(rob-(macierz[i][j] * wynik[j]));
        }

        if(fabs(macierz[i][i])<zero)
            return false;

        wynik[i]=(rob/macierz[i][i]);
    }
    return true;
}

//Kolumny i wiersze////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int znajdz_max_w_kolumnie(float **macierz,int kolumna,int wiersz_poczatkowy,int il_wierszy)//max element w kolumnie
{
    float maksimum=macierz[wiersz_poczatkowy][kolumna];
    int max_wiersz=wiersz_poczatkowy;
    for(int i=wiersz_poczatkowy+1; i<il_wierszy; i++)
    {
        if(macierz[i][kolumna]>maksimum)
        {
            maksimum=macierz[i][kolumna];
            max_wiersz=i;
        }
    }
    return max_wiersz;
}

void zamiana_wierszy(float **macierz,int n,int rzad_1,int rzad_2)//zamiana wierszy, max
{
    if(rzad_1!=rzad_2)
        for(int i=0; i<n; i++)
        {
            swap(macierz[rzad_1][i],macierz[rzad_2][i]);//funkcja 'swap' zamienia
        }
}

void macierz_trojkatna(float **macierz,int n)//utworzenie macierzy trojkatnej
{
    float wspolczynnik;
    int max_w_kolumnie;

    for(int k=0; k<n-1; k++)
    {
        max_w_kolumnie=znajdz_max_w_kolumnie(macierz,k,k,n);
        zamiana_wierszy(macierz,n+1,max_w_kolumnie,k);

        for(int i=k+1; i<n; i++)
        {
            wspolczynnik=macierz[i][k]/macierz[k][k];
            for(int j=k; j<=n; j++)
            {
                macierz[i][j]=((macierz[i][j])-((float)wspolczynnik*macierz[k][j]));
                if(macierz[i][j]<zero && macierz[i][j]>-zero)
                {
                    macierz[i][j]=0;
                }
            }
        }
    }

    for(int i=n-1; i>=0; i--)
    {
        for(int j=n; j>=i; j--)
        {
            macierz[i][j]=macierz[i][j]/macierz[i][i];
        }
    }
}

void wyniki(float **macierz,float *t,int n)//obliczenie niewiadomych oraz ich wyswietlenie
{
    t[n-1]=macierz[n-1][n];

    for(int i=n-2; i>=0; i--)
    {
        t[i]=macierz[i][n];
        for(int j=n-1; j>i; j--)
        {
            t[i]=((t[i])-(macierz[i][j]*t[j]));
        }
    }
}

//Pelny///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int znajdz_max_w_macierzy(float **macierz,int wiersz_poczatkowy,int il_wierszy)
{
    float maksimum=macierz[wiersz_poczatkowy][wiersz_poczatkowy];
    int max_wiersz=wiersz_poczatkowy;

    for(int i=wiersz_poczatkowy; i<il_wierszy; i++)
    {
        for(int j=wiersz_poczatkowy; j<il_wierszy; j++)
        {
            if(macierz[i][j]>maksimum)
            {
                maksimum=macierz[i][j];
                max_wiersz=i;
            }
        }
    }

    return max_wiersz;
}

void macierz_trojkatna_2(float **macierz,int n)
{
    float wspolczynnik;
    int max_w_kolumnie;
    for(int k=0; k<n-1; k++)
    {
        max_w_kolumnie=znajdz_max_w_macierzy(macierz,k,n);
        zamiana_wierszy(macierz,n+1,max_w_kolumnie,k);

        for(int i=k+1; i<n; i++)
        {
            wspolczynnik=macierz[i][k]/macierz[k][k];
            for(int j=k; j<=n; j++)
            {
                macierz[i][j]=((macierz[i][j])-((float)wspolczynnik*macierz[k][j]));
                if(macierz[i][j]<zero && macierz[i][j]>-zero)
                {
                    macierz[i][j]=0;
                }
            }
        }
    }

    for(int i=n-1; i>=0; i--)
    {
        for(int j=n; j>=i; j--)
        {
            macierz[i][j]=macierz[i][j]/macierz[i][i];
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{

    cout<<"1. Gauss podstawowy"<<endl;
    cout<<"2. Gauss kolumna"<<endl;
    cout<<"3. Gauss wiersz"<<endl;
    cout<<"4. Gauss pelny"<<endl;
    cout<<"5. Wyjscie"<<endl;

    char wybor;

    while(wybor!=5)
    {
        cout<<endl;
        cout<<"Podaj numer: "<<endl;
        cin>>wybor;

        switch(wybor)
        {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        case '1':
        {
            int n,i,j;
            float **macierz;
            float *wynik;

            cout<<"Podaj n (ilosc wierszy): "<<endl;
            cin>>n;

            macierz=new float * [n];
            wynik=new float [n];

            for(i=0; i<n; i++)
            {
                macierz[i]=new float[n+1];
            }

            for(i=0; i<n; i++)
            {
                for(j=0; j<=n; j++)
                {
                    cout<<"Podaj element macierzy ["<<i+1<<","<<j<<"]"<<endl;
                    cin>>macierz[i][j];
                }
            }

            cout<<endl;
            if(gauss_podstawowy(n,macierz,wynik))
            {
                for(i=0; i<n; i++)
                    cout<<"x"<<i+1<<" = "<<wynik[i]<<endl;
            }
            else
            {
                cout<<"Wystapilo dzielenie przez zero!"<<endl;
            }

            cin.get();
            break;
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        case '2':
        {

            int n;
            cout<<"Podaj n (ilosc wierszy): "<<endl;
            cin>>n;

            float **macierz=new float*[n];
            float *wynik=new float[n];
            for(int i=0; i<n; i++)
            {
                macierz[i]=new float[n+1];
            }

            for(int i=0; i<n; i++)
            {
                for(int j=0; j<=n; j++)
                {
                    cout<<"Podaj element macierzy ["<<i+1<<","<<j<<"]"<<endl;
                    cin>>macierz[i][j];
                }
            }

            macierz_trojkatna(macierz,n);
            wyniki(macierz,wynik,n);
            cout<<endl;
            for(int i=0; i<n; i++)
            {
                cout<<"x"<<i+1<<" = "<<wynik[i]<<endl;
            }

            cin.get();
            break;
        }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        case '3':
        {


            int n;
            cout<<"Podaj n (ilosc wierszy): "<<endl;
            cin>>n;

            float **macierz=new float*[n];
            float *wynik=new float[n];
            for(int i=0; i<n; i++)
            {
                macierz[i]=new float[n+1];
            }

            for(int i=0; i<n; i++)
            {
                for(int j=0; j<=n; j++)
                {
                    cout<<"Podaj element macierzy ["<<i+1<<","<<j<<"]"<<endl;
                    cin>>macierz[i][j];
                }
            }

            macierz_trojkatna(macierz,n);
            wyniki(macierz,wynik,n);
            cout<<endl;
            for(int i=0; i<n; i++)
            {
                cout<<"x"<<i+1<<" = "<<wynik[i]<<endl;
            }


            cin.get();
            break;
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        case '4':
        {
            int n;
            cout<<"Podaj n (ilosc wierszy): "<<endl;
            cin>>n;

            float **macierz=new float*[n];
            float *wynik=new float[n];
            for(int i=0; i<n; i++)
            {
                macierz[i]=new float[n+1];
            }

            for(int i=0; i<n; i++)
            {
                for(int j=0; j<=n; j++)
                {
                    cout<<"Podaj element macierzy ["<<i+1<<","<<j<<"]"<<endl;
                    cin>>macierz[i][j];
                }
            }

            macierz_trojkatna_2(macierz,n);
            wyniki(macierz,wynik,n);
            cout<<endl;
            for(int i=0; i<n; i++)
            {
                cout<<"x"<<i+1<<" = "<<wynik[i]<<endl;
            }

            cin.get();
            break;
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        case '5':
        {
            return 0;
            break;
        }

        }
    }
    return 0;
}
