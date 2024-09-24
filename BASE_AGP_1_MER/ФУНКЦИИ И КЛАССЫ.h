#include <random>
#define M_PI 3.14159265358979323846
using namespace std;
double HillFunc(double x)
{
    const int HH_COUNT = 14;
    double a = 0.0;
    double b = 1.0;
    //
    default_random_engine generator;
    uniform_real_distribution<double> distribution(a, b);
    double A[20] = {};
    double B[20] = {};
    //
    for (int i = 0; i < 20; ++i)
    {
        A[i] = -1.1 + distribution(generator) * 2.0;
        B[i] = -1.1 + distribution(generator) * 2.0;
    }
    //
    double res = A[0];
    for (int i = 1; i < HH_COUNT; i++)
    {
        res += A[i] * sin(i * 2 * M_PI * x) + B[i] * cos(i * 2 * M_PI * x);
    }
    return res;
}
double ShekelFunc(double x)
{
    const int SH_COUNT = 10;
    double a = 0.0;
    double b = 1.0;
    //
    default_random_engine generator;
    uniform_real_distribution<double> distribution(a, b);
    double A[20] = {};
    double B[20] = {};
    double K[20] = {};
    //
    for (int i = 0; i < 20; ++i)
    {
        A[i] = 10.0 * distribution(generator);
        B[i] = 1.0 + 0.2 * distribution(generator) * 2.0;
        K[i] = 5.0 + 20 * distribution(generator);
    }
    //
    double res = 0.0;
    for (int i = 1; i < SH_COUNT; i++)
    {
        res -= 1.0 / ((K[i] * pow((x - A[i]), 2.0) + B[i]));
    }
    return res;
}
double Characteristic(double m, double x1, double x2, double y1, double y2)
{
    return (m * (x2 - x1) + pow((y2 - y1), 2) / (m * (x2 - x1)) - 2 * (y2 + y1));
}
double Shag(double m, double x1, double x2, double y1, double y2)
{
    return ((x1 + x2) / 2) - ((y2 - y1) / (2 * m));
}
//
class Otrezok
{
protected:
    pair<double, double>* start;
    pair<double, double>* end;
    double M;
    double R;
    Otrezok* Next;
public:
    Otrezok(pair<double, double>* _start, pair<double, double>* _end)
    {
        start = _start;
        end = _end;
        M = abs((end->second - start->second) / (end->first - start->first));
        R = 0;
        Next = this;
    };
    void ChangeCharacteristic(double m)
    {
        R = Characteristic(m, start->first, end->first, start->second, end->second);
    }
    double GetCharacteristic()
    {
        return R;
    }
    void SetEnd(pair<double, double>* _end)
    {
        end = _end;
        M = abs((end->second - start->second) / (end->first - start->first));
    }
    pair<double, double>* GetEnd()
    {
        return end;
    }
    pair<double, double>* GetStart()
    {
        return start;
    }
    Otrezok* GetNext()
    {
        return Next;
    }
    void SetNext(Otrezok* _Next)
    {
        Next = _Next;
    }
    double GetM()
    {
        return M;
    }
    void SetM(double _M)
    {
        M = _M;
    }
};
//
class Otrezki
{
protected:
    Otrezok* Head;
    double Mmax;
    double m;
    double r;
    pair<double, double> x_Rmax;
    pair<double, double> y_Rmax;
public:
    Otrezki(Otrezok* _Head, double _r)
    {
        Head = _Head;
        r = _r;
        if (Head->GetM() != 0)
        {
            Mmax = Head->GetM();
            m = r * Mmax;
        }
        else
        {
            Mmax = 0;
            m = 1;
        }
        Head->ChangeCharacteristic(m);
        x_Rmax.first = Head->GetStart()->first;
        x_Rmax.second = Head->GetEnd()->first;
        y_Rmax.first = Head->GetStart()->second;
        y_Rmax.second = Head->GetEnd()->second;
    };
    void Add(pair<double, double>* tmp, double _a, double _b, double _epsilon_granichnoe)
    {
        Otrezok* curr = Head;
        while (curr->GetEnd()->first < tmp->first)
        {
            curr = curr->GetNext();
        }
        Otrezok* curr1 = new Otrezok(tmp, curr->GetEnd());
        curr->SetEnd(tmp);
        curr1->SetNext(curr->GetNext());
        curr->SetNext(curr1);
        //
        curr = Head;
        Otrezok* Otrezok_Rmax = curr;
        Mmax = curr->GetM();
        curr = curr->GetNext();
        while (curr != Head)
        {
            if (curr->GetM() > Mmax)
            {
                Mmax = curr->GetM();
                Otrezok_Rmax = curr;
            }
            curr = curr->GetNext();
        }
        Otrezok_Rmax->SetM(0.0);
        if (Mmax != 0)
        {
            m = r * Mmax;
        }
        else
        {
            m = 1;
        }
        //
        curr = Head;
        while (curr->GetNext() != Head)
        {
            curr = curr->GetNext();
            curr->ChangeCharacteristic(m);
        }
        curr = curr->GetNext();
        curr->ChangeCharacteristic(m);
        //
        double Rmax = -DBL_MAX;
        if (Shag(m, curr->GetStart()->first, curr->GetEnd()->first, curr->GetStart()->second, curr->GetEnd()->second) <= _b - _epsilon_granichnoe && Shag(m, curr->GetStart()->first, curr->GetEnd()->first, curr->GetStart()->second, curr->GetEnd()->second) >= _a + _epsilon_granichnoe)
        {
            Rmax = curr->GetCharacteristic();
        }
        curr = curr->GetNext();
        while (curr != Head)
        {
            if (curr->GetCharacteristic() > Rmax && Shag(m, curr->GetStart()->first, curr->GetEnd()->first, curr->GetStart()->second, curr->GetEnd()->second) <= _b - _epsilon_granichnoe && Shag(m, curr->GetStart()->first, curr->GetEnd()->first, curr->GetStart()->second, curr->GetEnd()->second) >= _a + _epsilon_granichnoe)
            {
                Rmax = curr->GetCharacteristic();
                x_Rmax.first = curr->GetStart()->first;
                x_Rmax.second = curr->GetEnd()->first;
                y_Rmax.first = curr->GetStart()->second;
                y_Rmax.second = curr->GetEnd()->second;
            }
            curr = curr->GetNext();
        }
    }
    double Getm()
    {
        return m;
    }
    pair <double, double> GetX_Rmax()
    {
        return x_Rmax;
    }
    pair <double, double> GetY_Rmax()
    {
        return y_Rmax;
    };
};
