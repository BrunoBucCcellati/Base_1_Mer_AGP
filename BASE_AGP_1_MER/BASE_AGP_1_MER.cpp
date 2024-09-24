#include <iostream>
#include "ФУНКЦИИ И КЛАССЫ.h"
unsigned short chislo_iteraciy = 0;
pair <double, double> Base_1_Mer_Algoritm(double a, double b, double epsilon, double r, double epsilon_granichnoe)
{
    pair <double, double> min_xy;
    //
    pair<double, double>* start = new pair<double, double>;
    start->first = a;
    pair<double, double>* end = new pair<double, double>;
    end->first = b;
    //
    start->second = ShekelFunc(a);
    end->second = ShekelFunc(b);
    //
    Otrezok* otrezok = new Otrezok(start, end);
    //
    Otrezki* interval = new Otrezki(otrezok, r);
    //
    pair<double, double> pred_i_sled_shag;
    pred_i_sled_shag.first = a;
    pred_i_sled_shag.second = b;
    //
    while (abs(pred_i_sled_shag.second - pred_i_sled_shag.first) > epsilon)
    {
        min_xy.first = pred_i_sled_shag.second;
        min_xy.second = ShekelFunc(min_xy.first);
        pred_i_sled_shag.first = pred_i_sled_shag.second;
        //
        pred_i_sled_shag.second = Shag(interval->Getm(), interval->GetX_Rmax().first, interval->GetX_Rmax().second, interval->GetY_Rmax().first, interval->GetY_Rmax().second);
        if (pred_i_sled_shag.second >= b)
        {
            srand(time(0));
            pred_i_sled_shag.second = (b - epsilon_granichnoe) + (double)rand() / RAND_MAX * (b - (b - epsilon_granichnoe));
        }
        if (pred_i_sled_shag.second <= a)
        {
            srand(time(0));
            pred_i_sled_shag.second = a + (double)rand() / RAND_MAX * ((a + epsilon_granichnoe) - a);
        }
        //
        pair<double, double>* promejutochnaya_tochka = new pair<double, double>;
        promejutochnaya_tochka->first = pred_i_sled_shag.second;
        promejutochnaya_tochka->second = ShekelFunc(pred_i_sled_shag.second);
        interval->Add(promejutochnaya_tochka, a, b, epsilon_granichnoe);
        chislo_iteraciy++;
    }
    return min_xy;
}
int main()
{
    pair <double, double> Extr = Base_1_Mer_Algoritm(0.0, 10.0, 0.0005, 2.0, pow(10.0, -6.0));
    cout << "Xmin = " << Extr.first << " " << "Ymin = " << Extr.second << endl;
    cout << chislo_iteraciy << endl;
}