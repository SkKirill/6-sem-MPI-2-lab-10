#include <iostream>
#include <vector>
#include <mpi.h>
#include <string>
#include <random>
#include <vector>
#include <math.h>
#include <iomanip>
#include <iostream>
using namespace std;

// seriesValue - функция, которая вовращает значение заданной по условию функции
// в нее передается точка в которой нужно расчитать (х) и погрешность, с которой расчитываем
double seriesValue(double x, double eps) {
    long double curRes = 1 / ((1 + x) * (1 + x)); // расчитываем точное значение, для сравнения с eps
    long double result = 1; // первое значение в ряду
    long double power = 1; // степень икса, на которую умнажается каждый элемент ряда
    long double current = 1;

    for (int i = 2; fabs(curRes - result) > eps; i++) {
        current *= -1;
        power *= x;// преобраовываем предыдущий элемент последовательности в новый
        result += current * (i * power); // добавляем значение к нашему реультату
    }
    return result; // возвращаем полученное значение
}

int main(int argc, char* argv[]) {

    int rank, nprocs, N; // rank - это номер процесса, nproc - это количество процессов, N - это количество точек на отрезке от a до b, в которых нужно посчитать функцию
    double a, b, eps; // a - начало отрезка, b - конец отрека, eps - значение точности расчетов

    vector<double> pointValues; // - вектор хранящий все точки нужные для расчета
    vector<double> valuesXForProc; // точки, которые потребуются для расчета конкретному процессу
    vector<double> results; // реультаты расчетов в точках каждого процесса
    vector<double> rootResults; // все результаты из всех процессов
    vector<double> exactValues; // точные значение расчитанные функциями C++

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // получаем значение номера процеса и помещаем в rank
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); // получаем значение, хранящее общее колличество процессов

    if (rank == 0) {
        // просим польователя ввести значения //
        cout << "Введите a:";
        cin >> a;
        cout << "Введите b:";
        cin >> b;
        cout << "Введите eps:";
        cin >> eps;
        cout << "Введите N:";
        cin >> N;

        double interval = (b - a) / N; // Вычисляем интервал между точками, чтобы знать все точки в которых нужно вычислять

        pointValues.push_back(a); // заносим в вектор точек начальную точку
        for (int i = 1; i < N - 1; i++)
        {
            pointValues.push_back(a + interval * i); // идем и добавляем все остальные точки в которых нужно расчитать функцию
        }
        pointValues.push_back(b); // заносим конечную точку

        /*
        после этих действий у нас будет вектор pointValues в котором есть все точки, которые потребуются для вычисления
        так же переменные
        a - начальная точка
        b - конечная точка
        eps - точность с которой нужно расчитать
        */
    };

    MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // отправляем из процесса 0 ранга значение переменной eps всем остальным в переменную eps
    MPI_Bcast(&N, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // отправляем из процесса 0 ранга значение переменной N всем остальным в переменную N

    int countPerProc; // countPerProc - количество точек расчитываемых в каждом процессе
    if ((N % nprocs) > 0) {
        countPerProc = (int)(N / nprocs) + 1; // если например процессов 7, а точек 12, то у тебя на процесс приходится 1.7точка.. В таком случае C++ округляет 1.7 в 1 и я добавляю еще 1 получаем 2
    }                                   // в такой ситуации последний процесс, тоесть 7 не будет считать ни в какой точке, но зато у тебя посчитаются во всех. 
    else {
        countPerProc = N / nprocs;
    }

    valuesXForProc.resize(countPerProc);

    if (pointValues.size() < countPerProc * nprocs) {
        pointValues.resize(countPerProc * nprocs); // в случае описанном на 100 строке нам понадобиться передать последнему процессу какието точки, 
        //для этого добавляем в pointValues в конец нулевые точки
    }

    MPI_Scatter(pointValues.data(), countPerProc, MPI_DOUBLE, valuesXForProc.data(), countPerProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // мы отправляем созданные точки из pointValues в количестве countPerProc в valuesXForProc в каждом процессе в количестве countPerProc, 
    //                          передаваемые типы MPI_DOUBLE, происходит в MPI_COMM_WORLD

    /*
    тут начинается расчет с каждого процесса
    */

    for (double point : valuesXForProc) // это тоже самое, что в C# foreach(double point in valuesXForProc)
    {
        results.push_back(seriesValue(point, eps)); // считаем значение в каждой точке и заносим в вектор results
    }

    rootResults.resize(countPerProc * nprocs); // создаем вектор в который будут помещаться все результаты

    MPI_Gather(results.data(), countPerProc, MPI_DOUBLE, rootResults.data(), countPerProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // поместили все результаты в rootResults в 0 процесс

    if (rank == 0) // 0 процессом выводим все полученные данные
    {
        for (double point : pointValues)
        {
            long double curRes = 1 / ((1 + point) * (1 + point));
            exactValues.push_back(curRes); // расчитываем для всех точек точное значение и помещаем в переменную exactValues
        }

        cout << setw(25) << "Точка" << setw(40) << "Значение функции" << setw(40) << "Точное значение" << endl;

        for (int i = 0; i < N; i++) // Так как при выводу мы идем до N, те точки, которые мы добавили на 109 строке выводиться не будут, хотя процесс расчитает и для них.
        {
            cout << setw(25) << pointValues[i] << setw(40) << rootResults[i] << setw(40) << exactValues[i] << endl;
        }
    }

    MPI_Finalize();
    return 0;
}