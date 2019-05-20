#include <iostream>
#include "matrix.h"

const int f_n = 20;
const int f_m = 8;

void cast_matrix(std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr.size() - 1; ++i)
        for (int j = 0; j < curr.size(); ++j)
            if (i != j)
                curr[i][j] = (i + 1.0 + j + 1.0) / (f_n + f_m);
            else
                curr[i][j] = f_n + f_m * f_m + (j + 1.0) / f_m + (i + 1.0) / f_n;
    
    for (int i = 0; i < curr.size(); ++i)
        curr[i][curr.size()] = 200 + 50 * (i + 1);;
}

int main()
{
    char type_of_input;
    std::cout << "Enter preferable type of input: file or terminal (F/T)\n";
    std::cin >> type_of_input;
    if (type_of_input == 'F')
        freopen("input.txt", "r", stdin);

    char type_of_matrix;
    std::cout << "Select type of matrix: functional or custom(F/C)";
    std::cin >> type_of_matrix;

    std::vector< std::vector<double> > curr;
    if (type_of_matrix == 'F')
    {
        curr.resize(f_n, std::vector<double> (f_n + 1));
        cast_matrix(curr);
    }
    else
    {
        std::cout << "Enter the size of matrix\n";
        int n;
        std::cin >> n;
        curr.resize(n, std::vector<double> (n + 1));
        std::cout << "Enter the values\n";

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n + 1; ++j)
                std::cin >> curr[i][j];
    }
    
    matrix m(curr);

    std::cout << "Choose the method: simple gauss or with the main elements(S/M)";
    char type_of_method;
    std::cin >> type_of_method;
    if (type_of_method == 'S')
        m.cast_gauss_equations();
    else
        m.cast_gauss_with_main_el();
    
    m.cast_norm();
    m.print();

    return 0;
}