#include <iostream>
#include "matrix.h"

void matrix::print()
{
    std::cout << "Solution:\n";
    for(int i = 0; i < res.size(); ++i)
        std::cout << res[i] << ' ';
    std::cout << '\n';

    std::cout << "Determinant: " << det << '\n';
    std::cout << "Reverse matrix :\n";
    for (int i = 0; i < reverse_cell.size(); ++i)
    {
        for (int j = 0; j < reverse_cell.size(); ++j)
            std::cout << reverse_cell[i][j] << ' ';
        std::cout << '\n';
    }

    std::cout << "Conditionality: " << cond << '\n';
}

void matrix::set_size(int n)
{
    cell.resize(n, std::vector<double> (n));
}

void matrix::string_exchange(int k, int l, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr[k].size(); ++i)
        std::swap(curr[k][i], curr[l][i]);
}

void matrix::column_exchange(int k, int l, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr.size(); ++i)
        std::swap(curr[i][k], curr[i][l]);
}

void matrix::string_multiplex_on_number(int k, double x, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr[k].size(); ++i)
        curr[k][i] *= x;
}

void matrix::column_multiplex_on_number(int k, double x, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr.size(); ++i)
        curr[i][k] *= x;
}

void matrix::string_devide_on_number(int k, double x, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr[k].size(); ++i)
        curr[k][i] /= x;
}

void matrix::column_devide_on_number(int k, double x, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr.size(); ++i)
        curr[i][k] /= x;
}

void matrix::string_add_to(int k, int l, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr[k].size(); ++i)
        curr[l][i] += curr[k][i];
}

void matrix::column_add_to(int k, int l, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr.size(); ++i)
        curr[i][l] += curr[i][k];
}

void matrix::string_subtract_from(int k, int l, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr[k].size(); ++i)
        curr[k][i] -= curr[l][i];
}
void matrix::column_subtract_from(int k, int l, std::vector< std::vector<double> > &curr)
{
    for (int i = 0; i < curr.size(); ++i)
        curr[i][k] -= curr[i][l];
}

void matrix::cast_gauss_equations()
{
    int its = 0;
    while (its < cell.size())
    {
        int it = its;
        while(it < cell.size() && cell[it][its] == 0)
            ++it;
        
        if (it != its)
        {
            det *= -1;
            string_exchange(it, its, cell);   
            string_exchange(it, its, reverse_cell);
            string_exchange(it, its, det_cell);
        }
        
        string_devide_on_number(its, cell[its][its], reverse_cell);
        string_devide_on_number(its, cell[its][its], cell);
        
        for (int i = 0; i < cell.size(); ++i)
        {
            if (i == its || cell[i][its] == 0)
                continue;

            double x = cell[i][its];
            string_multiplex_on_number(its, x, cell);
            string_subtract_from(i, its, cell);
            string_devide_on_number(its, x, cell);

            string_multiplex_on_number(its, x, reverse_cell);
            string_subtract_from(i, its, reverse_cell);
            string_devide_on_number(its, x, reverse_cell);
            
            x = det_cell[i][its];
            int backup = det_cell[its][its];
            string_devide_on_number(its, backup, det_cell);
            string_multiplex_on_number(its, x, det_cell);
            string_subtract_from(i, its, det_cell);
            string_devide_on_number(its, x, det_cell);
            string_multiplex_on_number(its, backup, det_cell);
            det *= det_cell[its][its];
        }
        ++its;
    }

    for (int i = 0; i < cell.size(); ++i)
        res.push_back(cell[i][cell.size()]);
}


void matrix::cast_gauss_with_main_el()
{
    double tmp;
    
    for (int i = 0; i < cell.size(); ++i)
    {
        double max = abs(cell[i][i]);
        int max_i = i;
        for (int j = i; j < cell.size(); ++j)
            if(max < abs(cell[j][i]))
            {
                max = abs(cell[j][i]);
                max_i = j;
            }
        
        if(cell[max_i][i] == 0)
        {
            std::cout << "Wrong input. Matrix has no solutions";
            return;
        }
        
        if(i != max_i)
        {
            string_exchange(i, max_i, cell);
            string_exchange(i, max_i, reverse_cell);
            det *= -1;
        }
        tmp = cell[i][i];
        det *= tmp;
        cell[i][cell.size()] /= tmp;
        for (int j = cell.size() - 1; j >= 0; --j)
        {
            cell[i][j] /= tmp;
            reverse_cell[i][j] /= tmp;
        }
        for (int j = i + 1; j < cell.size(); ++j)
        {
            tmp = cell[j][i];
            cell[j][cell.size()] -= tmp * cell[i][cell.size()];
            for (int k = cell.size() - 1; k >= 0; --k)
            {
                cell[j][k] -= tmp * cell[i][k];
                reverse_cell[j][k] -= tmp * reverse_cell[i][k];
            }
        }
    }
    
    res.resize(cell.size());
    res[cell.size() - 1] = cell[cell.size() - 1][cell.size()];
    for (int i = cell.size() - 2; i >= 0; --i)
    {
        res[i] = cell[i][cell.size()];
        for (int j = i + 1; j < cell.size(); ++j)
        {
            res[i] -= cell[i][j] * res[j];
            for (int k = 0; k < cell.size(); ++k)
                reverse_cell[i][k] -= cell[i][j] * reverse_cell[j][k];
        }
    }
}

void matrix::cast_norm()
{
    double curr_nrm;
    for (int i = 0; i < source_cell.size(); ++i)
    {
        curr_nrm = 0;
        for (int j = 0; j < source_cell.size() - 1; ++j)
            curr_nrm += abs(source_cell[i][j]);
    
        if (curr_nrm > source_nrm)
            source_nrm = curr_nrm;
    }
    for (int i = 0; i < reverse_cell.size(); ++i)
    {
        curr_nrm = 0;
        for (int j = 0; j < reverse_cell.size(); ++j)
            curr_nrm += abs(reverse_cell[i][j]);
        
        if (curr_nrm > reverse_nrm)
            reverse_nrm = curr_nrm;
    }
    cond = source_nrm * reverse_nrm;
}