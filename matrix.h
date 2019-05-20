#include <vector>

class matrix
{
    long double det = 1;
    long double source_nrm = 0;
    long double reverse_nrm = 0;
    long double cond;
    std::vector< std::vector<double> > cell;
    std::vector< std::vector<double> > source_cell;
    std::vector< std::vector<double> > reverse_cell;
    std::vector< std::vector<double> > det_cell;
    std::vector<double> res;

    public:
        void set_size(int);
        void print();
        void cast_gauss_equations();
        void cast_gauss_with_main_el();
        void cast_norm();

        matrix(std::vector< std::vector<double> > _cell)
        {
            source_cell = _cell;
            cell = _cell;
            det_cell = _cell;
            reverse_cell.resize(_cell.size(), std::vector<double> (_cell.size(), 0));
            for (int i = 0; i < reverse_cell.size(); ++i)
                reverse_cell[i][i] = 1;
        };

    private:
        void string_multiplex_on_number(int, double, std::vector< std::vector<double> > &);
        void column_multiplex_on_number(int, double, std::vector< std::vector<double> > &);
        void string_devide_on_number(int, double, std::vector< std::vector<double> > &);
        void column_devide_on_number(int, double, std::vector< std::vector<double> > &);
        void string_add_to(int, int, std::vector< std::vector<double> > &);
        void column_add_to(int, int, std::vector< std::vector<double> > &);
        void string_subtract_from(int, int, std::vector< std::vector<double> > &);
        void column_subtract_from(int, int, std::vector< std::vector<double> > &);
        void string_exchange(int, int, std::vector< std::vector<double> > &);
        void column_exchange(int, int, std::vector< std::vector<double> > &);
};