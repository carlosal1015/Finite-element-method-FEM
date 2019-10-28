#include<vector>

class Matrix
{
  public:
    void resize(int numRows_, int numCols_);
    void resize(int numRows_, int numCols_, double value);
    // access elements
    double& operator()(int i, int j);
    double  operator()(int i, int j) const;
    std::vector<double>& operator[](int i);
    const std::vector<double>& operator[](int i) const;
    // arithmetic functions
    Matrix& operator*=(double x);
    Matrix& operator+=(const Matrix& b);
    std::vector<double> solve(std::vector<double> b) const;
    // output
    void print() const;
    int rows() const
    {
        return numRows;
    }
    int cols() const
    {
        return numCols;
    }
    
    Matrix(int numRows_, int numCols_) :
            entries(numRows_), numRows(numRows), numCols(numCols_)
    {
        for (int i = 0; i < numRows_; ++i)
            entries[i].resize(numCols_);
    };
    
    Matrix(int dim) : Matrix(dim,dim)
    {};
    
    Matrix(int numRows_, int numCols_, double value)
    {
        resize(value,numRows,numCols);
    };
    
    Matrix(std::vector<std::vector<double> > a)
    {
        entries = a;
        numRows = a.size();
        if (numRows > 0)
            numCols = a[0].size();
        else
            numCols = 0;
    }
    
    Matrix(const Matrix& b)
    {
        numRows = b.numRows;
        numCols = b.numCols;
    }
    
  private:
    std::vector<std::vector<double> > entries;
    int numRows = 0;
    int numCols = 0;
};

std::vector<double> operator*(const Matrix& a,
                              const std::vector<double>& x);
Matrix operator*(const Matrix& a, double x);
Matrix operator*(double x, const Matrix& a);
Matrix operator+(const Matrix& a, const Matrix& b);
