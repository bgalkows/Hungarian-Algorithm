#include "minCostAlgorithm.cpp"


int main()
{
    std::vector< std::vector<double> > costMatrix = { {4, 1, 3}, {2, 0, 5}, {3, 2, 2} };

    minCostAlgorithm MCA;
    std::vector<int> assignment;

    double cost = MCA.Solve(costMatrix, assignment);

    for(unsigned int x = 0; x < costMatrix.size(); ++x)
    {
        std::cout << assignment[x] << '\t';
    }
    std::cout << "\ncost: " << cost << std::endl;

    return 0;
}


//      The output here should be:  1  0  2,   with a total cost of 5.