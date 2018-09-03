#ifndef MIN_COST
#define MIN_COST

#include <iostream>
#include <vector>

class minCostAlgorithm
{
public:
	minCostAlgorithm() = default;
	~minCostAlgorithm() = default;
	double Solve(std::vector< std::vector<double> >& costMatrix, std::vector<int>& Assignment);

private:
	void optimalAssignment(int *assignment, double *cost, double *costMatrix, int nOfRows, int nOfColumns);
	void buildAssignmentVector(int *assignment, bool *starMatrix, int nOfRows, int nOfColumns);
	void computeAssignmentCost(int *assignment, double *cost, double *costMatrix, int nOfRows);
	void step2a(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step2b(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step3(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step4(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col);
	void step5(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
};


#endif
