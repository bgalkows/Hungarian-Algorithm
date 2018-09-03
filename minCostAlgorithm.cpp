/*

    minCostAssignment.cpp - this is an implementation of the Hungarian algorithm,
    also known as the Kuhn-Munkres algorithm.  The goal is to find a column for each
    row in a given matrix such that the lowest possible sum is attained.


    Code is based off the Hungarian algorithm implementation by Markus Buehren and Cong Ma.

*/

#include <stdlib.h>
#include <cfloat> // for DBL_MAX
#include <cmath>  // for fabs()

#include "minCostAlgorithm.h"


//*****************************************************************************//
// A single function used to tie all others together and compute the minimum sum.
//*****************************************************************************//
double minCostAlgorithm::Solve(std::vector <std::vector<double> >& costMatrix, std::vector<int>& Assignment)
{
	unsigned int nRows = costMatrix.size();
	unsigned int nCols = costMatrix[0].size();

	double *costMatrixIn = new double[nRows * nCols];
	int *assignment = new int[nRows];
	double cost = 0.0;

	// Populate costMatrixIn with the contents of costMatrix.
    // The goal is to store costMatrix in a one-dimensional form.
	for (unsigned int i = 0; i < nRows; ++i)
		for (unsigned int j = 0; j < nCols; ++j)
			costMatrixIn[i + nRows * j] = costMatrix[i][j];
	
	// Call actual assignment solving function
	optimalAssignment(assignment, &cost, costMatrixIn, nRows, nCols);

	Assignment.clear();
	for (unsigned int r = 0; r < nRows; r++)
		Assignment.push_back(assignment[r]);

	delete[] costMatrixIn;
	delete[] assignment;
	return cost;
}


// Solve optimal solution using Hungarian algorithm.
void minCostAlgorithm::optimalAssignment(int *assignment, double *cost, double *costMatrixIn, int nOfRows, int nOfColumns)
{
	double *costMatrix, *costMatrixTemp, *costMatrixEnd, *columnEnd, value, minValue;
	bool *coveredColumns, *coveredRows, *starMatrix, *newStarMatrix, *primeMatrix;
	int nOfElements, minDim, row, col;

	// Initialize all of assignment to -1; starting state
	*cost = 0;
	for (row = 0; row<nOfRows; row++)
		assignment[row] = -1;

	// Construct a copy of cost matrix
    // Confirm that no elements are negative
	nOfElements = nOfRows * nOfColumns;
	costMatrix = (double *)malloc(nOfElements * sizeof(double));
	costMatrixEnd = costMatrix + nOfElements;

	for (row = 0; row<nOfElements; row++)
	{
		value = costMatrixIn[row];
		if (value < 0)
			std::cerr << "All matrix elements have to be non-negative." << std::endl;
		costMatrix[row] = value;
	}


	// Allocate memory for boolean arrays
	coveredColumns = (bool *)calloc(nOfColumns, sizeof(bool));
	coveredRows = (bool *)calloc(nOfRows, sizeof(bool));
	starMatrix = (bool *)calloc(nOfElements, sizeof(bool));
	primeMatrix = (bool *)calloc(nOfElements, sizeof(bool));
	newStarMatrix = (bool *)calloc(nOfElements, sizeof(bool)); // used in step4 

	// Preliminary setup
	if (nOfRows <= nOfColumns)
	{
		minDim = nOfRows;

		for (row = 0; row<nOfRows; row++)
		{
			// Locate minimum element in each row
			costMatrixTemp = costMatrix + row;
			minValue = *costMatrixTemp;
			costMatrixTemp += nOfRows;
			while (costMatrixTemp < costMatrixEnd)
			{
				value = *costMatrixTemp;
				if (value < minValue)
					minValue = value;
				costMatrixTemp += nOfRows;
			}

			// Subtract minimum element from each element in the row
			costMatrixTemp = costMatrix + row;
			while (costMatrixTemp < costMatrixEnd)
			{
				*costMatrixTemp -= minValue;
				costMatrixTemp += nOfRows;
			}
		}

		// For steps 1 and 2a
		for (row = 0; row<nOfRows; row++)
			for (col = 0; col<nOfColumns; col++)
				if (fabs(costMatrix[row + nOfRows*col]) < DBL_EPSILON)
					if (!coveredColumns[col])
					{
						starMatrix[row + nOfRows*col] = true;
						coveredColumns[col] = true;
						break;
					}
	}
	else //  Read: number of rows > number of columns
	{
		minDim = nOfColumns;

		for (col = 0; col<nOfColumns; col++)
		{
			// find the smallest element in the column 
			costMatrixTemp = costMatrix + nOfRows*col;
			columnEnd = costMatrixTemp + nOfRows;

			minValue = *costMatrixTemp++;
			while (costMatrixTemp < columnEnd)
			{
				value = *costMatrixTemp++;
				if (value < minValue)
					minValue = value;
			}

			// subtract the smallest element from each element of the column 
			costMatrixTemp = costMatrix + nOfRows*col;
			while (costMatrixTemp < columnEnd)
				*costMatrixTemp++ -= minValue;
		}

		// for steps 1 and 2a
		for (col = 0; col<nOfColumns; col++)
			for (row = 0; row<nOfRows; row++)
				if (fabs(costMatrix[row + nOfRows*col]) < DBL_EPSILON)
					if (!coveredRows[row])
					{
						starMatrix[row + nOfRows*col] = true;
						coveredColumns[col] = true;
						coveredRows[row] = true;
						break;
					}
		for (row = 0; row<nOfRows; row++)
			coveredRows[row] = false;

	}

	// Transition to step 2b
	step2b(assignment, costMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

	// Compute cost, remove invalid assignments
	computeAssignmentCost(assignment, cost, costMatrixIn, nOfRows);

	// Free memory
	free(costMatrix);
	free(coveredColumns);
	free(coveredRows);
	free(starMatrix);
	free(primeMatrix);
	free(newStarMatrix);

	return;
}


/********************************************************/
void minCostAlgorithm::buildAssignmentVector(int *assignment, bool *starMatrix, int nOfRows, int nOfColumns)
{
	int row, col;

	for (row = 0; row < nOfRows; ++row)
		for (col = 0; col < nOfColumns; ++col)
			if (starMatrix[row + nOfRows*col])
			{
#ifdef ONE_INDEXING
				assignment[row] = col + 1; // MATLAB-Indexing 
#else
				assignment[row] = col;
#endif
				break;
			}
}


/********************************************************/
void minCostAlgorithm::computeAssignmentCost(int *assignment, double *cost, double *costMatrix, int nOfRows)
{
	int row, col;

	for (row = 0; row<nOfRows; ++row)
	{
		col = assignment[row];
		if (col >= 0)
			*cost += costMatrix[row + nOfRows*col];
	}
}


/********************************************************/
void minCostAlgorithm::step2a(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool *starMatrixTemp, *columnEnd;
	int col;

	// Cover every column containing a starred zero 
	for (col = 0; col < nOfColumns; ++col)
	{
		starMatrixTemp = starMatrix + nOfRows*col;
		columnEnd = starMatrixTemp + nOfRows;
		while (starMatrixTemp < columnEnd){
			if (*starMatrixTemp++)
			{
				coveredColumns[col] = true;
				break;
			}
		}
	}

    // move to step 2b
    step2b(assignment, costMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

}


/*******************************************************/
void minCostAlgorithm::step2b(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{

	// Count number of covered columns
	int col, nOfCoveredColumns = 0;
	for (col = 0; col < nOfColumns; ++col)
		if (coveredColumns[col])
			nOfCoveredColumns++;

	if (nOfCoveredColumns == minDim)
	{
		// Finished, no more steps needed
		buildAssignmentVector(assignment, starMatrix, nOfRows, nOfColumns);
	}
	else
	{
		// move to step 3
		step3(assignment, costMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}

}




/********************************************************/
void minCostAlgorithm::step3(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool zerosFound;
	int row, col, starCol;

	zerosFound = true;
	while (zerosFound)
	{
		zerosFound = false;
		for (col = 0; col < nOfColumns; ++col)
			if (!coveredColumns[col])
				for (row = 0; row < nOfRows; ++row)
					if ((!coveredRows[row]) && (fabs(costMatrix[row + nOfRows*col]) < DBL_EPSILON))
					{
						// Prime zero
						primeMatrix[row + nOfRows*col] = true;

						// Find starred zero in current row
						for (starCol = 0; starCol<nOfColumns; ++starCol)
							if (starMatrix[row + nOfRows*starCol])
								break;

						if (starCol == nOfColumns) // No starred zero detected
						{
							// move to step 4
							step4(assignment, costMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim, row, col);
							return;
						}
						else
						{
							coveredRows[row] = true;
							coveredColumns[starCol] = false;
							zerosFound = true;
							break;
						}
					}
	}

	// move to step 5
	step5(assignment, costMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}


/********************************************************/
void minCostAlgorithm::step4(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col)
{
	int n, starRow, starCol, primeRow, primeCol;
	int nOfElements = nOfRows*nOfColumns;

	// Construct temporary copy of starMatrix
	for (n = 0; n < nOfElements; ++n)
		newStarMatrix[n] = starMatrix[n];

	// Star current zero
	newStarMatrix[row + nOfRows*col] = true;

	// Find starred zero in current column
	starCol = col;
	for (starRow = 0; starRow < nOfRows; starRow++)
		if (starMatrix[starRow + nOfRows*starCol])
			break;

	while (starRow<nOfRows)
	{
		// Unstar the starred zero
		newStarMatrix[starRow + nOfRows*starCol] = false;

		// Find primed zero in current row 
		primeRow = starRow;
		for (primeCol = 0; primeCol < nOfColumns; ++primeCol)
			if (primeMatrix[primeRow + nOfRows*primeCol])
				break;

		// Star the primed zero
		newStarMatrix[primeRow + nOfRows*primeCol] = true;

		// Find starred zero in current column
		starCol = primeCol;
		for (starRow = 0; starRow < nOfRows; ++starRow)
			if (starMatrix[starRow + nOfRows*starCol])
				break;
	}

	// Use temporary copy as the new starMatrix 
	// Delete all primes, uncover all rows 
	for (n = 0; n < nOfElements; ++n)
	{
		primeMatrix[n] = false;
		starMatrix[n] = newStarMatrix[n];
	}
	for (n = 0; n < nOfRows; ++n)
		coveredRows[n] = false;

	// move to step 2a 
	step2a(assignment, costMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}


/********************************************************/
void minCostAlgorithm::step5(int *assignment, double *costMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	double h, value;
	int row, col;

	// Find the lowest uncovered element - h
	h = DBL_MAX;
	for (row = 0; row < nOfRows; ++row)
		if (!coveredRows[row])
			for (col = 0; col < nOfColumns; ++col)
				if (!coveredColumns[col])
				{
					value = costMatrix[row + nOfRows*col];
					if (value < h)
						h = value;
				}

	// Add h to each covered row
	for (row = 0; row < nOfRows; ++row)
		if (coveredRows[row])
			for (col = 0; col < nOfColumns; ++col)
				costMatrix[row + nOfRows*col] += h;

	// Subtract h from each uncovered column 
	for (col = 0; col < nOfColumns; ++col)
		if (!coveredColumns[col])
			for (row = 0; row < nOfRows; ++row)
				costMatrix[row + nOfRows*col] -= h;

	// move to step 3 
	step3(assignment, costMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

