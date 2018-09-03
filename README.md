# Hungarian Algorithm

Here I have posted my implementation of the Hungarian algorithm, a powerful
solution to the linear sum assignment problem.  The goal is to find the optimal
pairing of rows to columns in a matrix in order to minimize the sum cost.

I needed the algorithm for a larger project involving convolutional neural
networks, the bulk of which cannot be shared due to an NDA.


To utilize the algorithm, declare an instance of the wrapper class and call
the Solve method while passing in a vector of ints and the 2D cost matrix.
The minimum cost will be returned in a double, and the vector parameter 
will be modified with the optimal assignments, with the index corresponding 
to each row and the value matching the cost matrix's column value.
