/**************************************************************************
INPUT: A 2D/3D wrapped phase image
Example: unwrapped = phaseUnwrapC(wrapped);

 Michael Mullen
 mulle399@umn.edu
 Created in Matlab R2016a
 version 3.0 (Feb 2020)
 uses 3D-SRNCP version of 2D-SRNCP
 an algorithm by Miguel Arevallilo Herraez, David R. Burton, 
 Michael J. Lalor, and Munther A. Gdeisat in Applied Optics, Vol. 
 41, No. 35, pp. 7437, 2002.
/*************************************************************************/

#include "mex.h"
#include <math.h>
#include "matrix.h"
#include <vector>
#include <algorithm>
#include <string>
#include <map> 
#include <iterator> 
#include <stdexcept>
#include <iostream>

//======================================================================
//                          BASIC FUNCTIONS HERE                        
//======================================================================
inline double M_PI = 4.0 * atan(1.0);

using namespace std;

//Function to unwrap two adjacent pixels
double pairwiseUnwrap(double phase1, double phase2){
    while(fabs(phase1-phase2) > M_PI){
        
        if(phase1 > phase2 + M_PI){
            phase1 = phase1 - 2*M_PI;
        }
        else if (phase1 < phase2 - M_PI){
            phase1 = phase1 + 2*M_PI;
        }
        else{
            break;
        }
        
    }
    double result = phase1 - phase2;
    return result;
};

//Information on edges
struct edgeInfo {
    char edgeType = 'n';      //horizontal, vertical, or page edge ('r', 'c', or 'p', respectively) 
    double R;           //reliability
    long int row;            //indices of pixel
    long int column;
    long int page;

    // for sorting, intentionally overloaded < operator with > operator to get decreasing sort.
    bool operator < (const edgeInfo& val) const
    {
        return (R > val.R);
    }
};

//Keeps track of which phase group a given pixel belongs to.
// Can I replace these with multimaps?
struct pixelGroup {
    long int group;
    long int nextIndex;
    long int firstIndex;
    long int lastIndex;
};

//makes code a little clearer when grabbing indices of adjacent voxels.
inline long int get3Dindex(const long int row, const long int col, const long int page, 
    const long int numCols, const long int numRows) {
    return col + (row)*numCols + page*numCols*numRows;
}

//calculate reliability of each edge.
void calcEdges(vector <pixelGroup>& groupArray, vector<edgeInfo>& edges, const long int numRows, 
    const long int numCols, const long int numPages, const vector<double>& Reliability){
    long int colIndex, rowIndex, pageIndex,count = 0;
    double rowEdge, columnEdge, pageEdge;

    for (long int index = 0; index < numRows * numCols * numPages; index++) {
        colIndex = index % numCols;
        rowIndex = ((index % (numRows * numCols)) - colIndex) / numCols;
        pageIndex = (index - rowIndex * numCols - colIndex) / (numRows * numCols);

        if (rowIndex < numRows - 1) {
            rowEdge = Reliability[get3Dindex(rowIndex, colIndex, pageIndex,numCols,numRows)]
                + Reliability[get3Dindex(rowIndex+1, colIndex, pageIndex, numCols, numRows)];
            edges[count].edgeType = 'r';
            edges[count].R = rowEdge;
            edges[count].row = rowIndex;
            edges[count].column = colIndex;
            edges[count].page = pageIndex;

            ++count;
        }

        if (colIndex < numCols - 1) {
            columnEdge = Reliability[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]
                + Reliability[get3Dindex(rowIndex, colIndex+1, pageIndex, numCols, numRows)];
            edges[count].edgeType = 'c';
            edges[count].R = columnEdge;
            edges[count].row = rowIndex;
            edges[count].column = colIndex;
            edges[count].page = pageIndex;

            ++count;
        }

        if (pageIndex < numPages - 1) {

            pageEdge = Reliability[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]
                + Reliability[get3Dindex(rowIndex, colIndex, pageIndex+1, numCols, numRows)];
            edges[count].edgeType = 'p';
            edges[count].R = columnEdge;
            edges[count].row = rowIndex;
            edges[count].column = colIndex;
            edges[count].page = pageIndex;

            ++count;
        }

        groupArray[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)].group = get3Dindex(rowIndex, colIndex, pageIndex,numCols,numRows);
        groupArray[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)].nextIndex = -1;
        groupArray[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)].lastIndex = get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows);
        groupArray[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)].firstIndex = get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows);
    }

}

//Reliability calculation
//Calculates second differenc in phase across all non-border pixels
//Includes diagonal phase differences
void calcReliability(vector<double>& Reliability, const double* unwrappedImage, const int& numRows, 
    const int& numCols, const long int numPages) {
    //horizontal, vertical, in-plane diagonals, page, through-page diagonals.
    double H, V, D1, D2;
    double P = 0, PD1 = 0, PD2 = 0, PD3 = 0, PD4 = 0;
    //End result stored in D
    double D;
    int colIndex, rowIndex, pageIndex;

    for (int index = 0; index < numRows * numCols * numPages; index++) {
        //protect against out of range errors hopefully. Had some past issues in that regard.
        try {
            colIndex = index % numCols;
            rowIndex = ((index % (numRows * numCols)) - colIndex) / numCols;
            pageIndex = (index - rowIndex * numCols - colIndex) / (numRows * numCols);

            if (rowIndex > 0 && rowIndex < numRows - 1 && colIndex > 0 && colIndex < numCols - 1) {
                H = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex-1, colIndex, pageIndex, numCols, numRows)],
                    unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                    pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex+1, colIndex, pageIndex, numCols, numRows)]);

                V = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex-1, pageIndex, numCols, numRows)],
                    unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                    pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex, colIndex+1, pageIndex, numCols, numRows)]);

                D1 = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex-1, colIndex-1, pageIndex, numCols, numRows)],
                    unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                    pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex+1, colIndex+1, pageIndex, numCols, numRows)]);

                D2 = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex-1, colIndex+1, pageIndex, numCols, numRows)],
                    unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                    pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex+1, colIndex-1, pageIndex, numCols, numRows)]);
                
                if (pageIndex > 0 && pageIndex < numPages - 1) {
                    P = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex - 1, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                        pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                            unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex + 1, numCols, numRows)]);
                    
                    PD1 = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex - 1, colIndex - 1, pageIndex - 1, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                        pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                            unwrappedImage[get3Dindex(rowIndex + 1, colIndex + 1, pageIndex + 1, numCols, numRows)]);
                    
                    PD2 = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex + 1, colIndex - 1, pageIndex - 1, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                        pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                            unwrappedImage[get3Dindex(rowIndex - 1, colIndex + 1, pageIndex + 1, numCols, numRows)]);
                    
                    PD3 = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex - 1, colIndex + 1, pageIndex - 1, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                        pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                            unwrappedImage[get3Dindex(rowIndex + 1, colIndex - 1, pageIndex + 1, numCols, numRows)]);
                    
                    PD4 = pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex + 1, colIndex + 1, pageIndex - 1, numCols, numRows)],
                        unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)]) -
                        pairwiseUnwrap(unwrappedImage[get3Dindex(rowIndex, colIndex, pageIndex, numCols, numRows)],
                            unwrappedImage[get3Dindex(rowIndex - 1, colIndex - 1, pageIndex + 1, numCols, numRows)]);
                    
                }
                
                D = sqrt(H * H + V * V + D1 * D1 + D2 * D2 + P * P + PD1 * PD1 + PD2 * PD2 + PD3 * PD3 + PD4 * PD4);
                Reliability[index] = 1 / D;

            }
            else {
                Reliability[index] = 0.0;
            }
        }
        catch (exception& e) {
            cerr << e.what() << endl;
            return;
        }
    }
}

//Need to clean up this function. This is the main calculation though. 
//Compares groups and unwraps them relative to each other.
void calcUnwrap(double* unwrappedImage, const long int numGroups, vector<edgeInfo>& edges, 
    vector <pixelGroup>& groupArray, const long int numCols, const long int numRows, const long int numPages) {

    long int adjGroup, adjIndex, currentCol, currentRow, currentPage, currentGroup;
    long int iters, numberOfJumps, jumpDirection, currentFirst, currentLast, nextFirst, nextLast;
    long int nextPixel;

    char edge;
    double phase1, phase2;
    
    int maxIters = 100; //use a modest number of iterations

    for (long int sortedIndex = 0; sortedIndex < numGroups; sortedIndex++) {

        //Retrieve parameters about the current edge type and its location
        currentRow = edges[sortedIndex].row;
        currentCol = edges[sortedIndex].column;
        currentPage = edges[sortedIndex].page;
        //the weird condition on numPages only is so that it still enters the loop in 2D, 
        //when currentPage = 0 always
        if (currentRow < numRows - 1 && currentCol < numCols - 1 && (currentPage < numPages - 1 || currentPage == 0)) {
            currentGroup = groupArray[get3Dindex(currentRow, currentCol, currentPage, numCols, numRows)].group;
            edge = edges[sortedIndex].edgeType;

            //Calculate index of adjacent pixel
            switch (edge) {
                case 'r': //row edge
                    adjIndex = get3Dindex(currentRow + 1, currentCol, currentPage, numCols, numRows);
                    break;
                case 'c': //column edge
                    adjIndex = get3Dindex(currentRow, currentCol + 1, currentPage, numCols, numRows);
                    break;
                case 'p': //page edge
                    adjIndex = get3Dindex(currentRow, currentCol, currentPage + 1, numCols, numRows);
                    break;
                default: {
                    mexPrintf("edgetype = %c \n", edge);
                    mexErrMsgTxt("ERROR: non-existent edge type");
                }
            }

            //group of adjacent pixel
            adjGroup = groupArray[adjIndex].group;
            //current phase of adjacent pixel
            phase1 = unwrappedImage[adjIndex];
            //current phase of current pixel
            phase2 = unwrappedImage[get3Dindex(currentRow, currentCol, currentPage, numCols, numRows)];

            //its is used in a while loop for phase unwrapping
            iters = 0;
            //m is the number of 2pi jumps
            numberOfJumps = 0;
            //jumpDirection indicates the direction of the 2 pi jumps (+1 or -1)
            jumpDirection = 0;
            //Only need to unwrap if the two different pixels are in different groups
            if (adjGroup != currentGroup) {

                //First condition is definition of a phase discontinuity and second condition just makes sure we don't get stuck in an infinite loop
                while (fabs(phase1 - phase2) > M_PI && iters < maxIters) {

                    //if phase needs to decrement by 2pi
                    if (phase1 > phase2 + M_PI) {
                        phase1 = phase1 - 2 * M_PI;
                        jumpDirection = -1;
                        ++numberOfJumps; //increment jump counter
                    }
                    //if phase needs to increment by 2pi
                    else if (phase1 < phase2 - M_PI) {
                        phase1 = phase1 + 2 * M_PI;
                        jumpDirection = 1;
                        ++numberOfJumps;
                    }

                    //if phase is just right
                    else {
                        iters = maxIters + 1;
                        break;
                    }

                    //Increase iteration count
                    ++iters;

                }

                //Merge the two groups by having last index of first group go to first index of second group
                currentFirst = groupArray[get3Dindex(currentRow, currentCol, currentPage, numCols, numRows)].firstIndex;
                currentLast = groupArray[get3Dindex(currentRow, currentCol, currentPage, numCols, numRows)].lastIndex;


                nextLast = groupArray[adjIndex].lastIndex;
                nextFirst = groupArray[adjIndex].firstIndex;

                groupArray[currentLast].nextIndex = nextFirst;
                groupArray[currentFirst].lastIndex = nextLast;

                nextPixel = groupArray[currentFirst].nextIndex;
                //next = -1 indicates end of group
                //Loop over all of group and add correct phase to adjGroup
                while (nextPixel != -1) {

                    if (groupArray[nextPixel].group == adjGroup) {
                        unwrappedImage[nextPixel] = unwrappedImage[nextPixel] + jumpDirection * 2 * M_PI * numberOfJumps;
                    }

                    groupArray[nextPixel].firstIndex = currentFirst;
                    groupArray[nextPixel].group = currentGroup;
                    groupArray[nextPixel].lastIndex = nextLast;

                    nextPixel = groupArray[nextPixel].nextIndex;

                }

            }
        }
    }
}


//========================================================================
//                         UNWRAPPER                                    
//========================================================================
void unwrap(const long int numRows, const long int numCols, const long int numPages, const double* inputArray,
        double* unwrappedImage)
{
    long int numGroups = (numCols - 1) * numRows * numPages + 
                         (numRows - 1)* numCols * numPages + 
                         (numPages - 1) * numRows * numCols; //Number of edges

    vector<edgeInfo> edges(numGroups);
    vector<double> Reliability(numRows * numCols * numPages);
    vector <pixelGroup> groupArray(numRows * numCols * numPages);
        
    //Defining reliability, using convention of paper mentioned in header
    //replce with call to function that returns reliability
    calcReliability(Reliability,unwrappedImage, numRows, numCols, numPages);
    
    //Define reliability of edges and store in struct edges
    //Do in a separate function
    calcEdges(groupArray, edges, numRows, numCols, numPages, Reliability);
   
    //====================== SORT EDGES HERE ==================
    std::sort(edges.begin(), edges.end());
    //========================================================

    //Do the unwrapping by looping over the groups
    //Export these tasks to a separate function
    calcUnwrap(unwrappedImage, numGroups, edges, groupArray, numCols, numRows, numPages);
    
}


//==========================================================================
//                          MEX FUNCTION HERE                            
//==========================================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //Ignoring the mwSize and such, classic C++ types seem more stable. 
    if (nlhs != 1){
        mexErrMsgTxt("ERROR:  Number of outputs must equal 1 \n");
    }
    
    if (nrhs != 1){
        mexErrMsgTxt("ERROR:  Number of inputs must equal 1. \n");
    }
    if(mxIsComplex(prhs[0]) == 1 || mxIsDouble(prhs[0]) != 1){
        mexErrMsgTxt("ERROR:  Input must be a real double array \n");

    }
    double *inputArray;             //Input phase image
    double *unwrappedImage;         //output array, unwrapped

    long int numRows = 0;
    long int numCols = 0;
    long int numPages = 0;
    long int ndims = 0; //will throw an out_of_range error soon if incorrect size.

    inputArray = mxGetPr(prhs[0]);
    
    // Allocate output array, note transposed size. 
    // since Matlab uses column major, C uses row major.
    const size_t* dims = mxGetDimensions(prhs[0]);
    ndims = (long int)mxGetNumberOfDimensions(prhs[0]);

    if (ndims < 2 || ndims > 3) {
        mexErrMsgTxt("ERROR:  ndims must equal 2 or 3 \n");

    }

    numRows = (long int)dims[0];
    numCols = (long int)dims[1];

    if (ndims == 3L) {
        numPages = (long int)dims[2];
    }
    else {
        numPages = 1; //used in multiplications later to determine size.
    }

    //Setup output image size.
    plhs[0] = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
    unwrappedImage = mxGetPr(plhs[0]);
    
    //Initialize unwrappedImage with the input image
    long int outputSize;
    outputSize = numRows * numCols * numPages * sizeof(double);
    memcpy(unwrappedImage, inputArray, outputSize);

    unwrap(numRows,numCols,numPages,inputArray,unwrappedImage);
    
}
