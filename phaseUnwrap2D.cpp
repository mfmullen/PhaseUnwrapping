/**************************************************************************
INPUT: A 2D wrapped phase image
Example: unwrapped = phaseUnwrapC(wrapped);

 Michael Mullen
 mulle399@umn.edu
 Created in Matlab R2016a
 version 3.0 (Feb 2020)
 uses 2D-SRNCP
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
#define M_PI 4.0*atan(1.0)

using namespace std;

//Function to unwrap two adjacent pixels
double pairwiseUnwrap(double phase1,double phase2){
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
    char edgeType;      //horizontal or vertical edge ('r' or 'c', respectively) 
    double R;           //reliability
    long int row;            //indices of pixel
    long int column;

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

//makes code a little clearer when grabbing indices.
inline long int get2Dindex(const long int row, const long int col, const long int numCols) {
    return col + (row)*numCols;
}

//calculate reliability of each edge.
void calcEdges(vector <pixelGroup>& groupArray, vector<edgeInfo>& edges, const long int numRows, 
    const long int numCols, const vector<double>& Reliability){
    long int colIndex, rowIndex, count = 0;
    double rowEdge, columnEdge;

    for (long int index = 0; index < numRows * numCols; index++) {
        colIndex = index % numCols;
        rowIndex = (index - colIndex) / numCols;

        if (rowIndex < numRows - 1) {
            rowEdge = Reliability[get2Dindex(rowIndex, colIndex, numCols)]
                + Reliability[get2Dindex(rowIndex+1, colIndex, numCols)];
            edges[count].edgeType = 'r';
            edges[count].R = rowEdge;
            edges[count].row = rowIndex;
            edges[count].column = colIndex;

            ++count;
        }

        if (colIndex < numCols - 1) {
            columnEdge = Reliability[get2Dindex(rowIndex, colIndex, numCols)]
                + Reliability[get2Dindex(rowIndex, colIndex + 1, numCols)];
            edges[count].edgeType = 'c';
            edges[count].R = columnEdge;
            edges[count].row = rowIndex;
            edges[count].column = colIndex;

            ++count;
        }
        groupArray[(rowIndex)*numCols + colIndex].group = get2Dindex(rowIndex, colIndex, numCols);
        groupArray[(rowIndex)*numCols + colIndex].nextIndex = -1;
        groupArray[(rowIndex)*numCols + colIndex].lastIndex = get2Dindex(rowIndex, colIndex, numCols);
        groupArray[(rowIndex)*numCols + colIndex].firstIndex = get2Dindex(rowIndex, colIndex, numCols);
    }
}

//Reliability calculation
//Calculated finite differences in phase across all non-border pixels
//Includes diagonal phase differences
void calcReliability(vector<double>& Reliability, const double* unwrappedImage, const int& numRows, const int& numCols) {
    double H, V, D1, D2, D;
    int colIndex, rowIndex;

    for (int index = 0; index < numRows * numCols; index++) {
        //protect against out of range errors hopefully. Had some past issues in that regard.
        try {
            colIndex = index % numCols;
            rowIndex = (index - colIndex) / numCols;
            if (rowIndex != 0 && rowIndex != numRows - 1 && colIndex != 0 && colIndex != numCols - 1) {
                H = pairwiseUnwrap(unwrappedImage[get2Dindex(rowIndex-1, colIndex, numCols)],
                    unwrappedImage[get2Dindex(rowIndex, colIndex, numCols)]) -
                    pairwiseUnwrap(unwrappedImage[get2Dindex(rowIndex, colIndex, numCols)],
                        unwrappedImage[get2Dindex(rowIndex+1, colIndex, numCols)]);

                V = pairwiseUnwrap(unwrappedImage[(rowIndex)*numCols + colIndex - 1],
                    unwrappedImage[(rowIndex)*numCols + colIndex]) -
                    pairwiseUnwrap(unwrappedImage[get2Dindex(rowIndex, colIndex, numCols)],
                        unwrappedImage[get2Dindex(rowIndex, colIndex+1, numCols)]);

                D1 = pairwiseUnwrap(unwrappedImage[get2Dindex(rowIndex-1, colIndex-1, numCols)],
                    unwrappedImage[get2Dindex(rowIndex, colIndex, numCols)]) -
                    pairwiseUnwrap(unwrappedImage[get2Dindex(rowIndex, colIndex, numCols)],
                        unwrappedImage[get2Dindex(rowIndex+1, colIndex+1, numCols)]);

                D2 = pairwiseUnwrap(unwrappedImage[get2Dindex(rowIndex-1, colIndex+1, numCols)],
                    unwrappedImage[get2Dindex(rowIndex, colIndex, numCols)]) -
                    pairwiseUnwrap(unwrappedImage[get2Dindex(rowIndex, colIndex, numCols)],
                        unwrappedImage[get2Dindex(rowIndex+1, colIndex-1, numCols)]);

                D = sqrt(H * H + V * V + D1 * D1 + D2 * D2);
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
    vector <pixelGroup>& groupArray, const long int numCols, const long int numRows) {

    long int adjGroup, adjIndex, currentCol, currentRow, currentGroup;
    long int iters, numberOfJumps, jumpDirection, currentFirst, currentLast, nextFirst, nextLast;
    long int nextPixel;

    char edge;
    double phase1, phase2;
    
    for (long int sortedIndex = 0; sortedIndex < numGroups; sortedIndex++) {

        //Retrieve parameters about the current edge type and its location
        currentRow = edges[sortedIndex].row;
        currentCol = edges[sortedIndex].column;
        currentGroup = groupArray[get2Dindex(currentRow, currentCol, numCols)].group;
        edge = edges[sortedIndex].edgeType;

        //Calculate index of adjacent pixel
        if (edge == 'r') {
            adjIndex = get2Dindex(currentRow+1, currentCol, numCols);
        }
        else {
            adjIndex = get2Dindex(currentRow, currentCol+1, numCols);
        }

        //group of adjacent pixel
        adjGroup = groupArray[adjIndex].group;
        //current phase of adjacent pixel
        phase1 = unwrappedImage[adjIndex];
        //current phase of current pixel
        phase2 = unwrappedImage[get2Dindex(currentRow, currentCol, numCols)];

        //its is used in a while loop for phase unwrapping
        iters = 0;
        //m is the number of 2pi jumps
        numberOfJumps = 0;
        //jumpDirection indicates the direction of the 2 pi jumps (+1 or -1)
        jumpDirection = 0;
        //Only need to unwrap if the two different pixels are in different groups
        if (adjGroup != currentGroup) {

            //First condition is definition of a phase discontinuity and second condition just makes sure we don't get stuck in an infinite loop
            while (fabs(phase1 - phase2) > M_PI && iters < 100) {

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
                    break;
                }

                //Increase iteration count
                ++iters;

            }

            //Merge the two groups by having last index of first group go to first index of second group
            currentFirst = groupArray[get2Dindex(currentRow, currentCol, numCols)].firstIndex;
            currentLast = groupArray[get2Dindex(currentRow, currentCol, numCols)].lastIndex;


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


//========================================================================
//                         UNWRAPPER                                    
//========================================================================
void unwrap(const long int numRows, const long int numCols, const double* inputArray,
        double* unwrappedImage)
{
    long int numGroups = 2*numRows*numCols - numRows - numCols; //Number of edges
    vector<edgeInfo> edges(numGroups);
    vector<double> Reliability(numRows * numCols);
    vector <pixelGroup> groupArray(numRows * numCols);
    
    long int size = numRows * numCols * sizeof(double);

    memcpy(unwrappedImage,inputArray,size);
    
    //Defining reliability, using convention of paper mentioned in header
    //replce with call to function that returns reliability
    calcReliability(Reliability,unwrappedImage, numRows, numCols);
    
    //Define reliability of edges and store in struct edges
    //Do in a separate function
    calcEdges(groupArray, edges, numRows, numCols, Reliability);
   
    //====================== SORT EDGES HERE ==================
    std::sort(edges.begin(), edges.end());
    //========================================================

    //Do the unwrapping by looping over the groups
    //Export these tasks to a separate function
    calcUnwrap(unwrappedImage, numGroups, edges, groupArray, numCols, numRows);
    
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

    long int numRows;
    long int numCols;
    
    inputArray = mxGetPr(prhs[0]);
    numCols = (long int)mxGetM(prhs[0]);
    numRows = (long int)mxGetN(prhs[0]);
    
    // Allocate output array, note transposed size. 
    // since Matlab uses column major, C uses row major.
    
    plhs[0] = mxCreateDoubleMatrix(numCols,numRows,mxREAL);
    unwrappedImage = mxGetPr(plhs[0]);
    
    unwrap(numRows,numCols,inputArray,unwrappedImage);
    
}
