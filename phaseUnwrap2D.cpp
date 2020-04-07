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

//Reliability calculation
//Calculated finite differences in phase across all non-border pixels
//Includes diagonal phase differences
void calcReliability(vector<double>& Reliability, const double* unwrappedImage, const int& numRows, const int& numCols) {
    double H, V, D1, D2, D;
    int colIndex, rowIndex;

    for (int index = 0; index < numRows * numCols; index++) {
        try {
            colIndex = index % numCols;
            rowIndex = (index - colIndex) / numCols;
            if (rowIndex != 0 && rowIndex != numRows - 1 && colIndex != 0 && colIndex != numCols - 1) {
                H = pairwiseUnwrap(unwrappedImage[(rowIndex - 1) * numCols + colIndex],
                    unwrappedImage[(rowIndex)*numCols + colIndex]) -
                    pairwiseUnwrap(unwrappedImage[(rowIndex)*numCols + colIndex],
                        unwrappedImage[(rowIndex + 1) * numCols + colIndex]);

                V = pairwiseUnwrap(unwrappedImage[(rowIndex)*numCols + colIndex + 1],
                    unwrappedImage[(rowIndex)*numCols + colIndex]) -
                    pairwiseUnwrap(unwrappedImage[(rowIndex)*numCols + colIndex],
                        unwrappedImage[(rowIndex)*numCols + colIndex + 1]);

                D1 = pairwiseUnwrap(unwrappedImage[(rowIndex - 1) * numCols + colIndex - 1],
                    unwrappedImage[(rowIndex)*numCols + colIndex]) -
                    pairwiseUnwrap(unwrappedImage[(rowIndex)*numCols + colIndex],
                        unwrappedImage[(rowIndex + 1) * numCols + colIndex + 1]);

                D2 = pairwiseUnwrap(unwrappedImage[(rowIndex - 1) * numCols + colIndex + 1],
                    unwrappedImage[(rowIndex)*numCols + colIndex]) -
                    pairwiseUnwrap(unwrappedImage[(rowIndex)*numCols + colIndex],
                        unwrappedImage[(rowIndex + 1) * numCols + colIndex - 1]);

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

    double columnEdge;
    double rowEdge;
    
    long int size = numRows * numCols * sizeof(double);
    long int currentRow,currentCol;
    long int group,adjIndex;
    long int currentGroup,adjGroup,next,nextLast,nextFirst,currentFirst,currentLast;
    char edge;
    long int rowIndex, colIndex;
    long int iters = 0;
    
    double phase1,phase2;
    
    long int jumpDirection = 0, numberOfJumps = 0;
    long int index,sortedIndex,count = 0;
    
    memcpy(unwrappedImage,inputArray,size);
    
    //Defining reliability, using convention of paper mentioned in header
    //replce with call to function that returns reliability
    calcReliability(Reliability,unwrappedImage, numRows, numCols);
    
    //Define reliability of edges and store in struct edges
    //Do in a separate function
    for (index = 0; index < numRows * numCols; index++){
        colIndex = index % numCols;
        rowIndex = (index-colIndex)/numCols;
        
        if (rowIndex < numRows - 1){
            rowEdge = Reliability[(rowIndex)*numCols + colIndex]
                    + Reliability[(rowIndex+1)*numCols + colIndex];
            edges[count].edgeType = 'r';
            edges[count].R = rowEdge;
            edges[count].row = rowIndex;
            edges[count].column = colIndex;
            
            ++count;
        }
        
        if (colIndex < numCols - 1){
            columnEdge = Reliability[(rowIndex)*numCols + colIndex]
                    + Reliability[(rowIndex)*numCols + colIndex+1];
            edges[count].edgeType = 'c';
            edges[count].R = columnEdge;
            edges[count].row = rowIndex;
            edges[count].column = colIndex;
            
            ++count;
        }
        groupArray[(rowIndex)*numCols + colIndex].group = (rowIndex)*numCols + colIndex;
        groupArray[(rowIndex)*numCols + colIndex].nextIndex = -1;
        groupArray[(rowIndex)*numCols + colIndex].lastIndex = (rowIndex)*numCols + colIndex;
        groupArray[(rowIndex)*numCols + colIndex].firstIndex = (rowIndex)*numCols + colIndex;


    }
    
    //multimap<double,double> myMap;
    
    //====================== SORT EDGES HERE ==================
    std::sort(edges.begin(), edges.end());
    //========================================================

     //Do the unwrapping by looping over the groups
    //Export these tasks to a separate function
    for (sortedIndex = 0; sortedIndex < numGroups; sortedIndex++){
        
        //Retrieve parameters about the current edge type and its location
        currentRow = edges[sortedIndex].row;
        currentCol = edges[sortedIndex].column;
        currentGroup = groupArray[(currentRow)*numCols + currentCol].group;
        edge = edges[sortedIndex].edgeType;
        
        //Calculate index of adjacent pixel
        if (edge == 'r'){
            adjIndex = (currentRow+1)*numCols + currentCol;
        }
        else{
            adjIndex = (currentRow)*numCols + currentCol + 1;
        }
        
        //group of adjacent pixel
        adjGroup = groupArray[adjIndex].group;
        //current phase of adjacent pixel
        phase1 = unwrappedImage[adjIndex];
        //current phase of current pixel
        phase2 = unwrappedImage[currentRow*numCols+currentCol];

        //its is used in a while loop for phase unwrapping
        iters = 0;
        //m is the number of 2pi jumps
        numberOfJumps = 0;
        //jumpDirection indicates the direction of the 2 pi jumps (+1 or -1)
        jumpDirection = 0;
        //Only need to unwrap if the two different pixels are in different groups
        if (adjGroup != currentGroup){
            
            //First condition is definition of a phase discontinuity and second condition just makes sure we don't get stuck in an infinite loop
            while (fabs(phase1 - phase2) > M_PI && iters < 100){   

                //if phase needs to decrement by 2pi
                if(phase1 > phase2+ M_PI){
                    phase1 = phase1 - 2*M_PI;
                    jumpDirection = -1;
                    ++numberOfJumps; //increment jump counter
                }
                //if phase needs to increment by 2pi
                else if(phase1 < phase2- M_PI){
                    phase1 = phase1 + 2*M_PI;
                    jumpDirection = 1;
                    ++numberOfJumps;
                }

                //if phase is just right
                else{
                    break;
                }

                //Increase iteration count
                ++iters;
                
            }
                        
            //Merge the two groups by having last index of first group go to first index of second group
            //Merge using multimaps?
            currentFirst = groupArray[(currentRow)*numCols + currentCol].firstIndex;
            currentLast = groupArray[(currentRow)*numCols + currentCol].lastIndex;


            nextLast = groupArray[adjIndex].lastIndex;
            nextFirst = groupArray[adjIndex].firstIndex;

            groupArray[currentLast].nextIndex = nextFirst;
            groupArray[currentFirst].lastIndex = nextLast;

            next = groupArray[currentFirst].nextIndex;
            //next = -1 indicates end of group
            //Loop over all of group and add correct phase to adjGroup
            while (next != -1) {

                if (groupArray[next].group == adjGroup) {
                    unwrappedImage[next] = unwrappedImage[next] + jumpDirection * 2 * M_PI * numberOfJumps;
                }

                groupArray[next].firstIndex = currentFirst;
                groupArray[next].group = currentGroup;
                groupArray[next].lastIndex = nextLast;

                next = groupArray[next].nextIndex;

            }
              
        }
        
    }
    
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
