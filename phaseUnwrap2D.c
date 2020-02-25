/**************************************************************************
INPUT: A 2D wrapped phase image
Example: unwrapped = phaseUnwrapC(wrapped);

 Michael Mullen
 mulle399@umn.edu
 Created in Matlab R2016a
 version 2.5 (October 2017)
 uses 2D-SRNCP
 an algorithm by Miguel Arevallilo Herraez, David R. Burton, 
 Michael J. Lalor, and Munther A. Gdeisat in Applied Optics, Vol. 
 41, No. 35, pp. 7437, 2002.
/*************************************************************************/

#include "mex.h"
#include "math.h"
#include "matrix.h"

/*************************************************************************/
/*                          BASIC FUNCTIONS HERE                         */
/*************************************************************************/
#define mexPI 3.14159265359

/*Function to unwrap two adjacent pixels*/
double pixelwiseUnwrap(double phase1,double phase2){
    while(fabs(phase1-phase2) > mexPI){
        
        if(phase1 > phase2 + mexPI){
            phase1 = phase1 - 2*mexPI;
        }
        else if (phase1 < phase2 - mexPI){
            phase1 = phase1 + 2*mexPI;
        }
        else{
            break;
        }
        
    }
    double result = phase1 - phase2;
    return result;
};

/*Information on edges*/
struct edgeInfo {
    char edgeType;      /*horizontal or vertical edge ('r' or 'c', respectively) */
    double R;           /*reliability*/
    mwSize row;            /*indices of pixel*/
    mwSize column;
};

struct pixelGroup {
    mwSize group;
    mwSize nextIndex;
    mwSize firstIndex;
    mwSize lastIndex;
};

/*Part of quick sort algorithm, sorting largest values to start of array*/
void quickSort(struct edgeInfo *edges, mwSize r) {
    
    if(r < 2){
        return;
    }
    
    mwSize upper, lower;
    struct edgeInfo t;
    double pivot;
    lower = 0;
    upper = r - 1;
    
    pivot = edges[lower].R;

        
    while(lower <= upper)
    {
        while( edges[lower].R > pivot){
            lower++;
        };
        while( edges[upper].R < pivot){
            upper--;
        };
        
        if( lower <= upper ){
            t.edgeType = edges[lower].edgeType;
            t.R = edges[lower].R;
            t.row = edges[lower].row;
            t.column = edges[lower].column;

            edges[lower].edgeType = edges[upper].edgeType;
            edges[lower].R = edges[upper].R;
            edges[lower].row = edges[upper].row;
            edges[lower].column = edges[upper].column;
            
            edges[upper].edgeType = t.edgeType;
            edges[upper].R = t.R;
            edges[upper].row = t.row;
            edges[upper].column = t.column;
                       
            lower++;
            upper--;
            
        }
       
    }
    quickSort(edges,upper+1);
    quickSort(&edges[lower],r-lower);
    
}

/*************************************************************************/
/*                          UNWRAPPER                                    */
/*************************************************************************/

void unwrap(mwSize numRows,mwSize numCols,double *inputArray,
        double *unwrappedImage)
{
    mwSize numGroups = 2*numRows*numCols - numRows - numCols; /*Number of edges*/
    
    struct edgeInfo *edges = mxMalloc(numGroups * sizeof(struct edgeInfo));      
    double *Reliability = mxMalloc(numRows * numCols * sizeof(double));
    struct pixelGroup *groupArray = mxMalloc(numRows * numCols * sizeof(struct pixelGroup));
    
    double columnEdge;
    double rowEdge;
    
    mwSize size = numRows * numCols * sizeof(double);
    
    mwSize currentRow,currentCol;
    mwSize group,adjIndex;
    mwSize currentGroup,adjGroup,next,nextLast,nextFirst,currentFirst,currentLast;
    char edge;
    mwSize rowIndex, colIndex;
    mwSize iters = 0;
    
    double H,V,D1,D2,D,phase1,phase2;
    
    mwSize flag = 0, m = 0;
    mwSize index,sortedIndex,count = 0;
    
    memcpy(unwrappedImage,inputArray,size);

    /*Defining reliability, using convention of paper mentioned in header*/
     for (index = 0; index < numRows * numCols; index++){
            colIndex = index % numCols;
            rowIndex = (index-colIndex)/numCols;
            if(rowIndex != 0 && rowIndex != numRows - 1 && colIndex != 0 && colIndex != numCols - 1){
            H = pixelwiseUnwrap(unwrappedImage[(rowIndex-1)*numCols + colIndex],
                    unwrappedImage[(rowIndex)*numCols + colIndex])-
                    pixelwiseUnwrap(unwrappedImage[(rowIndex+1)*numCols + colIndex],
                    unwrappedImage[(rowIndex)*numCols + colIndex]);
            
            V = pixelwiseUnwrap(unwrappedImage[(rowIndex)*numCols + colIndex + 1],
                    unwrappedImage[(rowIndex)*numCols + colIndex])-
                    pixelwiseUnwrap(unwrappedImage[(rowIndex)*numCols + colIndex+1],
                    unwrappedImage[(rowIndex)*numCols + colIndex]);
            
            D1 = pixelwiseUnwrap(unwrappedImage[(rowIndex-1)*numCols + colIndex-1],
                    unwrappedImage[(rowIndex)*numCols + colIndex])-
                    pixelwiseUnwrap(unwrappedImage[(rowIndex+1)*numCols + colIndex+1],
                    unwrappedImage[(rowIndex)*numCols + colIndex]);
            
            D2 = pixelwiseUnwrap(unwrappedImage[(rowIndex-1)*numCols + colIndex+1],
                    unwrappedImage[(rowIndex)*numCols + colIndex])-
                    pixelwiseUnwrap(unwrappedImage[(rowIndex+1)*numCols + colIndex-1],
                    unwrappedImage[(rowIndex)*numCols + colIndex]);
            
            D = sqrt(H*H + V*V + D1*D1 + D2*D2);
            Reliability[index] = 1/D;
            
            }
            else{
                Reliability[index] = 0.0;
            }
             
     }
    
    /*Define reliability of edges and store in struct edges*/
    
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
            
            count = count + 1;
        }
        
        if (colIndex < numCols - 1){
            columnEdge = Reliability[(rowIndex)*numCols + colIndex]
                    + Reliability[(rowIndex)*numCols + colIndex+1];
            edges[count].edgeType = 'c';
            edges[count].R = columnEdge;
            edges[count].row = rowIndex;
            edges[count].column = colIndex;
            
            count = count + 1;
        }
        groupArray[(rowIndex)*numCols + colIndex].group = (rowIndex)*numCols + colIndex;
        groupArray[(rowIndex)*numCols + colIndex].nextIndex = -1;
        groupArray[(rowIndex)*numCols + colIndex].lastIndex = (rowIndex)*numCols + colIndex;
        groupArray[(rowIndex)*numCols + colIndex].firstIndex = (rowIndex)*numCols + colIndex;


    }
    
    /*Done with these, free memory*/
    
    /************ SORT EDGES HERE ********************/
    quickSort(edges,numGroups); /*Initial pivot doesn't seem to matter much*/
    
    /*************************************************/
    
    for (sortedIndex = 0; sortedIndex < numGroups; sortedIndex++){
        iters = 0;
        
        m = 0;
        flag = 0;
        
        currentRow = edges[sortedIndex].row;
        currentCol = edges[sortedIndex].column;
        currentGroup = groupArray[(currentRow)*numCols + currentCol].group;
        
        edge = edges[sortedIndex].edgeType;
        
        if (edge == 'r'){
            adjIndex = (currentRow+1)*numCols + currentCol;
        }
        else{
            adjIndex = (currentRow)*numCols + currentCol + 1;
        }
        
        adjGroup = groupArray[adjIndex].group;
        phase1 = unwrappedImage[adjIndex];
        phase2 = unwrappedImage[currentRow*numCols+currentCol];
        
        if (adjGroup != currentGroup){
            
            while (fabs(phase1 - phase2) > mexPI && iters < 100){               
                if(phase1 > phase2+ mexPI){
                    phase1 = phase1 - 2*mexPI;
                    flag = -1;
                    m = m + 1;
                }
                else if(phase1 < phase2- mexPI){
                    phase1 = phase1 + 2*mexPI;
                    flag = 1;
                    m = m + 1;
                }
                else{
                    break;
                }
                iters = iters + 1;
                
            }
                        
            currentFirst = groupArray[(currentRow)*numCols + currentCol].firstIndex;
            currentLast = groupArray[(currentRow)*numCols + currentCol].lastIndex;
            /*mexPrintf("first = %d \t last = %d \n",currentFirst,currentLast);*/

            
            nextLast = groupArray[adjIndex].lastIndex;
            nextFirst = groupArray[adjIndex].firstIndex;
                    
            groupArray[currentLast].nextIndex = nextFirst;
            groupArray[currentFirst].lastIndex = nextLast;
            
            next = groupArray[currentFirst].nextIndex;
            
            while (next != -1){
                
                if(groupArray[next].group == adjGroup){
                     unwrappedImage[next] = unwrappedImage[next] + flag*2*mexPI*m;
                }
                
                groupArray[next].firstIndex = currentFirst;
                groupArray[next].group = currentGroup;
                groupArray[next].lastIndex = nextLast;
                
                next = groupArray[next].nextIndex;
            
            }
              
        }
        
    }
            
    mxFree(Reliability);    
    mxFree(edges);
    mxFree(groupArray);
}



/*************************************************************************/
/*                          MEX FUNCTION HERE                            /*
 * /*************************************************************************/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    if (nlhs != 1){
        mexErrMsgTxt("ERROR:  Number of outputs must equal 1 \n");
    }
    
    if (nrhs != 1){
        mexErrMsgTxt("ERROR:  Number of inputs must equal 1. \n");
    }
    if(mxIsComplex(prhs[0]) == 1 || mxIsDouble(prhs[0]) != 1){
        mexErrMsgTxt("ERROR:  Input must be a real double array \n");

    }
    double *inputArray;             /*Input phase image*/
    double *unwrappedImage;         /*output array, unwrapped*/
    mwSize numRows;
    mwSize numCols;
    
    inputArray = mxGetPr(prhs[0]);
    numCols = (mwSize)mxGetM(prhs[0]);
    numRows = (mwSize)mxGetN(prhs[0]);
    
    /* ===== Allocate output array, note transposed size. 
     *       since Matlab uses column major, C uses row major.           */
    
    plhs[0] = mxCreateDoubleMatrix(numCols,numRows,mxREAL);

    unwrappedImage = mxGetPr(plhs[0]);
    
    unwrap(numRows,numCols,inputArray,unwrappedImage);
    
}