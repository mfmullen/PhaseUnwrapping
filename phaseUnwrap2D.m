function unwrapped = phaseUnwrap2D(wrapped)
%INPUT: A 2D wrapped phase image
%Example: unwrapped = phaseUnwrap2D(wrapped);

% Michael Mullen
% mulle399@umn.edu
% Created in Matlab R2015b
% version 1.2 (October 2017)
% uses 2D-SRNCP
% an algorithm by Miguel Arevallilo Herra´ez, David R. Burton, 
% Michael J. Lalor, and Munther A. Gdeisat in Applied Optics, Vol. 
% 41, No. 35, pp. 7437, 2002.

if length(size(wrapped)) ~= 2
    error('Input array is not 2D');
end

if ~isreal(wrapped)
    error('Input array must be real-valued.');
end

%************* DETERMINE PIXEL RELIABILITIES ************************
%initialize matrices to determine second differences
H = zeros(size(wrapped));
V = H;
D1 = H;
D2 = H;

pixelGroup = H;
unwrapped = wrapped;        %initial phase map
maxIters = numel(wrapped);  %maximum number of iterations for unwrappin groups
                            %most likely do not need

numRows = size(H,1);        %determine number of rows in image
numCols = size(H,2);        %determine number of columns in image

%determine second differences, ignoring edges (object almost never goes
%there anyway)

for row = 2:numRows - 1
    for col = 2:numCols - 1
        
        H(row,col) = pixelwiseUnwrap(wrapped(row-1,col),wrapped(row,col)) - ...
            pixelwiseUnwrap(wrapped(row+1,col),wrapped(row,col));
        
        V(row,col) = pixelwiseUnwrap(wrapped(row,col-1),wrapped(row,col)) - ...
            pixelwiseUnwrap(wrapped(row,col+1),wrapped(row,col));
        
        D1(row,col) = pixelwiseUnwrap(wrapped(row-1,col-1),wrapped(row,col)) - ...
            pixelwiseUnwrap(wrapped(row+1,col+1),wrapped(row,col));
        
        D2(row,col) = pixelwiseUnwrap(wrapped(row-1,col+1),wrapped(row,col)) - ...
            pixelwiseUnwrap(wrapped(row+1,col-1),wrapped(row,col));
                
    end
end
D = sqrt(H.^2 + V.^2 + D1.^2 + D2.^2);
        
Reliability = 1./D;
indices = find(isinf(Reliability));
Reliability(indices) = zeros(size(indices));

%************* ASSIGN EDGE RELIABILITIES *****************
%also initialize each pixel to its own group

edges = struct();            %keeping edge information in structure array
                             %tracks reliability, edge type, and pixel

columnEdges = Reliability(:,1:end-1) + Reliability(:,2:end);
rowEdges = Reliability(1:end-1,:) + Reliability(2:end,:);

count = 1;
for row = 1:numRows
    for col = 1:numCols
        
        if row < numRows
            edges(count).edgeType = 'row';
            edges(count).R = rowEdges(row,col);
            edges(count).pixel = [row col];
            count = count + 1;
        end
        
        if col < numCols
            edges(count).edgeType = 'column';
            edges(count).R = columnEdges(row,col);
            edges(count).pixel = [row col];
            count = count + 1;
            
        end
        
        pixelGroup(row,col) = (row-1) * numCols + col;
                    
    end
end

%************* SORT EDGES ***********************
[~, indices] = sort([edges.R],'descend');
Sorted = edges(indices);

%************* BEGIN UNWRAPPING PROCEDURE ***********************
sortedIndex = 1;

while sortedIndex <= length(Sorted)
    
    currentPixel = Sorted(sortedIndex).pixel;
    currentRow = currentPixel(1);
    currentCol = currentPixel(2);
    currentGroup = pixelGroup(currentRow,currentCol);
    edge = Sorted(sortedIndex).edgeType;
    
    %unwrap two pixels
    if strcmp(edge,'row')
        
        adjGroup = pixelGroup(currentRow+1,currentCol);
        
        if adjGroup ~= currentGroup
            [pixelGroup, unwrapped] = groupUnwrap(unwrapped,currentRow,currentCol,edge,currentGroup,adjGroup,maxIters,pixelGroup);            
        end
        
    else
        
        adjGroup = pixelGroup(currentRow,currentCol+1);
       
        if adjGroup ~= currentGroup
            [pixelGroup,unwrapped] = groupUnwrap(unwrapped,currentRow,currentCol,edge,currentGroup,adjGroup,maxIters,pixelGroup);            
        end
        
    end
    
    sortedIndex = sortedIndex + 1;
    
end

end

function [newPixelGroups, newPhaseImage] = groupUnwrap(phaseImage,currentRow,currentCol,edge,group1,group2,maxIters,pixelGroups)

m = 0;          %internal variable to determine # of 2 pi steps
iters = 1;      %put maximum number of iterations, probably not necessary
flag = 0;       %determine whether to add or subtract multiples of 2 pi
newPhaseImage = phaseImage;     %initialize output matrix
newPixelGroups = pixelGroups;   %initialize other output matrix

if strcmp(edge,'row') == 1
    while abs(phaseImage(currentRow+1,currentCol) - phaseImage(currentRow,currentCol)) > pi && iters < maxIters
        if phaseImage(currentRow+1,currentCol) > phaseImage(currentRow,currentCol) + pi
            phaseImage(currentRow+1,currentCol) = phaseImage(currentRow+1,currentCol) - 2*pi;
            flag = -1;
            m = m + 1;
            
        elseif phaseImage(currentRow+1,currentCol) < phaseImage(currentRow,currentCol) - pi
            phaseImage(currentRow+1,currentCol) = phaseImage(currentRow+1,currentCol) + 2*pi;
            flag = 1;
            m = m + 1;
            
        else
            
            break;
        end
        iters = iters + 1;
    end
    
elseif strcmp(edge,'column') == 1
    while abs(phaseImage(currentRow,currentCol+1) - phaseImage(currentRow,currentCol)) > pi && iters < maxIters
        if phaseImage(currentRow,currentCol+1) > phaseImage(currentRow,currentCol) + pi
            flag = -1;
            m = m + 1;           
            phaseImage(currentRow,currentCol+1) = phaseImage(currentRow,currentCol+1) - 2*pi;
        
        elseif phaseImage(currentRow,currentCol+1) < phaseImage(currentRow,currentCol) - pi
            phaseImage(currentRow,currentCol+1) = phaseImage(currentRow,currentCol+1) + 2*pi;
            flag = 1;
            m = m + 1;
            
        else
            
            break;
        end
        iters = iters + 1;
        
    end
    
else
    error(['Unknown edge type ',edge,' not recognized']);
end

%determine smaller group
k1 = find(pixelGroups==group1);
k2 = find(pixelGroups==group2);

if length(k1) ~= length(k2)
    minLength = min(length(k1),length(k2));
else
   minLength = length(k1); 
end

%merge groups
if minLength == length(k1)
    k = k1;
    newPixelGroups(k) = group2 * ones(size(pixelGroups(k)));
else %minLength == length(k2), group 2 is smaller
    k = k2;
    newPixelGroups(k) = group1 * ones(size(pixelGroups(k)));
end

%add phase to the SMALLER group
newPhaseImage(k2) = newPhaseImage(k2) + flag*2*pi*m;

end

function result = pixelwiseUnwrap(phase1,phase2)
%initialize parameters so we don't get stuck in a loop. Probably don't need
iters = 1;
maxIters = 1e4;

%determine unwrapped phase difference
while abs(phase1 - phase2) > pi && iters < maxIters
    if phase1 > phase2 + pi
        phase1 = phase1 - 2*pi;
                
    elseif phase1 < phase2 - pi
        phase1 = phase1 + 2*pi;
                
    else
        break;
    end
    iters = iters + 1;
end

%output result
result = phase1 - phase2;

end