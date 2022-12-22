
def makeMat(strainNamesL,distD):
    '''Construct a matrix from distance dictionary.'''
    M=[]
    for strainI in strainNamesL:
        rowL=[]
        for strainJ in strainNamesL:
            rowL.append(format(distD[(strainI,strainJ)],".3f"))
        M.append(rowL)
    return M
    
def makeColList(strainNamesL,M):
    '''Take (row based) matrix and list of column names, and make a list
       of matrix columns (including strain names first).'''
    colL=[]
    for col in range(len(strainNamesL)):
        thisCol=[strainNamesL[col]]
        for row in range(len(M)):
            thisCol.append(M[row][col])
        colL.append(thisCol)
    return colL

def fixedWidthCol(column):
    '''Take a column and pad each element with spaces so all are the same length.'''
    mx=0
    for s in column:
        if len(s)>mx:
            mx = len(s)
    newColumn=[]
    for s in column:
        newColumn.append(s+((mx+1)-len(s))*' ')
    return newColumn

def printColRange(strainNamesL,colsToPrintL):
    '''Print the columns specified in the the matrix M.'''

    # slice out rows we want
    strainNamesFixedL=fixedWidthCol(['']+strainNamesL)
    rowStrL=strainNamesFixedL
    for col in colsToPrintL:
        fwCol=fixedWidthCol(col)
        for i in range(len(fwCol)):
            row=fwCol[i]
            rowStrL[i]+=' '+row
    # now print
    for rowStr in rowStrL:
        print(rowStr)
    print()

def printDists(strainNamesL,distD,colsPerSection):
    '''Take a distance dictionary and a list of strain names, and print out nicely.'''

    M=makeMat(strainNamesL,distD)
    colL=makeColList(strainNamesL,M)
   
    for stCol in range(0,len(strainNamesL),colsPerSection):
        colsToPrintL=colL[stCol:stCol+colsPerSection]
        printColRange(strainNamesL,colsToPrintL)
