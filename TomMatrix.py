###############################
##~~Written by Tom Blanchet~~##
#~~~~~Written  2012-2013~~~~~~#
###############################
#import itertools
#import Tkinter
import math
import sys
import random
import time
import mpmath
#import copy

MULTIPLIERTHRESHHOLD = 350000 ##Yes, actually...

def timer(func, *args, **keywords):
    start = time.clock()
    retVal = func(*args, **keywords)
    stop = time.clock()
    return retVal, stop-start

class EigenError(StandardError):
    pass
class MatrixError(StandardError):
    pass
class matrix:
    ##Notes:
    ##  Need to Do:
    ##      Tune the mulitplier threshhold
    ##      Efficient algorithms for sparse Matricies
    ##      Make algorithms work with keyvalues
    ##      "exact" Eigenvalues (When possible)
    ##      Make QRDecomposition more efficient
    ##      Make Schur more efficient

    """This is a matrix class that deals with operations and manipulations of
        real matrices.

        A few ways to make a matrix:
        matrix(<number>, rows, cols) makes a diagonal matrix with the number along the diagonal.
            E.G.
            >>> print matrix(2.5, 3, 5)
            [2.5, 0.0, 0.0, 0.0, 0.0]
            [0.0, 2.5, 0.0, 0.0, 0.0]
            [0.0, 0.0, 2.5, 0.0, 0.0]

        matrix(<list>, rows, cols) fills a matrix with the elements from the list.
            The list length must be < rows*cols, and will fill by row, then column
            E.G.
            >>> print matrix([1, 2, 3], 2, 2)
            [1, 2]
            [3, 0]
            >>> matrix([1, 2, 3, 4, 5], 2, 2)
            Traceback (most recent call last):
              File "<pyshell#189>", line 1, in <module>
                matrix([1, 2, 3, 4, 5], 2, 2)
              File "C:/Python27/TomMatrix.py", line 67, in __init__
                raise BufferError("data must be able to fit in the matrix")
            BufferError: data must be able to fit in the matrix

        matrix(<dict>, rows, cols) fills the matrix based on keys from a dictionary.
            The keys of the dictionary MUST be 2-tuples to denote rows and columns.
            The keys may contain strs, ints, or longs. To view str rows or columns,
            use matrix.KeyPrint(). If Rows and Cols = 0, this matrix will only
            contain the keys that are given in the dictionary. Note that the size
            of the matrix does NOT inclue key valued entries.
            (Currently, this feature is not fully supported)

        matrix(None, rows, cols) is the same as matrix(0, rows, cols)

        By default, matrix() makes a 2x2 matrix with all zero entries. 

        For an instanciated matrix m:
            Get or set an entry: m[row, col] or m[index]
            Get a submatrix: m[rowStart, colStart: rowStop, colStop]
            
            Get a row: m.getRow(row)
            Pop a row: m.popRow(row)
            Add a row: m.addRow(row)
            Get a column: m.getCol(col)
            Pop a column: m.popCol(col)
            Add a column: m.addCol(col)

            Make a copy of m: m.Copy()
            Transpose m: m.transpose() or m.t()
            Make a copy of the ranspose of m: m.tCopy() or m.tC()

            To Print: m.Print() or print m (the latter for a unicode free print)
                      m.KeyPrint() prints the matrix and key valued entries.

        For 2 matricies m and n:
            m+n is the sum of m and n
            m*n is the matrix product of m and n
            
        If m is a square matrix (check with m.isSquare()):
            m.inv() = m**-1 calculates the inverse of m
            n/m = n*(m.inv())
            m**k = m*m*m*...*m k times. (k must be an integer)
            m.det() is the determinant of m

        See below for other operations the matrix can do.
        """
    
    def __init__(self, data=None, rows=2, cols=2, makeMpf = False):
        self._numRows = rows
        self._numCols = cols
        self._printer = printer()
        self._keyPrinter = keyValuedPrinter()
        self._data = dict()
        self.isMpf = False
##        if data == None:
##            data = 0
##            self._data = [0]*rows*cols
        if any([type(data)==i for i in [int, float, long, complex, mpmath.mpf]]):
            for i in range(min(rows, cols)):
                self._data[i, i] = data
##            self._data = [data*(i==j) for i in range(rows) for j in range(cols)]
        elif type(data) == dict:
            self._data = data
            self._checkKeys()
        elif type(data) == list:
            if len(data)>rows*cols:
                raise BufferError("data must be able to fit in the matrix")
            for i in xrange(len(data)):
                if all([not isinstance(data[i], t) for t in [int, float, long, complex, mpmath.mpf]]):
                    raise TypeError("Entries must be an int, float, long, complex, or mpmath.mpf.")
                self._data[i//self._numCols, i%self._numCols] = data[i]
        self._printer(self, doPrint = False)
        self._keyPrinter(self, doPrint=False)
        if makeMpf:
            self.makeMpf()

    def _checkKeys(self):
        for key in self._data.iterkeys():
            if type(key)!=tuple and len(key)!=2:
                raise MatrixError("Keys must be 2 tuples.")
            elif not any(type(key[0])==i for i in (int, long, str)):
                raise MatrixError("Keys must only contain ints, longs, or strings.")
            elif not any(type(key[1])==i for i in (int, long, str)):
                raise MatrixError("Keys must only contain ints, longs, or strings.")

    def Print(self):
        self._printer()

    def KeyPrint(self):
        self._keyPrinter()
        
    def __len__(self):
        return len(self._data)
    
    def copy(self):
        return matrix(self._data.copy(), self.dims()[0], self.dims()[1])

    def getNumCols(self):
        return self._numCols

    def getNumRows(self):
        return self._numRows

    def __eq__(self, other):
        if isinstance(other, matrix):
            return all([self._data == other._data,
                        self._numRows == other._numRows, self._numCols == self._numCols])
        else: return False

    def __ne__(self, other):
        return not self == other

    def dims(self):
        return self._numRows, self._numCols
    
    def __getitem__(self, key):
        if type(key) == int:
            try:
                return self._data[key//self._numRows, key%self._numRows]
            except KeyError:
                if key//self._numRows<self._numCols:
                    return float(0)
                else:
                    raise IndexError("Index out of range")
        if type(key) == tuple:
            if len(key) == 2:
                ##Currently, slices do not parse correctly.
                ##This code corrects that.
                if type(key[0]) == slice:
                    j = (key[0].stop, key[1])##hacky hacky hacky
                    NewMat = matrix(rows = j[0], cols=j[1])
                    for row in range(j[0]):
                        for col in range(j[1]):
                            NewMat[row, col] = self[row, col]
                    return NewMat ###Return
                elif type(key[1]) == slice and key[1].stop == None:
                    i = (key[0], key[1].start)##hacky hacky hacky
                    nrows = self._numRows-i[0]
                    ncols = self._numCols-i[1]
                    NewMat = matrix(0, nrows, ncols)
                    for row in range(nrows):
                        for col in range(ncols):
                            NewMat[row, col] = self[i[0]+row, i[1]+col]
                    return NewMat ###Return
                elif self._data.has_key(key):
                    return self._data[key]
                elif (key[0]<self._numRows)&(key[1]<self._numCols)&(key[0]>=0)&(key[1]>=0):
                    return 0
                elif type(key[0])==str or type(key[1])==str:
                    return 0
                else:
                    raise IndexError("Index out of Range")

        if type(key) == slice:
            return self[key.start[0], key.start[1]:key.stop[0], key.stop[1]]

        if (len(key) == 3) and (type(key[1]) == slice):
            i = (key[0], key[1].start)##hacky hacky hacky
            j = (key[1].stop, key[2])
            if type(i) == tuple:
                if type(j) == tuple:
                    nrows = j[0]-i[0]
                    ncols = j[1]-i[1]
                    NewMat = matrix(0, nrows, ncols)
                    for row in range(nrows):
                        for col in range(ncols):
                            NewMat[row, col] = self[i[0]+row, i[1]+col]
                    return NewMat ###Return

        raise TypeError("Key must be in a proper form.")

                        

    def __setitem__(self, key, value):
        if self.isMpf:
            value = mpmath.mpf(value)
##        if type(value) !=  self._mtype:
##            raise TypeError("must be the same type as the matrix")
        if type(key) == int:
            self._data[key//self._numRows, key%self._numRows] = value
            return
        if type(key) == tuple:
            if len(key) == 2:
                if (key[0]<self._numRows)&(key[1]<self._numCols)&(key[0]>=0)&(key[1]>=0):
##                    self._data[self._numCols*key[0] + key[1]] = value
                    self._data[key[0],key[1]] = value
                    return
                elif type(key[0]) == str and any(type(key[1])==i for i in (str, int, long)):
                    self._data[key[0], key[1]] = value
                    return
                elif (type(key[1]) == str) and any(type(key[0])==i for i in (str, int, long)):
                    self._data[key[0], key[1]] = value
                    return
                else:
                    raise IndexError("list index out of range")
        raise TypeError("must be an int or two ints separated by a comma")

    def makeMpf(self):
        for k in self._data.iterkeys():
            self._data[k] = mpmath.mpf(self._data[k])
        self.isMpf = True

##    def makeFloat(self):
##        for k in self._data.iterkeys():
##            self._data[k] = float(self._data[k])
##        self.isMpf = False
        
    def addCol(self, position, data=None):
        """Adds a new column to the matrix"""
        self.t()
        self.addRow(position, data)
        self.t()

    def addRow(self, position, data=None):
        """Adds a new row to the matrix"""
        if position<0:
            position += self._numRows+1
        if data!=None:
            if len(data)>self._numCols:
                raise BufferError("The data must fit in the matrix")
        self._numRows+=1
        for col in xrange(self._numCols):
            for row in xrange(self._numRows, position-1, -1):
                if self._data.has_key((row, col)):
                    self._data[row+1, col] = self._data.pop((row, col))
        if data == None:
            return
        else:
            for i in xrange(len(data)):
                self._data[position, i]=data[i]
                
##        if position<0:
##            position+=self._numRows+1
##        if data == None:
##            data = [0]*self._numCols
##        if len(data)!=self._numCols:
##            raise BufferError("The data must fit in the matrix")
##        if position>self._numRows or position<0:
##            raise IndexError("Position is out of range")
##        dataPos = self._numCols*position
##        listA = []
##        listB = []
##        for i in range(len(self._data)):
##            if i<dataPos:
##                listA.append(self._data[i])
##            else:
##                listB.append(self._data[i])
##        self._data = listA+data+listB
##        self._numRows+=1

    def popCol(self, position):
        """Returns and removes the column at the position."""
        self.t()
        ret = self.popRow(position)
        self.t()
        return ret
    
    def popRow(self, position):
        """Returns and removes the row at the position."""
        
        if position<0:
            position+=self._numRows
        if position<0 or position>=self._numRows:
            raise IndexError("position out of Range")
        ret = []
        for col in xrange(self._numCols):
            if self._data.has_key((position, col)):
                ret.append(self._data.pop((position, col)))
            else:
                ret.append(0)
            for row in xrange(position+1, self._numRows):
                if self._data.has_key((row, col)):
                    self[row-1, col] = self._data.pop((row, col))
        self._numRows-=1
        return ret
##        if position>self._numRows or position<0:
##            raise IndexError("Position is out of range")
##        dataPos = self._numCols*position
##        listA = []
##        listB = []
##        ret = []
##        for i in range(len(self._data)):
##            if i<dataPos:
##                listA.append(self._data[i])
##            elif i>=dataPos+self._numCols:
##                listB.append(self._data[i])
##            else:
##                ret.append(self._data[i])
##        self._data = listA+listB
##        self._numRows-=1
##        return ret

    def swapCol(self, pos1, pos2):
        """Swaps the columns at pos1 and pos2."""
        self.t()
        self.swapRow(pos1, pos2)
        self.t()

    def swapRow(self, pos1, pos2):
        if pos1<0:
            pos1+= self._numRows
        if pos2<0:
            pos2+= self._numRows
        if pos1>self._numRows or pos1<0:
            raise IndexError("Pos1 is out of range")
        if pos2>self._numRows or pos2<0:
            raise IndexError("Pos2 is out of range")
        row1 = self.getRow(pos1)
        row2 = self.getRow(pos2)
        for k in range(self._numCols):
            self[pos2, k] = row1[k]
            self[pos1, k] = row2[k]
        
    def __str__(self):
        ret = ""
        return ret.join(str(self.getRow(i))+"\n" for i in range(self._numRows))[:-1]
    
    def __call__(self):
        return dict(data=self._data, rows=self._numRows, cols=self._numCols, asStr=self.__str__())

    def dotProduct(self, other):
        """Takes the Dot Product
            This function treats each matrix like a vector."""
        if not isinstance(other, matrix):
            raise NotImplemented
        if (self._numRows != other._numRows)|(self._numCols != other._numCols):
            raise MatrixError("must be the same size of matrix")
        return sum(self[i]*other[i] for i in range(len(self)))

    def iterkeys(self):
        return self._data.iterkeys()

    def getKeysInRows(self):
        ret = dict()
        for key in self.iterkeys():
            if ret.has_key(key[0]):
                ret[key[0]].append(key)
            else:
                ret[key[0]]=[key]
        return ret

    def getKeysInCols(self):
        ret = dict()
        for key in self.iterkeys():
            if ret.has_key(key[1]):
                ret[key[1]].append(key)
            else:
                ret[key[1]]=[key]
        return ret

    def getKeysInRowsAndCols(self):
        retRows = dict()
        retCols = dict()
        for key in self.iterkeys():
            if retRows.has_key(key[0]):
                retRows[key[0]].append(key)
            else:
                retRows[key[0]]=[key]
            if retCols.has_key(key[1]):
                retCols[key[1]].append(key)
            else:
                retCols[key[1]]=[key]
        return retRows, retCols

    def has_key(self):
        return self._data.has_key
    
    def __add__(self, other):
        if not isinstance(other, matrix):
            raise NotImplemented
        if (self._numRows != other._numRows)|(self._numCols != other._numCols):
            raise MatrixError("must be the same size of matrix")
##            ret = matrix(rows=max(self._numRows, other._numRows),\
##                     cols=max(self._numCols, other._numCols))
##            for j in range(ret._numCols):
##                for i in range(ret._numRows):
##                    ret[i, j] += self[i, j]
##                    ret[i, j] += other[i, j]
        else:
            ret = self.copy()
            for key in other._data.iterkeys():
                ret[key] += other[key]
        return ret

        
    def __sub__(self, other):
        if not isinstance(other, matrix):
            raise NotImplemented
        if (self._numRows != other._numRows)|(self._numCols != other._numCols):
            raise MatrixError("must be the same size of matrix")
##            ret = matrix(rows=max(self._numRows, other._numRows),\
##                     cols=max(self._numCols, other._numCols))
##            for j in range(ret._numCols):
##                for i in range(ret._numRows):
##                    ret[i, j] += self[i, j]
##                    ret[i, j] -= other[i, j]
        else:
            ret = self.copy()
            for key in other._data.iterkeys():
                ret[key] -= other[key]
        return ret

    def transpose(self):
        """Transposes the matrix."""
        newRows = self._numCols
        newCols = self._numRows
        new_data = dict()
        for key in self._data.iterkeys():
            new_data[key[1], key[0]]= self._data[key]
        self._numRows = newRows
        self._numCols = newCols
        self._data = new_data
        
    def tCopy(self):
        """Returns the Transpose of the Matrix"""
        ret = self.copy()
        ret.t()
        return ret
    
    t = transpose ##alias
    tC = tCopy
        
    def isSquare(self):
        return self._numRows == self._numCols

    def makeSquare(self):
        dif = self._numRows-self._numCols
        if dif>0:
            for i in range(dif):
                self.addCol(-1)
        elif dif<0:
            for i in range(-dif):
                self.addRow(-1)

    def getRow(self, rowNum):
        return [self[rowNum, i] for i in range(self._numCols)]

    
    def getCol(self, colNum):
        return [self[i, colNum] for i in range(self._numRows)]

    def oldMul(self, other):
        if any([type(other)==i for i in [int, float, long, complex]]):
            return matrix([j*other for j in self._data], self._numRows, self._numCols)
        if not isinstance(other, matrix):
            raise NotImplemented
        if self._numCols != other._numRows:
            raise MatrixError(
                "in m*n, the number of rows of n must match the number of columns of m")
        new_matrix = matrix(rows = self._numRows, cols = other._numCols)
        for i in xrange(new_matrix._numRows):
            for j in xrange(new_matrix._numCols):
                new_matrix[i, j] = sum(self[i, k]*other[k, j] for k in\
                                        xrange(min(self._numCols, other._numRows)))
        
        return new_matrix
    
    def __mul__(self, other):
        if any([type(other)==i for i in [int, float, long, complex, mpmath.mpf]]):
            return matrix([j*other for j in self], self._numRows, self._numCols)
        if not isinstance(other, matrix):
            raise NotImplemented
        if self._numCols != other._numRows:
            raise MatrixError(
                "in m*n, the number of rows of n must match the number of columns of m")
        mindim = min(self._numRows, self._numCols, other._numRows, other._numCols)
        if mindim*(mpmath.mp.prec**self.isMpf)<=MULTIPLIERTHRESHHOLD or mindim<=1: 
            new_matrix = matrix(rows = self._numRows, cols = other._numCols)
##            for j in xrange(new_matrix._numCols):    
##                for i in xrange(new_matrix._numRows):
##                    new_matrix[i, j] = fastDot([(self[i, k],other[k, j]) for k in\
##                                           xrange(self._numCols)])
            SRows, SCols = self.getKeysInRowsAndCols()
            ORows, OCols = other.getKeysInRowsAndCols()
            SColsAndORows = SCols.viewkeys()&ORows.viewkeys()
            if len(SColsAndORows)!=0:
                if self.isMpf and other.isMpf:
                    for SRow in SRows.iterkeys():
                        for OCol in OCols.iterkeys():
                            new_matrix[SRow, OCol] = mpmath.fdot((self[SRow, k], other[k, OCol]) for k in SColsAndORows)
                else:
                    for SRow in SRows.iterkeys():
                        for OCol in OCols.iterkeys():
                            new_matrix[SRow, OCol] = fastDot((self[SRow, k], other[k, OCol]) for k in SColsAndORows)
##                if iterORows.has_key(keyS[1]):
##                    for keyO in iterORows[keyS[1]]:
##                        if toDotDict.has_key((keyS[0], keyO[1])):
##                            toDotDict[keyS[0], keyO[1]].append((self[keyS],other[keyO]))
##                        else:
##                            toDotDict[keyS[0], keyO[1]] = [(self[keyS],other[keyO])]
##            if self.isMpf and other.isMpf:
##                for key in toDotDict.iterkeys():
##                        new_matrix[key] = mpmath.fdot(toDotDict[key])
##            else:
##                for key in toDotDict.iterkeys():
##                        new_matrix[key] = fastDot(toDotDict[key])
            return new_matrix
        else: #Strassen
            selfC = self.copy()
            otherC = other.copy()
            selfC.makeSquare()
            otherC.makeSquare()
            dif = otherC.getNumRows()-selfC.getNumRows()
            if dif>0:
                for i in range(dif):
                    selfC.addRow(-1)
                    selfC.addCol(-1)
            elif dif<0:
                for i in range(-dif):
                    otherC.addRow(-1)
                    otherC.addCol(-1)
            if self.getNumRows()%2==1:
                selfC.addRow(-1)
                selfC.addCol(-1)
                otherC.addRow(-1)
                otherC.addCol(-1)
                
            A_segs = (selfC._numRows//2, selfC._numCols//2)
            B_segs = (otherC._numRows//2, otherC._numCols//2)
            A = [[selfC[:A_segs[0], A_segs[1]],selfC[0, A_segs[1]:A_segs[0],selfC._numRows] ],\
                 [selfC[A_segs[0], 0: selfC._numRows, A_segs[1]], selfC[A_segs[0], A_segs[1]:]]]
            B = [[otherC[:B_segs[0], B_segs[1]],otherC[0, B_segs[1]:B_segs[0],otherC._numRows] ],\
                 [otherC[B_segs[0], 0: otherC._numRows, B_segs[1]], otherC[B_segs[0], B_segs[1]:]]]

            M_1 = (A[0][0]+A[1][1])*(B[0][0]+B[1][1])
            M_2 = (A[1][0]+A[1][1])*B[0][0]
            M_3 = A[0][0]*(B[0][1]-B[1][1])
            M_4 = A[1][1]*(B[1][0]-B[0][0])
            M_5 = (A[0][0]+A[0][1])*B[1][1]
            M_6 = (A[1][0]-A[0][0])*(B[0][0]+B[0][1])
            M_7 = (A[0][1]-A[1][1])*(B[1][0]+B[1][1])
            blockMat = [[M_1+M_4-M_5+M_7, M_3+M_5],\
                        [M_2+M_4        , M_1-M_2+M_3+M_6]]
            
            new_matrix = matrix(None, self._numRows, other._numCols)
            for i in range(new_matrix._numRows):
                for j in range(new_matrix._numCols):
                    isRow0 = i<blockMat[0][0]._numRows
                    isCol0 = j<blockMat[0][0]._numCols
                    new_matrix[i, j] =\
                                  blockMat[not isRow0][not isCol0]\
                                  [i-((not isRow0)*blockMat[0][0]._numRows),\
                                   j-((not isCol0)*blockMat[0][0]._numCols)]
            return new_matrix

    def __rmul__(self, other):
        if any([type(other)==i for i in [int, float, long, complex]]):
            return self*other
        else:
            raise NotImplemented

    def HadamardProd(self, other):
        """Returns the Hadamard Product of the matrices"""
        if not isinstance(other, matrix):
            raise NotImplemented
        if self.dims()!=other.dims():
            raise MatrixError(
                "Matrices must be the same size.")
        new_matrix = matrix(None, self._numRows, self._numCols)
        for i in range(new_matrix._numRows):
            for j in range(new_matrix._numCols):
                new_matrix[i, j] = self[i, j]*other[i, j]
        return new_matrix

    def KroneckerProd(self, other):
        """Returns the Kronecker Product of the matrices"""
        rows = self._numRows * other._numRows
        cols = self._numCols * other._numCols
        outMat = matrix(rows=rows, cols=cols)
        for j in range(cols):
            for i in range(rows):
                outMat[i, j] = self[i/other._numRows, j/other._numCols]*other[i%other._numRows,j%other._numCols]
        return outMat

    def __div__(self, other):
        if isinstance(other, matrix):
            return self*other.inv()
        elif any([type(other)==i for i in [int, float, long, complex]]):
            return matrix([j/other for j in self._data], self._numRows, self._numCols)
        else:
            raise NotImplemented
        
    def __rdiv__(self, other):
        if any([type(other)==i for i in [int, float, long, complex]]):
            return other*self.inv()
        else:
            raise NotImplemented
        
    def __mod__(self, other):
        if not ((type(other)==int) or (type(other)==long)):
            raise NotImplemented
        else:
            i = 0
            stopCon = len(self._data)
            while(i<stopCon):
                self._data[i] %= other
                i += 1
    def __pow__(self, other):
        if not self.isSquare():
            raise MatrixError(
                "must be a square matrix")
        k = bin(abs(other))
        n = self
        ptr = 1
        ret = matrix(rows=self._numRows, cols=self._numCols)##Will make the unit matrix
        for i in range(self._numRows): 
            for j in range(self._numCols):
                if i==j:
                    ret[i,j] = 1
        while(k[-ptr] != 'b'):
            if(n == 1):
                    return ret 	#exit
            if(k[-ptr] == '1'):
                    ret = (ret*n)
            n = (n*n)
            ptr = ptr + 1
        if other<0:
            return ret.inv()
        else:
            return ret

    def inv(self):
        """Calculates the inverse of the matrix."""
        ###This is a somewhat novel algortihm.
        ###It is based on Back Substitution.
        if not self.isSquare():
            raise MatrixError("The matix must be square to find the inverse.")
        Q, R, sign = self.QRDecomposition()
        invR = matrix(0, self._numCols, self._numRows)
        rowindex = range(self._numRows)
        rowindex.reverse()
        colindex = range(self._numCols)
        colindex.reverse()
        for k in colindex: ##for each column
            for j in rowindex: #for each row
                if j<=k:
                    if R[j, j] == 0:
                        raise MatrixError("The matrix is not invertable.")
                    invR[j, k] = ((j==k)-sum(R[j, i]*invR[i, k] for i in range(j+1, self._numRows)))/float(R[j,j])
        return invR*Q.tC()
        
##    def LUPDecomp(self):
##        """Returns the LUP Decomposition of the matrix.
##
##            A fourth value is returned as well, this is 1
##            if there were an even number of row exchanges in P,
##            and -1 if not."""

    def Img(self, threshhold = 0.1**5):
        """Returns the image of the matrix
            The threshhold is for determining at what point a
            vector is a zero vector, to account for floating
            point errors. A good estimate of the threshhold makes
            for faster calculation."""
        ##This uses RRQR to determine an orthonomal basis for the image
        CurA = self.copy()
        stayFlag = True
        Q = matrix(1, self._numRows, self._numRows)
##        qList = []
        count = 0
        reduced = False
        for j in range(min(self._numCols, self._numRows)):
            stayFlag = True
            skipFlag = False
            while(stayFlag):
                if count == CurA._numCols:
                    Q.t()
                    retList = []
                    for i in range(Q._numCols - CurA._numRows):
                        retList.append(Q.popCol(0))
                    return retList ###Exit
                stayFlag=False
                alpha = (sum(CurA[i, 0]**2 for i in range(CurA._numRows))**.5)
                u = vector([0]*j+[CurA[i, 0]+(i==0)*alpha for i in range(CurA._numRows)])
                magU = sum(uIn**2 for uIn in u)**.5
                if magU<=threshhold:
                    skipFlag = True
                if skipFlag == False:
                    Q.HouseholderRef(u)
                    CurA.HouseholderRef(u[j:])
                skipFlag = False
##                v = u/magU
##                q = matrix(1, self._numCols, self._numCols)-v.outerProduct(v*2)
##                qList.append(q)
##                qCut = matrix([q[row, col] for col in range(j, q._numRows) for row in range(j, q._numRows)],
##                              CurA._numRows, CurA._numRows)
##                CurA = qCut*CurA
                if abs(CurA[0,0]) < threshhold:
                    if CurA._numCols == 1:
                        count += 1
                        stayFlag = True
                        continue
                    PushCol = CurA.popCol(0)
                    CurA.addCol(-1, PushCol)
                    count += 1
                    stayFlag = True
            if j < min(self._numCols, self._numRows)-1:
                CurA = matrix([CurA[j, i]  for j in range(1, CurA._numRows) for i in range(1, CurA._numCols)],
                              CurA._numRows-1, CurA._numCols-1)
        Q.t()        
        ret = [Q.popCol(0) for i in range(Q._numCols)]
        if abs(CurA[0,0]) < threshhold:
            ret = ret[:-1]
        return ret
        
    
##        AList = [vector(self.getCol(i)) for i in range(self._numCols)]
##        UList = [AList[0]]
##        EList = []
##        ##If the varience of the entries are above the threshhold,
##        ##then we concider the vector to be "non-zero" (varience
##        ##assuming the expected value of the vector is 0)
##        if (UList[0].norm()/float(self._numRows))**.5>threshhold:
##            EList = [AList[0]/AList[0].norm()]
##        for a in AList[1:]:
##            UList.append(a - sum(e.proj(a) for e in EList))
##            if (UList[-1].norm()/float(self._numRows))**.5>threshhold:
##                EList.append(UList[-1]/(UList[-1]).norm())
##            else:
##                UList = UList[:-1]
##        return [list(k) for k in EList]

    def Ker(self, threshhold = 0.1**5):
        """Returns the kernel of the matrix.
            The threshhold is for determining at what point a
            vector is a zero vector, to account for floating
            point errors. A good estimate of the threshhold makes
            for faster calculation."""
        
        ##This uses RRQR to determine an orthonomal basis for the kernel
        CurA = self.copy()
        CurA.t()
        stayFlag = True
        Q = matrix(1, self._numCols, self._numCols)
        if self.isMpf:
            Q.makeMpf()
            CurA.makeMpf()
##        qList = []
        count = 0
        reduced = False
        for j in range(min(self._numCols, self._numRows)):
            stayFlag = True
            skipFlag = False
            while(stayFlag):
                if count == CurA._numCols:
                    Q.t()
                    retList = []
                    for i in range(CurA._numRows):
                        retList.append(Q.popCol(-1))
                    return retList ###Exit
                stayFlag=False
                alpha = (sum(CurA[i, 0]**2 for i in range(CurA._numRows))**.5)
                u = vector([0]*j+[CurA[i, 0]+(i==0)*alpha for i in range(CurA._numRows)])
                magU = sum(uIn**2 for uIn in u)**.5
                if magU<=threshhold:
                    skipFlag = True
                if skipFlag == False:
                    Q.HouseholderRef(u)
                    CurA.HouseholderRef(u[j:])
                skipFlag = False
##                v = u/magU
##                q = matrix(1, self._numCols, self._numCols)-v.outerProduct(v*2)
##                qList.append(q)
##                qCut = matrix([q[row, col] for col in range(j, q._numRows) for row in range(j, q._numRows)],
##                              CurA._numRows, CurA._numRows)
##                CurA = qCut*CurA
                if abs(CurA[0,0]) < threshhold:
                    if CurA._numCols == 1:
                        count += 1
                        stayFlag = True
                        continue
                    PushCol = CurA.popCol(0)
                    CurA.addCol(-1, PushCol)
                    count += 1
                    stayFlag = True
            if j < min(self._numCols, self._numRows)-1:
                CurA = matrix([CurA[j, i]  for j in range(1, CurA._numRows) for i in range(1, CurA._numCols)],
                              CurA._numRows-1, CurA._numCols-1)
        Q.t()
        retList = []
        if abs(CurA[0,0]) < threshhold:
            retList.append(Q.popCol(-1))
        for i in range(CurA._numRows-1):
            retList.append(Q.popCol(-1))
        return retList ###Exit

##        ###this approximates, designed for non-rational form matricies.
##        ###Finding an orthonormal basis that spans the rows of the matrix.
##        ###Using the Gram-Schmidt process
##        AList = [vector(self.getRow(i)) for i in range(self._numRows)]
##        UList = [AList[0]]
##        EList = []
##        ##If the varience of the entries are above the threshhold,
##        ##then we concider the vector to be "non-zero" (varience
##        ##assuming the expected value of the vector is 0)
##        if (UList[0].norm()/float(self._numCols))**.5>threshhold: 
##            EList = [AList[0]/AList[0].norm()]
##        for a in AList[1:]:
##            UList.append(a - sum(e.proj(a) for e in EList))
##            if (UList[-1].norm()/float(self._numCols))**.5>threshhold:
##                EList.append(UList[-1]/(UList[-1]).norm())
##            else:
##                UList = UList[:-1]
##        IdMatrix = matrix(1000, self._numCols, self._numCols)
##        IList = [vector(IdMatrix.getCol(i)) for i in range(self._numCols)]
##        ##AList contains the space that is orthogonal to the rows of the matrix
##        AList = [i - sum(e.proj(i) for e in EList) for i in IList]
##        UList = [AList[0]]
##        EList = []
##        if (UList[0].norm()/float(self._numCols))**.5>threshhold:
##            EList = [AList[0]/AList[0].norm()]
##        for a in AList[1:]:
##            UList.append(a - sum(e.proj(a) for e in EList))
##            if (UList[-1].norm()/float(self._numCols))**.5>threshhold:
##                EList.append(UList[-1]/(UList[-1]).norm())
##            else:
##                UList = UList[:-1]
####        rmVal = []
####        for ai in range(len(AList)):
####            if AList[ai].norm()>threshhold:
####                AList[ai] = AList[ai]/AList[ai].norm()
####            else:
####                rmVal.append(AList[ai])
####        for val in rmVal:
####            AList.remove(val)
#
#        return [list(k) for k in EList]

    def trace(self):
        """Returns the trace of the matrix."""
        return sum(self[i, i] for i in range(min(self._numRows, self._numCols)))

    def det(self):
        """Finds the determinant."""
        if not self.isSquare():
            raise MatrixError("The matix must be square to find the determinant.")
        if self._numRows == 1:
            return self[0,0]
        QRdecomp = self.QRDecomposition()
        return QRdecomp[2]*product(*[QRdecomp[1][i,i] for i in range(self._numRows)])

    def rank(self, threshhold = .1**5):
        """Finds the rank of the matrix."""
        
        CurA = self.copy()
        stayFlag = True
        count = 0
        for j in range(min(self._numCols, self._numRows)):
            stayFlag = True
            skipFlag = False
            while(stayFlag):
                if count == CurA._numCols:
                    return j ####Here is the exit for small ranks!
                stayFlag=False
                alpha = (sum(CurA[i, 0]**2 for i in range(CurA._numRows))**.5)
                u = vector([0]*j+[CurA[i, 0]+(i==0)*alpha for i in range(CurA._numRows)])
                magU = sum(uIn**2 for uIn in u)**.5
                if magU<=threshhold:
                    skipFlag = True
                if skipFlag == False:
                    CurA.HouseholderRef(u[j:])
                skipFlag = True
##                v = u/magU
##                q = matrix(1, self._numRows, self._numRows)-v.outerProduct(v*2)
##                qCut = matrix([q[row, col] for col in range(j, q._numRows) for row in range(j, q._numRows)],
##                              CurA._numRows, CurA._numRows)
##                CurA = qCut*CurA
                if abs(CurA[0,0]) < threshhold:
                    if CurA._numCols == 1:
                        stayFlag = True
                        count += 1
                        continue
                    PushCol = CurA.popCol(0)
                    CurA.addCol(-1, PushCol)
                    count += 1
                    stayFlag = True
            if j < min(self._numCols, self._numRows)-1:
                CurA = matrix([CurA[j, i]  for j in range(1, CurA._numRows) for i in range(1, CurA._numCols)],
                              CurA._numRows-1, CurA._numCols-1)
        if abs(CurA[0,0]) < threshhold:
            return min(self._numCols, self._numRows)-1
        return min(self._numCols, self._numRows)

    def nullity(self):
        """Returns the nullity of the matrix"""
        return self._numCols - self.rank()

    def REF(self, threshhold = .1**5):
        """Returns a row echelon form of the matrix"""
        selfC = self.copy()
        pivotRow = 0
        pivotCol = 0
        while(pivotRow<self._numRows and pivotCol<self._numCols):
            if abs(selfC[pivotRow,pivotCol])<threshhold:
                ##can I swap?
                swapped = False
                for row in range(pivotRow+1, min(self._numCols, self._numRows)):
                    if selfC[row, pivotCol]>=threshhold:
                        selfC.swapRow(pivotRow, row)
                        swapped=True
                        break
                ##If I can't, move to the next column
                if not swapped:
                    pivotCol += 1
                    if pivotCol>=self._numRows:
                        break
            for col in range(pivotCol+1, self._numCols):
                selfC[pivotRow, col] = selfC[pivotRow, col]/float(selfC[pivotRow, pivotCol])
            selfC[pivotRow, pivotCol] = 1.0
            for row in range(pivotRow+1, self._numRows):
                for col in range(pivotCol+1, self._numCols):
                    selfC[row, col] = selfC[row, col]-(selfC[pivotRow, col]*selfC[row, pivotCol])
                selfC[row, pivotCol] = 0.0
            pivotRow+=1
            pivotCol+=1
        return selfC

    def RREF(self, threshhold = .1**5):
        """Returns the reduced row echelon form of the matrix."""
        mat = self.REF(threshhold)
        pivotRow = 0
        pivotCol = 0
        while(pivotRow<mat._numRows and pivotCol<mat._numCols):
            while(mat[pivotRow, pivotCol] != 1.0):
                pivotCol += 1
                if pivotCol>=mat._numCols:
                    return mat ##Exit
            for row in range(pivotRow):
                for col in range(pivotCol+1, mat._numCols):
                    mat[row, col] = mat[row, col] - mat[pivotRow, col]*mat[row, pivotCol]
                mat[row, pivotCol] = 0.0
            pivotRow += 1
        return mat
        
    def Eigen(self, itr = 100, confidence= 0.1**10, threshhold = 0.1**10, kerThreshhold = 0.1**10):
        """Finds the eigendecomposition (eigenvalues, eigenvectors) of the matrix.

            iterMax is the max iterations of the Schur decomposition algorithm, while the confidence is
            the confidence of zeros below the diagonal of the Schur form of the matrix.
            The threshhold is the threshhold for determining at what point two eigenvectors/eigenvalues
            are the same (to account for floating point errors). kerThreshhold is
            the threshhold variable in the kernel operation for finding eigenvectors.
            This algorithm only finds eigenvectors of geometric multiplicity 1.
            (this means solutions to ((M-lambda*I)**n)*v = 0 for n=1 and not for n > 1, where
            M is the matrix, v is a solution vector, lambda is the eigenvalue, and I is the identity matrix.)"""
        if not self.isSquare():
            raise MatrixError("The matix must be square to find the eigendecomposition.")
        kerThreshhold = kerThreshhold*(sum(abs(self[i, i]) for i in range(self._numRows)))/float(self.getNumRows())
##        selfHessen= self.Hessenburg()[0]
##        selfSchur = selfHessen.Schur(itr, confidence)
        if confidence == None:
            confidence = sys.maxint
        ##Finds the Charactoristic polynomial, then its roots
        selfCharPoly = self.CharPoly()
        selfCharPoly.reverse()
        eigenvalues = mpmath.polyroots(selfCharPoly)
##        eigenvalues = [selfSchur[0][i,i] for i in range(self._numRows)]
##        boxCheck = [abs(selfSchur[0][i+1, i])>confidence for i in range(self._numRows-1)]
##        ##This block is checking for 2x2 matrices along the diagonal
##        ##Then, it swaps those for the complex eigenvalues in "eigenvalues"
##        ##UhFlag is for when there are two "non-zero" entries in a row in the
##        ##subdiagonal after the specified iterations. This should never happen.
##        UhFlag = False
##        for bCi in range(len(boxCheck)):
##            if boxCheck[bCi] == True:
##                if UhFlag == True:
##                    raise EigenError("Failure to execute. matrix.Schur is broken.")
##                swapMat = selfSchur[0][bCi, bCi:bCi+2, bCi+2]
##                sMtrace = swapMat.trace()
##                sMdet = swapMat.det()
##                if (sMtrace**2-4*sMdet)<0:
##                    eigenvalues[bCi] = complex(sMtrace/2.0, (abs(sMtrace**2-4*sMdet)**.5)/2.0)
##                    eigenvalues[bCi+1] = complex(sMtrace/2.0, -(abs(sMtrace**2-4*sMdet)**.5/2.0))
##                else:
##                    eigenvalues[bCi] = sMtrace/2.0 + (sMtrace**2-4*sMdet)**.5/2.0
##                    eigenvalues[bCi+1] = sMtrace/2.0 - (sMtrace**2-4*sMdet)**.5/2.0
##                UhFlag = True
##            else:
##                UhFlag = False
        
        ##Now I'm clearing duplicates from "eigenvalues"
        rmvals = []
        for ev1 in range(len(eigenvalues)):
                for ev2 in range(ev1+1, len(eigenvalues)):
                    if abs(eigenvalues[ev1]-eigenvalues[ev2])<=threshhold:
                        rmvals.append(eigenvalues[ev1])
                        break
        for val in rmvals:
            eigenvalues.remove(val)
                        
        ##This block finds complex eigenvectors through block matricies.
        ##That is, for an eigenvalue = a+bi, it finds the kernel of the
        ##Following block matrix (where self = M, identity = I):
        ##[M-aI  bI  ]
        ##[-bI   M-aI]
        ##The kernel vector(s) will be [--v_1--,--v_2--], where v_1 + v_2i
        ##is an eigenvector
        eigenvectors = []
        for eival in eigenvalues:
            if any([type(eival) == i for i in (complex, mpmath.mpc)]):
                matMinusReal = self  - matrix(eival.real, self._numRows, self._numCols)
                imagI = matrix(eival.imag, self._numRows, self._numCols)
                I = matrix(1, 2, 2)
                skewI = matrix([0, 1, -1, 0], 2, 2)
                kernMat = I.KroneckerProd(matMinusReal) + skewI.KroneckerProd(imagI)
                if self.isMpf:
                    kernMat.makeMpf()
                eigenVec = kernMat.Ker(kerThreshhold)
                for iEiVec in range(len(eigenVec)):
                    eigenVec[iEiVec] = [complex(eigenVec[iEiVec][i],
                                                eigenVec[iEiVec][i+self._numRows]) for i in range(self._numRows)]
                ##This block is to remove vectors that are simply complex multiples of eachother
                ##This is done by assuring the vectors are close enough by log ratio
                ##as dictated by the threshhold parameter
                rmIndecies = set()
                for iEiVec1 in range(len(eigenVec)):
                    for iEiVec2 in range(iEiVec1+1, len(eigenVec)):
                        quotientVec = []
                        for i in range(len(eigenVec[iEiVec1])):
                            ##These 4 if statements are for the special
                            ##cases near 0, to account for floating point errors
                            ##This does use the same threshhold metric
                            ##as used for the log ratio. This may be
                            ##changed in the future.
                            if abs(eigenVec[iEiVec2][i]) <= threshhold:
                                if abs(eigenVec[iEiVec1][i])>threshhold:
                                    break
                                else:
                                    continue
                            if abs(eigenVec[iEiVec1][i]) <= threshhold:
                                if abs(eigenVec[iEiVec2][i])>threshhold:
                                    break
                                else:
                                    continue
                            quotientVec += [eigenVec[iEiVec1][i]/eigenVec[iEiVec2][i]]
                        if len(quotientVec) == 0:
                            avgQuo = 0
                        else:
                            avgQuo = sum(quotientVec)/float(len(quotientVec))
                        if all(min(abs(avgQuo-QuoCoeff), abs(avgQuo+QuoCoeff))<=threshhold for QuoCoeff in quotientVec):
                            rmIndecies.add(iEiVec2)
                rmVec = [eigenVec[i] for i in rmIndecies]
                for v in rmVec:
                    eigenVec.remove(v)
                eigenvectors.append(eigenVec)
            else:
                matMinusEigen = self - matrix(eival, self._numRows, self._numCols)
                if self.isMpf:
                    matMinusEigen.makeMpf()
                eigenvectors.append(matMinusEigen.Ker(kerThreshhold))
        return eigenvalues, eigenvectors

    def CharPoly(self):
        """Returns the characteristic polynomial of the matrix.
        [a_0, a_1, a_2, ..., a_(n-1), 1.0] implies
        a_0 + a_1*lambda + a_2*lambda**2+...+a_(n-1)*lambda**(n-1)+ lambda**n
        is the characteristic polynomial."""
        ##This algorithm uses Neuton Identities on the trace of the powers of the matrix to find the
        ##coefficients of the Charactoristic Polynomial. These traces will be the power sums
        ##of the eigenvalues

        
        if not self.isSquare():
            raise MatrixError("The matix must be square to find the charactoristic polynomial.")
        selfMpf = self.copy()
        selfMpf.makeMpf()
        CurA = selfMpf.copy()
        prevPrec = mpmath.mp.prec
        mpmath.mp.prec = max(self.getNumRows()*5, mpmath.mp.prec)
        pSumList = []
##        CurA = CurA.Hessenburg(retQ=False)
        for i in range(self._numRows):
            pSumList.append(CurA.trace())
            CurA=CurA*selfMpf
        pSumList.append(CurA.trace())
        eList = []
        for i in range(len(pSumList)-1):
            sign = 1
            accu = 0
            for eIndex in range(len(eList)):
                accu+=sign*pSumList[eIndex]*eList[len(eList)-1-eIndex]
                sign*=-1
            accu+=sign*pSumList[len(eList)]
            eList.append(accu/mpmath.mpf(i+1))
        for i in range(len(eList)): 
            eList[i]*= (-1)**(i+1)
        eList.reverse()
        eList.append(mpmath.mpf(1.0))
        mpmath.mp.prec = prevPrec
        return eList
                
            
    def SatPoly(self, threshhold = 0.1**5, kerThreshhold = 0.1**5):
        """Finds the coefficients of the minimal (monic) polynomial the matrix satisfies.
            Coefficients are ordered from least significant to most.
            [a_0, a_1, a_2, ..., a_(n-1), a_n] implies
            a_0 + a_1*lambda + a_2*lambda**2+...+a_(n-1)*lambda**(n-1)+ a_n*lambda**n
            is the satisfactory polynomial."""
        if not self.isSquare():
            raise MatrixError("The matix must be square to find the satisfactory polynomial.")
        matrixPowArray = [matrix(1, self._numRows, self._numCols)]
        m = self.copy()
        for i in range(self._numRows):
            matrixPowArray.append(m)
            m = m*self
        powMatrix = matrix(0, self._numRows+1, self._numRows**2)
        ##This unpacks the array into powMatrix,
        ##where each row corresponds to a matrix as a vector
        for matI in range(self._numRows+1):
            for i in range(self._numRows**2):
                powMatrix[matI, i] = matrixPowArray[matI][i]
##        imgPowMat = [vector(i) for i in powMatrix.Img()]
##        idVec = vector([0]*self._numRows + [1])
##        ##The returned value is the a kernel vector of
##        ##the transpose of powMatrix with a non-zero
##        ##last coefficient. (This will exist)
##        return list(idVec - sum(v.proj(idVec) for v in imgPowMat))
        powMatrix.t()
        kernel = powMatrix.Ker(kerThreshhold)
        if len(kernel) == 0:
            raise MatrixError("Unable to find a satisfactory polynomial.")
        kernel[0].reverse()
        kernMat = matrix(kernel[0], 1, len(kernel[0]))
        for i in range(1, len(kernel)):
            kernel[i].reverse()
            kernMat.addRow(-1, kernel[i])
        retMat = kernMat.REF(threshhold)
        retPoly = retMat.getRow(retMat._numRows-1)
        retPoly.reverse()
        return retPoly
        

    def SVD(self):
        """Finds the singular values, and associated U and V. (M = USigmaV*)"""
        A = [(self*self.tC()).Hessenburg(), (self.tC()*self).Hessenburg()]
        #the two matrices in A will be symmetric tridiagonal.
        U = [A[0][1], A[1][1]]#Both of these will be square
        A = [A[0][0], A[1][0]]#just getting the matricies
        
        singularVals = [None, None]
        
        for i in range(2):
            ##Recursively finding the charactoristic polynomial
            poly_a = [1]
            poly_b = [A[i][0,0], -1]
            alpha = [A[i][k, k] for k in range(1, A[i].dims()[0])]
            beta = [A[i][k, k+1]**2 for k in range(A[i].dims()[0]-1)]
            for row in range(A[i].dims()[0]-1):
                poly_c = poly_b
                poly_a = poly_a + [0]
                poly_b = [alpha[row]*poly_b[0] - beta[row]*poly_a[0]]+\
                        [alpha[row]*poly_b[k]-poly_b[k-1] - beta[row]*poly_a[k]\
                        for k in range(1, len(poly_b))]
                poly_a = poly_c
            
            charPoly = poly_b
            
            charPoly.reverse()
            singularVals = mpmath.polyroots(charPoly)
            kernVecs = []
            for sV in singularVals:
                lamMat = matrix(sV, A[i].dims()[0], A[i].dims()[1], True)
                kernMat = A[i]-lamMat
                kernVecs.append(kernMat.Ker())
            Udat = [kernVecs[i][j][k] for i in range(len(kernVecs))\
                for j in range(len(kernVecs[i])) for k in range(len(kernVecs[i][j]))]
            U[i] = matrix(Udat, self.dims()[i], self.dims()[i])
            
        sigMat = U[0].tC()*self*U[1]
        sigVals = [sigMat[k, k] for k in range(min(self.dims()))]
        return sigVals, U[0], U[1]
        
#        for i in range(2):
#            temp = A[i].Schur(10)
#            U[i] = U[i]*temp[1]
#        sigma = U[1].tC()*self*U[0]
#        singVals = [sigma[i,i] for i in range(min(sigma.dims()))]
#        return singVals, U[1], U[0].tC()
        
        
#        A = (self.tC()*self).Eigen()
#        B = (self*self.tC()).Eigen()
#        Vdat = [A[1][i][j][k] for i in range(len(A[1]))\
#                for j in range(len(A[1][i])) for k in range(len(A[1][i][j]))]
#        Udat = [B[1][i][j][k] for i in range(len(B[1]))\
#                for j in range(len(B[1][i])) for k in range(len(B[1][i][j]))]
#        V = matrix(Vdat, self._numCols, self._numCols)
#        U = matrix(Udat, self._numRows, self._numRows)
#        if len(A[0])<=len(B[0]):
#            singularVals = [i**.5 for i in A[0]]
#        else:
#            singularVals = [i**.5 for i in B[0]]
#        return singularVals, U.tC(), V.tC()

    def pseudoinverse(self, threshold = 10**-5):
        sing, U, V = self.SVD()
        for i in range(len(sing)):
            if sing[i]>threshold:
                sing[i] = 1/sing[i]
            else:
                sing[i] = 0
        sigma = makeDiagMatrix(sing, (V.getNumCols(), U.getNumCols()))
        return V*sigma*U.tC()

    def orthogonalProjector(self, option = 'rows'):
        """Returns an orthogonal projector matrix.
            Options:'rows', 'cols', 'ker', 'tker'"""
        
        pseudoInv = self.pseudoinverse()
        if option == 'rows':
            ret = self*pseudoInv
        elif option == 'cols':
            ret = pseudoInv*self
        elif option == 'tker':
            P = self*pseudoInv
            I = makeDiagMatrix([1]*P.getNumRows())
            ret = I-P
        elif option == 'ker':
            Q = pseudoInv*self
            I = makeDiagMatrix([1]*Q.getNumRows())
            ret = I-Q
        return ret
    
    def QRDecomposition(self):
        """Returns the QR Decomposition of the matrix (with det of Q)"""

##        if not self.isSquare():
##            raise StandardError("The matix must be square to QR decompose.")
        detQ = 1
        CurA = self.copy()
        retA = self.copy()
#        QList = []
        Q = matrix(1, self._numRows, self._numRows)
        for j in range(min(self._numCols, self._numRows)):
            alpha = (sum(CurA[i, 0]**2 for i in range(CurA._numRows))**.5)
            u = vector([0]*j+[CurA[i, 0]+(i==0)*alpha for i in range(CurA._numRows)])
            magU = sum(uIn**2 for uIn in u)**.5
            if magU==0:
                if j<min(self._numCols, self._numRows)-1:
                    CurA = matrix([CurA[j, i]  for j in range(1, CurA._numRows) for i in range(1, CurA._numCols)],
                                  CurA._numRows-1, CurA._numCols-1)
                continue
            Q.HouseholderRef(u)
            detQ *= -1
            CurA.HouseholderRef(u[j:])
            retA.HouseholderRef(u)
##            v = u/magU
##            q = matrix(1, self._numRows, self._numRows)-v.outerProduct(v*2)
##            qCut = matrix([q[row, col] for col in range(j, q._numRows) for row in range(j, q._numRows)],
##                          CurA._numRows, CurA._numRows)
##            QList.append(q)
##            CurA = qCut*CurA
            if j<min(self._numCols, self._numRows)-1:
                CurA = matrix([CurA[j, i]  for j in range(1, CurA._numRows) for i in range(1, CurA._numCols)],
                              CurA._numRows-1, CurA._numCols-1)
##        QList.reverse()
##        Q = product(*[q for q in QList])
        return Q.tC(), retA, detQ
                    
##        AList = [vector(self.getCol(i)) for i in range(self._numCols)]
##        UList = [AList[0]]
##        try:
##            EList = [UList[0]/(UList[0]).norm()]
##        except ZeroDivisionError:
##            raise StandardError("Must be an invertable matrix.")
##        for a in AList[1:]:
##            UList.append(a - sum(e.proj(a) for e in EList))
##            try:
##                EList.append(UList[-1]/(UList[-1]).norm())
##            except ZeroDivisionError:
##                raise StandardError("Must be an invertable matrix.")
##        Qmatrix = matrix(rows = len(AList), cols = len(AList))##transforming into Q matrix
##        for e_index in range(len(EList)):
##            for row in range(len(EList)):
##                Qmatrix[row, e_index] = (EList[e_index])[row]
##        Rmatrix = matrix(rows = len(AList), cols = len(AList))##transforming into R matrix
##        for e_index_row in range(len(AList)):
##            for a_index_col in range(len(AList)):
##                if a_index_col >= e_index_row:
##                    Rmatrix[e_index_row, a_index_col] = innerProduct(AList[a_index_col], EList[e_index_row])
##        if (vector(((Qmatrix*Rmatrix)-self)._data)).norm() > confidence:
##            if haltIfInaccurate:
##                raise Warning("Not an accurate QR decomposition. (the matrix may not me invertable)")
##        return Qmatrix, Rmatrix

    def Schur(self, itr = 200, confidence=0.1**10):##Crude Implementation
        """Returns the Schur form of the matrix, with the associated orthogonal matrix.
        itr is the number of iterations per loop, while confidence is the confidence
        that there are no two non-zero entries in a row in the sub-diagonal.
        (If confidence == None, this will be ignored)
        It is recommended that the matrix is put into upper Hessenburg form before
        executing this function.
        
        The algorithm implemented is the QR-step algorithm.
        This algorithm works with most matricies, however it will find strange results with
        matricies with complex eigenvalues (specifically, 2x2 blocks along the diagonal),
        and this algorithm will not converge for matricies with eigenbasis V s.t.
        there is a leading principle minor of V**-1 that is zero."""
        prevPrec = mpmath.mp.prec
        if confidence == None:
            confidence = sys.maxint
        if (self.isSquare() == False):
            raise MatrixError("The matix must be square.")
        mpmath.mp.prec = max(int(mpmath.log10(confidence)**2), mpmath.mp.prec)
        A = self.copy()
        A.makeMpf()
##        prevA = matrix(self._numRows, self._numCols)
        Q = matrix(1, self.dims()[0], self.dims()[1])
        Q.makeMpf()
        #for i in range(iterations):
        orgItr = itr
        confNotSat = True
        while(confNotSat):
            itr = orgItr
            while itr > 0:
                ##the error bound is simply looking to see how close the matrix is to an upper triangular
                QR = A.QRDecomposition()
#                prevA = A.copy()
                A = QR[1]*QR[0]
                Q = Q*QR[0]
                itr-=1
            boxCheck = [abs(A[i+1, i])>confidence for i in range(self._numRows-1)]
            confNotSat = any([boxCheck[i]*boxCheck[i+1] for i in range(len(boxCheck)-1)])
        mpmath.mp.prec = prevPrec
        return A, Q

    def Hessenburg(self, retQ = True):
        """Returns the upper Hessenburg form of the matrix (with Q)
        This procedure is taken from the book: Numerical Analysis, Burden and Faires, 8th Edition
        via wikipedia. However, the code itself is original."""
        if (self.isSquare() == False):
            raise MatrixError("The matix must be square.")
        A = self.copy()
        PList = []
        for k in range(self._numRows-1):
            a = ((A[k+1, k]<0)*2-1)*(sum([A[j, k]**2 for j in range(k+1, A._numRows)])**.5)
            r = ((a**2 - A[k+1, k]*a)/2.0)**.5
            if r == 0:
                continue
            v = vector([0]*A._numRows)
            v[k+1] = (A[k+1, k] - a)/(2*r)
            for i in range(k+2, A._numRows):
                v[i] = A[i, k]/(2.0*r)
            PList.append(matrix(1, A._numRows, A._numCols) - 2*vector.outerProduct(v, v))
            A = PList[-1]*A*PList[-1]
        ###Now to clear everything below the subdiagonal
        for i in range(A._numRows):
            for j in range(A._numCols):
                if j+1<i:
                    if A._data.has_key((i, j)):
                        del A._data[i, j]
        if retQ:
            return A, product(*PList)
        else:
            return A

    def HouseholderRef(self, v):
        """Efficiently performs a Householder Reflection on the matrix (in place) with a vector v."""
        if len(v) != self._numRows:
            raise MatrixError("v must have the same number of entries as the matrix does rows.")
        isMpf = self.isMpf
        v = vector(v)
        v = v/v.norm()
##        self.t()
        if not isMpf:
            for col in range(self._numCols):
                curCol = self.getCol(col)
                twoTimesInProd = 2*sum(v[i]*curCol[i] for i in range(self._numRows))
                for row in range(self._numRows):
                    self[row, col] = self[row, col]-twoTimesInProd*v[row]
        else:
            for col in range(self._numCols):
                curCol = self.getCol(col)
                twoTimesInProd = 2*mpmath.fdot((v[i],curCol[i]) for i in range(self._numRows))
                for row in range(self._numRows):
                    self[row, col] = self[row, col]-twoTimesInProd*v[row]
##        self.t()
        
class modMatrix(matrix):
    def __init__(self, data=None, rows=2, cols=2, mod=2):
        self._numRows = rows
        self._numCols = cols
        self._mod = mod
        self._printer = printer()
        if data == None:
            self._data = [0]*rows*cols
        elif len(data)>rows*cols:
            return BufferError("data must be able to fit in the matrix")
        else:
            self._data = data + [0]*(rows*cols-len(data))
        self._printer(self, doPrint = False)
        self%=mod
        
    def copy(self):
        return modMatrix([i for i in self._data], self.dims()[0], self.dims()[1], self._mod)
        
    def __eq__(self, other):
        if isinstance(other, modMatrix):
            return all([self._data == other._data, self._numRows == other._numRows, self._numCols == self._numCols, self._mod == other._mod])
        else: return False

    def __add__(self, other):
        if not isinstance(other, modMatrix):
            raise NotImplemented
        if self._mod != other._mod:
            raise NotImplemented
        if (self._numRows != other._numRows)|(self._numCols != other._numCols):
            raise MatrixError("must be the same size of matrix")
        return modMatrix(rows=self._numRows, cols=self._numCols, data=[self[i]+other[i] for i in range(len(self))], mod=self._mod)

    def __sub__(self, other):
        if not isinstance(other, modMatrix):
            raise NotImplemented
        if self._mod != other._mod:
            raise NotImplemented
        if (self._numRows != other._numRows)|(self._numCols != other._numCols):
            raise MatrixError("must be the same size of matrix")
        return modMatrix(rows=self._numRows, cols=self._numCols, data=[self[i]-other[i] for i in range(len(self))], mod=self._mod)
        
    def __mul__(self, other):
        if not isinstance(other, matrix):
            raise NotImplemented
        if self._numCols != other._numRows:
            raise MatrixError(
                "in m*n, the number of rows of n must match the number of columns of m")
        new_matrix = modMatrix(None, self._numRows, other._numCols, mod=self._mod)
        for i in range(new_matrix._numRows):
            for j in range(new_matrix._numCols):
                new_matrix[i, j] = sum(self[i, k]*other[k, j] for k in range(self._numCols))%self._mod
        return new_matrix
    
    def __pow__(self, other):
        if not self.isSquare():
            raise MatrixError(
                "must be a square matrix")
        k = bin(abs(other))
        n = self
        ptr = 1
        ret = modMatrix(rows=self._numRows, cols=self._numCols, mod=self._mod)##Will make the unit matrix
        for i in range(self._numRows): 
            for j in range(self._numCols):
                if i==j:
                    ret[i,j] = 1
        while(k[-ptr] != 'b'):
            if(n == 1):
                    return ret 	#exit
            if(k[-ptr] == '1'):
                    ret = (ret*n)
            n = (n*n)
            ptr = ptr + 1
        if other<0:
            return ret.inv()
        else:
            return ret

    def det(self):
        return matrix.det(self)%self._mod

class printer:
    ##NOTES:
    ##Add options
    ##Make better default
    _print = None
    _reference = None ##THIS IS THE ACTUAL OBJECT, TAKE CARE OF IT!!
    def __init__(self, *args, **tags):
        self._data = dict(*args, **tags)

    def __call__(self, newJob = None, doPrint = True):
        if newJob != None:
            self._reference = newJob
        if isinstance(self._reference, matrix):
            self._format_matrixPrint()
        else:
            self._print = newJob.__str__()
        if doPrint:
            print self._print

    def _format_matrixPrint(self):
        matrix = self._reference
        colWidths = [max(len(str(matrix[i, j])) for i in range(matrix.getNumRows())) for j in range(matrix.getNumCols())] #will contain width of each col
        ##self._print = u"\u2308"+ " "*(sum(colWidths)+(len(colWidths)-1)*3) + u"\u2309" + "\n"
        self._print = ""
        space = " "
        NumRows = matrix.getNumRows()
        NumCols = matrix.getNumCols()
        for row in range(NumRows):
            if NumRows == 1:
                self._print = "[" + space.join(str(matrix[row, j])+space*(colWidths[j]-len(str(matrix[row, j]))) for j in range(NumCols)) + " ]"
                break
            if row == 0:
                self._print += u"\u2308" + space.join((matrix[row, j]).__str__()+space*(colWidths[j]-len((matrix[row, j]).__str__())) for j in range(NumCols)) + u"\u2309" + "\n"
            if row != 0 and row != NumRows-1:
                self._print += u"\u2223" + u"\u2009"+ u"\u200A"+ space.join((matrix[row, j]).__str__()+space*(colWidths[j]-len((matrix[row, j]).__str__())) for j in range(NumCols)) + u"\u2009" + u"\u200A" + u"\u2223" + "\n"
            if row == NumRows-1:
                self._print += u"\u230A" + space.join((matrix[row, j]).__str__()+space*(colWidths[j]-len((matrix[row, j]).__str__())) for j in range(NumCols)) + u"\u230B"
        ##self._print += u"\u230A" + space*(sum(colWidths)+(len(colWidths)-1)*3) + u"\u230B"

    def getPrintStr(self):
        return self._print

class keyValuedPrinter(printer):
    
    def _format_matrixPrint(self):
        matrix = self._reference
        keysInRows, keysInCols = matrix.getKeysInRowsAndCols()
        iterRows, iterCols = [i for i in keysInRows.keys()], [i for i in keysInCols.keys()]
        for i in range(matrix.getNumRows()):
            try:
                k = iterRows.index(i)
            except:
                iterRows.append(i)
        for i in range(matrix.getNumCols()):
            try:
                k = iterCols.index(i)
            except:
                iterCols.append(i)

        iterRows.sort()
        iterCols.sort()
        colWidths = [max([len(str(matrix[i, j])) for i in iterRows]+[len(str(j))]) for j in iterCols] #will contain width of each col
        try:
            maxRowKeyWidth = max(len(str(i)) for i in iterRows)
        except ValueError:
            maxRowKeyWidth = 0
        ##self._print = u"\u2308"+ " "*(sum(colWidths)+(len(colWidths)-1)*3) + u"\u2309" + "\n"
        self._print = ""
        space = " "
#        NumRows = matrix.getNumRows()
#        NumCols = matrix.getNumCols()
        self._print = space*(maxRowKeyWidth+1) + space.join(str(iterCols[j])+space*(colWidths[j]-len(str(iterCols[j]))) for j in range(len(iterCols))) + space + "\n"
        for row in range(len(iterRows)):
            if len(iterRows) == 1:
                self._print += str(iterRows[row]) + "[" + space.join(str(matrix[iterRows[row], iterCols[j]])+space*(colWidths[j]-len(str(matrix[iterRows[row], iterCols[j]])))\
                                                          for j in range(len(iterCols))) + " ]"
                break
            if row == 0:
                self._print += str(iterRows[row])+space*(maxRowKeyWidth-len(str(iterRows[row])))+ u"\u2308" +\
                               space.join((matrix[iterRows[row], iterCols[j]]).__str__()+space*(colWidths[j]-len((matrix[iterRows[row], iterCols[j]]).__str__()))\
                                                      for j in range(len(iterCols))) + u"\u2309" + "\n"
            if row != 0 and row != len(iterRows)-1:
                self._print += str(iterRows[row])+space*(maxRowKeyWidth-len(str(iterRows[row])))+ u"\u2223" + u"\u2009"+ u"\u200A"+ space.join((matrix[iterRows[row], iterCols[j]]).__str__()+\
                                                                            space*(colWidths[j]-len((matrix[iterRows[row], iterCols[j]]).__str__()))\
                                                                            for j in range(len(iterCols))) + u"\u2009" + u"\u200A" + u"\u2223" + "\n"
            if row == len(iterRows)-1:
                self._print += str(iterRows[row])+space*(maxRowKeyWidth-len(str(iterRows[row])))+ u"\u230A" + space.join((matrix[iterRows[row], iterCols[j]]).__str__()+\
                                                                          space*(colWidths[j]-len((matrix[iterRows[row], iterCols[j]]).__str__()))\
                                                                          for j in range(len(iterCols))) + u"\u230B"
        ##self._print += u"\u230A" + space*(sum(colWidths)+(len(colWidths)-1)*3) + u"\u230B"


    
class matrixEditor:
    ##NOTES:
    ##Implement GUI shit
    ##May be useful with printer
    def __init__(self, matrix):
        self._matrix = matrix
        
    def refresh(self):
        raise NotImplemented
    
    def loadMatrix(self, Matrix):
        self._matrix = matrix
        self.refresh()

class vector(list):
    ##NOTES:
    ##Finish methods
    ##Make sure to return vectors where normally lists are returned.
    def __init__(self, inList):
        self.extend(i for i in inList)

    def __mul__(self, other):
        if isinstance(other, vector):
            return vector(self[i]*other[i] for i in range(min(len(self), len(other))))
        elif any([isinstance(other, i) for i in (int, float, long, complex, mpmath.mpf)]):
            return vector(self[i]*other for i in range(len(self)))
        else:
            raise NotImplemented

    def __sub__(self, other):
        if isinstance(other, vector):
            return vector(self[i]-other[i] for i in range(min(len(self), len(other))))
        elif any([isinstance(other, i) for i in (int, float, long, complex, mpmath.mpf)]):
            return vector(self[i]-other for i in range(len(self)))
        else:
            raise NotImplemented

    def __add__(self, other):
        if isinstance(other, vector):
            return vector(self[i]+other[i] for i in range(min(len(self), len(other))))
        elif any([isinstance(other, i) for i in (int, float, long, complex, mpmath.mpf)]):
            return vector(self[i]+other for i in range(len(self)))
        else:
            raise NotImplemented

    def __radd__(self, other):
        if isinstance(other, vector):
            return vector(other[i]+self[i] for i in range(min(len(self), len(other))))
        elif any([isinstance(other, i) for i in (int, float, long, complex, mpmath.mpf)]):
            return vector(other+self[i] for i in range(len(self)))
        else:
            raise NotImplemented

    def __div__(self, other):

        if isinstance(other, vector):
            return vector(self[i]/other[i] for i in range(min(len(self), len(other))))
        elif any([isinstance(other, i) for i in (int, float, long, complex, mpmath.mpf)]):
            return vector(self[i]/other for i in range(len(self)))
        else:
            raise NotImplemented

    def norm(self):
        return vector.innerProduct(self,self)**.5

    def proj(self, other):
        return  self*(vector.innerProduct(self, other)/float(vector.innerProduct(self, self)))
        
    def innerProduct(self, other):
        if isinstance(other, vector):
            return sum(self*other.conjugate())
        else:
             raise NotImplemented
            
    def conjugate(self):
        return vector([i.conjugate() for i in self])
    
    def outerProduct(self, other):
        if isinstance(other, vector):
            mdata = [i*j for i in self for j in other]
            return matrix(mdata, len(self), len(other))

def makeDiagMatrix(diagonal, dims=None):
    if dims == None:
        retMat = matrix(rows = len(diagonal), cols = len(diagonal))
    else:
        if min(dims) != len(diagonal):
            raise IndexError("the min dim must be the same size as the diagonal!")
        retMat = matrix(rows = dims[0], cols = dims[1])
    for i in range(len(diagonal)):
        retMat[i, i] = diagonal[i]
    return retMat

def product(*numList):
    """Returns the product of all the arguements."""
    if len(numList)==0:
        return 1
    ret = [n for n in numList]
    while(len(ret) != 1):
        iterList = (2*i for i in range((len(ret))/2))
        if len(ret)%2 == 0:
            ret = [ret[i]*ret[i+1] for i in iterList]
        else:
            temp = ret[-1]
            ret = [ret[i]*ret[i+1] for i in iterList]
            ret.append(temp)
    return ret[0]

def randomMatrix(rows, cols, low=0, high=1):
    """Returns a matrix with random coefficients from [low, high)."""
    data = [random.random()*(high-low) + low for i in range(rows*cols)]
    return matrix(data, rows, cols)

def randomEigenMatrix(eigenvalues, low=-1, high=1):
    """Returns a random square matrix with the specified (complex) eigenvalues.
    low and high are the bounds on the uniform coefficients of the eigenvectors.
    (note, with low=-1 and high = 1, this effectively covers every
    possible eigenvector, although not with equal probability.)"""
    eiVals = []
    for ev in eigenvalues:
        if type(ev)==complex:
            if eiVals.count((ev.real, abs(ev.imag)))==0:
                eiVals.append((ev.real, abs(ev.imag)))
        elif any([type(ev)==i for i in [float, int, long]]):
            eiVals.append(ev)
    size = sum((type(i)==tuple)+1 for i in eiVals)
    eigenMat = matrix(0, size, size)
    eiI = 0
    diaI = 0
    while(eiI<len(eiVals)):
        if any([type(eiVals[eiI])==i for i in [float, int, long]]):
            eigenMat[diaI, diaI] = eiVals[eiI]
            diaI+=1
            eiI+=1
        elif type(eiVals[eiI])==tuple:
            eigenMat[diaI, diaI] = eiVals[eiI][0]
            eigenMat[diaI+1, diaI+1] = eiVals[eiI][0]
            eigenMat[diaI+1, diaI] = eiVals[eiI][1]
            eigenMat[diaI, diaI+1] = -eiVals[eiI][1]
            diaI+=2
            eiI+=1
    randMat = randomMatrix(size, size, low, high)
    while(randMat.det()== 0): ##makes sure it's invertable
        randMat = randomMatrix(size, size, low, high)
    return randMat*eigenMat*randMat.inv()

def randomUnitVector(size):
    """Returns a random, uniformly distributed unit vector."""
    randVec=[(random.random()*2)-1]
    for i in range(size-1):
        mag = sum([j**2 for j in randVec])
        randVec.append(((random.random()*2)-1)*(1-mag)**.5)
    retvec = vector(randVec)
    if retvec.norm()==0:
        return randomUnitVector(size)
    return retvec/retvec.norm()

def randomUnitaryMatrix(size):
    """Returns a random Unitarty Matrix."""
    if size>30:
        print "Currently, the algorthm is unstable for size over"
        print "30."
    unitVectors = [list(randomUnitVector(i+1)) for i in range(size)]
#    unitaryVecs = []
    i = 1
    curMat = matrix(unitVectors[0], 1, 1)
    while(i<size):
        ###I will first find an orthonormal basis
        ###for the kernel of the next vector= v,
        ###then they will be put in the columns of
        ###a matrix. I multiply "curMat"
        ###by this matrix, and then put v
        ###as a new column of the resulting matrix.

        kernel = (matrix(unitVectors[i], 1, i+1)).Ker()
        kernDat = [k for vec in kernel for k in vec]
        kernMat = matrix(kernDat, len(kernel),len(kernel[0]))
        kernMat.t()
        curMat = kernMat*curMat
        curMat.addCol(0, unitVectors[i])
        i+=1
    return curMat 

def changeMpfDecimalPrecision(newDps):
    mpmath.mp.dps = newDps

def fastSum(numList, makeMpf=True, usePySum = False):
    """An edited form of mpmath.mpf_sum."""
    if usePySum:
        return sum(numList)
    man = 0
    exp = 0
    prec = 0
    rnd = mpmath.libmp.round_down
    absolute = False
    max_extra_prec = mpmath.mp.prec*2 or 1000000  # XXX
    special = None
#    flagMpf = False
    for x in numList:
        if any([type(x) == i for i in (float, int, long)]):
            ##may throw an exception in python 2.5
            xman, xexp,  = math.frexp(x)
            xsign, xman = xman<0, abs(xman)
            xbc = 53
        elif type(x) == mpmath.mpf:
            xsign, xman, xexp, xbc = x._mpf_
        else:
            raise TypeError("Wrong type to fastSum.")
        if xman:
            if xsign and not absolute:
                xman = -xman
            delta = xexp - exp
            if xexp >= exp:
                # x much larger than existing sum?
                # first: quick test
                if (delta > max_extra_prec) and \
                    ((not man) or delta-mpmath.libmp.bitcount(abs(man)) > max_extra_prec):
                    man = xman
                    exp = xexp
                else:
                    man += (xman*(2**delta))
            else:
                delta = -delta
                # x much smaller than existing sum?
                if delta-xbc > max_extra_prec:
                    if not man:
                        man, exp = xman, xexp
                else:
                    man = (man*(2**delta)) + xman
                    exp = xexp
        elif xexp:
            if absolute:
                x = mpmath.libmp.mpf_abs(x)
            special = mpmath.libmp.mpf_add(special or mpmath.libmp.fzero,\
                                           mpmath.libmp.from_man_exp(man, exp, xbc), 1)
    # Will be inf or nan
    if special:
        return mpmath.mpf(special)
    if not makeMpf:
        return man*(2**exp)
    else:
        return mpmath.mpf(mpmath.libmp.from_man_exp(int(man*(1<<53)), exp-53, prec, rnd))

def fastDot(numPairsList, makeMpf = False):
    if not makeMpf:
        return sum(k[0]*k[1] for k in numPairsList)

    man = 0
    exp = 0
    prec = mpmath.mp.prec
    rnd = mpmath.libmp.round_down
    absolute = False
    max_extra_prec = mpmath.mp.prec*2 or 1000000  # XXX
    special = None
    flagMpf = False
    for x1, x2 in numPairsList:
        if any([type(x1) == i for i in (float, int, long)]):
            ##may throw an exception in python 2.5
            x1man, x1exp,  = math.frexp(x1)
            x1sign, x1man = x1man<0, abs(x1man)
            x1bc=53
        elif type(x1) == mpmath.mpf:
            x1sign, x1man, x1exp, x1bc = x1._mpf_
        else:
            raise TypeError("Wrong type to fastDot.")
        if any([type(x2) == i for i in (float, int, long)]):
            ##may throw an exception in python 2.5
            x2man, x2exp  = math.frexp(x2)
            x2sign, x2man = x2man<0, abs(x2man)
            x2bc=53
        elif type(x2) == mpmath.mpf:
            x2sign, x2man, x2exp, x2bc = x2._mpf_
        else:
            raise TypeError("Wrong type to fastDot.")
        if x1sign:
                x1man = -x1man
        if x2sign:
                x2man = -x2man
##        xsign, xman, xexp, xbc = \
##                mpmath.libmp.mpf_mul(mpmath.libmp.from_man_exp(x1man, x1exp, x1bc),\
##                                    mpmath.libmp.from_man_exp(x2man, x2exp, x2bc))
        xsign = x1sign^x2sign
        xman = x1man*x2man
        xexp = x1exp+x2exp
        xbc = max(x1bc, x2bc)
        if xman:
            if xsign and not absolute:
                xman = -xman
            delta = xexp - exp
            if xexp >= exp:
                # x much larger than existing sum?
                # first: quick test
                if (delta > max_extra_prec) and \
                    ((not man) or delta-mpmath.libmp.bitcount(abs(man)) > max_extra_prec):
                    man = xman
                    exp = xexp
                else:
                    man += (xman*(2**delta))
            else:
                delta = -delta
                # x much smaller than existing sum?
                if delta-xbc > max_extra_prec:
                    if not man:
                        man, exp = xman, xexp
                else:
                    man = (man*(2**delta)) + xman
                    exp = xexp
        elif xexp:
            if absolute:
                x = mpmath.libmp.mpf_abs(x)
            special = mpmath.libmp.mpf_add(special or mpmath.libmp.fzero,\
                                            mpmath.libmp.from_man_exp(man, exp, xbc), 1)
    # Will be inf or nan
    if special:
        return mpmath.mpf(special)
    return mpmath.mpf(mpmath.libmp.from_man_exp(long(man*(1<<53)), exp-53, prec, rnd))
    
