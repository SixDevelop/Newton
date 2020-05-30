using System;

namespace Newton{
    public class Matrix{
        private Random rnd = new Random();
        public double eps;
        public double[,] _mat;              
        public double[,] matL,matU;      
        public double[,] matP,matQ;      
        public double[,] Q, R;
        public int colSize,rowSize;           
        public bool isLU;                       
        public bool isQR; 
        public bool isSing;
        public int size;                 
        private int _rank;                      
        public double norm;
        public double _determinant;   
        public int jacobyInter;
        public int seidelIter;
        public bool isSquare;
        public double[,] revMatrix;
        public double condNumber;
        public int numOfOperations;

        public Matrix(double[,] matrix){

            colSize = matrix.GetLength(1);
            rowSize = matrix.GetLength(0);

            _mat = new double[rowSize, colSize];

            for(int i = 0; i < rowSize; i++)
            for(int j = 0; j < colSize; j++)
                _mat[i,j] = matrix[i,j];
            isLU = false;

            matL = new double[rowSize, rowSize];
            matU = new double[rowSize, colSize];
            matP = new double[rowSize, rowSize];
            matQ = new double[colSize, colSize];

            isSquare = false;

            if(colSize > rowSize)
                _rank = rowSize;
            else if(colSize < rowSize) 
                _rank = colSize;
            else {
                isSquare = true;
                size = colSize;
                Reverse();
                norm = Norm(_mat);
                condNumber = Norm(revMatrix) * Norm(_mat);
                _rank = size;
            }
            
            for(int i =0; i < rowSize; i++)
                matP[i,i] = 1;
            for(int j = 0; j < colSize; j++)
                matQ[j,j] = 1; 
            norm = Norm(_mat);           
            eps = 1e-14 * norm;
        }  

        public int GetRank() {
            if(isLU)
                return _rank;
            else {
                LUDecompose();
                return _rank;
            }
        }
        public double GetDeterminant() {
            if(isLU)
                return _determinant;
            else {
                LUDecompose();
                return _determinant;
            }
        }
        
        public void Reverse() {
            if(!isLU)
                LUDecompose();
            
            double[] vec = new double[size];
            revMatrix = new double[size, size];
            for(int k = 0; k < size; k++) {
                vec[k] = 1;
                vec = multiplyVec(matP, vec, size);

                for(int j = 0; j < size; j++) {
                    for(int i = j + 1; i < size; i++)
                        vec[i] -= matL[i,j] * vec[j];
                }
                for(int j = size - 1; j > -1; j--) {
                    vec[j] /= matU[j,j];
                    for(int i = 0; i < j; i++)
                        vec[i] -= matU[i,j] * vec[j];
                }
                for(int i = 0; i < size; i++)
                    revMatrix[i,k] = vec[i];
                for(int i = 0; i < size; i++)
                    vec[i] = 0;
            }

        }
        public double Norm(double[,] matrix) {
            double current;
            double _norm = 0;
            for(int i = 0; i < matrix.GetLength(0); i++) {
                current = 0;
                for(int j = 0; j < matrix.GetLength(1); j++) {
                    current += Math.Abs(matrix[i,j]);
                }
                if(current > _norm)
                    _norm = current;
            }
            return _norm;
        } 
        public void LUDecompose() {
            if(isSquare) {
                numOfOperations = 0;
                isLU = true;
                isSing = true;
                matL = new double[size, size];
                matU = new double[size, size];
                matP = new double[size, size];
                matQ = new double[size, size];

                for(int i = 0; i < size; i++) {
                    for(int j = 0; j < size; j++) {
                        matU[i,j] = _mat[i,j];
                    }
                    matP[i,i] = 1;
                    matQ[i,i] = 1;
                }

                int line1,line2;
                double leadElem;
                int numSwap = 0;

                for(int i = 0; i < size; i++) {
                    if(matU[i,i] == 0) {   
                        isSing = false;
                        for(int j = i + 1; j < size; j++) {
                            if(matU[j,i] != 0) {
                                line1 = i;
                                line2 = j;
                                matU.SwapRows(line1, line2);
                                matL.SwapRows(line1, line2);
                                matP.SwapRows(line1, line2);
                                i = i - 1;
                                numSwap++;
                                isSing = true;
                                break;
                            }
                        }
                        if(i >= _rank)
                            break;
                        if(isSing == false) {
                            if(_rank - 1 > i) {
                                matU.SwapCols(i, _rank - 1);
                                matQ.SwapCols(i, _rank - 1);
                                numSwap++;
                            }
                            _rank--;
                            i = i - 1;
                        }
                }
                else {
                    isSing = true;
                    leadElem = matU[i,i];
                    line1 = i;
                    line2 = i;
                    for(int j = i; j < size; j++) {
                        if(Math.Abs(matU[j,i]) > leadElem && matU[j,i] != 0) {
                            leadElem = Math.Abs(matU[j,i]);
                            line2 = j;
                        }
                    }
                    if(line1 != line2) {
                        matU.SwapRows(line1, line2);
                        matL.SwapRows(line1, line2);
                        matP.SwapRows(line1, line2);
                        numSwap++;
                    }
                    double coeff;        
                    for(int k = i + 1; k < size; k++) {
                        if(matU[i,i] != 0) {
                            coeff = matU[k, i] / matU[i, i];
                            for(int j = i; j < size; j++) {
                                matU[k, j] = matU[k, j] - matU[i, j] * coeff;
                                if(Math.Abs(matU[k,j]) < eps)
                                    matU[k, j] = 0;

                                numOfOperations += 3;
                            }
                            matL[k, i] = coeff;
                        }
                    }
                    
                }
            }
            for(int i = 0; i < size; i++)
                matL[i,i] = 1;
                _determinant = 1;
            for(int i = 0; i < size; i++)
                _determinant *= matU[i,i];
            if(numSwap % 2 != 0)
                _determinant *= (-1);
            }
            else {
                isLU = true;
                isSing = true;
                for(int i = 0; i < rowSize; i++) {
                    for(int j = 0; j < colSize; j++) {
                        matU[i,j] = _mat[i,j];
                    }
                }

                double leadEl = 0;
                int line1,line2;
                for(int i = 0; i < rowSize; i++) {
                    if(i < colSize) {
                        if(matU[i,i] == 0) {
                            isSing = false;
                            if(i == rowSize - 1) {
                                for(int j = i + 1; j < colSize; j++) {
                                    if(matU[i,j] != 0) {
                                        isSing = true;
                                        break;
                                    }
                                }
                            }
                            for(int j = i + 1; j < rowSize; j++) {
                                if(matU[j,i] != 0) {
                                    line1 = i;
                                    line2 = j;
                                    matU.SwapRows(line1, line2);
                                    matL.SwapRows(line1, line2);
                                    matP.SwapRows(line1, line2);
                                    i = i - 1;
                                    isSing = true;
                                    break;
                                }
                            }
                            if(i >= _rank)
                                break;
                            if(isSing == false) {
                                if(_rank - 1 > i) {
                                    matU.SwapCols(i, _rank - 1);
                                    matQ.SwapCols(i, _rank - 1);
                                }
                                _rank--;
                                i = i - 1;
                            }
                        }
                        else {
                        leadEl = matU[i,i];
                        line1 = i;
                        line2 = i;
                        for(int j = i; j < rowSize; j++) {
                            if(Math.Abs(matU[j,i]) > leadEl && matU[j,i] != 0) {
                                leadEl = matU[j,i];
                                line2 = j;
                            }

                        }
                        if(line1 != line2) {
                            matU.SwapRows(line1, line2);
                            matP.SwapRows(line1, line2);
                            matL.SwapRows(line1, line2);
                        }

                        double coeff;
                        for(int k = i + 1; k < rowSize; k++) {
                            if(matU[i,i] != 0) {
                                coeff = matU[k,i] / matU[i,i];
                                for(int j = i; j < colSize; j++) {
                                    matU[k,j] -= matU[i,j] * coeff;
                                    if(Math.Abs(matU[k,j]) < eps)
                                        matU[k,j] = 0;
                                }
                                matL[k,i] = coeff;
                            }
                        }
                        }
                    }
                }
                    for(int i = 0; i < rowSize; i++)
                        matL[i,i] = 1;
            }
        }

        public static void Output(double[,] matrix,int sizeRows,int sizeCols) {
            for(int i = 0; i < sizeRows; i++) {
                for(int j = 0; j < sizeCols; j++)
                    Console.Out.Write(Math.Round(matrix[i,j],5) + " ");
                Console.Out.Write("\n");
            }
            Console.Out.Write("\n" + "\n");
        }

        public static double[,] Multiply(double[,] matr1,double[,] matr2,int sizeRows1,int c,int sizeCols2) {
            double temp = 0;
            double[,] matrix = new double[sizeRows1,sizeCols2];
            for(int i = 0; i < sizeRows1; i++) {
                for(int j = 0; j < sizeCols2; j++) {
                    for(int k = 0; k < c; k++) {
                        temp += matr1[i,k] * matr2[k,j];
                    }
                    matrix[i,j] = temp;
                    temp = 0;
                }
            }
            return matrix;
        }
        
        public static double[] multiplyVec(double[,] matrix,double[] vec,int vecSize) {
            double[] result = new double[matrix.GetLength(0)];
            for(int i = 0; i < matrix.GetLength(0); i++) {
                for(int j = 0; j < matrix.GetLength(1); j++) {
                     result[i] += vec[j] *   matrix[i,j];
                }
            }
            return result;
        }
        
        public double[,] Solution(double[] vec) {
            double[,] resultSolution;
            if(!isLU)
                LUDecompose();

            double[] X = new double[size];
            double[] Y = new double[size];
            int fsrNum = size - _rank;
                if(fsrNum !=0)
                    resultSolution = new double[fsrNum,size];
                else {
                    resultSolution = new double[1,size];
                    fsrNum = 1;
                }


            vec = multiplyVec(matP,vec,size);        
            isSing = true;

            for(int i = 0; i < size; i++){
                if(vec[i] != 0){
                    isSing = true;
                    break;
                }
            }
            for(int i = 0;i < fsrNum;i++){
                for(int j = size - 1; j >_rank - 1; j--) {
                    X[j] = rnd.Next(0,2) ;
                }
                // look for y :
                // Ly = b
                for(int j = 0; j < size; j++){
                    Y[j] = vec[j];
                        for (int k = 0; k < j; k++)
                            Y[j] -= matL[j,k] * Y[k];
                }
                // look for x :
                // Ux = y
                for(int j = _rank - 1; j > -1; j--){
                    X[j] = Y[j] / matU[j,j];
                    for(int k = j+1;k < size;k++)
                        X[j] -= matU[j,k] * X[k] / matU[j,j];
                }
                X = multiplyVec(matQ,X,size);            
                for(int j = 0; j < size; j++)
                    resultSolution[i,j] = X[j];
                bool isSame = true;
                if(i>0)
                for(int k = i; k >0; k--)
                for(int j = 0; j < size; j++) {
                    if(resultSolution[i,j] != resultSolution[k-1,j])
                        isSame = false;
                }
                
                if(isSame && i > 0)
                    i--;

            }
                // find matrix condition number
                double vecNorm = 0;
                double solutionNorm = 0;
                for(int i = 0; i < vec.Length - 1; i++){
                    vecNorm += vec[i] * vec[i];
                    solutionNorm += X[i] * X[i];
                }
                vecNorm = Math.Sqrt(vecNorm);
                solutionNorm = Math.Sqrt(solutionNorm);
            return resultSolution;
        }   
        public void QRDecompose() {
            if(isSquare) {
                isQR = true;
                Q = new double[size,size];
                R = new double[size,size];
                double[,] rotateMat = new double[size,size];
                double sin,cos;                         
                for(int i = 0; i < size; i++) {
                    rotateMat[i,i] = 1;
                    Q[i,i] = 1;
                    for(int j = 0; j < size; j++)
                        R[i,j] = _mat[i,j];
                }
                for(int j = 0; j < size - 1; j++) {
                    for(int i = j + 1; i < size; i++) {
                        if(R[i,j] != 0) {
                            cos = R[j,j] / (Math.Sqrt( Math.Pow (R[j,j], 2) 
                                                                    + Math.Pow (R[i,j], 2)));
                            sin = -R[i,j] / (Math.Sqrt( Math.Pow (R[j,j], 2) 
                                                                    + Math.Pow (R[i,j], 2)));
                            
                            rotateMat[j,j] = cos;    rotateMat[j,i] = -sin;
                            rotateMat[i,j] = sin;    rotateMat[i,i] = cos;

                            R = Multiply(rotateMat,R,size,size,size);
                            Q = Multiply(rotateMat,Q,size,size,size);
                            
                            rotateMat[j,j] = 1;    rotateMat[j,i] = 0;
                            rotateMat[i,j] = 0;    rotateMat[i,i] = 1;

                        }
                    }
                }
                for(int i = 0; i < size; i++)
                for(int j = 0; j < size; j++) {
                    if(Math.Abs(Q[i,j]) < eps)
                        Q[i,j] = 0;
                    if(Math.Abs(R[i,j]) < eps)
                        R[i,j] = 0;
                }
            }
            else throw new Exception("Матрица не квадратная");
        }
        public double[] QRSolution(double[] b) {
            if(!isQR) {
                QRDecompose();
            }
            double[] X = new double[size];
            b = multiplyVec(Q,b,size);
            for(int j = size - 1; j > -1; j--) {
                X[j] = b[j] / R[j,j];

                for(int i = j+1; i < size; i++)
                    X[j] -= R[j,i] * X[i] / R[j,j];
            }

            return X;
        }
        public double[] methodJacoby(double[] b,double[] x0) {
            double[,] matrixB = new double[size,size];          
            double[] Xk = new double[size];                     
            double[] Xk1 = new double[size];                    
            jacobyInter = 0;                                    

            Xk = x0;
            for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++)
                matrixB[i,j] = _mat[i,j];
            }
            for(int i = 0; i < size; i++) {
                matrixB[i,i] = 0;                               
                for(int j = 0;j < size;j++)
                    matrixB[i,j] /= _mat[i,i];
            }
            for(int i = 0; i < size; i++) {
                b[i] /= _mat[i,i];
                Xk[i] = b[i];
            }
            double dif;
            if(Norm(matrixB) < 0.5)
                dif = (1 - Norm(matrixB)) /Norm(matrixB) ;
            else dif = eps;
            
            do {
                double normDif = 0;
                jacobyInter++;
                for(int i = 0; i < size; i++)
                    Xk1[i] = Xk[i];
                Xk = multiplyVec(matrixB,Xk,size);

                for(int i = 0; i < size; i++) {
                    Xk[i] = -Xk[i] + b[i];
                }

                for(int j = 0; j < size; j++) {
                    normDif += (Xk[j] - Xk1[j]) * (Xk[j] - Xk1[j]);
                }
                normDif = Math.Sqrt(normDif);

                if(normDif < dif)
                    break;

            }while(true);


            return Xk;
        }   
       
        public double[] methodSeidel(double[] b,double[] x0) {

            double[,] matrixB = new double[size,size];          
            double[] Xk = new double[size];                     
            double[] Xk1 = new double[size];                    
            seidelIter = 0;                                    
            double[] Xp = new double[size];

            for(int i = 0; i < size; i++)
                Xk[i] = x0[i];

            for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++)
                matrixB[i,j] = _mat[i,j];
            }
            for(int i = 0; i < size; i++) {
                matrixB[i,i] = 0;                               
                for(int j = 0; j < size; j++)
                    matrixB[i,j] /= _mat[i,i];
            }
            for(int i = 0; i < size; i++) {
                b[i] /= _mat[i,i];
                Xk[i] = b[i];
            }

            double dif;
            if(Norm(matrixB) < 0.5)
                dif = (1 - Norm(matrixB)) /Norm(matrixB) ;
            else dif = eps;

            do {
                seidelIter++;

                for(int i = 0; i < size; i++)
                    Xp[i] = Xk[i];

                double normDif = 0;
                
                for(int i = 0; i < size; i++) {
                    Xk1[i] = b[i];
                    for(int j = 0; j < size; j++) {
                        Xk1[i] -= Xk[j] * matrixB[i,j];
                    }
                    Xk[i] = Xk1[i];
                }

                for(int j = 0; j < size; j++) {
                    normDif += (Xk[j] - Xp[j]) * (Xk[j] - Xp[j]);
                }
                normDif = Math.Sqrt(normDif);

                if(normDif < dif)
                    break;


            }while(true);
            
                return Xk;
        }
    }
}