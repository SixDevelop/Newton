using System;
using System.Diagnostics;

namespace Newton{
    public class Newton{
        private double[,] F;                     
        private double[] F1;
        private double[] x,xk1;                 
        private double[] currX;              
        private Matrix J;               
        private int size;                       
        private double eps;
        private int amountOfOperations;
        private double norm;                    
        private int iter;                       
        public Stopwatch sw;                    
        public Newton() { 
            size = 10;
            F = new double[1,size];
            eps = 1e-3;
            x = new double[size];
            xk1 = new double[size];
            currX = new double[size];
            F1 = new double[size];
            iter = 0;
            norm = 0;
            sw = new Stopwatch();
        }

        public double[] Method() {
            Console.Out.Write("Newton with search reverse matrix :" + "\n");
            sw.Start();


            x.InitX();
            iter = 0;
            amountOfOperations = 0;
            while(true) {
                iter++;

                F1.InitialF(x);
                for(int i = 0; i < size; i++)
                    F[0,i] = F1[i];
                norm = 0;
                J = new Matrix(InitializationExtension.InitialJ(x));
                J.LUDecompose();
                amountOfOperations += J.numOfOperations;

                currX = Matrix.multiplyVec(J.revMatrix,F1,size);
                
                for(int i = 0; i < size; i++) {
                    xk1[i] = x[i] - currX[i];
                
                    if (Math.Abs(xk1[i] - x[i]) > norm)
                        norm = Math.Abs(xk1[i] - x[i]);

                    x[i] = xk1[i];

                    amountOfOperations += 5;
                }

                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Number of iterations:" + iter + "\n");
            Console.Out.Write("Number of arithmetic operations:" + amountOfOperations + "\n");
            Console.Out.Write("Spent time:" + sw.Elapsed + "\n");
            return xk1;
        }
        public double[] ModifiedMethod() {
           Console.Out.Write("Modified Newton:" + "\n");
           sw.Start();
           x.InitX();
           iter = 0;
           amountOfOperations = 0;
           J = new Matrix(InitializationExtension.InitialJ(x));

           J.LUDecompose();
           amountOfOperations += J.numOfOperations;

            while(true) {
               iter++;
               norm = 0;
               F1.InitialF(x);
               for(int i = 0;i < size;i++)
                F[0,i] = F1[i];
               var a = J.Solution(F1);
                for(int i = 0; i < size; i++)
                    currX[i] = a[0,i];
                for(int i = 0;i < size;i++) {
                   xk1[i] = x[i] - currX[i];

                    if (Math.Abs(xk1[i] - x[i]) > norm)
                        norm = Math.Abs(xk1[i] - x[i]);

                    x[i] = xk1[i];

                    amountOfOperations += 5;
                }
                
                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Number of iterations:" + iter + "\n");
            Console.Out.Write("Number of arithmetic operations:" + amountOfOperations + "\n");
            Console.Out.Write("Spent time:" + sw.Elapsed + "\n");

            return xk1;
         }   
        public double[] TransitionToModMeth() {
            Console.Out.Write("Classic to modified Newton:" + "\n");
            int step = 8;
            Console.WriteLine("Step k = " + step);
            x.InitX();

            double[] xFix = new double[size];
            xFix.InitX();

            double[] x0 = new double[size];
            x0.InitX();
            iter = 0;
            amountOfOperations = 0;
            sw.Start();
            while(true){
                iter++;
                norm = 0;

                F1.InitialF(x);
                for(int i = 0; i < size; i++)
                    F[0,i] = F1[i];

                if(step > 0){
                    J = new Matrix(InitializationExtension.InitialJ(x));
                    J.LUDecompose();
                    amountOfOperations += J.numOfOperations;

                    for(int i = 0; i < size; i++)
                        xFix[i] = x[i];
                    currX = Matrix.multiplyVec(J.revMatrix,F1,size);
                    step--;
                }
                else {
                    if(step == 0) {
                        J = new Matrix(InitializationExtension.InitialJ(xFix));
                        J.LUDecompose();
                    }
                    amountOfOperations += J.numOfOperations;
                    currX = Matrix.multiplyVec(J.revMatrix,F1,size);

                    step--;
                }
                for(int i = 0; i < size; i++){
                    xk1[i] = x[i] - currX[i];

                    if (Math.Abs(xk1[i] - x[i]) > norm)
                        norm = Math.Abs(xk1[i] - x[i]);

                    x[i] = xk1[i];

                    amountOfOperations += 5;
                }
                if(norm < eps)
                    break;

            }
            sw.Stop();
            Console.Out.Write("Number of iterations:" + iter + "\n");
            Console.Out.Write("Number of arithmetic operations:" + amountOfOperations + "\n");
            Console.Out.Write("Spent time:" + sw.Elapsed + "\n");
            return xk1;
        }

        public double[] MethodWithRewriteRevMat(){
            int step = 7;
            Console.Out.Write("Searching reverse matrix every k = " + step + " iterations :" + "\n");
            
            x.InitX();
            J = new Matrix(InitializationExtension.InitialJ(x));

            J.LUDecompose();
            amountOfOperations += J.numOfOperations;

            iter = 0;
            amountOfOperations = 0;
            sw.Start();
            while(true) {
                iter++;
                norm = 0;
                F1.InitialF(x);
                for(int i = 0; i < size; i++)
                    F[0,i] = F1[i];
                if(iter % step == 1){
                J = new Matrix(InitializationExtension.InitialJ(x));

                    J.LUDecompose();
                    amountOfOperations += J.numOfOperations;

                    J.Reverse();
                }
                var a = J.Solution(F1);
                for(int i = 0; i < size; i++)
                    currX[i] = a[0,i];

                for(int i = 0; i < size; i++){
                    xk1[i] = x[i] - currX[i];

                    if (Math.Abs(xk1[i] - x[i]) > norm)
                        norm = Math.Abs(xk1[i] - x[i]);

                    x[i] = xk1[i];

                    amountOfOperations += 5;
                }

                if(norm < eps)
                    break;
            }
            sw.Stop();
            Console.Out.Write("Number of iterations:" + iter + "\n");
            Console.Out.Write("Number of arithmetic operations:" + amountOfOperations + "\n");
            Console.Out.Write("Spent time:" + sw.Elapsed + "\n");
            return xk1;
        }
    }
}