import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Backend {
    private int mod;

    public void setMod(int mod){
        this.mod=mod;
    }

    public int getLastDeg(ArrayList<Integer> coeffs){
        int lastDeg = coeffs.size()-1;
        while(lastDeg>=0 && coeffs.get(lastDeg)==0){
            lastDeg-=1;
        }
        return lastDeg;
    }

    private int modulo(int num, int mod){
        if(num>0){
            return (num+mod)%mod;
        }
        while(num<=0){
            num+=mod;
        }
        num+=mod;
        return num%mod;
    }

    public int modInverse(int a, int m){
        a=a&m;
        for(int x=1;x<m;x++){
            if((a*x)%m==1){
                return x;
            }
        }
        return 1;
    }

    public int inverseMatrix(int num){
        for(int i=1;i<this.mod;i++){
            if((num*i)%this.mod==1){
                return i;
            }
        }
        return 1;
    }

    public void gaussJordanElimination(List<List<Integer>> matrix){
        //приведение матрицы к ступенчатому виду

        int n=matrix.size();
        int m=matrix.get(0).size();

        for(int i=0; i<n; i++){
            //находим ведущий элемент
            if(matrix.get(i).get(i)==0){
                for(int k=i+1;k<n;k++){
                    if(matrix.get(k).get(i)!=0){
                        List<Integer> temp = matrix.get(i);
                        matrix.set(i,matrix.get(k));
                        matrix.set(k,temp);
                        break;
                    }
                }
            }
            if(matrix.get(i).get(i)==0) continue;//если ведущий элемент равен нулю - пропуск шага

            //нормализовать ведущий элемент до 1
            int inv = modInverse(matrix.get(i).get(i),this.mod);
            for(int j=0;j<m;j++){
                matrix.get(i).set(j,modulo(matrix.get(i).get(j)*inv,this.mod));
            }

            //обнулить элементы в текущем столбце над и под ведущим элементом
            for(int k=0;k<n;k++){
                if(k!=i && matrix.get(k).get(i)!=0){
                    int factor = matrix.get(k).get(i);
                    for(int j=0;j<m;j++){
                        matrix.get(k).set(j,modulo(matrix.get(k).get(j)-(factor*matrix.get(i).get(j)),this.mod));
                    }
                }
            }
        }
    }

    public List<List<Integer>> findBasisOfSolutionSpace(List<List<Integer>> matrix){
        //функуция для нахождения базиса пространства решений однородной системы
        int n=matrix.size();
        int m=matrix.get(0).size();
        List<List<Integer>> basis = new ArrayList<>();

        gaussJordanElimination(matrix);


        ArrayList<Boolean> isFreeVariable = new ArrayList<>();
        for(int i=0;i<m;i++){
            isFreeVariable.add(true);
        }

        //определение свободных переменных
        for(int i=0;i<n;i++){
            int leadingEntry=-1;
            for(int j=0;j<m;j++){
                if(matrix.get(i).get(j)!=0){
                    leadingEntry=j;
                    break;
                }
            }
            if(leadingEntry!=-1){
                isFreeVariable.set(leadingEntry,false);
            }
        }

        //формирование базиса пространства решений
        for(int j=0;j<m;j++){
            if(isFreeVariable.get(j)){
                ArrayList<Integer> solution = new ArrayList<>();
                for(int abc=0;abc<m;abc++){
                    solution.add(0);
                }
                solution.set(j,1);
                for(int i=0;i<n;i++){
                    int leadingEntry=-1;
                    for(int k=0;k<m;k++){
                        if(matrix.get(i).get(k)!=0){
                            leadingEntry=k;
                            break;
                        }
                    }
                    if(leadingEntry!=-1){
                        solution.set(leadingEntry,modulo(-1*matrix.get(i).get(j),this.mod));
                    }
                }
                basis.add(solution);
            }
        }
        return basis;
    }

    public List<List<Integer>> solveFundamentalSystem(List<List<Integer>> matrix){
        //функция для решения фундаментальной системы решений

        //применяем модуль ко всем элементам
        for(int i=0;i<matrix.size();i++){
            for(int j=0;j<matrix.get(i).size();j++){
                matrix.get(i).set(j,modulo(matrix.get(i).get(j),this.mod));
            }
        }

        //возвращаем базис пространства решений
        return findBasisOfSolutionSpace(matrix);
    }

    public void gaussElimination(List<List<Integer>> matrix){
        int rows = matrix.size();
        int cols = matrix.get(0).size();
        int lead=0;

        for(int r=0;r<rows;r++){
            if(cols<=lead){
                break;
            }
            int i=r;
            while(matrix.get(i).get(lead)==0){
                ++i;
                if(rows==i){
                    i=r;
                    ++lead;
                    if(cols == lead){
                        return;
                    }
                }
            }
            List<Integer> temp = matrix.get(i);
            matrix.set(i,matrix.get(r));
            matrix.set(r,temp);

            int lv = matrix.get(r).get(lead);
            for(int j=0;j<cols;j++){
                matrix.get(r).set(j,modulo(matrix.get(r).get(j)*inverseMatrix(lv),this.mod));
            }
            for(i=0;i<rows;i++){
                if(i!=r){
                    lv = matrix.get(i).get(lead);
                    for(int j=0;j<cols;j++){
                        int temp2 = matrix.get(i).get(j);
                        matrix.get(i).set(j,modulo(temp2-(matrix.get(r).get(j)*lv),this.mod));
                    }
                }

            }
            ++lead;
        }
    }

    public void printMatrix(List<List<Integer>> matrix){
        Log logger = new Log();
        for(int i=0;i<matrix.size();i++){
            for(int n=0;n<matrix.get(i).size();n++){
                logger.writeLog(matrix.get(i).get(n)+" ",false);
            }
            logger.writeLog("\n",false);
        }
    }

    public List<List<Integer>> buildSylvesterMatrix(ArrayList<Integer> coeffs1, ArrayList<Integer> coeffs3, int c){
        //функция строения матрицы сильвестра
        int n1 = getLastDeg(coeffs1);
        int n3 = getLastDeg(coeffs3);

        int sizeMatrix=n1+n3;//размер матрицы

        List<List<Integer>> sylvesterMatrix = new ArrayList<>();
        for(int i=0;i<sizeMatrix;i++){
            ArrayList<Integer> temp = new ArrayList<>();
            for(int m=0;m<sizeMatrix;m++){
                temp.add(0);
            }
            sylvesterMatrix.add(temp);
        }

        int filledRows = 0;//Переменная для отслеживания количества заполненных строк

        //Заполнение коэффициентами многочлена coeffs1
        for(int j=0;j<n3;++j){
            for(int i=0;i<=n1;++i) {
                sylvesterMatrix.get(j).set((i+j),coeffs1.get(n1-i));
            }
            filledRows++;
        }

        //Заполнение коэффициентами многочлена coeffs3
        for(int j=filledRows;j<sizeMatrix;j++){
            for(int i=0;i<=n3;i++){
                sylvesterMatrix.get(j).set((i+j-n3),coeffs3.get(n3-i));
            }
        }

        //Заполнение c
        for(int j=filledRows;(j<(filledRows+n1) && j<sizeMatrix);j++){
            sylvesterMatrix.get(j).set((n3+j-filledRows),modulo(coeffs3.get(0)-c,this.mod));
        }
        return sylvesterMatrix;
    }

    public int determinant(List<List<Integer>> matrix){
        int n = matrix.size();
        //System.out.println("n: "+n);
        int det = 1;

        for(int i = 0; i < n; i++) {
            if(matrix.get(i).get(i) == 0) {
                boolean swapped = false;
                for (int j = (i + 1); j < n; j++) {
                    //System.out.println(matrix.get(j).get(i));
                    if (matrix.get(j).get(i) != 0) {
                        ArrayList<Integer> temp = (ArrayList<Integer>) matrix.get(i);
                        matrix.set(i,matrix.get(j));
                        matrix.set(j,temp);
                        det = (this.mod - det) % this.mod;
                        swapped = true;
                        break;
                    }
                }
                if (!swapped) return 0;
            }

            det = (det * matrix.get(i).get(i)) % this.mod;
            int inv = 1;
            int b = matrix.get(i).get(i);
            for (int j = 0; j < this.mod - 2; j++) {
                inv = (inv * b) % this.mod;
            }

            for (int j = i + 1; j < n; j++) {
                int factor = (matrix.get(j).get(i) * inv) % this.mod;
                for (int k = i; k < n; k++) {
                    matrix.get(j).set(k,modulo((matrix.get(j).get(k)-(factor*matrix.get(i).get(k))+this.mod),this.mod));
                }
            }
        }
        return det;
    }

    public ArrayList<Polynomial> addNonZeroPolynomialsInOrder(List<List<Polynomial>> allAdjustedPolys, ArrayList<Polynomial> resultPoly){
        int m = allAdjustedPolys.size();
        int numPolys = allAdjustedPolys.get(0).size();

        for(int j=0;j<numPolys;j++){
            for(int i=0;i<m;i++){
                if(!allAdjustedPolys.get(i).get(j).isZero()){
                    resultPoly.add(allAdjustedPolys.get(i).get(j));
                }
            }
        }
        return resultPoly;
    }

    public void run(int[] coefficients){
        Log logger = new Log();
        ArrayList<Integer> coeffs1 = new ArrayList<>();
        for(int i=0;i<coefficients.length;i++){
            coeffs1.add(coefficients[i]);
        }

        logger.writeLog("User mod input: "+this.mod,true);

        Polynomial poly1 = new Polynomial(coeffs1, this.mod);
        logger.writeLog("Your polinomial: ",true);
        poly1.printPolynomial();

        logger.writeLog("With module "+this.mod+"\n",false);

        int n=getLastDeg(coeffs1);
        logger.writeLog("Last polynomial degree: "+n,true);

        List<List<Integer>> coeffs2 = new ArrayList<>();


        //Создание коэффициентов для деления
        //Например: если MOD = 5; то при i = 1: coeffs_i = {0, 0, 0, 0, 0, 1} или x^5 (здесь обратный порядок)
        for(int i=0;i<n;i++){
            int[] temp = new int[this.mod*i+1];
            ArrayList<Integer> coeffs_i = new ArrayList<>();
            for(int m=0;m<temp.length;m++){
                coeffs_i.add(temp[m]);
            }
            coeffs_i.set(coeffs_i.size()-1,1);
            coeffs2.add(coeffs_i);
        }

        List<List<Integer>> matrixQ = new ArrayList<>();

        for(int i=0;i<coeffs2.size();i++){
            Polynomial poly2 = new Polynomial((ArrayList<Integer>)coeffs2.get(i),this.mod);
            //printPolynomial(poly2);
            Polynomial[] divisionResult = poly2.divide(poly1);
            matrixQ.add(divisionResult[1].getCoefficients());
        }

        for(int i=0;i<matrixQ.size();i++){
            int expectedSize = coeffs1.size()-1;
            while(matrixQ.get(i).size()<expectedSize){
                matrixQ.get(i).add(0);
            }
        }

        logger.writeLog("Matrix Q:\n",true);
        printMatrix(matrixQ);

        for(int i=0;i<n;i++){
            for(int j=0;j<i;j++){
                //System.out.println("i: "+i+" ; j: "+j);
                int c = matrixQ.get(i).get(j);
                matrixQ.get(i).set(j, matrixQ.get(j).get(i));
                matrixQ.get(j).set(i,c);
            }
        }

        logger.writeLog("Transposed matrix QT: \n",false);
        printMatrix(matrixQ);

        //создание единичной матрицы
        List<List<Integer>> matrixE = new ArrayList<>();
        int[] temp1 = new int[n];
        for(int i=0;i<n;i++){
            ArrayList<Integer> temp2 = new ArrayList<>();
            for(int m=0;m<n;m++){
                temp2.add(temp1[m]);
            }
            matrixE.add(temp2);
        }
        for(int i=0;i<n;i++){
            matrixE.get(i).set(i,1);
        }

        logger.writeLog("Matrix E:\n",false);
        printMatrix(matrixE);


        //вычитаем из матрицы QT единичную
        List<List<Integer>> matrixQTMinusE = new ArrayList<>();
        temp1 = new int[n];
        for(int i=0;i<n;i++){
            ArrayList<Integer> temp2 = new ArrayList<>();
            for(int m=0;m<n;m++){
                temp2.add(temp1[m]);
            }
            matrixQTMinusE.add(temp2);
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                matrixQTMinusE.get(i).set(j,(matrixQ.get(i).get(j)-matrixE.get(i).get(j)));
            }
        }

        for(int i=0;i<n;i++){//modulo for every elem
            for(int j=0;j<n;j++){
                matrixQTMinusE.get(i).set(j,modulo(matrixQTMinusE.get(i).get(j),this.mod));
            }
        }

        logger.writeLog("Matrix QT-E:\n",false);
        printMatrix(matrixQTMinusE);

        //приведение матрицы к ступенчатому виду
        gaussElimination(matrixQTMinusE);

        //модуль для каждого элемента ступенчатой матрицы
        for(int i=0;i<n;i++){//modulo for every elem
            for(int j=0;j<n;j++){
                matrixQTMinusE.get(i).set(j,modulo(matrixQTMinusE.get(i).get(j),this.mod));
            }
        }

        logger.writeLog("Matrix QT-E in step view:\n",false);
        printMatrix(matrixQTMinusE);


        //считаем ранг ступенчатой матрицы
        int rank=0;
        for(int i=0;i<matrixQTMinusE.size();i++){
            boolean allZeros = true;
            for(int m=0;m<matrixQTMinusE.get(i).size();m++){
                if(matrixQTMinusE.get(i).get(m)!=0){
                    allZeros=false;
                    break;
                }
            }
            if(!allZeros){
                ++rank;
            }
        }

        List<List<Integer>> resultMatrix = new ArrayList<>();


        //идем по строкам ступенчатой матрицы
        for(int i=0;i<matrixQTMinusE.size();i++){
            boolean allZeros=true;
            for(int m=0;m<matrixQTMinusE.get(i).size();m++){
                if(matrixQTMinusE.get(i).get(m)!=0){
                    allZeros=false;
                    break;
                }
            }
            //если строка ненулевая, добавляем ее в вывод
            if(!allZeros){
                resultMatrix.add(matrixQTMinusE.get(i));
            }
        }

        //итоговая матрица
        logger.writeLog("Final matrix:\n",false);
        printMatrix(resultMatrix);

        logger.writeLog("Matrix rank: "+rank,true);

        //количество множителей
        int k = n - rank;
        logger.writeLog("Amount of multipliers: "+k,true);
        if(k<=1){
            System.out.println("No soluion");
            return;
        }


        //Создаем вектор для хранения ФСР (фундаментальной системы решений)
        List<List<Integer>> FSRMatrix = solveFundamentalSystem(resultMatrix);

        //Выводим фундаментальную систему решений
        logger.writeLog("Basis of the solution space:\n",true);
        printMatrix(FSRMatrix);

        //Удаляем первую строку
        FSRMatrix.removeFirst();

        //Arraylist for polinomial storage
        ArrayList<Polynomial> polys = new ArrayList<>();

        //Преобразование каждой строки матрицы в многочлен и сохранение в виде ArrayList
        for(int i=0;i<FSRMatrix.size();i++){
            Polynomial poly = new Polynomial((ArrayList<Integer>)FSRMatrix.get(i),this.mod);
            polys.add(poly);
        }

        //Итоговые многочлены:
        logger.writeLog("Polynomials: ",true);
        for(int i=0;i<polys.size();i++){
            polys.get(i).printPolynomial();
        }
        logger.writeLog(" ",true);

        ArrayList<Integer> coeffs3 = (ArrayList<Integer>) FSRMatrix.get(0);
        Polynomial poly3 = new Polynomial(coeffs3,this.mod);


        //ArrayList для хранения результатов
        ArrayList<Integer> resultants = new ArrayList<>();
        for(int i=0;i<this.mod;i++){
            resultants.add(0);
        }

        //Создаем матрицу Сильвестра и находим его определители
        for(int c=0;c<this.mod;c++){
            List<List<Integer>> sylvesterMatrix = buildSylvesterMatrix(coeffs1,coeffs3,c);
            resultants.set(c,modulo(determinant(sylvesterMatrix),this.mod));
        }

        logger.writeLog("Results for different \"c\" values:",true);
        for(int c=0;c<this.mod;c++){
            logger.writeLog("c = "+c+": "+resultants.get(c),true);
        }

        ArrayList<Integer> root = new ArrayList<>();

        //Записываем корни равные нулевому значению
        for(int c=0;c<this.mod;c++){
            if(resultants.get(c)==0){
                root.add(c);
            }
        }

        logger.writeLog("Polynomial roots (where polynomial = 0):\n",true);
        for(int i=0;i<root.size();i++){
            logger.writeLog(root.get(i)+" ",false);
        }
        logger.writeLog("\n",false);

        //хранение для хранения всех результатов adjustedPolynomial
        List<List<Polynomial>> allAdjustedPolys = new ArrayList<>();
        for(int i=0;i<k;i++){
            allAdjustedPolys.add(new ArrayList<>());
        }

        //перебор значений для -root
        for(int i=0;i<root.size();i++){
            logger.writeLog("Subtracting "+root.get(i)+"\n",true);
            ArrayList<Integer> temp = new ArrayList<>();
            temp.add(root.get(i));
            Polynomial adjustedPoly = poly3.substract(new Polynomial(temp,this.mod));
            allAdjustedPolys.get(i).add(adjustedPoly);
            adjustedPoly.printPolynomial();
            logger.writeLog("\n",false);
        }

        ArrayList<Polynomial> resultPoly = new ArrayList<>();

        addNonZeroPolynomialsInOrder(allAdjustedPolys, resultPoly);

        logger.writeLog("Final polynomials:\n",true);
        for(int i=0;i<resultPoly.size();i++){
            resultPoly.get(i).printPolynomial();
        }
        logger.writeLog(" ",true);

        //определение размеров матрицы
        int maxDeg = 0;
        for(int i=0;i<resultPoly.size();i++){
            int deg = getLastDeg(resultPoly.get(i).getCoefficients());
            maxDeg = Math.max(maxDeg,deg);
        }
        //создание матрицы коэффициентов
        List<List<Integer>> resultPolyMatrix = new ArrayList<>();
        for(int i=0;i<resultPoly.size();i++){
            ArrayList<Integer> temp = new ArrayList<>();
            for(int m=0;m<maxDeg+1;m++){
                temp.add(0);
            }
            resultPolyMatrix.add(temp);
        }

        //заполнение матрицы коэффициентов
        for(int i=0;i<resultPoly.size();i++){
            ArrayList<Integer> coeffs = resultPoly.get(i).getCoefficients();
            for(int j=0;j<coeffs.size();j++){
                resultPolyMatrix.get(i).set(j,coeffs.get(j));
            }
        }

        logger.writeLog("Coefficient matrix resultPolyMatrix:\n",true);
        printMatrix(resultPolyMatrix);
        logger.writeLog(" ",true);

        logger.writeLog("GCD:\n",true);
        System.out.println("\n----- ANSWER -----");
        poly1.printPolynomialAnswer();
        System.out.print(" = ");
        for(int i=0;i<resultPolyMatrix.size();i++){
            poly3 = new Polynomial((ArrayList<Integer>) resultPolyMatrix.get(i),this.mod);
            Polynomial gcdResult = poly1.gcd(poly1,poly3);
            gcdResult.printPolynomial();
            gcdResult.printPolynomialAnswer();
            if(i!=resultPolyMatrix.size()-1) System.out.print("*");
        }
        System.out.print(" mod "+this.mod);

        //Сгенерированные примеры:
        //{1, 20, 5, 14, 18, 6, 7} MOD 23
        //{1, 9, 8, 15, 13, 7, 6, 5, 16, 10, 18} MOD 23
        //{1, 18, 16, 11, 3, 22, 5, 19, 8, 0, 18, 6, 21, 21, 2, 10} MOD 23


    }

}
