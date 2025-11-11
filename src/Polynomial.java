import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Polynomial {
    //polynomial operations
    private ArrayList<Integer> coefficients = new ArrayList<>();//коэффициенты хранятся в обратном порядке.
    private int mod;//change it
    public Polynomial(int mod){
        this.coefficients.add(0);
        this.mod=mod;
    }
    public Polynomial(ArrayList<Integer> coeffs,int mod){
        this.mod=mod;
        for(int i=0;i<coeffs.size();i++){
            this.coefficients.add(modulo(coeffs.get(i),this.mod));
        }
        normalize();
    }
    public Polynomial(List<List<Integer>> matrix, int mod){
        this.mod=mod;
        for(int i=0;i< matrix.size();i++){
            for(int n=0;n<matrix.get(i).size();n++){
                this.coefficients.add(modulo(matrix.get(i).get(n),this.mod));
            }
        }
    }

    public ArrayList<Integer> getCoefficients() {
        return this.coefficients;
    }

    public int getMod(){
        return this.mod;
    }

    public void setMod(int mod){
        this.mod=mod;
    }

    private void normalize(){//удаление ведущих нулей
        while(!this.coefficients.isEmpty() && this.coefficients.getLast()==0){
            this.coefficients.removeLast();
        }
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

    private int inverse(int num){//ищет обратный элемент по модулю
        for(int i = 1; i<this.mod;i++){
            if(modulo(num*i,this.mod)==1){
                return i;
            }
        }
        return 1;//если num==1, то обратный ==1
    }

    public int countNonZeroCoefficients(){
        int count=0;
        for(int i=0;i<this.coefficients.size();i++){
            if(this.coefficients.get(i)!=0){
                count+=1;
            }
        }
        return count;
    }

    public Polynomial add(Polynomial pol){
        ArrayList<Integer> result = new ArrayList<>();
        int maxSize = Math.max(this.coefficients.size(),pol.coefficients.size());

        for(int i=0;i<maxSize;i++){
            int coeff1 = 0;
            int coeff2 = 0;
            if(i<this.coefficients.size()) coeff1 = this.coefficients.get(i);
            if(i<pol.coefficients.size()) coeff2 = pol.coefficients.get(i);
            result.add(modulo((coeff1+coeff2),this.mod));
        }
        return new Polynomial(result,this.mod);
    }

    public Polynomial substract(Polynomial pol){
        ArrayList<Integer> result = new ArrayList<>();
        int maxSize = Math.max(this.coefficients.size(),pol.coefficients.size());

        for(int i=0;i<maxSize;i++){
            int coeff1 = 0;
            int coeff2 = 0;
            if(i<this.coefficients.size()) coeff1 = this.coefficients.get(i);
            if(i<pol.coefficients.size()) coeff2 = pol.coefficients.get(i);
            result.add(modulo((coeff1-coeff2),this.mod));
        }
        return new Polynomial(result,this.mod);
    }

    public Polynomial multiplyByScalar(int scalar){
        ArrayList<Integer> result = new ArrayList<>();
        for(int i=0;i<this.coefficients.size();i++){
            result.add(modulo(this.coefficients.get(i)*scalar,this.mod));
        }
        return new Polynomial(result,this.mod);
    }

    public Polynomial multiply(Polynomial pol){
        int[] initialFill = new int[this.coefficients.size() + pol.coefficients.size() - 1];//filled with zeros by default
        ArrayList<Integer> result = new ArrayList<>();
        for(int i=0;i<initialFill.length;i++){
            result.add(initialFill[i]);
        }

        for(int i=0;i<this.coefficients.size();i++){
            for(int n=0;n<pol.coefficients.size();n++){
                int temp=result.get(i+n);
                result.set(i+n, temp + modulo(this.coefficients.get(i)*pol.coefficients.get(n),this.mod));

                temp=result.get(i+n);
                result.set(i+n,modulo(temp,this.mod));
            }
        }

        return new Polynomial(result,this.mod);
    }

    public Polynomial[] divide(Polynomial divisior){
        //q=quotient; r=remainder
        int[] initialFill = new int[this.coefficients.size()];//filled with zeros by default
        ArrayList<Integer> qCoeffs = new ArrayList<>();
        for(int i=0;i<initialFill.length;i++){
            qCoeffs.add(initialFill[i]);
        }
        ArrayList<Integer> rCoeffs = new ArrayList<>(this.coefficients);

        while(rCoeffs.size() >= divisior.coefficients.size()){
            int degDiff = rCoeffs.size() - divisior.coefficients.size();
            int q = modulo(rCoeffs.getLast() * inverse(divisior.coefficients.getLast()),this.mod);

            qCoeffs.set(degDiff,q);

            for(int i=0; i<divisior.coefficients.size();i++){
                int temp=rCoeffs.get(i+degDiff);
                rCoeffs.set(i+degDiff,modulo(temp-divisior.coefficients.get(i)*q,this.mod));
            }

            while(!rCoeffs.isEmpty() && rCoeffs.getLast()==0){
                rCoeffs.removeLast();
            }//maybe modify "normalize()" ?
        }
        return new Polynomial[] {new Polynomial(qCoeffs,this.mod) , new Polynomial(rCoeffs,this.mod)};
    }

    public Polynomial gcd(Polynomial a, Polynomial b){
        Polynomial x = new Polynomial(a.coefficients,a.mod); //creating copy
        Polynomial y = new Polynomial(b.coefficients,b.mod); //creating copy
        while(!y.getCoefficients().isEmpty()){
            Polynomial temp = x;
            x=y;
            y=temp.divide(y)[1];//second
        }
        //Normalizing polynomial in order to make greatest deg equal 1
        ArrayList<Integer> coeffs = x.getCoefficients();
        if(coeffs.isEmpty()) return x; // Return zero polynomial as gcd
        int leadingCoeff = coeffs.getLast();
        int inverseLC=inverse(leadingCoeff);
        for(int i=0;i<coeffs.size();i++){
            coeffs.set(i,modulo(coeffs.get(i)*inverseLC,this.mod));
        }
        return new Polynomial(coeffs,this.mod);
    }

    public boolean isZero(){
        for(int i=0;i<this.coefficients.size();i++){
            if(coefficients.get(i)!=0){
                return false;
            }
        }
        return true;
    }

    public void printPolynomial(){
        Log logger = new Log();
        boolean firstTerm = true;

        for(int i=this.coefficients.size()-1;i>=0;i--){
            if(this.coefficients.get(i)!=0){
                if(!firstTerm){
                    logger.writeLog(" + ",false);
                }else{
                    firstTerm=false;
                }
                logger.writeLog(this.coefficients.get(i)+"X^"+i,false);
            }
        }
        if(firstTerm){
            logger.writeLog("0",false);//If polynomial is null (containg zeros)
        }
        logger.writeLog("\n",false);
    }

    public void printPolynomialAnswer(){
        boolean firstTerm = true;
        System.out.print("(");
        for(int i=this.coefficients.size()-1;i>=0;i--){
            if(this.coefficients.get(i)!=0){
                if(!firstTerm){
                    System.out.print(" + ");
                }else{
                    firstTerm=false;
                }
                System.out.print(this.coefficients.get(i)+"X^"+i);
            }
        }
        if(firstTerm){
            System.out.print("0");//If polynomial is null (containg zeros)
        }
        System.out.print(")");
    }
}
