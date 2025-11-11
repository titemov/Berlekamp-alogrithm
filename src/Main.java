import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class Main {

    public static int[] reverseTheArray(int[] array){
        int len = array.length;
        int[] result= new int[len];

        for(int i=0;i<len;i++){
            result[i]=array[len-1-i];
        }

        return result;
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        Log logger = new Log();
        logger.initialEntry();
        try {
            System.out.println("Enter coefficients using \"space\" as separator (for example: \"1 0 4 1\"):");
            //fix inverse
            String str = scanner.nextLine();
            String temp[] = str.split(" ");
            int[] temp2 = new int[temp.length];
            for(int i=0;i< temp.length;i++){
                temp2[i]=Integer.parseInt(temp[i]);
            }
            int[] coeffs = reverseTheArray(temp2);
            logger.writeLog("User input:"+Arrays.toString(str.split(" ")),true);
            //System.out.println(Arrays.toString(str.split(" ")));
            //System.out.println(Arrays.toString(coeffs));
            System.out.println("Enter prime number (mod):");
            int mod = scanner.nextInt();

            Backend backend = new Backend();
            backend.setMod(mod);
            backend.run(coeffs);
        }catch (Exception e){
            System.out.println("ERROR! "+e);
        }
    }
}