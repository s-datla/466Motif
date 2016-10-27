package bioinfo.motif;
import java.util.*;
import java.util.Random;

public class Main {
    public static double[][] randomMotif(int ML, double ICPC){
		double[][] motif = new double[4][ML];
		double[] columnValues = new double[ML];
		double p;
		if (ICPC == 1){
			p=0.8105;
		}
		else if (ICPC == 1.5){
			p=0.9245;
		}
		else if (ICPC == 2){
			p=1;
		}
		else{
			return null;
		}
		for (int i = 0; i<ML; i++){
			int nucleotideToEmphasize = (int) (Math.random()*4);
			//System.out.println(nucleotideToEmphasize);
			for (int j = 0; j<4; j++){
				motif[j][i] = (1-p);
			}
			motif[nucleotideToEmphasize][i]=p;
		}
		return motif;
	}
	
    private int icpc = 0;
    private int motifLength = 0;
    private int sequenceLength = 20;
    private int sequenceCount = 0;

    public static void main(String[] args) {
        new Main().run();
    }

    void run(){
        takeInput();
        for (int i = 0; i < this.sequenceCount; i++){
            generateSequences();
        }
    }

    void takeInput(){

        Scanner input = new Scanner(System.in);
        while(true){
            System.out.println("Please enter the Information Content Per Column (ICPC): ");
            try {
                this.icpc = input.nextInt();
            } catch (Exception e){
                System.out.println("Invalid Entry, please try again");
                input.nextLine();
                continue;
            }
            System.out.println("Please enter the Motif Length: ");
            try {
                this.motifLength= input.nextInt();
            } catch (Exception e){
                System.out.println("Invalid Entry, please try again");
                input.nextLine();
                continue;
            }
            System.out.println("Please enter the Sequence Length: ");
            try {
                this.sequenceLength= input.nextInt();
            } catch (Exception e){
                System.out.println("Invalid Entry, please try again");
                input.nextLine();
                continue;
            }
            System.out.println("Please enter the Sequence Count: ");
            try {
                this.sequenceCount= input.nextInt();
            } catch (Exception e){
                System.out.println("Invalid Entry, please try again");
                input.nextLine();
                continue;
            }
            break;
        }


    }

    void generateSequences(){
        Random newr = new Random();
        String sequence = "";
        String nucleotides = "ACGT";
        for(int i =0; i < this.sequenceLength; i++){
            int nextval = newr.nextInt(4);
            sequence += nucleotides.charAt(nextval);
        }
        System.out.println(sequence);
    }


}
