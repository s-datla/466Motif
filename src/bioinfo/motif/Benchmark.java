package bioinfo.motif;
import java.util.*;
import java.util.Random;
import java.io.*;

public class Benchmark {
	Iterator iter;
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
				motif[j][i] = (1-p)/3;
			}
			motif[nucleotideToEmphasize][i]=p;
		}
		return motif;
	}

	public static String bindingSite(double[][] motif){
		String s = "";
		for (int i = 0; i<motif[0].length; i++){
			double p = Math.random();
			if (p<motif[0][i]){
				s=s+"A";
			}
			else if (p<motif[0][i]+motif[1][i]){
				s=s+"C";
			}
			else if (p<motif[0][i]+motif[1][i]+motif[2][i]){
				s=s+"G";
			}
			else{
				s=s+"T";
			}
		}
		return s;
	}

    private double icpc = 0;
    private int motifLength = 0;
    private int sequenceLength = 20;
    private int sequenceCount = 0;
    private String dataSet ="";

    public static void main(String[] args) {
        new Benchmark().run();
    }

    void run(){
        //delete the files from previous run
        deleteOldFiles();
        double[][] randomMotif = randomMotif(this.motifLength, this.icpc);
        String bindingSite = bindingSite(randomMotif);
        takeInput();
        List<String> s = generateSequences();
        plantMotif(s, randomMotif);
    }

    public void createDataSet(double icpc, int ml, int sl, int sc, String dataSet) {
        this.icpc = icpc;
        this.motifLength = ml;
        this.sequenceLength = sl;
        this.sequenceCount = sc;
        this.dataSet = dataSet;

        //delete the files from previous run
        deleteOldFiles();
        double[][] randomMotif = randomMotif(this.motifLength, this.icpc);
        String bindingSite = bindingSite(randomMotif);
        List<String> s = generateSequences();
        plantMotif(s,randomMotif);
        printMotif(randomMotif);
        printMotifLength();
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

    List<String> generateSequences(){
        Random newr = new Random();
        String sequence = "";
        List<String> sequences = new ArrayList<String>();
        String nucleotides = "ACGT";
        for (int j=0; j < this.sequenceCount; j++ ) {
            sequence = "";
            for (int i = 0; i < this.sequenceLength; i++) {
                int nextval = newr.nextInt(4);
                sequence += nucleotides.charAt(nextval);
            }
            sequences.add(sequence);
        }
        return sequences;
    }

    void plantMotif(List<String> sequences, double[][] randomMotif){
    	
        Random newr = new Random();
        for (int j=0; j<sequences.size(); j++) {
        	String motif = bindingSite(randomMotif);
            int newPosition = newr.nextInt(sequences.get(j).length() - motif.length());
            System.out.println(newPosition);
            char[] chars = sequences.get(j).toCharArray();
            char[] motifc = motif.toCharArray();
            for (int i = 0; i < motif.length(); i++) {
                chars[newPosition + i] = motifc[i];
            }
            String finalseq = String.valueOf(chars);
            printFASTA(finalseq);
            printLocation(String.valueOf(newPosition));
            System.out.println(sequences.get(j) + "\n" + finalseq + "\n");
        }
    }

    void printFASTA(String s){
        try {
            File f = new File(dataSet + "/sequences.fa");
            if (!(f.exists() && !f.isDirectory())) {
                    f.createNewFile();
            }
            FileWriter fw = new FileWriter(f, true);
            System.out.println("REACHED CHECKPOINT 1");
            fw.write(s);
            fw.write("\n");
            fw.flush();
            fw.close();
        } catch (Exception e){
            e.printStackTrace();
        }

    }

    void printLocation(String location){
        try {
            File f = new File(dataSet + "/sites.txt");
            if (!(f.exists() && !f.isDirectory())) {
                f.createNewFile();
            }
            FileWriter fw = new FileWriter(f, true);
            fw.write(location);
            fw.write("\n");
            fw.flush();
            fw.close();
        } catch (Exception e){
            e.printStackTrace();
        }

    }

    void printMotif(double[][] motif){
        try {
            File f = new File(dataSet + "/motif.txt");
            if (!(f.exists() && !f.isDirectory())) {
                f.createNewFile();
            }
            FileWriter fw = new FileWriter(f);
            fw.write(">MOTIF  " + this.motifLength + "\n" );
            for (int i=0; i<this.motifLength; i++) {
                for (int j=0; j<4 ; j++) {
                    fw.write(motif[j][i]+ " ");

                }
                fw.write("\n");
            }
            fw.flush();
            fw.close();
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    void printMotifLength(){
        try {
            File f = new File(dataSet + "/motiflength.txt");
            if (!(f.exists() && !f.isDirectory())) {
                f.createNewFile();
            }
            FileWriter fw = new FileWriter(f);
            fw.write(""+this.motifLength);
            fw.flush();
            fw.close();
        } catch (Exception e){
            e.printStackTrace();
        }

    }

    private void deleteOldFiles() {
        File f = new File(dataSet+"/sites.txt");
        if (f.exists() && !f.isDirectory()) {
            f.delete();
        }
        f = new File(dataSet + "/sequences.fa");
        if (f.exists() && !f.isDirectory()) {
            f.delete();
        }

    }
}
