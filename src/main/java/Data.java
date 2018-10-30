import java.util.Random;

public class Data {

    public int[][] randomTable(){
        int[][] aaChanges = new int[64][20];
        Random rand = new Random();
        for (int i = 0; i < aaChanges.length; i++) {
            for (int j = 0; j < aaChanges[i].length; j++) {
                aaChanges[i][j] = rand.nextInt(500);
            }
        }
        return aaChanges;
    }



    public void codeTableFieldOverTime(int codon, int aa){

    }

}
