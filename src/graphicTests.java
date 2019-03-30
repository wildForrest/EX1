import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;

public class graphicTests {
	public static void main(String args[]){
		hashTest();
	   
	  

}
	 public static void readAndWriteAnImage(){
		 int width = 340;    //width of the image
		    int height = 940;   //height of the image
		    BufferedImage image = null;
		    File f = null;

		    //read image
		    try{
		      f = new File("strawberry.jpg"); //image file path
		      image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		      image = ImageIO.read(f);
		      System.out.println("Reading complete.");
		    }catch(IOException e){
		      System.out.println("Error: "+e);
		    }
		    
		  //write image
		    try{
		      f = new File("strawberry2.jpg");  //output file path
		      ImageIO.write(image, "jpg", f);
		      System.out.println("Writing complete.");
		    }catch(IOException e){
		      System.out.println("Error: "+e);
		    }
		    
	    }
	 public static void hashTest(){
		 HashMap<Integer,HashMap<Integer,Double>> d = new HashMap<>();
		 dictionaryPut(d,0,0,5.5);
		 System.out.println(dictionaryContainsKey(d,1,0));
	 }
	 
	 public static void dictionaryPut(HashMap<Integer,HashMap<Integer,Double>> d,
			 int i, int j,double val){
		 if(d.get(i)==null){
			 d.put(i, new HashMap<Integer,Double>());
		 }
		 d.get(i).put(j, val);
	 }
	 public static double dictionaryGet(HashMap<Integer,HashMap<Integer,Double>> d,
			 int i, int j){
		 return d.get(i).get(j);
	 }
	 public static boolean dictionaryContainsKey(HashMap<Integer,HashMap<Integer,Double>> d,
			 int i, int j){
		 if(d.containsKey(i)){
			 if(d.get(i).containsKey(j)){
				 return true;
			 }
		 }
		 return false;
		 
		 
	 }
	
}
