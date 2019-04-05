import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;

public class graphicTests {
	public static void main(String args[]){
		imageRepresentationTest();
	   
	  

}
	public static void imageRepresentationTest(){
		try{
		    BufferedImage img = null;
			img = new BufferedImage(100, 200, BufferedImage.TYPE_INT_ARGB);
			System.out.println("image Height:"+img.getHeight()+"image Width:"+img.getWidth());

			int  clr1   = img.getRGB(99, 200);
			
		      //File f = new File("images\\Empty.jpg");  //output file path
		      //ImageIO.write(img, "jpg", f);
		    }
		catch(Exception e){
		      System.out.println("Error: "+e);
		    }
	}
	public static void setPixel(BufferedImage img,int i, int j, int RGB,int A){
		/**author: Roee**/
		int p;
		p = (A<<24) | (RGB<<16) | (RGB<<8) | RGB;
        img.setRGB(i, j, p);	
	}

	public static double[] extractRGB(BufferedImage img,int x, int y){
		 int  clr1   = img.getRGB(x, y);
		 int a1 = (clr1>>24)&0xff;
		 int  r1   = (clr1 & 0x00ff0000) >> 16;
	     int  g1 = (clr1 & 0x0000ff00) >> 8;
	     int  b1  =  clr1 & 0x000000ff;
	     double[] rgb={r1,g1,b1,a1};
	     return rgb;
	     
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
